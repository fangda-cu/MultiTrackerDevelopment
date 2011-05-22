/**
 * \file BARodStepper.cc
 *
 * \author smith@cs.columbia.edu
 * \date 02/16/2010
 */

//#define KEEP_ONLY_SOME_RODS

#include <typeinfo>
#include "BARodStepper.hh"
#include "../../Threads/MultithreadedStepper.hh"
#include "../../Core/Timer.hh"
#include "../../Collisions/Collision.hh"

#include <iostream.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

namespace BASim
{

BARodStepper::BARodStepper(std::vector<ElasticRod*>& rods, std::vector<TriangleMesh*>& trimeshes,
        std::vector<ScriptingController*>& scripting_controllers, std::vector<RodTimeStepper*>& steppers, const double& dt,
        const double time, const int num_threads, const PerformanceTuningParameters perf_param) :
            m_num_dof(0),
            m_rods(rods),
            m_number_of_rods(m_rods.size()),
            m_triangle_meshes(trimeshes),
            m_scripting_controllers(scripting_controllers),
            m_steppers(steppers),
            m_dt(dt),
            m_base_indices(),
            m_base_triangle_indices(),
            m_edges(),
            m_faces(),
            m_vertex_radii(),
            m_edge_radii(),
            m_face_radii(),
            m_masses(),
            m_collision_immune(),
            m_respns_enbld(true),
            m_nan_enc(false),
            m_inf_enc(false),
            m_lt0_enc(false),
            m_gt0_enc(false),
            m_collision_detector(NULL),
            m_obj_start(-1),
            m_t(time),
            m_rod_labels(),
            m_simulationFailed(false),
            m_stopOnRodError(false),
            m_perf_param(perf_param),
            m_level(0),
            m_geodata(m_xn, m_vnphalf, m_vertex_radii, m_masses, m_collision_immune, m_obj_start,
                    m_perf_param.m_implicit_thickness, m_perf_param.m_implicit_stiffness)
{
    // For debugging purposes
#ifdef KEEP_ONLY_SOME_RODS
    std::set<int> keep_only;
    keep_only.insert(9);

    std::vector<ElasticRod*>::iterator rod = m_rods.begin();
    std::vector<RodTimeStepper*>::iterator stepper = m_steppers.begin();
    for (int i = 0; i < m_number_of_rods; i++)
    {
        if (keep_only.find(i) == keep_only.end())
        {
            for (int j = 0; j < (*rod)->nv(); j++)
            (*rod)->setVertex(j, 0 * ((*rod)->getVertex(j)));
            m_rods.erase(rod);
            m_steppers.erase(stepper);
        }
        else
        {
            rod++;
            stepper++;
        }
    }
    m_number_of_rods = m_rods.size();
    std::cerr << "Number of rods remaining: " << m_number_of_rods << std::endl;
#endif

    for (std::vector<RodTimeStepper*>::iterator stepper = m_steppers.begin();
    	 stepper != m_steppers.end();
     	 ++stepper) 
    {
       (*stepper)->setMaxIterations(perf_param.m_maximum_number_of_solver_iterations);
    }


#ifdef DEBUG
    for( int i = 0; i < (int) m_number_of_rods; ++i ) assert( m_rods[i] != NULL );
    for( int i = 0; i < (int) m_triangle_meshes.size(); ++i ) assert( m_triangle_meshes[i] != NULL );
    for( int i = 0; i < (int) m_steppers.size(); ++i ) assert( m_steppers[i] != NULL );
#endif
    assert(m_dt > 0.0);

    if (num_threads > 0)
    {
        m_num_threads = num_threads;
        std::cerr << "User-set number of threads = " << m_num_threads << std::endl;
    }
    else
    {
        m_num_threads = sysconf(_SC_NPROCESSORS_ONLN);
        std::cerr << "Default-set number of threads = " << m_num_threads << std::endl;
    }
#ifdef HAVE_OPENMP
    omp_set_num_threads(m_num_threads);
#endif
    // Update internal state, prepare for execution
    prepareForExecution();

#ifdef DEBUG
    // Number of degrees of freedom is non-negative multiple of 3 (3 coords per vertex)
    assert( m_num_dof >= 0 );
    assert( m_num_dof%3 == 0 );

    // Base indices refers to rods
    assert( m_base_indices.size() == m_number_of_rods );
    // Each base index is a non-negative multiple of 3
    for( int i = 0; i < (int) m_base_indices.size(); ++i ) assert( m_base_indices[i] >= 0 );
    for( int i = 0; i < (int) m_base_indices.size(); ++i ) assert( m_base_indices[i]%3 == 0 );
    // Each base index must be greater than last
    for( int i = 0; i < (int) m_base_indices.size()-1; ++i ) assert( m_base_indices[i] < m_base_indices[i+1] );

    // Base tirangle indices refers to triangles
    assert( m_base_triangle_indices.size() == m_triangle_meshes.size() );
    // Each base index is a non-negative multiple of 3
    for( int i = 0; i < (int) m_base_triangle_indices.size(); ++i ) assert( m_base_triangle_indices[i] >= 0 );
    for( int i = 0; i < (int) m_base_triangle_indices.size(); ++i ) assert( m_base_triangle_indices[i]%3 == 0 );
    // Each base index must be greater than last
    for( int i = 0; i < (int) m_base_triangle_indices.size()-1; ++i ) assert( m_base_triangle_indices[i] < m_base_triangle_indices[i+1] );

    // Check that we computed the proper start location of tirangle objects in the globale DOF array
    if( m_base_triangle_indices.size() > 0 ) assert( m_obj_start == (int) m_base_triangle_indices.front()/3 );

    // All edges and faces should contain valid vertices
    for( int i = 0; i < (int) m_edges.size(); ++i ) assert( m_edges[i].first >= 0 );
    for( int i = 0; i < (int) m_edges.size(); ++i ) assert( m_edges[i].second >= 0 );
    for( int i = 0; i < (int) m_edges.size(); ++i ) assert( m_edges[i].first < m_num_dof/3 );
    for( int i = 0; i < (int) m_edges.size(); ++i ) assert( m_edges[i].second < m_num_dof/3 );

    for( int i = 0; i < (int) m_faces.size(); ++i ) assert( m_faces[i].idx[0] >= 0 );
    for( int i = 0; i < (int) m_faces.size(); ++i ) assert( m_faces[i].idx[1] >= 0 );
    for( int i = 0; i < (int) m_faces.size(); ++i ) assert( m_faces[i].idx[2] >= 0 );
    for( int i = 0; i < (int) m_faces.size(); ++i ) assert( m_faces[i].idx[0] < m_num_dof/3 );
    for( int i = 0; i < (int) m_faces.size(); ++i ) assert( m_faces[i].idx[1] < m_num_dof/3 );
    for( int i = 0; i < (int) m_faces.size(); ++i ) assert( m_faces[i].idx[2] < m_num_dof/3 );

    // In our case, all face vertices should belong to triangles
    for( int i = 0; i < (int) m_faces.size(); ++i ) assert( m_faces[i].idx[0] >= m_obj_start );
    for( int i = 0; i < (int) m_faces.size(); ++i ) assert( m_faces[i].idx[1] >= m_obj_start );
    for( int i = 0; i < (int) m_faces.size(); ++i ) assert( m_faces[i].idx[2] >= m_obj_start );

    // Vertex radii must equal the number of verts!
    assert( (int) m_vertex_radii.size() == m_num_dof/3 );
    // Edge radii must equal the number of edges!
    assert( m_edge_radii.size() == m_edges.size() );
    // Face radii must equal the number of faces!
    assert( m_face_radii.size() == m_faces.size() );
    // All radii must be greater or equal to 0
    for( int i = 0; i < (int) m_vertex_radii.size(); ++i ) assert( m_vertex_radii[i] >= 0 );
    for( int i = 0; i < (int) m_edge_radii.size(); ++i ) assert( m_edge_radii[i] >= 0 );
    for( int i = 0; i < (int) m_face_radii.size(); ++i ) assert( m_face_radii[i] >= 0 );

    // In our case, face radii must be 0. All other radii must be positive.
    for( int i = 0; i < (int) m_face_radii.size(); ++i ) assert( m_face_radii[i] == 0 );
    for( int i = 0; i < (int) m_edge_radii.size(); ++i ) assert( m_edge_radii[i] < std::numeric_limits<double>::infinity() );

    // Number of masses must equal number of verts
    assert( (int) m_masses.size() == m_num_dof/3 );
    // Masses must be positive
    for( int i = 0; i < (int) m_masses.size(); ++i ) assert( m_masses[i] >= 0.0 );

    // Check that rod masses are positive doubles, that face masses are infs
    // TODO: Scripted verts get infinite mass, clean up this check later
    //for( int i = 0; i < m_obj_start; ++i ) assert( m_masses[i] < std::numeric_limits<double>::infinity() );
    for( int i = m_obj_start; i < m_num_dof/3; ++i ) assert( m_masses[i] == std::numeric_limits<double>::infinity() );

    // For each edge, ensure that both vertices are either rod edges or face edges
    for( int i = 0; i < (int) m_edges.size(); ++i )
    assert( (m_edges[i].first<m_obj_start && m_edges[i].second<m_obj_start) || (m_edges[i].first>=m_obj_start && m_edges[i].second>=m_obj_start) );

    // For each triangle, ensure that all vertices do indeed belong to a triangulated object
    for( int i = 0; i < (int) m_faces.size(); ++i )
    assert( m_faces[i].idx[0]>=m_obj_start && m_faces[i].idx[1]>=m_obj_start && m_faces[i].idx[2]>=m_obj_start );

    // TODO: Furhter, check that they all belong to same rod or triangle obj
#endif

#ifdef TIMING_ON
    IntStatTracker::getIntTracker("CONVERGENCE_FAILURES_PROPAGATED_TO_BARodStepper");
#endif

    /**
     *  Prepare backup structures.
     */
    m_startForces = new VecXd*[m_number_of_rods];
    for (int i = 0; i < m_number_of_rods; i++)
        m_startForces[i] = new VecXd(m_rods[i]->ndof());

    m_preCollisionForces = new VecXd*[m_number_of_rods];
    for (int i = 0; i < m_number_of_rods; i++)
        m_preCollisionForces[i] = new VecXd(m_rods[i]->ndof());

    m_endForces = new VecXd*[m_number_of_rods];
    for (int i = 0; i < m_number_of_rods; i++)
        m_endForces[i] = new VecXd(m_rods[i]->ndof());

    m_rodbackups.resize(m_number_of_rods);
    int i = 0;
    for (std::vector<ElasticRod*>::const_iterator rod = m_rods.begin(); rod != m_rods.end(); rod++)
        m_rodbackups[i++].resize(**rod);

    m_objbackups.resize(m_triangle_meshes.size());
    i = 0;
    for (std::vector<TriangleMesh*>::const_iterator mesh = m_triangle_meshes.begin(); mesh != m_triangle_meshes.end(); mesh++)
        m_objbackups[i++].resize(**mesh);

    // Initially all rods passed from Maya will be simulated
    for (int i = 0; i < m_number_of_rods; i++)
        m_simulated_rods.push_back(i);

}

BARodStepper::~BARodStepper()
{
    delete m_collision_detector;

    for (int i = 0; i < m_number_of_rods; i++)
    {
        delete m_startForces[i];
        delete m_preCollisionForces[i];
        delete m_endForces[i];
    }
    delete[] m_startForces;
    delete[] m_preCollisionForces;
    delete[] m_endForces;
}

// TODO: Check the indices here
void BARodStepper::prepareForExecution()
{
    delete m_collision_detector;
    m_collision_detector = NULL;

    for (int i = 0; i < m_number_of_rods; ++i)
    {
        assert(m_rods[i] != NULL);
	m_rods[i]->globalRodIndex = i;
    }

    //std::cout << "About to extract rod information" << std::endl;
    for (int i = 0; i < m_number_of_rods; ++i)
    {
        assert(m_rods[i] != NULL);
#ifdef DEBUG
        // Sanity checks for collision detection purposes
        ensureCircularCrossSection( *m_rods[i] );
        ensureNoCollisionsByDefault( *m_rods[i] );
#endif

        // Extract edges from the new rod
        for (int j = 0; j < m_rods[i]->nv() - 1; ++j)
        {
            m_edges.push_back(std::pair<int, int>(getNumVerts() + j, getNumVerts() + j + 1));
            assert(m_rods[i]->getRadiusScale() * m_rods[i]->radiusA(j) > 0.0);
            m_edge_radii.push_back(m_rods[i]->getRadiusScale() * m_rods[i]->radiusA(j));
        }
        assert(m_edges.size() == m_edge_radii.size());

        for (int j = 0; j < m_rods[i]->nv() - 1; ++j)
        {
            assert(m_rods[i]->getRadiusScale() * m_rods[i]->radiusA(j) > 0.0);
            // Radii associated with edges ... what to do if at a vertex with two radii? Average?
            m_vertex_radii.push_back(m_rods[i]->getRadiusScale() * m_rods[i]->radiusA(j));
        }
        m_vertex_radii.push_back(m_vertex_radii.back()); // < TODO: What the $^#! is this call?

        // Update vector that tracks the rod DOF in the system
        m_base_indices.push_back(getNumDof());

        // Extract masses from the new rod
        for (ElasticRod::vertex_iter itr = m_rods[i]->vertices_begin(); itr != m_rods[i]->vertices_end(); ++itr)
        {
            assert(m_rods[i]->getVertexMass(*itr, -1) > 0.0);
            m_masses.push_back(m_rods[i]->getVertexMass(*itr, -1));
        }

        // Update total number of DOF in the system
        m_num_dof += 3 * m_rods[i]->nv();

        assert((int) m_masses.size() == getNumVerts());
        assert((int) m_vertex_radii.size() == getNumVerts());
    }
    assert(m_number_of_rods == m_base_indices.size());
    //std::cout << "Extracted rod information: " << m_num_dof / 3 << " vertices" << std::endl;

    m_obj_start = m_base_indices.back() / 3 + m_rods.back()->nv();

    //std::cout << "About to extract tri mesh information" << std::endl;
    for (int i = 0; i < (int) m_triangle_meshes.size(); ++i)
    {
        assert(m_triangle_meshes[i] != NULL);

        // Extract faces from the tri_mesh
        for (TriangleMesh::face_iter fit = m_triangle_meshes[i]->faces_begin(); fit != m_triangle_meshes[i]->faces_end(); ++fit)
        {
            TriangularFace triface;
            int j = 0;
            for (TriangleMesh::FaceVertexIter fvit = m_triangle_meshes[i]->fv_iter(*fit); fvit; ++fvit, ++j)
            {
                assert(j >= 0);
                assert(j < 3);
                triface.idx[j] = getNumVerts() + (*fvit).idx();
            }
            m_faces.push_back(triface);
            m_face_radii.push_back(0.0);
        }
        assert(m_faces.size() == m_face_radii.size());
        //std::cout << "Finished extracting face stuff: " << i << std::endl;

        // Extract the vertex radii from the tri_mesh
        for (int j = 0; j < m_triangle_meshes[i]->nv(); ++j)
            m_vertex_radii.push_back(0.0);

        // Extract masses from the tri_mesh (just infinity for this object)
        for (int j = 0; j < m_triangle_meshes[i]->nv(); ++j)
            m_masses.push_back(std::numeric_limits<double>::infinity());

        // Add the mesh to the system
        //m_triangle_meshes.push_back(m_triangle_meshes[i]);

        // Update vector that tracks the ScriptedTriangleMesh DOF in the system
        m_base_triangle_indices.push_back(getNumDof());

        // Update total number of DOF in the system
        m_num_dof += 3 * m_triangle_meshes[i]->nv();

        assert((int) m_masses.size() == getNumVerts());
        assert((int) m_vertex_radii.size() == getNumVerts());
    }
    assert(m_base_triangle_indices.size() == m_triangle_meshes.size());
    //std::cout << "Extracted tri mesh information" << std::endl;

    // Resize the internal storage
    m_xn.resize(getNumDof());
    m_xn.setConstant(std::numeric_limits<double>::signaling_NaN());
    m_xnp1.resize(getNumDof());
    m_xnp1.setConstant(std::numeric_limits<double>::signaling_NaN());
    m_xdebug.resize(getNumDof());
    m_xdebug.setConstant(std::numeric_limits<double>::signaling_NaN());
    m_vnphalf.resize(getNumDof());
    m_vnphalf.setConstant(std::numeric_limits<double>::signaling_NaN());

    //std::cout << "About to extract positions" << std::endl;
    // Load positions for initial construction of the BVH
    RodSelectionType selected_rods;
    for (int i = 0; i < m_number_of_rods; i++)
        selected_rods.push_back(i);
    extractPositions(m_xn, selected_rods);
    extractVelocities(m_vnphalf, selected_rods);
    //std::cout << "Extracted positions" << std::endl;

    //  std::cout << "About to create CollisionDetector" << std::endl;
    m_collision_detector = new CollisionDetector(m_geodata, m_edges, m_faces, m_dt, m_perf_param.m_skipRodRodCollisions,
            m_num_threads);
    // std::cout << "Created CollisionDetector" << std::endl;

    // for (std::vector<double>::const_iterator i = m_masses.begin(); i != m_masses.end(); i++)
    //       m_collision_immune.push_back(*i == std::numeric_limits<double>::infinity());
    m_collision_immune.resize(getNumVerts());

    if (m_perf_param.m_enable_penalty_response)
        enableImplicitPenaltyImpulses();

}

bool BARodStepper::execute()
{
    std::cerr << "Executing time step " << m_t << std::endl;

    m_collision_detector->buildBVH();
    std::cerr << "BVH has been (re)built" << std::endl;

    if (!m_collision_disabled_rods.empty())
    {
        std::cerr << "The following rods were collision-disabled in the previous time step: ";
        for (std::set<int>::const_iterator rod = m_collision_disabled_rods.begin(); rod != m_collision_disabled_rods.end(); rod++)
            std::cerr << *rod << " ";
        std::cerr << std::endl;
    }
    m_collision_disabled_rods.clear();

    Timer::getTimer("BARodStepper::execute").start();
    bool do_adaptive = true;
    bool result;

    // std::cout << "BARodStepper::execute: listing scripted vertices... " << std::endl;
    int k = 0;
    for (int i = 0; i < m_number_of_rods; ++i)
    {
        // Extract masses from the new rod
        for (ElasticRod::vertex_iter itr = m_rods[i]->vertices_begin(); itr != m_rods[i]->vertices_end(); ++itr)
        {
            if (m_rods[i]->getBoundaryCondition()->isVertexScripted((*itr).idx()))
            {
                //       std::cout << "BARodStepper::execute: rod " << i << " vertex " << (*itr).idx() << " prescribed."
                //               << std::endl;
                m_masses[k++] = std::numeric_limits<double>::infinity();
            }
            else
            {
                assert(m_rods[i]->getVertexMass(*itr, -1) > 0.0);
                m_masses[k++] = m_rods[i]->getVertexMass(*itr, -1);
            }
        }
    }
    assert(k = m_masses.size());

    // Prepare the list initially containing all rods.
    RodSelectionType selected_rods = m_simulated_rods;

    if (do_adaptive)
    {
        assert(m_level == 0);
        result = adaptiveExecute(m_dt, selected_rods);
    }
    else
        result = nonAdaptiveExecute(m_dt, selected_rods);

    Timer::getTimer("BARodStepper::execute").stop();
    // Timer::report();

    std::cerr << std::endl;

    return result;
}

bool BARodStepper::nonAdaptiveExecute(double dt, RodSelectionType& selected_rods)
{
    setDt(dt);
    setTime(m_t + dt);
    //for (int i = 0; i < m_scripting_controllers.size(); ++i)
    //  m_scripting_controllers[i]->setTime(m_t);
    return step(selected_rods);
}

bool BARodStepper::adaptiveExecute(double dt, RodSelectionType& selected_rods)
{
    std::cout << "BARodStepper::adaptiveExecute starting at level " << m_level << " with m_t = " << m_t << " and dt = " << dt
            << std::endl;

    // Backup all selected rods
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        m_rodbackups[*rod].backupRod(*m_rods[*rod]);

    // Backup all objects
    for (int i = 0; i < (int) m_triangle_meshes.size(); ++i)
        m_objbackups[i].backupMesh(*m_triangle_meshes[i]);

    // Backup the current simulation time
    double time = m_t;

    // Set the desired timestep
    // std::cout << "Setting dt to: " << dt << std::endl;
    setDt(dt);
    // Advance the current time
    setTime(m_t + dt);

    // Attempt a full time step
    step(selected_rods);
    if (m_simulationFailed)
    {
        std::cout << "t = " << m_t << ": simulation failed and is now stopped" << std::endl;
        return true;
    }
    if (selected_rods.empty())
    {
        std::cout << "t = " << m_t << ": adaptiveExecute has simulated (or killed) all rods" << std::endl;
        // Success!
        return true;
    }

    std::cout << "t = " << m_t << ": adaptiveExecute left " << selected_rods.size() << " rods for substepping" << std::endl;

    // Otherwise do two half time steps
    // std::cout << "Adaptive stepping in Bridson stepper" << std::endl;
    // std::cout << "Number of rods remaining: " << selected_rods.size() << std::endl;

    m_level++;

    // Restore all rods that remained selected after the step
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        m_rodbackups[*rod].restoreRod(*m_rods[*rod]);

    // Restore all objects
    for (int i = 0; i < (int) m_triangle_meshes.size(); ++i)
        m_objbackups[i].restoreMesh(*m_triangle_meshes[i]);

    // Restore the time
    setTime(time);

    // Back up rod selection for time step 2
    RodSelectionType selected_rods_2 = selected_rods;

    std::cout << "t = " << m_t << " selected_rods: adaptiveExecute substepping (part 1) " << selected_rods.size() << " rods"
            << std::endl;

    // Otherwise attempt two steps of half length
    bool first_success = adaptiveExecute(0.5 * dt, selected_rods);
    if (!first_success)
    {
        setDt(dt);
        return false;
    }

    std::cout << "t = " << m_t << " selected_rods: adaptiveExecute substepping (part 2) " << selected_rods_2.size() << " rods"
            << std::endl;

    // Remove from the rod selection any one that might have been removed from the simulation list on the first time step
    for (RodSelectionType::iterator rod = selected_rods_2.begin(); rod != selected_rods_2.end(); rod++)
        if (find(m_simulated_rods.begin(), m_simulated_rods.end(), *rod) == m_simulated_rods.end())
        {
            std::cerr << "Erasing from second time step rod number " << *rod << std::endl;
            selected_rods_2.erase(rod--);
        }

    bool second_success = adaptiveExecute(0.5 * dt, selected_rods_2);

    if (!second_success)
    {
        setDt(dt);
        return false;
    }

    // std::cout << "Finished two adaptive steps" << std::endl;
    setDt(dt);
    m_level--;

    return first_success && second_success;
}

bool BARodStepper::step(RodSelectionType& selected_rods)
{
    if (m_simulationFailed)
        return true;

    std::cout << "t = " << m_t << ": BARodStepper::step() begins with " << selected_rods.size() << " rods" << std::endl;

    assert(m_edges.size() == m_edge_radii.size());
    assert((int) m_masses.size() == m_xn.size() / 3);
    assert(m_xn.size() == m_xnp1.size());
    assert(m_xn.size() == m_vnphalf.size());
    assert(m_rod_labels.size() == 0 || m_rod_labels.size() == m_number_of_rods);

    // Sanity check to ensure rods are not "internally colliding" because radius is bigger than edge length
#ifdef DEBUG
    for( int i = 0; i < (int) m_number_of_rods; ++i ) ensureNoCollisionsByDefault( *m_rods[i] );
#endif
    // Sanity check to ensure different parts of sim have same time/timetep
#ifdef DEBUG
    for( int i = 0; i < (int) m_scripting_controllers.size(); ++i ) assert( m_scripting_controllers[i]->getTime() == m_t );
    for( int i = 0; i < (int) m_scripting_controllers.size(); ++i ) assert( m_scripting_controllers[i]->getDt() == m_dt );
    for( int i = 0; i < (int) m_steppers.size(); ++i ) assert( m_steppers[i]->getTimeStep() == m_dt );
    for( int i = 0; i < (int) m_number_of_rods; ++i ) assert( m_rods[i]->getTimeStep() == m_dt );
#endif

    // std::cerr << "This step will treat " << selected_rods.size() << " remaining rod" << (selected_rods.size() > 1 ? "s" : "")
    //        << std::endl;

    // Prepare list of steppers to be executed.
    std::vector<RodTimeStepper*> selected_steppers;
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        selected_steppers.push_back(m_steppers[*rod]);

    if (m_perf_param.m_enable_explosion_detection)
        computeForces(m_startForces, selected_rods);

    // Save the pre-timestep positions
    extractPositions(m_xn, selected_rods);

    //  std::cout << "Pre-timestep positions: " << m_xn << std::endl;

    // Track whether or not this solve succeeds entirely
    bool dependable_solve = true;
    // Determine which vertex are to be considered collision-immune for this step
    computeImmunity(selected_rods);

    // Step rods forward ignoring collisions
    START_TIMER("BARodStepperDynamics");

    // Step scripted objects forward, set boundary conditions
    for (std::vector<ScriptingController*>::const_iterator scripting_controller = m_scripting_controllers.begin(); scripting_controller
            != m_scripting_controllers.end(); scripting_controller++)
        (*scripting_controller)->execute();

    // Jungseock's implicit penalty
    std::list<Collision*> penalty_collisions;
    if (m_perf_param.m_enable_penalty_response)
    {
        // Clear existing penalties
        for (std::vector<RodPenaltyForce*>::const_iterator penalty_force = m_implicit_pnlty_forces.begin(); penalty_force
                != m_implicit_pnlty_forces.end(); penalty_force++)
            (*penalty_force)->clearPenaltyForces();

        executeImplicitPenaltyResponse(penalty_collisions, selected_rods);
    }

#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < selected_steppers.size(); i++)
    {
        dependable_solve = dependable_solve && selected_steppers[i]->execute();
    }

    std::cerr << "Dynamic step is " << (dependable_solve ? "" : "not ") << "entirely dependable!" << std::endl;


    /*
     // Launch num_threads threads which will execute all elements of m_steppers.
     MultithreadedStepper<std::list<RodTimeStepper*> > multithreaded_stepper(selected_steppers, 4);//m_num_threads); // FIXME
     if (!multithreaded_stepper.Execute()) // if at least one of the steppers has not solved
     {
     dependable_solve = false;
     m_simulationFailed = true;
     std::cout << "Dynamic step is not entirely dependable!" << std::endl;
     }
     */

    // Clean up penalty collisions list
    for (std::list<Collision*>::iterator i = penalty_collisions.begin(); i != penalty_collisions.end(); i++)
        delete *i;

    STOP_TIMER("BARodStepperDynamics");

    // If we do rod-rod collisions (meaning no selective adaptivity) and global dependability failed, we might as well stop here.
    if (!m_perf_param.m_skipRodRodCollisions && !dependable_solve)
    {
        std::cout << "t = " << m_t << " selected_rods: step() failed (due to rod-rod) for " << selected_rods.size() << " rods"
                << std::endl;
        return false;
    }

    // Post time step position
    extractPositions(m_xnp1, selected_rods);


    if (m_perf_param.m_enable_explosion_detection)
        computeForces(m_preCollisionForces, selected_rods);

    // Average velocity over the timestep just completed
    m_vnphalf = (m_xnp1 - m_xn) / m_dt;


    for (int i = 0; i < selected_steppers.size(); i++)
    {
	applyInextensibilityVelocityFilter(selected_steppers[i]->getRod()->globalRodIndex);
    }

    for (int rodidx = 0; rodidx < m_rods.size(); ++rodidx)
    {
        std::cout << "Edge lengths after step() calls inextensibility(): ";
        int rodbase = m_base_indices[rodidx];
 
        for (int i = 0; i < m_rods[rodidx]->nv() - 1; ++i)
        {
	    Vec3d x0 = m_xn.segment<3> (rodbase + 3 * i);
	    Vec3d v0 = m_vnphalf.segment<3> (rodbase + 3 * i);
	    Vec3d x1 = m_xn.segment<3> (rodbase + 3 * i + 3);
	    Vec3d v1 = m_vnphalf.segment<3> (rodbase + 3 * i + 3);
	    Vec3d x0N = x0 + m_dt * v0;
	    Vec3d x1N = x1 + m_dt * v1;
	    Vec3d eN = (x1N - x0N);
	    std::cout << " " << eN.norm();
	}  
	std::cout << std::endl;
    }

    // Mark invalid rods as entirely collision-immune, so we don't waste time on colliding them.

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        if (!m_steppers[*rod]->HasSolved())
        {
            // std::cerr << "Rod number " << *rod << " failed to solve" << std::endl;
            for (int j = 0; j < m_rods[*rod]->nv(); ++j)
                m_collision_immune[m_base_indices[*rod] / 3 + j] = true;
        }

    //if( m_pnlty_enbld ) executePenaltyResponse();
    START_TIMER("BARodStepperResponse");
    std::cerr << "Starting collision response" << std::endl;
    bool all_collisions_succeeded = true;
    std::vector<bool> failed_collisions_rods(m_number_of_rods);
    if (m_perf_param.m_maximum_number_of_collisions_iterations > 0)
    {
        if (!executeIterativeInelasticImpulseResponse(failed_collisions_rods))
        {
            std::cout << "Some collision responses failed!" << std::endl;
            all_collisions_succeeded = false;
        }
    }
    std::cerr << "Finished collision response" << std::endl;
    STOP_TIMER("BARodStepperResponse");

    // Store the response part for visualization
    m_vnresp = m_vnphalf - (m_xnp1 - m_xn) / m_dt;

    // Compute final positions from corrected velocities
    m_xnp1 = m_xn + m_dt * m_vnphalf;

    // Ensure boundary conditions respected by corrected positions


#ifdef DEBUG
    // For each selected rod
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
    {
        RodBoundaryCondition* boundary = m_rods[*rod]->getBoundaryCondition();
        int rodbase = m_base_indices[*rod];

        // For each vertex of the current rod
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
        {
            // If that vertex has a prescribed position
            if (boundary->isVertexScripted(j))
            {
                //std::cout << "BridsonTimeStepper is calling RodBoundaryCondition at m_t = " << m_t << std::endl;
                Vec3d desiredposition = boundary->getDesiredVertexPosition(j, m_t);
                Vec3d actualvalue = m_xnp1.segment<3> (rodbase + 3 * j);
                assert(approxEq(desiredposition, actualvalue, 1.0e-6));
            }
        }
    }
#endif

    // Copy new positions and velocities back to rods
    restorePositions(m_xnp1, selected_rods);
    restoreVelocities(m_vnphalf, selected_rods);
    // Also copy response velocity to rods (for visualisation purposes only)
    restoreResponses(m_vnresp, selected_rods);

    // Update frames and such in the rod (Is this correct? Will this do some extra stuff?)
    for (RodSelectionType::const_iterator selected_rod = selected_rods.begin(); selected_rod != selected_rods.end(); selected_rod++)
        m_rods[*selected_rod]->updateProperties();

    // Sanity check to ensure rod's internal state is consistent
#ifdef DEBUG
    for( int i = 0; i < (int) m_number_of_rods; ++i ) m_rods[i]->verifyProperties();
#endif

    // Post time step position
    //extractPositions(m_xdebug);

    //std::cout << "Post-timestep positions, again: " << m_xdebug << std::endl;

    // Explosion detection
    if (m_perf_param.m_enable_explosion_detection)
        computeForces(m_endForces, selected_rods);
    bool explosions_detected = false;
    std::vector<bool> exploding_rods(m_number_of_rods);
    if (m_perf_param.m_enable_explosion_detection)
        explosions_detected = checkExplosions(exploding_rods, failed_collisions_rods, selected_rods);

    // Decide whether to substep or kill some rods
    for (RodSelectionType::iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
    {
        if (m_steppers[*rod]->HasSolved() && !exploding_rods[*rod] && !failed_collisions_rods[*rod])
        {
            selected_rods.erase(rod--);
            continue;
        }
        if (!m_steppers[*rod]->HasSolved())
        {
            if (m_level >= m_perf_param.m_max_number_of_substeps_for_solver)
            {
                switch (m_perf_param.m_in_case_of_solver_failure)
                {
                case PerformanceTuningParameters::IgnoreError:
                    // Nothing
                    break;
                case PerformanceTuningParameters::KillTheRod:
                    std::cerr << "Solver Failure: ";
                    killTheRod(*rod);
                    break;
                case PerformanceTuningParameters::HaltSimulation:
                    m_simulationFailed = true;
                    break;
                }
                selected_rods.erase(rod--);
                continue;
            }
        }
        if (exploding_rods[*rod])
        {
            if (m_level >= m_perf_param.m_max_number_of_substeps_for_explosion)
            {
                switch (m_perf_param.m_in_case_of_explosion_failure)
                {
                case PerformanceTuningParameters::IgnoreError:
                    // Nothing
                    break;
                case PerformanceTuningParameters::KillTheRod:
                    std::cerr << "Explosion Failure: ";
                    killTheRod(*rod);
                    break;
                case PerformanceTuningParameters::HaltSimulation:
                    m_simulationFailed = true;
                    break;
                }
                selected_rods.erase(rod--);
                continue;
            }
        }
        if (failed_collisions_rods[*rod])
        {
            if (m_level >= m_perf_param.m_max_number_of_substeps_for_collision)
            {
                switch (m_perf_param.m_in_case_of_collision_failure)
                {
                case PerformanceTuningParameters::IgnoreError:
                    // Nothing
                    break;
                case PerformanceTuningParameters::KillTheRod:
                    std::cerr << "Collision Failure: ";
                    killTheRod(*rod);
                    break;
                case PerformanceTuningParameters::HaltSimulation:
                    m_simulationFailed = true;
                    break;
                }
                selected_rods.erase(rod--);
                continue;
            }
        }
    }

    for (int rodidx = 0; rodidx < m_rods.size(); ++rodidx)
    {
        std::cout << "Edge lengths after step(): ";
        int rodbase = m_base_indices[rodidx];
 
        for (int i = 0; i < m_rods[rodidx]->nv() - 1; ++i)
        {
	    Vec3d x0 = m_xn.segment<3> (rodbase + 3 * i);
	    Vec3d v0 = m_vnphalf.segment<3> (rodbase + 3 * i);
	    Vec3d x1 = m_xn.segment<3> (rodbase + 3 * i + 3);
	    Vec3d v1 = m_vnphalf.segment<3> (rodbase + 3 * i + 3);
	    Vec3d x0N = x0 + m_dt * v0;
	    Vec3d x1N = x1 + m_dt * v1;
	    Vec3d eN = (x1N - x0N);
	    std::cout << " " << eN.norm();
	}  
	std::cout << std::endl;
    }

    bool all_rods_are_ok = dependable_solve && all_collisions_succeeded && !explosions_detected;
    std::cout << "This step is " << (all_rods_are_ok ? "" : "\033[31;1mNOT\033[m ") << "dependable." << std::endl;

    return all_rods_are_ok;
}

/**
 * Extracting/Restoring
 */
void BARodStepper::extractPositions(VecXd& positions, const RodSelectionType& selected_rods) const
{
    assert(m_number_of_rods == m_base_indices.size());
    assert(getNumDof() == positions.size());

    if (getNumDof() == 0)
        return;

#ifdef DEBUG
    positions.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
        {
            assert(m_base_indices[*rod] + 3 * j + 2 < positions.size());
            positions.segment<3> (m_base_indices[*rod] + 3 * j) = m_rods[*rod]->getVertex(j);
        }

    assert(m_triangle_meshes.size() == m_base_triangle_indices.size());

    //    std::cerr << "positions.size() = " << positions.size() << std::endl;
    //    std::cerr << "m_base_triangle_indices.size() = " << m_base_triangle_indices.size() << std::endl;

    for (int i = 0; i < (int) m_triangle_meshes.size(); ++i)
    {
        int j = 0;
        for (TriangleMesh::vertex_iter vit = m_triangle_meshes[i]->vertices_begin(); vit
                != m_triangle_meshes[i]->vertices_end(); ++vit, ++j)
        {
            assert(m_base_triangle_indices[i] + 3 * j + 2 < positions.size());
            positions.segment<3> (m_base_triangle_indices[i] + 3 * j) = m_triangle_meshes[i]->getVertex(*vit);
        }
    }

    //    assert((positions.cwise() == positions).all());

    // Ensure boundary conditions loaded properly
#ifdef DEBUG
    // For each rod
    for( int i = 0; i < (int) m_number_of_rods; ++i )
    {
        RodBoundaryCondition* boundary = m_rods[i]->getBoundaryCondition();
        int rodbase = m_base_indices[i];

        // For each vertex of the current rod, if that vertex has a prescribed position
        for( int j = 0; j < m_rods[i]->nv(); ++j ) if( boundary->isVertexScripted(j) )
        {
            std::cout << "BridsonTimeStepper is calling RodBoundaryCondition at m_t = " << m_t << std::endl;
            Vec3d desiredposition = boundary->getDesiredVertexPosition(j, m_t);
            Vec3d actualvalue = positions.segment<3>(rodbase+3*j);
            assert( approxEq(desiredposition, actualvalue, 1.0e-6) );
        }
    }
#endif
}

void BARodStepper::extractVelocities(VecXd& velocities, const RodSelectionType& selected_rods) const
{
    assert(m_number_of_rods == m_base_indices.size());
    assert(getNumDof() == velocities.size());

    if (getNumDof() == 0)
        return;

#ifdef DEBUG
    velocities.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
        {
            assert(m_base_indices[*rod] + 3 * j + 2 < velocities.size());
            velocities.segment<3> (m_base_indices[*rod] + 3 * j) = m_rods[*rod]->getVelocity(j);
        }

    assert(m_triangle_meshes.size() == m_base_triangle_indices.size());
    for (int i = 0; i < (int) m_triangle_meshes.size(); ++i)
    {
        int j = 0;
        for (TriangleMesh::vertex_iter vit = m_triangle_meshes[i]->vertices_begin(); vit
                != m_triangle_meshes[i]->vertices_end(); ++vit, ++j)
        {
            assert(m_base_triangle_indices[i] + 3 * j + 2 < velocities.size());
            velocities.segment<3> (m_base_triangle_indices[i] + 3 * j) = Vec3d::Zero();
        }
    }

    //    assert((velocities.cwise() == velocities).all());
}

void BARodStepper::restorePositions(const VecXd& positions, const RodSelectionType& selected_rods)
{
    assert(m_number_of_rods == m_base_indices.size());

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
            if (!m_rods[*rod]->getBoundaryCondition()->isVertexScripted(j))
                m_rods[*rod]->setVertex(j, positions.segment<3> (m_base_indices[*rod] + 3 * j));
}

void BARodStepper::restoreVelocities(const VecXd& velocities, const RodSelectionType& selected_rods)
{
    assert(m_number_of_rods == m_base_indices.size());

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
            if (!m_rods[*rod]->getBoundaryCondition()->isVertexScripted(j))
                m_rods[*rod]->setVelocity(j, velocities.segment<3> (m_base_indices[*rod] + 3 * j));
}

void BARodStepper::restoreResponses(const VecXd& responses, const RodSelectionType& selected_rods)
{
    assert(m_number_of_rods == m_base_indices.size());

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
            if (!m_rods[*rod]->getBoundaryCondition()->isVertexScripted(j))
                m_rods[*rod]->setResponse(j, responses.segment<3> (m_base_indices[*rod] + 3 * j));
}

/**
 * Enabling/disabling procedures
 */
void BARodStepper::enableImplicitPenaltyImpulses()
{
    m_perf_param.m_enable_penalty_response = true;
    for (int i = 0; i < (int) m_number_of_rods; i++)
    {
        RodPenaltyForce *pnlty = new RodPenaltyForce();
        m_implicit_pnlty_forces.push_back(pnlty);
        m_steppers[i]->addExternalForce(pnlty);
    }

    std::cerr << "Implicit penalty response is now enabled" << std::endl;

}

void BARodStepper::disableImplicitPenaltyImpulses()
{
    m_perf_param.m_enable_penalty_response = false;

    for (int i = 0; i < (int) m_number_of_rods; i++)
    {
        std::vector<RodExternalForce*>& forces = m_steppers[i]->getExternalForces();
        for (int j = 0; j < (int) forces.size(); j++)
        {
            RodPenaltyForce* rod_penalty_force = dynamic_cast<RodPenaltyForce*> (forces[j]);
            if (rod_penalty_force)
            {
                forces.erase(forces.begin() + j);
                break;
            }
        }
    }
    m_implicit_pnlty_forces.clear();
}

void BARodStepper::enableResponse()
{
    m_respns_enbld = true;
}

void BARodStepper::disableResponse()
{
    m_respns_enbld = false;
}

/*
 void BARodStepper::enableIterativeInelasticImpulses()
 {
 m_itrv_inlstc_enbld = true;
 }

 void BARodStepper::disableIterativeInelasticImpulses()
 {
 m_itrv_inlstc_enbld = false;
 }
 */

void BARodStepper::setNumInelasticIterations(const int& num_itr)
{
    assert(num_itr >= 0);
    m_perf_param.m_maximum_number_of_collisions_iterations = num_itr;
}

void BARodStepper::computeImmunity(const RodSelectionType& selected_rods)
{
    // Initially, the rods not treated in this step are collision-immune
    for (std::vector<bool>::iterator i = m_collision_immune.begin(); i != m_collision_immune.begin() + m_obj_start; i++)
        *i = true;
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        if (m_collision_disabled_rods.find(*rod) == m_collision_disabled_rods.end()) // If the rod is not in the disabled set
            for (int j = 0; j < m_rods[*rod]->nv(); ++j)
                m_collision_immune[m_base_indices[*rod] / 3 + j] = false;

    // Find the initial vertex-face intersections and mark them collision-immune
    std::list<Collision*> collisions;
    m_collision_detector->getCollisions(collisions, EdgeFace);
    while (!collisions.empty())
    {
        EdgeFaceIntersection* intersection = dynamic_cast<EdgeFaceIntersection*> (collisions.front());
        collisions.pop_front();

        if (intersection)
        {
            m_collision_immune[intersection->v0] = true;
            // std::cerr << "BARodStepper::step: Vertex " << intersection->v0 << " has been marked collision-immune\n";
        }
    }
}

/**
 * Utilities
 */
double BARodStepper::computeTotalForceNorm() const
{
    double totalforcenormsqr = 0.0;
    for (int i = 0; i < (int) m_number_of_rods; ++i)
    {
        VecXd force(m_rods[i]->ndof());
        force.setZero();
        m_rods[i]->computeForces(force);
        totalforcenormsqr += force.squaredNorm();
    }
    return sqrt(totalforcenormsqr);
}

void BARodStepper::setDt(double dt)
{
    assert(dt > 0.0);
    m_dt = dt;

    // Set the timestep for the rod controllers
    for (int i = 0; i < (int) m_steppers.size(); ++i)
        m_steppers[i]->setTimeStep(dt);

    // Set the timestep for the scripted object controllers
    for (int i = 0; i < (int) m_scripting_controllers.size(); ++i)
        m_scripting_controllers[i]->setDt(dt);
}

void BARodStepper::setTime(double time)
{
    m_t = time;
    // std::cout << "settingTime in BARodStepper to be " << m_t << std::endl;

    // Set the time for the rod controllers
    for (int i = 0; i < (int) m_steppers.size(); ++i)
        m_steppers[i]->setTime(m_t);

    // Set the time for the scripted object controllers
    for (int i = 0; i < (int) m_scripting_controllers.size(); ++i)
        m_scripting_controllers[i]->setTime(m_t);
}

double BARodStepper::getDt() const
{
    return m_dt;
}

double BARodStepper::getTime() const
{
    //std::cout << "BARodStepper::getTime() = " << m_t << std::endl;
    return m_t;
}

void BARodStepper::setRodLabels(const std::vector<std::string>& rod_labels)
{
    assert(rod_labels.size() == m_number_of_rods);
    m_rod_labels = rod_labels;
}

int BARodStepper::getContainingRod(int vert_idx) const
{
    assert(vert_idx >= 0);
    assert(vert_idx < getNumVerts());
    int rod = -1;

    int accm = 0;
    for (int i = 0; i < (int) m_number_of_rods; ++i)
    {
        accm += m_rods[i]->nv();
        if (vert_idx < accm)
        {
            rod = i;
            break;
        }
    }
    return rod;
}

bool BARodStepper::isRodVertex(int vert) const
{
    assert(vert >= 0);
    assert(vert < getNumVerts());

    // Is a vertex if index is less than start of object vertices in global array
    return vert < m_obj_start;
}

bool BARodStepper::vertexAndFaceShareVertex(const int& v, const int& f0, const int& f1, const int& f2) const
{
    return v == f0 || v == f1 || v == f2;
}

bool BARodStepper::vertexAndFaceShareVertex(const int& vertex, const int& face) const
{
    assert(face < (int) m_faces.size());
    assert(vertex >= 0);
    assert(vertex < getNumVerts());
    assert(m_faces[face].idx[0] >= 0);
    assert(m_faces[face].idx[0] < getNumVerts());
    assert(m_faces[face].idx[1] >= 0);
    assert(m_faces[face].idx[1] < getNumVerts());
    assert(m_faces[face].idx[2] >= 0);
    assert(m_faces[face].idx[2] < getNumVerts());

    return vertexAndFaceShareVertex(vertex, m_faces[face].idx[0], m_faces[face].idx[1], m_faces[face].idx[2]);
}

bool BARodStepper::isProperCollisionTime(double time)
{
    if (time != time)
    {
        if (!m_nan_enc)
            std::cerr
                    << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered NaN collision time from root finder. Supressing further messages of this type."
                    << std::endl;
        m_nan_enc = true;
        return false;
    }
    if (time == std::numeric_limits<double>::infinity())
    {
        if (!m_inf_enc)
            std::cerr
                    << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered INF collision time from root finder. Supressing further messages of this type."
                    << std::endl;
        m_inf_enc = true;
        return false;
    }
    if (time < 0.0)
    {
        if (!m_lt0_enc)
            std::cerr << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered scaled collision time " << time
                    << " less than 0.0. Supressing further messages of this type." << std::endl;
        m_lt0_enc = true;
        return false;
    }
    if (time > 1.0)
    {
        if (!m_gt0_enc)
            std::cerr << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered scaled collision time " << time
                    << " greater than 1.0. Supressing further messages of this type." << std::endl;
        m_gt0_enc = true;
        return false;
    }
    return true;
}

int BARodStepper::getNumDof() const
{
    assert(m_num_dof >= 0);
    return m_num_dof;
}

int BARodStepper::getNumVerts() const
{
    assert(m_num_dof % 3 == 0);
    return m_num_dof / 3;
}

// Ensures each rod edge has circular cross section.
void BARodStepper::ensureCircularCrossSection(const ElasticRod& rod) const
{
    // Ensure circular cross section
    for (int i = 0; i < (int) rod.ne(); ++i)
    {
        if (rod.getRadiusScale() * rod.radiusA(i) != rod.getRadiusScale() * rod.radiusB(i))
        {
            std::cerr
                    << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Contact currently not supported for non-circular cross sections. Assuming circular cross section."
                    << std::endl;
        }
    }
}

// Ensures each internal rod edge has length less than sum of neighbors' radii.
void BARodStepper::ensureNoCollisionsByDefault(const ElasticRod& rod) const
{
    // Ensure "non-attached" edges are not colliding by default
    for (int i = 1; i < (int) rod.ne() - 1; ++i)
    {
        double edgelen = rod.getEdge(i).norm();
        double radsum = rod.getRadiusScale() * rod.radiusA(i - 1) + rod.getRadiusScale() * rod.radiusA(i + 1);
        if (edgelen <= radsum)
        {
            std::cerr
                    << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Detected edges that collide by default. Instabilities may result."
                    << std::endl;
        }
    }
}

void BARodStepper::killTheRod(int rod) // TODO: remove the rod properly in Maya
{
    std::cerr << "Rod number " << rod << " is killed" << std::endl;

    for (int j = 0; j < m_rods[rod]->nv(); j++)
        m_rods[rod]->setVertex(j, 0 * m_rods[rod]->getVertex(j));
    m_simulated_rods.erase(find(m_simulated_rods.begin(), m_simulated_rods.end(), rod));
}

void BARodStepper::computeForces(VecXd ** Forces, const RodSelectionType& selected_rods)
{
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
    {
        Forces[*rod]->setZero();
        m_rods[*rod]->computeForces(*Forces[*rod]);
    }
}

}

