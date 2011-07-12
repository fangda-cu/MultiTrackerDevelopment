/**
 * \file BARodStepper.cc
 *
 * \author smith@cs.columbia.edu
 * \date 02/16/2010
 */

//#define KEEP_ONLY_SOME_RODS

#include <typeinfo>
#include "BAGroomingStepper.hh"
#include "../../Threads/MultithreadedStepper.hh"
#include "../../Core/Timer.hh"
#include "../../Collisions/Collision.hh"

#include "RodLevelSetForce.hh"
#include "RodClumpingForce.hh"

#include <iostream>
#include <fstream>
#include <sstream>

#include <weta/PantaRay/knn/kNN.hh>

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

namespace BASim
{
template<class T> std::ostream& operator<<(std::ostream& os, std::vector<T>& v)
{
    os << "{";
    for (size_t i = 0; i < v.size(); ++i)
    {
        os << v[i];
        if (i < v.size() - 1)
            os << ",";
    }
    os << "}";
    return os;
}

BAGroomingStepper::BAGroomingStepper(std::vector<ElasticRod*>& rods, std::vector<TriangleMesh*>& trimeshes,
        std::vector<ScriptingController*>& scripting_controllers, std::vector<GroomingTimeStepper*>& steppers,
        const double& dt, const double time, const int num_threads, const PerformanceTuningParameters perf_param,
        std::vector<LevelSet*>& levelSets, const int rods_per_clump) :
            m_num_dof(0),
            m_rods(rods),
            m_number_of_rods(m_rods.size()),
            m_triangle_meshes(trimeshes),
            m_scripting_controllers(scripting_controllers),
            m_steppers(steppers),
            m_dt(dt),
            m_base_dof_indices(),
            m_base_vtx_indices(),
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
            m_useKineticDamping(false),
            m_stopOnRodError(false),
            m_perf_param(perf_param),
            m_level(0),
            m_geodata(m_xn, m_vnphalf, m_vertex_radii, m_masses, m_collision_immune, m_obj_start,
                    m_perf_param.m_implicit_thickness, m_perf_param.m_implicit_stiffness), m_level_sets(levelSets),
            m_rods_per_clump(rods_per_clump)

//m_timers(NULL)
{
    m_levelset_forces.resize(m_level_sets.size());

    for (size_t i = 0; i < m_level_sets.size(); ++i)
    {
        m_levelset_forces[i] = new RodLevelSetForce(m_level_sets[i]);

        for (std::vector<GroomingTimeStepper*>::iterator stepper = m_steppers.begin(); stepper != m_steppers.end(); ++stepper)
        {
            (*stepper)->addExternalForce(m_levelset_forces[i]);
        }
    }

    for (std::vector<GroomingTimeStepper*>::iterator stepper = m_steppers.begin(); stepper != m_steppers.end(); ++stepper)
    {
        (*stepper)->setMaxIterations(m_perf_param.m_solver.m_max_iterations);
    }

    g_log = new TextLog(std::cerr, MsgInfo::kInfo, true);
    InfoStream(g_log, "") << "Started logging BAGroomingStepper\n";

#ifndef NDEBUG
    for (int i = 0; i < (int) m_number_of_rods; ++i)
        assert(m_rods[i] != NULL);
    for (int i = 0; i < (int) m_triangle_meshes.size(); ++i)
        assert(m_triangle_meshes[i] != NULL);

    // Do not check if any level sets are null as that may be valid as not all triangle meshes
    // may have a level set.
    //for( int i = 0; i < (int) m_level_sets.size(); ++i ) assert( m_level_sets[i] != NULL );
    if (m_level_sets.size() > 0)
    {
        // If there are any level sets then there has to be as one level set for every triangle mesh
        assert(m_level_sets.size() == m_triangle_meshes.size());
    }
    for (int i = 0; i < (int) m_steppers.size(); ++i)
        assert(m_steppers[i] != NULL);
#endif
    assert(m_dt > 0.0);

    if (num_threads > 0)
    {
        m_num_threads = num_threads;
        CopiousStream(g_log, "") << "User-set number of threads = " << m_num_threads << "\n";
    }
    else
    {
        m_num_threads = sysconf(_SC_NPROCESSORS_ONLN);
        CopiousStream(g_log, "") << "Default-set number of threads = " << m_num_threads << "\n";
    }
#ifdef HAVE_OPENMP
    omp_set_num_threads(m_num_threads);
#endif
    // Update internal state, prepare for execution
    prepareForExecution();

#ifndef NDEBUG
    // Number of degrees of freedom is non-negative multiple of 3 (3 coords per vertex)
    assert(m_num_dof >= 0);
    assert(m_num_dof % 3 == 0);

    // Base indices refers to rods
    assert(m_base_dof_indices.size() == m_number_of_rods);
    // Each base index is a non-negative multiple of 3
    for (int i = 0; i < (int) m_base_dof_indices.size(); ++i)
        assert(m_base_dof_indices[i] >= 0);
    for (int i = 0; i < (int) m_base_dof_indices.size(); ++i)
        assert(m_base_dof_indices[i] % 3 == 0);
    // Each base index must be greater than last
    for (int i = 0; i < (int) m_base_dof_indices.size() - 1; ++i)
        assert(m_base_dof_indices[i] < m_base_dof_indices[i + 1]);

    // Base tirangle indices refers to triangles
    assert(m_base_triangle_indices.size() == m_triangle_meshes.size());
    // Each base index is a non-negative multiple of 3
    for (int i = 0; i < (int) m_base_triangle_indices.size(); ++i)
        assert(m_base_triangle_indices[i] >= 0);
    for (int i = 0; i < (int) m_base_triangle_indices.size(); ++i)
        assert(m_base_triangle_indices[i] % 3 == 0);
    // Each base index must be greater than last
    for (int i = 0; i < (int) m_base_triangle_indices.size() - 1; ++i)
        assert(m_base_triangle_indices[i] < m_base_triangle_indices[i + 1]);

    // Check that we computed the proper start location of tirangle objects in the globale DOF array
    if (m_base_triangle_indices.size() > 0)
        assert(m_obj_start == (int) m_base_triangle_indices.front() / 3);

    // All edges and faces should contain valid vertices
    for (int i = 0; i < (int) m_edges.size(); ++i)
        assert(m_edges[i].first >= 0);
    for (int i = 0; i < (int) m_edges.size(); ++i)
        assert(m_edges[i].second >= 0);
    for (int i = 0; i < (int) m_edges.size(); ++i)
        assert(m_edges[i].first < m_num_dof / 3);
    for (int i = 0; i < (int) m_edges.size(); ++i)
        assert(m_edges[i].second < m_num_dof / 3);

    for (int i = 0; i < (int) m_faces.size(); ++i)
        assert(m_faces[i].idx[0] >= 0);
    for (int i = 0; i < (int) m_faces.size(); ++i)
        assert(m_faces[i].idx[1] >= 0);
    for (int i = 0; i < (int) m_faces.size(); ++i)
        assert(m_faces[i].idx[2] >= 0);
    for (int i = 0; i < (int) m_faces.size(); ++i)
        assert(m_faces[i].idx[0] < m_num_dof / 3);
    for (int i = 0; i < (int) m_faces.size(); ++i)
        assert(m_faces[i].idx[1] < m_num_dof / 3);
    for (int i = 0; i < (int) m_faces.size(); ++i)
        assert(m_faces[i].idx[2] < m_num_dof / 3);

    // In our case, all face vertices should belong to triangles
    for (int i = 0; i < (int) m_faces.size(); ++i)
        assert(m_faces[i].idx[0] >= m_obj_start);
    for (int i = 0; i < (int) m_faces.size(); ++i)
        assert(m_faces[i].idx[1] >= m_obj_start);
    for (int i = 0; i < (int) m_faces.size(); ++i)
        assert(m_faces[i].idx[2] >= m_obj_start);

    // Vertex radii must equal the number of verts!
    assert((int) m_vertex_radii.size() == m_num_dof / 3);
    // Edge radii must equal the number of edges!
    assert(m_edge_radii.size() == m_edges.size());
    // Face radii must equal the number of faces!
    assert(m_face_radii.size() == m_faces.size());
    // All radii must be greater or equal to 0
    for (int i = 0; i < (int) m_vertex_radii.size(); ++i)
        assert(m_vertex_radii[i] >= 0);
    for (int i = 0; i < (int) m_edge_radii.size(); ++i)
        assert(m_edge_radii[i] >= 0);
    for (int i = 0; i < (int) m_face_radii.size(); ++i)
        assert(m_face_radii[i] >= 0);

    // In our case, face radii must be 0. All other radii must be positive.
    for (int i = 0; i < (int) m_face_radii.size(); ++i)
        assert(m_face_radii[i] == 0);
    for (int i = 0; i < (int) m_edge_radii.size(); ++i)
        assert(m_edge_radii[i] < std::numeric_limits<double>::infinity());

    // Number of masses must equal number of verts
    assert((int) m_masses.size() == m_num_dof / 3);
    // Masses must be positive
    for (int i = 0; i < (int) m_masses.size(); ++i)
        assert(m_masses[i] >= 0.0);

    // Check that rod masses are positive doubles, that face masses are infs
    // TODO: Scripted verts get infinite mass, clean up this check later
    //for( int i = 0; i < m_obj_start; ++i ) assert( m_masses[i] < std::numeric_limits<double>::infinity() );
    for (int i = m_obj_start; i < m_num_dof / 3; ++i)
        assert(m_masses[i] == std::numeric_limits<double>::infinity());

    // For each edge, ensure that both vertices are either rod edges or face edges
    for (int i = 0; i < (int) m_edges.size(); ++i)
        assert(
                (m_edges[i].first < m_obj_start && m_edges[i].second < m_obj_start) || (m_edges[i].first >= m_obj_start
                        && m_edges[i].second >= m_obj_start));

    // For each triangle, ensure that all vertices do indeed belong to a triangulated object
    for (int i = 0; i < (int) m_faces.size(); ++i)
        assert(m_faces[i].idx[0] >= m_obj_start && m_faces[i].idx[1] >= m_obj_start && m_faces[i].idx[2] >= m_obj_start);

    // TODO: Furhter, check that they all belong to same rod or triangle obj
#endif

#ifdef TIMING_ON
    IntStatTracker::getIntTracker("CONVERGENCE_FAILURES_PROPAGATED_TO_BAGroomingStepper");
#endif

    /**
     *  Prepare backup structures.
     */
    m_startForces.resize(m_number_of_rods);
    for (int i = 0; i < m_number_of_rods; i++)
        m_startForces[i] = new VecXd(m_rods[i]->ndof());

    m_preDynamicForces.resize(m_number_of_rods);
    for (int i = 0; i < m_number_of_rods; i++)
        m_preDynamicForces[i] = new VecXd(m_rods[i]->ndof());

    m_preCollisionForces.resize(m_number_of_rods);
    for (int i = 0; i < m_number_of_rods; i++)
        m_preCollisionForces[i] = new VecXd(m_rods[i]->ndof());

    m_endForces.resize(m_number_of_rods);
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

    // For debugging purposes
#ifdef KEEP_ONLY_SOME_RODS
    WarningStream(g_log, "", MsgInfo::kOncePerMessage)
    << "WARNING: KEEP_ONLY_SOME_RODS: Simulating only a specified subset of rods!\n***********************************************************\n";
    std::set<int> keep_only;

    //keep_only.insert(33);
    //keep_only.insert(0);
    //keep_only.insert(6);
    keep_only.insert(2);
    //keep_only.insert(10);
    // keep_only.insert(189);
    //keep_only.insert(710);

    // Only the rods in the keep_only set are kept, the others are killed.
    for (int i = 0; i < m_number_of_rods; i++)
    {
        if (keep_only.find(i) == keep_only.end())
        {
            for (int j = 0; j < m_rods[i]->nv(); j++)
            {
                m_rods[i]->setVertex(j, 0 * m_rods[i]->getVertex(j));
            }
        }
        else
        {
            m_simulated_rods.push_back(i);
        }
    }
#else
    // Initially all rods passed from Maya will be simulated
    for (int i = 0; i < m_number_of_rods; i++)
        m_simulated_rods.push_back(i);
#endif
    m_killed_rods.clear();
    InfoStream(g_log, "") << "STEPPER DT " << m_dt << "\n";

    activateClumpingForce();

    CopiousStream(g_log, "") << "Finished BAGroomingStepper constructor\n";
}

BAGroomingStepper::~BAGroomingStepper()
{
    delete m_collision_detector;

    for (int i = 0; i < m_number_of_rods; i++)
    {
        delete m_startForces[i];
        delete m_preDynamicForces[i];
        delete m_preCollisionForces[i];
        delete m_endForces[i];
    }

    delete m_clumpingForce;
    delete g_log;
}

// TODO: Check the indices here
void BAGroomingStepper::prepareForExecution()
{
    delete m_collision_detector;
    m_collision_detector = NULL;

    for (int i = 0; i < m_number_of_rods; ++i)
    {
        assert(m_rods[i] != NULL);
        m_rods[i]->globalRodIndex = i;
        // std::cerr << "Address of rod Nr " << i << ": " << m_rods[i] << '\n';
    }

    CopiousStream(g_log, "") << "About to extract rod information\n";

    for (int i = 0; i < m_number_of_rods; ++i)
    {
        assert(m_rods[i] != NULL);
#ifndef NDEBUG
        // Sanity checks for collision detection purposes
        ensureCircularCrossSection(*m_rods[i]);
        ensureNoCollisionsByDefault(*m_rods[i]);
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
        m_base_dof_indices.push_back(getNumDof());
        m_base_vtx_indices.push_back(getNumDof() / 3);

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
    assert(m_number_of_rods == m_base_dof_indices.size());
    CopiousStream(g_log, "") << "Extracted rod information: " << m_num_dof / 3 << " vertices\n";

    m_obj_start = m_base_vtx_indices.back() + m_rods.back()->nv();

    CopiousStream(g_log, "") << "About to extract tri mesh information\n";
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
        CopiousStream(g_log, "") << "Finished extracting face stuff: " << i << "\n";

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
    CopiousStream(g_log, "") << "Extracted tri mesh information\n";

    // Resize the internal storage
    m_xn.resize(getNumDof());
    m_xn.setConstant(std::numeric_limits<double>::signaling_NaN());
    m_xnp1.resize(getNumDof());
    m_xnp1.setConstant(std::numeric_limits<double>::signaling_NaN());
    m_xdebug.resize(getNumDof());
    m_xdebug.setConstant(std::numeric_limits<double>::signaling_NaN());
    m_vnphalf.resize(getNumDof());
    m_vnphalf.setConstant(std::numeric_limits<double>::signaling_NaN());

    CopiousStream(g_log, "") << "About to extract positions\n";
    // Load positions for initial construction of the BVH
    RodSelectionType selected_rods;
    for (int i = 0; i < m_number_of_rods; i++)
        selected_rods.push_back(i);
    extractPositions(m_xn, selected_rods, 0.0);
    extractVelocities(m_vnphalf, selected_rods);
    CopiousStream(g_log, "") << "Extracted positions\n";

    CopiousStream(g_log, "") << "About to create collision detector\n";
    m_collision_detector = new CollisionDetectorType(m_geodata, m_edges, m_faces, m_dt, m_perf_param.m_skipRodRodCollisions,
            m_num_threads);
    CopiousStream(g_log, "") << "Created collision detector\n";

    m_collision_immune.resize(getNumVerts());

    if (m_perf_param.m_enable_penalty_response)
        enableImplicitPenaltyImpulses();

    // DEBUG
    m_total_solver_killed = m_total_collision_killed = m_total_explosion_killed = m_total_stretching_killed = 0;

    m_initialLengths.resize(m_number_of_rods);
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 1; j < m_rods[*rod]->nv(); j++)
            m_initialLengths[*rod] += (m_rods[*rod]->getVertex(j) - m_rods[*rod]->getVertex(j - 1)).norm();
}

bool BAGroomingStepper::execute()
{
    START_TIMER("BAGroomingStepper::execute")

    // Step scripted objects forward, set boundary conditions
    for (std::vector<ScriptingController*>::const_iterator scripting_controller = m_scripting_controllers.begin(); scripting_controller
            != m_scripting_controllers.end(); scripting_controller++)
        (*scripting_controller)->execute();

    for (size_t i = 0; i < m_levelset_forces.size(); ++i)
    {
        assert(m_levelset_forces[i]);
        m_levelset_forces[i]->setStiffness(m_perf_param.m_implicit_stiffness);
        m_levelset_forces[i]->setThickness(m_perf_param.m_implicit_thickness);
        m_levelset_forces[i]->setSubsampling(m_perf_param.m_levelset_subsampling);
    }

    // Prepare the list initially containing all rods.
    RodSelectionType selected_rods = m_simulated_rods;

    step(selected_rods);

    STOP_TIMER("BAGroomingStepper::execute")

    return true;
}

void BAGroomingStepper::step(RodSelectionType& selected_rods)
{
    START_TIMER("BAGroomingStepper::step")

    START_TIMER("BAGroomingStepper::step/setup");

    // Prepare list of steppers to be executed.
    std::vector<GroomingTimeStepper*> selected_steppers;
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        selected_steppers.push_back(m_steppers[*rod]);

    // Save the pre-timestep positions
    extractPositions(m_xn, selected_rods, m_t - m_dt);

    STOP_TIMER("BAGroomingStepper::step/setup");

    START_TIMER("BAGroomingStepper::step/steppers");

    TraceStream(g_log, "BAGroomingStepper") << "Stepping forward " << selected_steppers.size() << " rod(s).\n";

    bool dependable_solve = true;
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < selected_steppers.size(); i++)
    {
        GroomingTimeStepper* const stepper = selected_steppers[i];

        bool result = stepper->execute();
        if (!result)
            TraceStream(g_log, "") << stepper->getDiffEqSolver().getName() << " solver for rod "
                    << stepper->getRod()->globalRodIndex << " failed to converge\n";
    }

    STOP_TIMER("BAGroomingStepper::step/steppers");

    STOP_TIMER("BAGroomingStepper::step");
}

/*
 * Extracting/Restoring
 */
void BAGroomingStepper::extractPositions(VecXd& positions, const RodSelectionType& selected_rods, const double time) const
{
    assert(m_number_of_rods == m_base_dof_indices.size());
    assert(getNumDof() == positions.size());

    if (getNumDof() == 0)
        return;

#ifndef NDEBUG
    positions.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
        {
            assert(m_base_dof_indices[*rod] + 3 * j + 2 < positions.size());
            positions.segment<3> (m_base_dof_indices[*rod] + 3 * j) = m_rods[*rod]->getVertex(j);
        }

    assert(m_triangle_meshes.size() == m_base_triangle_indices.size());

    //    std::cerr << "positions.size() = " << positions.size() << '\n';
    //    std::cerr << "m_base_triangle_indices.size() = " << m_base_triangle_indices.size() << '\n';

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
#ifndef NDEBUG
    // For each rod in the selected list
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
    {
        RodBoundaryCondition* boundary = m_rods[*rod]->getBoundaryCondition();
        int rodbase = m_base_dof_indices[*rod];

        // For each vertex of the current rod, if that vertex has a prescribed position
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
            if (boundary->isVertexScripted(j))
            {
                //  std::cout << "BridsonTimeStepper is calling RodBoundaryCondition at m_t = " << m_t << '\n';
                Vec3d desiredposition = boundary->getDesiredVertexPosition(j, time);
                Vec3d actualvalue = positions.segment<3> (rodbase + 3 * j);
                assert(approxEq(desiredposition, actualvalue, 1.0e-6));
            }
    }
#endif
}

void BAGroomingStepper::extractVelocities(VecXd& velocities, const RodSelectionType& selected_rods) const
{
    assert(m_number_of_rods == m_base_dof_indices.size());
    assert(getNumDof() == velocities.size());

    if (getNumDof() == 0)
        return;

#ifndef NDEBUG
    velocities.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
        {
            assert(m_base_dof_indices[*rod] + 3 * j + 2 < velocities.size());
            velocities.segment<3> (m_base_dof_indices[*rod] + 3 * j) = m_rods[*rod]->getVelocity(j);
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

void BAGroomingStepper::restorePositions(const VecXd& positions, const RodSelectionType& selected_rods)
{
    assert(m_number_of_rods == m_base_dof_indices.size());

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
            if (!m_rods[*rod]->getBoundaryCondition()->isVertexScripted(j))
                m_rods[*rod]->setVertex(j, positions.segment<3> (m_base_dof_indices[*rod] + 3 * j));
}

void BAGroomingStepper::restoreVelocities(const VecXd& velocities, const RodSelectionType& selected_rods)
{
    assert(m_number_of_rods == m_base_dof_indices.size());

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
            if (!m_rods[*rod]->getBoundaryCondition()->isVertexScripted(j))
                m_rods[*rod]->setVelocity(j, velocities.segment<3> (m_base_dof_indices[*rod] + 3 * j));
}

void BAGroomingStepper::restoreResponses(const VecXd& responses, const RodSelectionType& selected_rods)
{
    assert(m_number_of_rods == m_base_dof_indices.size());

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
            if (!m_rods[*rod]->getBoundaryCondition()->isVertexScripted(j))
                m_rods[*rod]->setResponse(j, responses.segment<3> (m_base_dof_indices[*rod] + 3 * j));
}

/**
 * Enabling/disabling procedures
 */
void BAGroomingStepper::enableImplicitPenaltyImpulses()
{
    // m_perf_param.m_enable_penalty_response = true;
    // for (int i = 0; i < (int) m_number_of_rods; i++)
    // {
    //     RodPenaltyForce *pnlty = new RodPenaltyForce();
    //     m_implicit_pnlty_forces.push_back(pnlty);
    //     m_steppers[i]->addExternalForce(pnlty);
    // }

    // TraceStream(g_log, "") << "Implicit penalty response is now enabled\n";

}

void BAGroomingStepper::disableImplicitPenaltyImpulses()
{
    // m_perf_param.m_enable_penalty_response = false;

    // for (int i = 0; i < (int) m_number_of_rods; i++)
    // {
    //     std::vector<RodExternalForce*>& forces = m_steppers[i]->getExternalForces();
    //     for (int j = 0; j < (int) forces.size(); j++)
    //     {
    //         RodPenaltyForce* rod_penalty_force = dynamic_cast<RodPenaltyForce*> (forces[j]);
    //         if (rod_penalty_force)
    //         {
    //             forces.erase(forces.begin() + j);
    //             break;
    //         }
    //     }
    // }
    // m_implicit_pnlty_forces.clear();
}

void BAGroomingStepper::enableResponse()
{
    m_respns_enbld = true;
}

void BAGroomingStepper::disableResponse()
{
    m_respns_enbld = false;
}

void BAGroomingStepper::setNumInelasticIterations(const int& num_itr)
{
    assert(num_itr >= 0);
    m_perf_param.m_collision.m_max_iterations = num_itr;
}

void BAGroomingStepper::computeImmunity(const RodSelectionType& selected_rods)
{
    // Initially, the rods not treated in this step are collision-immune
    for (std::vector<bool>::iterator i = m_collision_immune.begin(); i != m_collision_immune.begin() + m_obj_start; i++)
        *i = true;
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        if (m_collision_disabled_rods.find(*rod) == m_collision_disabled_rods.end()) // If the rod is not in the disabled set
            for (int j = 0; j < m_rods[*rod]->nv(); ++j)
                m_collision_immune[m_base_vtx_indices[*rod] + j] = false;

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
            // std::cerr << "BAGroomingStepper::step: Vertex " << intersection->v0 << " has been marked collision-immune\n";
        }
    }
}

/**
 * Utilities
 */
double BAGroomingStepper::computeTotalForceNorm() const
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

void BAGroomingStepper::setDt(double dt)
{
    InfoStream(g_log, "") << "STEPPER DT " << dt << "\n";
    assert(dt > 0.0);
    m_dt = dt;

    // Set the timestep for the rod controllers
    for (int i = 0; i < (int) m_steppers.size(); ++i)
        m_steppers[i]->setTimeStep(dt);

    // Set the timestep for the scripted object controllers
    for (int i = 0; i < (int) m_scripting_controllers.size(); ++i)
        m_scripting_controllers[i]->setDt(dt);
}

void BAGroomingStepper::setTime(double time)
{
    m_t = time;
    // std::cout << "settingTime in BAGroomingStepper to be " << m_t << '\n';

    // Set the time for the rod controllers
    for (int i = 0; i < (int) m_steppers.size(); ++i)
        m_steppers[i]->setTime(m_t);

    // Set the time for the scripted object controllers
    for (int i = 0; i < (int) m_scripting_controllers.size(); ++i)
        m_scripting_controllers[i]->setTime(m_t);
}

double BAGroomingStepper::getDt() const
{
    return m_dt;
}

double BAGroomingStepper::getTime() const
{
    //std::cout << "BAGroomingStepper::getTime() = " << m_t << '\n';
    return m_t;
}

void BAGroomingStepper::skipRodRodCollisions(bool skipRodRodCollisions)
{
    TraceStream(g_log, "") << "Switching rod-rod collisions " << (skipRodRodCollisions ? "OFF" : "ON") << '\n';
    m_perf_param.m_skipRodRodCollisions = skipRodRodCollisions;

    if (m_collision_detector)
        m_collision_detector->setSkipRodRodCollisions(skipRodRodCollisions);
}

void BAGroomingStepper::setRodLabels(const std::vector<std::string>& rod_labels)
{
    assert(rod_labels.size() == m_number_of_rods);
    m_rod_labels = rod_labels;
}

int BAGroomingStepper::getContainingRod(int vert_idx) const
{
    assert(vert_idx >= 0);
    assert(vert_idx < getNumVerts());

    return upper_bound(m_base_vtx_indices.begin(), m_base_vtx_indices.end(), vert_idx) - m_base_vtx_indices.begin() - 1;
}

bool BAGroomingStepper::isRodVertex(int vert) const
{
    assert(vert >= 0);
    assert(vert < getNumVerts());

    // Is a vertex if index is less than start of object vertices in global array
    return vert < m_obj_start;
}

bool BAGroomingStepper::vertexAndFaceShareVertex(const int& v, const int& f0, const int& f1, const int& f2) const
{
    return v == f0 || v == f1 || v == f2;
}

bool BAGroomingStepper::vertexAndFaceShareVertex(const int& vertex, const int& face) const
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

bool BAGroomingStepper::isProperCollisionTime(double time)
{
    if (time != time)
    {
        if (!m_nan_enc)
            DebugStream(g_log, "")
                    << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered NaN collision time from root finder. Supressing further messages of this type."
                    << '\n';
        m_nan_enc = true;
        return false;
    }
    if (time == std::numeric_limits<double>::infinity())
    {
        if (!m_inf_enc)
            DebugStream(g_log, "")
                    << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered INF collision time from root finder. Supressing further messages of this type."
                    << '\n';
        m_inf_enc = true;
        return false;
    }
    if (time < 0.0)
    {
        if (!m_lt0_enc)
            DebugStream(g_log, "") << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered scaled collision time " << time
                    << " less than 0.0. Supressing further messages of this type.\n";
        m_lt0_enc = true;
        return false;
    }
    if (time > 1.0)
    {
        if (!m_gt0_enc)
            DebugStream(g_log, "") << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Encountered scaled collision time " << time
                    << " greater than 1.0. Supressing further messages of this type.\n";
        m_gt0_enc = true;
        return false;
    }
    return true;
}

int BAGroomingStepper::getNumDof() const
{
    assert(m_num_dof >= 0);
    return m_num_dof;
}

int BAGroomingStepper::getNumVerts() const
{
    assert(m_num_dof % 3 == 0);
    return m_num_dof / 3;
}

// Ensures each rod edge has circular cross section.
void BAGroomingStepper::ensureCircularCrossSection(const ElasticRod& rod) const
{
    // Ensure circular cross section
    for (int i = 0; i < (int) rod.ne(); ++i)
    {
        if (rod.getRadiusScale() * rod.radiusA(i) != rod.getRadiusScale() * rod.radiusB(i))
        {
            DebugStream(g_log, "")
                    << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Contact currently not supported for non-circular cross sections. Assuming circular cross section."
                    << '\n';
        }
    }
}

// Ensures each internal rod edge has length less than sum of neighbors' radii.
void BAGroomingStepper::ensureNoCollisionsByDefault(const ElasticRod& rod) const
{
    // Ensure "non-attached" edges are not colliding by default
    for (int i = 1; i < (int) rod.ne() - 1; ++i)
    {
        double edgelen = rod.getEdge(i).norm();
        double radsum = rod.getRadiusScale() * rod.radiusA(i - 1) + rod.getRadiusScale() * rod.radiusA(i + 1);
        if (edgelen <= radsum)
        {
            DebugStream(g_log, "")
                    << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Detected edges that collide by default. Instabilities may result."
                    << '\n';
        }
    }
}

void BAGroomingStepper::killTheRod(int rod) // TODO: remove the rod properly in Maya
{
    m_killed_rods.push_back(rod);
    for (int j = 0; j < m_rods[rod]->nv(); j++)
        m_rods[rod]->setVertex(j, 0 * m_rods[rod]->getVertex(j));
    m_simulated_rods.erase(find(m_simulated_rods.begin(), m_simulated_rods.end(), rod));

#ifndef KEEP_ONLY_SOME_RODS
    assert(m_simulated_rods.size() + m_killed_rods.size() == m_number_of_rods);
#endif
}

void BAGroomingStepper::computeForces(std::vector<VecXd*> Forces, const RodSelectionType& selected_rods)
{
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
    {
        Forces[*rod]->setZero();
        //m_rods[*rod]->computeForces(*Forces[*rod]);
        m_steppers[*rod]->evaluatePDot(*Forces[*rod]);
    }
}

void BAGroomingStepper::addRod(ElasticRod* rod, GroomingTimeStepper* stepper)
{
    // Add the rod
    m_rods.push_back(rod);
    m_simulated_rods.push_back(m_number_of_rods);
    m_rods.back()->globalRodIndex = m_number_of_rods++;
    // Add the stepper
    m_steppers.push_back(stepper);
    m_steppers.back()->setMaxIterations(m_perf_param.m_solver.m_max_iterations);

    assert(m_rods.size() == m_steppers.size());
    assert(m_number_of_rods == m_rods.size());

    // Extend the force vectors used to detect explosions
    m_startForces.push_back(new VecXd(m_rods.back()->ndof()));
    m_preDynamicForces.push_back(new VecXd(m_rods.back()->ndof()));
    m_preCollisionForces.push_back(new VecXd(m_rods.back()->ndof()));
    m_endForces.push_back(new VecXd(m_rods.back()->ndof()));

    m_rodbackups.resize(m_number_of_rods); // growing by 1
    m_rodbackups.back().resize(*m_rods.back());

    // Extract edges from the new rod
    for (int j = 0; j < m_rods.back()->nv() - 1; ++j)
    {
        m_edges.push_back(std::pair<int, int>(getNumVerts() + j, getNumVerts() + j + 1));
        assert(m_rods.back()->getRadiusScale() * m_rods.back()->radiusA(j) > 0.0);
        m_edge_radii.push_back(m_rods.back()->getRadiusScale() * m_rods.back()->radiusA(j));
    }
    assert(m_edges.size() == m_edge_radii.size());

    for (int j = 0; j < m_rods.back()->nv() - 1; ++j)
    {
        assert(m_rods.back()->getRadiusScale() * m_rods.back()->radiusA(j) > 0.0);
        // Radii associated with edges ... what to do if at a vertex with two radii? Average?
        m_vertex_radii.push_back(m_rods.back()->getRadiusScale() * m_rods.back()->radiusA(j));
    }
    m_vertex_radii.push_back(m_vertex_radii.back()); // < TODO: What the $^#! is this call?

    // Update vector that tracks the rod DOF in the system
    m_base_dof_indices.push_back(getNumDof());
    m_base_vtx_indices.push_back(getNumDof() / 3);

    // Extract masses from the new rod
    for (ElasticRod::vertex_iter itr = m_rods.back()->vertices_begin(); itr != m_rods.back()->vertices_end(); ++itr)
    {
        assert(m_rods.back()->getVertexMass(*itr, -1) > 0.0);
        m_masses.push_back(m_rods.back()->getVertexMass(*itr, -1));
    }

    // Update total number of DOF in the system
    m_num_dof += 3 * m_rods.back()->nv();

    assert((int) m_masses.size() == getNumVerts());
    assert((int) m_vertex_radii.size() == getNumVerts());

    // Resize the internal storage
    m_xn.resize(getNumDof());
    m_xnp1.resize(getNumDof());
    m_xdebug.resize(getNumDof());
    m_vnphalf.resize(getNumDof());

    m_collision_immune.resize(getNumVerts());

    if (m_perf_param.m_enable_penalty_response)
        enableImplicitPenaltyImpulses(); // TODO: probably a memory leak to fix here.

    // Update the rods in the collision detector
    m_collision_detector->rebuildRodElements(m_edges);

    // Initial length of the recently added rod
    m_initialLengths.push_back(0);
    for (int j = 1; j < m_rods.back()->nv(); j++)
        m_initialLengths.back() += (m_rods.back()->getVertex(j) - m_rods.back()->getVertex(j - 1)).norm();

}

void BAGroomingStepper::removeRod(int rodIdx)
{
    killTheRod(rodIdx);
}

}

/*
 * BAGroomingStepper_CollisionResponse.cc
 *
 *  Created on: 19/05/2011
 *
 */

#include <typeinfo>
#include "BAGroomingStepper.hh"
#include "../../Threads/MultithreadedStepper.hh"
#include "../../Core/Timer.hh"
#include "../../Collisions/Collision.hh"

#include <iostream.h>
#ifdef HAVE_OPENMP
#include <omp.h>
#endif

namespace BASim
{

/**
 * Inelastic response
 */
void BAGroomingStepper::exertVertexImpulse(const Vec3d& I, const double& m, const int& idx, VecXd& v)
{
    assert(m > 0.0);
    assert(idx >= 0);
    assert(idx < getNumVerts());
    assert(v.size() == getNumDof());

    v.segment<3> (3 * idx) += I / m;
}

void BAGroomingStepper::exertEdgeImpulse(const Vec3d& I, const double& m0, const double& m1, const double& alpha,
        const int& idx0, const int& idx1, VecXd& v)
{
    assert(m0 > 0.0);
    assert(m1 > 0.0);
    assert(alpha >= 0.0);
    assert(alpha <= 1.0);
    assert(idx0 >= 0);
    assert(idx0 < getNumVerts());
    assert(idx1 >= 0);
    assert(idx1 < getNumVerts());
    assert(v.size() == getNumDof());

    v.segment<3> (3 * idx0) += ((1 - alpha) / m0) * I;
    v.segment<3> (3 * idx1) += (alpha / m1) * I;
}

void BAGroomingStepper::exertFaceImpulse(const Vec3d& I, const double& m0, const double& m1, const double& m2, const double& u,
        const double& v, const double& w, const int& idx0, const int& idx1, const int& idx2, VecXd& vel)
{
    assert(m0 > 0.0);
    assert(m1 > 0.0);
    assert(m2 > 0.0);
    //assert( u >= 0.0 ); assert( u <= 1.0 );
    //assert( v >= 0.0 ); assert( v <= 1.0 );
    //assert( w >= 0.0 ); assert( w <= 1.0 );
    assert(approxEq(u + v + w, 1.0));
    assert(idx0 >= 0);
    assert(idx0 < getNumVerts());
    assert(idx1 >= 0);
    assert(idx1 < getNumVerts());
    assert(idx2 >= 0);
    assert(idx2 < getNumVerts());
    assert(vel.size() == getNumDof());

    vel.segment<3> (3 * idx0) += u * I / m0;
    vel.segment<3> (3 * idx1) += v * I / m1;
    vel.segment<3> (3 * idx2) += w * I / m2;
}

void BAGroomingStepper::exertInelasticImpulse(EdgeEdgeCTCollision& eecol)
{
    assert(eecol.e0_v0 >= 0);
    assert(eecol.e0_v0 < getNumVerts());
    assert(eecol.e0_v1 >= 0);
    assert(eecol.e0_v1 < getNumVerts());
    assert(eecol.e1_v0 >= 0);
    assert(eecol.e1_v0 < getNumVerts());
    assert(eecol.e1_v1 >= 0);
    assert(eecol.e1_v1 < getNumVerts());

    // Compute the relative velocity of the edges at the collision point
    //   Vec3d relvel = m_geodata.computeRelativeVelocity(clssn.e0_v0, clssn.e0_v1, clssn.e1_v0, clssn.e1_v1, clssn.s, clssn.t);
    //   double magrelvel = relvel.dot(clssn.n);
    //    double magrelvel = eecol.computeRelativeVelocity(m_geodata);

    // TraceStream(g_log, "") << "BAGroomingStepper:exertInelasticImpulse: pre-impulse e-e relative velocity = "
    //         << eecol.computeRelativeVelocity() << '\n';

    if (eecol.GetCachedRelativeVelocity() >= 0.0)
    {
        DebugStream(g_log, "")
                << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Relative velocity computed \
                  incorrectly before applying edge-edge inelastic impulse (bug in normal \
                  computation?). Magnitude of relative velocity: "
                << eecol.computeRelativeVelocity() << ". Rod: " << getContainingRod(eecol.e0_v0) << '\n';
    }
    assert(eecol.GetCachedRelativeVelocity() < 0.0);

    // Add some extra "kick" to relative velocity to account for FPA errors
    //eecol.ApplyRelativeVelocityKick();
    Vec3d I = eecol.computeInelasticImpulse();
    //computeEdgeEdgeInelasticImpulse(m_masses[eecol.e0_v0], m_masses[eecol.e0_v1], m_masses[eecol.e1_v0],
    //m_masses[eecol.e1_v1], eecol.s, eecol.t, eecol.GetRelativeVelocity(), eecol.n);

    // TraceStream(g_log, "") << "BAGroomingStepper::exertInelasticEdgeEdgeImpulse: (" << eecol.e0_v0 << "-" << eecol.e0_v1 << ", " << eecol.e1_v0
    //         << "-" << eecol.e1_v1 << ") s=" << eecol.s << " t=" << eecol.t << " I=" << I << '\n';

    exertEdgeImpulse(-I, m_masses[eecol.e0_v0], m_masses[eecol.e0_v1], eecol.s, eecol.e0_v0, eecol.e0_v1, m_vnphalf);
    exertEdgeImpulse(I, m_masses[eecol.e1_v0], m_masses[eecol.e1_v1], eecol.t, eecol.e1_v0, eecol.e1_v1, m_vnphalf);

    // TraceStream(g_log, "") << "BAGroomingStepper:exertInelasticImpulse: post-impulse e-e relative velocity = "
    //         << eecol.computeRelativeVelocity() << '\n';

    assert(eecol.computeRelativeVelocity() >= 0);
}

void BAGroomingStepper::exertInelasticImpulse(VertexFaceCTCollision& vfcol)
{
    assert(vfcol.v0 >= 0);
    assert(vfcol.v0 < getNumVerts());
    assert(vfcol.f0 >= 0);
    assert(vfcol.f0 < getNumVerts());
    assert(vfcol.f1 >= 0);
    assert(vfcol.f1 < getNumVerts());
    assert(vfcol.f2 >= 0);
    assert(vfcol.f2 < getNumVerts());

    // Compute the relative velocity of the edges at the collision point
    //   double magrelvel = vfcol.computeRelativeVelocity(m_geodata);
    assert(vfcol.GetCachedRelativeVelocity() < 0.0);

    // TraceStream(g_log, "") << "BAGroomingStepper:exertInelasticImpulse: pre-impulse v-f relative velocity = "
    //         << vfcol.GetCachedRelativeVelocity() << '\n';

    // Add some extra "kick" to relative velocity to account for FPA errors
    //vfcol.ApplyRelativeVelocityKick();
    Vec3d I = vfcol.computeInelasticImpulse();
    // computeVertexFaceInelasticImpulse(m_masses[vfcol.v0], m_masses[vfcol.f0], m_masses[vfcol.f1], m_masses[vfcol.f2],
    // vfcol.u, vfcol.v, vfcol.w, eecol.GetRelativeVelocity(), vfcol.n);

    // TraceStream(g_log, "") << "BAGroomingStepper::exertInelasticImpulse: (" << vfcol.v0 << ", " << vfcol.f0 << "-" << vfcol.f1 << "-"
    //         << vfcol.f2 << ") (u,v,w)=(" << vfcol.u << "," << vfcol.v << "," << vfcol.w << ") I=" << I << '\n';

    exertFaceImpulse(-I, m_masses[vfcol.f0], m_masses[vfcol.f1], m_masses[vfcol.f2], vfcol.u, vfcol.v, vfcol.w, vfcol.f0,
            vfcol.f1, vfcol.f2, m_vnphalf);
    exertVertexImpulse(I, m_masses[vfcol.v0], vfcol.v0, m_vnphalf);

    // TraceStream(g_log, "") << "BAGroomingStepper:exertInelasticImpulse: post-impulse v-f relative velocity = "
    //         << vfcol.computeRelativeVelocity() << '\n';

    assert(vfcol.computeRelativeVelocity() >= 0);
}

/**
 * Compliant inelastic response
 */
bool BAGroomingStepper::executeIterativeInelasticImpulseResponse(std::vector<bool>& failed_collisions_rods,
        std::vector<bool>& stretching_rods)
{
    bool all_rods_collisions_ok = true;
    return all_rods_collisions_ok; /*
     // Check whether the solver left some rods stretched
     checkLengths(stretching_rods);

     // Detect continuous time collisions
     std::list<Collision*> collisions_list;
     TraceStream(g_log, "") << "Detecting collisions...\n";
     m_collision_detector->getCollisions(collisions_list, ContinuousTime);
     TraceStream(g_log, "") << "Initial potential collisions: " << m_collision_detector->m_potential_collisions << "\n";

     // Iterativly apply inelastic impulses
     for (int itr = 0; !collisions_list.empty() && itr < m_perf_param.m_collision.m_max_iterations; ++itr)
     {
     TraceStream(g_log, "") << "CTcollision response iteration " << itr << '\n';
     TraceStream(g_log, "") << "Detected " << collisions_list.size() << " continuous time collisions (potential: "
     << m_collision_detector->m_potential_collisions << ")\n";

     // Just sort the collision times to maintain some rough sense of causality
     collisions_list.sort(CompareTimes);

     if (m_perf_param.m_skipRodRodCollisions)
     {
     // Keep only one collision per rod
     std::set<int> already_collided_rods;
     for (std::list<Collision*>::iterator col_it = collisions_list.begin(); col_it != collisions_list.end(); col_it++)
     {
     CTCollision* collision = dynamic_cast<CTCollision*> (*col_it);
     if (collision)
     {
     // Only keep collisions involving a rod we see for the first time
     int colliding_rod = getContainingRod(collision->GetRodVertex());
     if (already_collided_rods.find(colliding_rod) == already_collided_rods.end())
     already_collided_rods.insert(colliding_rod);
     else
     {
     delete collision;
     collisions_list.erase(col_it--);
     }
     }
     }
     TraceStream(g_log, "") << "of which " << collisions_list.size() << " are on different rods\n";
     // Now apply response to the remaining collisions
     for (std::list<Collision*>::iterator col_it = collisions_list.begin(); col_it != collisions_list.end(); col_it++)
     {
     CTCollision* collision = dynamic_cast<CTCollision*> (*col_it);
     if (collision)
     exertCompliantInelasticImpulse(collision);
     // So, which rod was it?
     int colliding_rod = getContainingRod(collision->GetRodVertex());

     /*
     // Check for explosion here

     m_rods[colliding_rod]->computeForces(*m_endForces[colliding_rod]);
     TraceStream(g_log, "") << "Rod " << colliding_rod << " Force norms: initial: "
     << (*(m_startForces[colliding_rod])).norm() << " pre-collisions: "
     << (*(m_preCollisionForces[colliding_rod])).norm() << " post-collisions: "
     << (*(m_endForces[colliding_rod])).norm() << "\n";
     for (int j = 0; j < m_rods[colliding_rod]->ndof(); ++j)
     {
     const double s = (*(m_startForces[colliding_rod]))[j];
     const double p = (*(m_preCollisionForces[colliding_rod]))[j];
     const double e = (*(m_endForces[colliding_rod]))[j];
     const double rate = fabs(s - e) / (fabs(s) + m_perf_param.m_explosion_damping);
     if (isnan(rate) || rate > m_perf_param.m_explosion_threshold)
     {
     //                        explosions_detected = true;
     //                        exploding_rods[colliding_rod] = true;
     TraceStream(g_log, "") << "Rod number " << colliding_rod
     << " had an explosion during collision response: s = " << s << " p = " << p << " e = " << e
     << " rate = " << rate << " \n";
     // If the rod just had an explosion, give up trying resolving its collisions
     failed_collisions_rods[colliding_rod] = true;
     for (int v = 0; v < m_rods[colliding_rod]->nv(); ++v)
     m_collision_immune[m_base_dof_indices[colliding_rod] / 3 + v] = true;
     break;
     }
     }
     */
    /*
     // Test for stretching and if so, mark the rod immune
     stretching_rods[colliding_rod] = !checkLength(colliding_rod);
     if (stretching_rods[colliding_rod])
     {
     // Declare the rod collision-immune for the rest of the time step
     for (int j = 0; j < m_rods[colliding_rod]->nv(); ++j)
     m_collision_immune[m_base_vtx_indices[colliding_rod] + j] = false;
     }

     delete collision;
     collisions_list.erase(col_it--);
     }
     }
     else
     while (!collisions_list.empty())
     {
     CTCollision* collision = dynamic_cast<CTCollision*> (collisions_list.front());
     collisions_list.pop_front();
     if (collision)
     exertCompliantInelasticImpulse(collision);
     delete collision;
     m_collision_detector->updateContinuousTimeCollisions();
     }

     // Detect remaining collisions (including at the end of the last iteration, so we know what failed)
     TraceStream(g_log, "") << "Detecting collisions...\n";
     m_collision_detector->getCollisions(collisions_list, ContinuousTime, false); // No need to update the mesh bvh bounding boxes.
     }

     if (!collisions_list.empty())
     {
     all_rods_collisions_ok = false;

     TraceStream(g_log, "") << "Remains " << collisions_list.size() << " unresolved collisions (potential: "
     << m_collision_detector->m_potential_collisions << ")\n";

     // Just in case we haven't emptied the collisions but exited when itr == m_num_inlstc_itrns
     for (std::list<Collision*>::iterator col = collisions_list.begin(); col != collisions_list.end(); col++)
     {
     // Let's see which collisions are remaining.
     EdgeEdgeCTCollision* eecol = dynamic_cast<EdgeEdgeCTCollision*> (*col);
     VertexFaceCTCollision* vfcol = dynamic_cast<VertexFaceCTCollision*> (*col);
     if (vfcol)
     {
     assert(isRodVertex(vfcol->v0));
     failed_collisions_rods[getContainingRod(vfcol->v0)] = true;
     }
     if (eecol)
     {
     if (isRodVertex(eecol->e0_v0) && isRodVertex(eecol->e0_v1))
     {
     assert(getContainingRod(eecol->e0_v0) == getContainingRod(eecol->e0_v1));
     failed_collisions_rods[getContainingRod(eecol->e0_v0)] = true;
     }
     if (isRodVertex(eecol->e1_v0) && isRodVertex(eecol->e1_v1))
     {
     assert(getContainingRod(eecol->e1_v0) == getContainingRod(eecol->e1_v1));
     failed_collisions_rods[getContainingRod(eecol->e1_v0)] = true;
     }
     }
     delete *col;
     }
     }

     //#ifdef TIMING_ON
     //if( itr >= 2 ) IntStatTracker::getIntTracker("STEPS_WITH_MULTIPLE_IMPULSE_ITERATIONS") += 1;
     //#endif

     // checkLengths(stretching_rods);

     return all_rods_collisions_ok;*/
}

void BAGroomingStepper::computeCompliantLHS(MatrixBase* lhs, int rodidx)
{
    assert(lhs != NULL);
    assert(rodidx >= 0);
    assert(rodidx < m_number_of_rods);

    //double emphasizeStaticEquilibium = 0;
    //InfoStream(g_log,"") << "BAGroomingStepper::computeCompliantLHS: emphasizing static equilibrium = " << emphasizeStaticEquilibium;

    // lhs = -h^2*dF/dx
    lhs->setZero();
    //TraceStream(g_log, "") << "WARNING: COMPLIANCE is disabled!" << '\n';
    //m_rods[rodidx]->computeJacobian(0, -m_dt * m_dt - emphasizeStaticEquilibium, *lhs);
    m_rods[rodidx]->computeJacobian(0, -m_dt * m_dt, *lhs);
    lhs->finalize();

    // lhs = M - h^2*dF/dx
    for (int i = 0; i < m_rods[rodidx]->ndof(); ++i)
        lhs->add(i, i, m_rods[rodidx]->getMass(i));
    lhs->finalize();

    assert(lhs->isApproxSymmetric(1.0e-6));

    // Set rows/cols of twist DOFs to identity
    // TODO: Handle this better
    EPropHandle<int> edgeIdx;
    m_rods[rodidx]->property_handle(edgeIdx, "edge index");
    std::vector<int> dofs_to_clear;
    for (ElasticRod::edge_iter eit = m_rods[rodidx]->edges_begin(); eit != m_rods[rodidx]->edges_end(); ++eit)
    {
        int twist_dof = m_rods[rodidx]->property(edgeIdx)[*eit];
        dofs_to_clear.push_back(twist_dof);
    }

    lhs->zeroRows(dofs_to_clear, 1.0);
    lhs->finalize();
    lhs->zeroCols(dofs_to_clear, 1.0);
    lhs->finalize();

    assert(lhs->isApproxSymmetric(1.0e-6));

    // TEMP
    //  ElasticRod::RodForces& forces = m_rods[rodidx]->getForces();
    //  ElasticRod::RodForces::iterator fIt;
    //  for (fIt = forces.begin(); fIt != forces.end(); ++fIt)
    //    //(*fIt)->globalJacobian(baseidx, scale, J);
    //    TraceStream(g_log, "") << "Force: " << (*fIt)->getName() << '\n';
    // END TEMP
}

void BAGroomingStepper::exertCompliantInelasticImpulse(const CTCollision* cllsn)
{
    assert(cllsn->isAnalysed());
    const EdgeEdgeCTCollision* eecol = dynamic_cast<const EdgeEdgeCTCollision*> (cllsn);
    const VertexFaceCTCollision* vfcol = dynamic_cast<const VertexFaceCTCollision*> (cllsn);
    //TraceStream(g_log, "") << "BAGroomingStepper:exertCompliantInelasticImpulse: pre-impulse e-e relative velocity = " << cllsn->GetRelativeVelocity() << '\n';

    if (eecol)
    {
        exertCompliantInelasticEdgeEdgeImpulse(*eecol);
        //exertInelasticImpulse(cllsn.getEdgeEdge());
    }
    else if (vfcol)
    {
        exertCompliantInelasticVertexFaceImpulse(*vfcol);
        //exertCompliantInelasticVertexFaceImpulse(cllsn.getVertexFace());
        //exertInelasticImpulse(cllsn.getVertexFace());
    }
    //TraceStream(g_log, "") << "BAGroomingStepper:exertCompliantInelasticImpulse: post-impulse e-e relative velocity = " << cllsn->GetRelativeVelocity() << '\n';
}

void BAGroomingStepper::exertCompliantInelasticVertexFaceImpulse(const VertexFaceCTCollision& vfcol)
{
    //   TraceStream(g_log, "") << "Vertex-face compliant inelastic impulse" << '\n';
    //   TraceStream(g_log, "") << vfcol << '\n';

    // TraceStream(g_log, "") << "BAGroomingStepper::exertCompliantInelasticImpulse: (" << vfcol.v0 << ", " << vfcol.f0 << "-" << vfcol.f1 << "-"
    //         << vfcol.f2 << ") (u,v,w)=(" << vfcol.u << "," << vfcol.v << "," << vfcol.w << ")" << '\n';

    // For now, assume vertex is free and entire face is fixed
    assert(!m_geodata.isVertexFixed(vfcol.v0));
    assert(YATriangle(vfcol.f0, vfcol.f1, vfcol.f2).IsFixed(m_geodata));

    // Determine which rod the free vertex belongs to
    int rodidx = getContainingRod(vfcol.v0);

    // Ensure the free vertex belongs to a rod
    assert(rodidx >= 0);
    assert(rodidx < (int) m_number_of_rods);

    // If the rod has not solved properly, no need to compute its collision response
    if (!m_steppers[rodidx]->HasSolved())
    {// NB: this is redundant because we declared the rod collision-immune already.
        DebugStream(g_log, "") << "WARNING: attempt to do vertex-face collision with non-dependable rod\n";
        return;
    }

    // Determine which vertex of the rod the free vertex is
    assert(m_base_dof_indices[rodidx] % 3 == 0);
    int v0 = vfcol.v0 - m_base_dof_indices[rodidx] / 3;
    assert(v0 >= 0);
    assert(v0 < m_rods[rodidx]->nv());

    TraceStream(g_log, "") << "CTcollision: vertex " << v0 << " of rod " << rodidx << '\n';

    // Compute the relative velocity of the collision
    //double magrelvel = vfcol.computeRelativeVelocity(m_geodata);
    assert(vfcol.computeRelativeVelocity() < 0.0);

    // Get storage for lhs of linear system, get a solver
    LinearSystemSolver* lss = m_solver_collection.getLinearSystemSolver(m_rods[rodidx]->ndof());
    MatrixBase* lhs = lss->m_lhs;
    assert(lhs != NULL);
    LinearSolverBase* solver = lss->m_solver;
    assert(solver != NULL);

    // Compute M - h^2*dF/dx
    computeCompliantLHS(lhs, rodidx);

    // Compute the 'base' index of the vertex
    int base0 = m_rods[rodidx]->vertIdx(v0, 0);
    assert(base0 % 4 == 0);

    int ndof = m_rods[rodidx]->ndof();

    int rodbase = m_base_dof_indices[rodidx];

    // Compute the 'global' normal
    std::vector<VecXd> normal;

    normal.push_back(VecXd(ndof));
    normal.back().setZero();
    normal.back().segment<3> (base0) = vfcol.GetNormal();

    // Compute three normals per scripted vertex
    const std::vector<int>& scriptedverts = m_rods[rodidx]->getBoundaryCondition()->scriptedVertices();
    for (size_t i = 0; i < scriptedverts.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            normal.push_back(VecXd(ndof));
            normal.back().setZero();
            // Zero the velocity of this particular dof of this particular vertex
            normal.back()(m_rods[rodidx]->vertIdx(scriptedverts[i], j)) = 1.0;
        }
    }
    int numconstraints = (int) (normal.size());
    int nvdof = 3 * m_rods[rodidx]->nv();

    //TraceStream(g_log, "") << "nc: " << numconstraints << '\n';
    //for( size_t i = 0; i < n.size(); ++i ) TraceStream(g_log, "") << n[i] << '\n';

    // Determine the desired values for each constraint
    VecXd desired_values(numconstraints);
#ifndef NDEBUG
    desired_values.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    // First constraint, the collision constraint, is 0
    desired_values(0) = 0.0;
    // The 'scripted vertex' constraints are the respective components of the velocity
    int curdof = 1;
    for (size_t i = 0; i < scriptedverts.size(); ++i)
    {
        assert(scriptedverts[i] >= 0);
        assert((int) scriptedverts[i] < m_rods[rodidx]->nv());
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 0);
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 1);
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 2);
    }
    assert(curdof == numconstraints);
    assert(desired_values.size() == numconstraints);

    //TraceStream(g_log, "") << "Constraint values: " << desired_values << '\n';

    // Currently, all fixed vertex constraints normalized
#ifndef NDEBUG
    for (int i = 0; i < numconstraints; ++i)
        assert(approxEq(normal[i].norm(), 1.0, 1.0e-9));
#endif

    std::vector<VecXd> ntilde;
    for (int i = 0; i < numconstraints; ++i)
    {
        ntilde.push_back(VecXd(ndof));
        ntilde.back().setZero();
        int status = solver->solve(ntilde.back(), normal[i]);
        if (status < 0)
            DebugStream(g_log, "")
                    << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. exertCompliantInelasticVertexFaceImpulse. Time: "
                    << m_t << '\n';
    }

    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (ElasticRod::edge_iter eit = m_rods[rodidx]->edges_begin(); eit != m_rods[rodidx]->edges_end(); ++eit)
    {
        for (int i = 0; i < numconstraints; ++i)
            assert(ntilde[i](m_rods[rodidx]->edgeIdx(*eit)) == 0.0);
    }
#endif

    // Vectors restriced to just vertex DoFs. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn;
    for (int i = 0; i < numconstraints; ++i)
        posnn.push_back(VecXd(nvdof));
#ifndef NDEBUG
    for (int i = 0; i < numconstraints; ++i)
        posnn[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    std::vector<VecXd> posnntilde;
    for (int i = 0; i < numconstraints; ++i)
        posnntilde.push_back(VecXd(nvdof));
#ifndef NDEBUG
    for (int i = 0; i < numconstraints; ++i)
        posnntilde[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (int i = 0; i < numconstraints; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rodidx]->vertices_begin(); vit != m_rods[rodidx]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rodidx]->vertIdx(*vit, 0);
            posnn[i].segment<3> (3 * theidx) = normal[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    for (int i = 0; i < numconstraints; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rodidx]->vertices_begin(); vit != m_rods[rodidx]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rodidx]->vertIdx(*vit, 0);
            posnntilde[i].segment<3> (3 * theidx) = ntilde[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    // Compute the Lagrange multipliers
    // TODO: This matrix is symmetric, exploit that fact to avoid computations
    MatXd lglhs(numconstraints, numconstraints);
#ifndef NDEBUG
    lglhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    for (int i = 0; i < numconstraints; ++i)
        for (int j = 0; j < numconstraints; ++j)
            lglhs(i, j) = posnn[i].dot(posnntilde[j]);
    assert(approxSymmetric(lglhs, 1.0e-6));

    Eigen::VectorXd lgrhs(numconstraints);
#ifndef NDEBUG
    lgrhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    lgrhs(0) = vfcol.computeRelativeVelocity(); //posnN0.dot(m_vnphalf.segment(m_base_dof_indices[rodidx],posnN0.size()));
    assert(lgrhs(0) < 0.0);
    for (int i = 1; i < numconstraints; ++i)
        lgrhs(i) = posnn[i].dot(m_vnphalf.segment(rodbase, nvdof));
    lgrhs *= -1.0;

    //TraceStream(g_log, "") << lgrhs << '\n';

    lgrhs += desired_values;

    //TraceStream(g_log, "") << lgrhs << '\n';
    //TraceStream(g_log, "") << desired_values << '\n';

    Eigen::VectorXd alpha(numconstraints); // = lglhs.inverse()*lgrhs;
    alpha = lglhs.lu().solve(lgrhs);

    assert(alpha(0) >= 0.0);

    double magrelvel = vfcol.computeRelativeVelocity();

    // TraceStream(g_log, "") << "BAGroomingStepper::exertCompliantInelasticVertexFaceImpulse: pre-impulse view of collision: " << vfcol
    //         << '\n';
    //TraceStream(g_log, "") << "BAGroomingStepper::exertCompliantInelasticVertexFaceImpulse: pre-impulse velocities: " << m_vnphalf.segment(rodbase, nvdof) << '\n';

    for (int i = 0; i < numconstraints; ++i)
        m_vnphalf.segment(rodbase, nvdof) += alpha(i) * posnntilde[i];

    // TraceStream(g_log, "") << "BAGroomingStepper::exertCompliantInelasticVertexFaceImpulse: post-impulse view of collision: " << vfcol
    //         << '\n';
    //TraceStream(g_log, "") << "BAGroomingStepper::exertCompliantInelasticVertexFaceImpulse: post-impulse velocities: " << m_vnphalf.segment(rodbase, nvdof) << '\n';

    // Ensure that the impulse eliminated the realtive velocity
    //#ifndef NDEBUG
    double postmagrelvel = vfcol.computeRelativeVelocity();
    // TraceStream(g_log, "") << "BAGroomingStepper::exertCompliantInelasticVertexFaceImpulse: relative velocity pre-impulse = " << magrelvel
    //         << " post-impulse = " << postmagrelvel << '\n';
    // Ensure the inelastic impulse decreased the realtive velocity
    assert(fabs(postmagrelvel) <= fabs(magrelvel));
    // Should add some 'extra kick' to ensure collision gets killed, but for now just be content with small velocity
    assert(fabs(postmagrelvel) < 1.0e-9);
    //assert( postmagrelvel < 0.0 );
    //#endif

    applyInextensibilityVelocityFilter(rodidx);

    // Ensure the 'scripted vertices' achieved the desired velocity
#ifndef NDEBUG
    curdof = 1;
    for (size_t i = 0; i < scriptedverts.size(); ++i)
    {
        assert(scriptedverts[i] >= 0);
        assert((int) scriptedverts[i] < m_rods[rodidx]->nv());
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+0) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rodbase + 3 * scriptedverts[i] + 0), desired_values(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+1) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rodbase + 3 * scriptedverts[i] + 1), desired_values(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+2) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rodbase + 3 * scriptedverts[i] + 2), desired_values(curdof), 1.0e-6));
        ++curdof;
    }
#endif
}

void BAGroomingStepper::exertCompliantInelasticEdgeEdgeImpulse(const EdgeEdgeCTCollision& eecol)
{
    //   TraceStream(g_log, "") << "Edge-edge compliant inelastic impulse" << '\n';
    //   TraceStream(g_log, "") << eecol << '\n';

    // Determine if either edge is totally fixed
    bool rod0fixed = YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(m_geodata);
    bool rod1fixed = YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata);

    // Skip collisions in which both edges are fixed
    if (rod0fixed && rod1fixed)
        return;
    // Neither rod is totally fixed

    else if (!rod0fixed && !rod1fixed)
        exertCompliantInelasticEdgeEdgeImpulseBothFree(eecol);
    // One rod is totally fixed

    else
        exertCompliantInelasticEdgeEdgeImpulseOneFixed(eecol);
}

void BAGroomingStepper::exertCompliantInelasticEdgeEdgeImpulseBothFree(const EdgeEdgeCTCollision& eecol)
{
    // Ensure both edges have two free vertices
    //assert( isEntireEdgeFree(eecol.e0_v0,eecol.e0_v1) );
    //assert( isEntireEdgeFree(eecol.e1_v0,eecol.e1_v1) );

    assert(!YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(m_geodata) && !YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata));

    // TraceStream(g_log, "") << "BAGroomingStepper::exertCompliantInelasticEdgeEdgeImpulseBothFree: (x[" << eecol.e0_v0 << "]="
    //         << m_geodata.GetPoint(eecol.e0_v0) << " - x[" << eecol.e0_v1 << "]=" << m_geodata.GetPoint(eecol.e0_v1) << ",   x["
    //         << eecol.e1_v0 << "]=" << m_geodata.GetPoint(eecol.e1_v0) << " - x[" << eecol.e1_v1 << "]=" << m_geodata.GetPoint(
    //         eecol.e1_v1) << ")" << " s=" << eecol.s << " t=" << eecol.t << '\n';

    // Determine which rod each edge belongs to
    int rod0 = getContainingRod(eecol.e0_v0);
    assert(rod0 == getContainingRod(eecol.e0_v1));
    int rod1 = getContainingRod(eecol.e1_v0);
    assert(rod1 == getContainingRod(eecol.e1_v1));

    // Don't do self-collisions, for now
    assert(rod0 != rod1);

    if (!m_steppers[rod0]->HasSolved() || !m_steppers[rod1]->HasSolved())
    { // This should never happen because rod-rod collisions preclude selective adaptivity.
        DebugStream(g_log, "") << "WARNING: attempt to do rod-rod collision on non-dependable rods";
        return;
    }

    TraceStream(g_log, "") << "CTcollision: rod " << rod0 << " vs. rod " << rod1 << '\n';

    // Compute the relative velocity, which must be negative
    //    Vec3d relvel = m_geodata.computeRelativeVelocity(eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t);
    //    double magrelvel = relvel.dot(eecol.n);
    //   double magrelvel = eecol.computeRelativeVelocity(m_geodata);
    assert(eecol.computeRelativeVelocity() < 0.0);

    // Determine the total number of degrees of freedom in each rod
    int rod0ndof = m_rods[rod0]->ndof();
    int rod1ndof = m_rods[rod1]->ndof();

    // Determine the number of vertices in each rod
    int rod0nv = m_rods[rod0]->nv();
    int rod1nv = m_rods[rod1]->nv();

    // Determine the number of vertices in each rod
    int rod0nvdof = 3 * rod0nv;
    int rod1nvdof = 3 * rod1nv;

    //TraceStream(g_log, "") << "Ndof 0: " << rod0ndof << '\n';
    //TraceStream(g_log, "") << "Ndof 1: " << rod1ndof << '\n';
    //TraceStream(g_log, "") << "NV 0: " << rod0nv << '\n';
    //TraceStream(g_log, "") << "NV 1: " << rod1nv << '\n';
    //TraceStream(g_log, "") << "NVdof 0: " << rod0nvdof << '\n';
    //TraceStream(g_log, "") << "NVdof 1: " << rod1nvdof << '\n';

    // Determine where in the 'global' vertex dof pool each rod begins
    int rod0base = m_base_dof_indices[rod0];
    assert(rod0base % 3 == 0);
    int rod1base = m_base_dof_indices[rod1];
    assert(rod1base % 3 == 0);

    //TraceStream(g_log, "") << "rod0base: " << rod0base << '\n';
    //TraceStream(g_log, "") << "rod1base: " << rod1base << '\n';

    // Vertex number in 'position N'
    int i0 = eecol.e0_v0 - rod0base / 3;
    assert(i0 >= 0);
    assert(i0 < m_rods[rod0]->nv());
    int i1 = eecol.e0_v1 - rod0base / 3;
    assert(i1 >= 0);
    assert(i1 < m_rods[rod0]->nv());
    int j0 = eecol.e1_v0 - rod1base / 3;
    assert(j0 >= 0);
    assert(j0 < m_rods[rod1]->nv());
    int j1 = eecol.e1_v1 - rod1base / 3;
    assert(j1 >= 0);
    assert(j1 < m_rods[rod1]->nv());

    // Base indices within each rod's storage
    int base_i0 = m_rods[rod0]->vertIdx(i0, 0);
    assert(base_i0 >= 0);
    assert(base_i0 < rod0ndof);
    int base_i1 = m_rods[rod0]->vertIdx(i1, 0);
    assert(base_i1 >= 0);
    assert(base_i1 < rod0ndof);
    int base_j0 = m_rods[rod1]->vertIdx(j0, 0);
    assert(base_j0 >= 0);
    assert(base_j0 < rod1ndof);
    int base_j1 = m_rods[rod1]->vertIdx(j1, 0);
    assert(base_j1 >= 0);
    assert(base_j1 < rod1ndof);

    // Constraint normals for rod 0
    std::vector<VecXd> n0;
    // Collision constraint
    n0.push_back(VecXd(rod0ndof));
    n0.back().setZero();
    n0.back().segment<3> (base_i0) = -(1.0 - eecol.s) * eecol.GetNormal();
    n0.back().segment<3> (base_i1) = -eecol.s * eecol.GetNormal();
    // 'Scripted vertex' constraints
    const std::vector<int>& scriptedverts0 = m_rods[rod0]->getBoundaryCondition()->scriptedVertices();
    for (size_t i = 0; i < scriptedverts0.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            n0.push_back(VecXd(rod0ndof));
            n0.back().setZero();
            // Zero the velocity of this particular dof of this particular vertex
            n0.back()(m_rods[rod0]->vertIdx(scriptedverts0[i], j)) = 1.0;
        }
    }
    // Determine the desired value of each constraint
    int nc0 = (int) (n0.size());
    VecXd cval0(nc0);
#ifndef NDEBUG
    cval0.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    // First constraint, the collision constraint, is 0
    cval0(0) = 0.0;
    // The 'scripted vertex' constraints are the respective components of the velocity
    int curdof = 1;
    for (size_t i = 0; i < scriptedverts0.size(); ++i)
    {
        assert(scriptedverts0[i] >= 0);
        assert((int) scriptedverts0[i] < m_rods[rod0]->nv());
        cval0(curdof++) = m_vnphalf(rod0base + 3 * scriptedverts0[i] + 0);
        cval0(curdof++) = m_vnphalf(rod0base + 3 * scriptedverts0[i] + 1);
        cval0(curdof++) = m_vnphalf(rod0base + 3 * scriptedverts0[i] + 2);
    }
    assert(curdof == nc0);
    assert(cval0.size() == nc0);
    // Ensure collision constraint adds up to actual normal
#ifndef NDEBUG
    Vec3d testn0 = n0[0].segment<3> (base_i0) + n0[0].segment<3> (base_i1);
    Vec3d actln0 = -eecol.GetNormal();
    assert(approxEq(testn0, actln0, 1.0e-6));
#endif
    // Currently, all fixed vertex constraints normalized
#ifndef NDEBUG
    for (int i = 1; i < nc0; ++i)
        assert(approxEq(n0[i].norm(), 1.0, 1.0e-9));
#endif
    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (int i = 0; i < nc0; ++i)
        for (ElasticRod::edge_iter eit = m_rods[rod0]->edges_begin(); eit != m_rods[rod0]->edges_end(); ++eit)
            assert(n0[i](m_rods[rod0]->edgeIdx(*eit)) == 0.0);
#endif

    //for( int i = 0; i < nc0; ++i ) TraceStream(g_log, "") << "n0:    " << n0[i] << '\n';
    //for( int i = 0; i < nc0; ++i ) TraceStream(g_log, "") << "cval0: " << cval0[i] << '\n';

    // Constraint normals for rod 1
    std::vector<VecXd> n1;
    // Collision constraint
    n1.push_back(VecXd(rod1ndof));
    n1.back().setZero();
    n1.back().segment<3> (base_j0) = (1.0 - eecol.t) * eecol.GetNormal();
    n1.back().segment<3> (base_j1) = eecol.t * eecol.GetNormal();
    // 'Scripted vertex' constraints
    const std::vector<int>& scriptedverts1 = m_rods[rod1]->getBoundaryCondition()->scriptedVertices();
    for (size_t i = 0; i < scriptedverts1.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            n1.push_back(VecXd(rod1ndof));
            n1.back().setZero();
            // Zero the velocity of this particular dof of this particular vertex
            n1.back()(m_rods[rod1]->vertIdx(scriptedverts1[i], j)) = 1.0;
        }
    }
    // Determine the desired value of each constraint
    int nc1 = (int) (n1.size());
    VecXd cval1(nc1);
#ifndef NDEBUG
    cval1.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    // First constraint, the collision constraint, is 0
    cval1(0) = 0.0;
    // The 'scripted vertex' constraints are the respective components of the velocity
    curdof = 1;
    for (size_t i = 0; i < scriptedverts1.size(); ++i)
    {
        assert(scriptedverts1[i] >= 0);
        assert((int) scriptedverts1[i] < m_rods[rod1]->nv());
        cval1(curdof++) = m_vnphalf(rod1base + 3 * scriptedverts1[i] + 0);
        cval1(curdof++) = m_vnphalf(rod1base + 3 * scriptedverts1[i] + 1);
        cval1(curdof++) = m_vnphalf(rod1base + 3 * scriptedverts1[i] + 2);
    }
    assert(curdof == nc1);
    assert(cval1.size() == nc1);
    // Ensure collision constraint adds up to actual normal
#ifndef NDEBUG
    Vec3d testn1 = n1[0].segment<3> (base_j0) + n1[0].segment<3> (base_j1);
    Vec3d actln1 = eecol.GetNormal();
    assert(approxEq(testn1, actln1, 1.0e-6));
#endif
    // Currently, all fixed vertex constraints normalized
#ifndef NDEBUG
    for (int i = 1; i < nc1; ++i)
        assert(approxEq(n1[i].norm(), 1.0, 1.0e-9));
#endif
    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (int i = 0; i < nc1; ++i)
        for (ElasticRod::edge_iter eit = m_rods[rod1]->edges_begin(); eit != m_rods[rod1]->edges_end(); ++eit)
            assert(n1[i](m_rods[rod1]->edgeIdx(*eit)) == 0.0);
#endif

    //for( int i = 0; i < nc1; ++i ) TraceStream(g_log, "") << "n1:    " << n1[i] << '\n';
    //for( int i = 0; i < nc1; ++i ) TraceStream(g_log, "") << "cval1: " << cval1[i] << '\n';

    // Compute ntilde for each rod 0 constraint
    LinearSystemSolver* lss = m_solver_collection.getLinearSystemSolver(rod0ndof);
    MatrixBase* lhs0 = lss->m_lhs;
    assert(lhs0 != NULL);
    assert(lhs0->rows() == rod0ndof);
    assert(lhs0->cols() == rod0ndof);
    LinearSolverBase* solver0 = lss->m_solver;
    assert(solver0 != NULL);
    computeCompliantLHS(lhs0, rod0);
    std::vector<VecXd> ntilde0;
    for (int i = 0; i < nc0; ++i)
    {
        ntilde0.push_back(VecXd(rod0ndof));
        ntilde0.back().setZero();
        int status = solver0->solve(ntilde0.back(), n0[i]);
        if (status < 0)
            DebugStream(g_log, "") << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. "
                    << '\n';
    }
    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (int i = 0; i < nc0; ++i)
        for (ElasticRod::edge_iter eit = m_rods[rod0]->edges_begin(); eit != m_rods[rod0]->edges_end(); ++eit)
            assert(ntilde0[i](m_rods[rod0]->edgeIdx(*eit)) == 0.0);
#endif

    //for( int i = 0; i < nc0; ++i ) TraceStream(g_log, "") << "ntilde0: " << ntilde0[i] << '\n';

    // Compute ntilde for each rod 1 constraint
    lss = m_solver_collection.getLinearSystemSolver(rod1ndof);
    MatrixBase* lhs1 = lss->m_lhs;
    assert(lhs1->rows() == rod1ndof);
    assert(lhs1->cols() == rod1ndof);
    assert(lhs1 != NULL);
    LinearSolverBase* solver1 = lss->m_solver;
    assert(solver1 != NULL);
    computeCompliantLHS(lhs1, rod1);
    std::vector<VecXd> ntilde1;
    for (int i = 0; i < nc1; ++i)
    {
        ntilde1.push_back(VecXd(rod1ndof));
        ntilde1.back().setZero();
        int status = solver1->solve(ntilde1.back(), n1[i]);
        if (status < 0)
            DebugStream(g_log, "") << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. "
                    << '\n';
    }
    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (int i = 0; i < nc1; ++i)
        for (ElasticRod::edge_iter eit = m_rods[rod1]->edges_begin(); eit != m_rods[rod1]->edges_end(); ++eit)
            assert(ntilde1[i](m_rods[rod1]->edgeIdx(*eit)) == 0.0);
#endif

    //for( int i = 0; i < nc1; ++i ) TraceStream(g_log, "") << "ntilde1: " << ntilde1[i] << '\n';

    // Vectors restriced to just vertex DoFs for rod 0. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn0;
    for (int i = 0; i < nc0; ++i)
        posnn0.push_back(VecXd(rod0nvdof));
#ifndef NDEBUG
    for (int i = 0; i < nc0; ++i)
        posnn0[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    std::vector<VecXd> posnntilde0;
    for (int i = 0; i < nc0; ++i)
        posnntilde0.push_back(VecXd(rod0nvdof));
#ifndef NDEBUG
    for (int i = 0; i < nc0; ++i)
        posnntilde0[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (int i = 0; i < nc0; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rod0]->vertices_begin(); vit != m_rods[rod0]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rod0]->vertIdx(*vit, 0);
            posnn0[i].segment<3> (3 * theidx) = n0[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    //for( int i = 0; i < nc0; ++i ) TraceStream(g_log, "") << "posnn0: " << posnn0[i] << '\n';

    for (int i = 0; i < nc0; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rod0]->vertices_begin(); vit != m_rods[rod0]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rod0]->vertIdx(*vit, 0);
            posnntilde0[i].segment<3> (3 * theidx) = ntilde0[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    //for( int i = 0; i < nc0; ++i ) TraceStream(g_log, "") << "posnntilde0: " << posnntilde0[i] << '\n';


    // Vectors restriced to just vertex DoFs for rod 1. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn1;
    for (int i = 0; i < nc1; ++i)
        posnn1.push_back(VecXd(rod1nvdof));
#ifndef NDEBUG
    for (int i = 0; i < nc1; ++i)
        posnn1[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    std::vector<VecXd> posnntilde1;
    for (int i = 0; i < nc1; ++i)
        posnntilde1.push_back(VecXd(rod1nvdof));
#ifndef NDEBUG
    for (int i = 0; i < nc1; ++i)
        posnntilde1[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (int i = 0; i < nc1; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rod1]->vertices_begin(); vit != m_rods[rod1]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rod1]->vertIdx(*vit, 0);
            posnn1[i].segment<3> (3 * theidx) = n1[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    //for( int i = 0; i < nc1; ++i ) TraceStream(g_log, "") << "posnn1: " << posnn1[i] << '\n';

    for (int i = 0; i < nc1; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rod1]->vertices_begin(); vit != m_rods[rod1]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rod1]->vertIdx(*vit, 0);
            posnntilde1[i].segment<3> (3 * theidx) = ntilde1[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    //for( int i = 0; i < nc1; ++i ) TraceStream(g_log, "") << "posnntilde1: " << posnntilde1[i] << '\n';


    // Compute the Lagrange multipliers for all constraints
    // TODO: This matrix is symmetric, exploit that fact to avoid computations
    int numalpha = nc0 + nc1 - 1;
    MatXd lglhs(numalpha, numalpha);
    //#ifndef NDEBUG
    //  lglhs.setConstant(std::numeric_limits<double>::signaling_NaN());
    //#endif
    lglhs.setZero();

    //TraceStream(g_log, "") << "System size: " << numalpha << '\n';

    // TODO: Make sure indexing is correct here
    // Entry (0,0)
    lglhs(0, 0) = posnn0[0].dot(posnntilde0[0]) + posnn1[0].dot(posnntilde1[0]);
    // Entires (0,1::numalpha)
    int curel = 1;
    for (int i = 1; i < nc0; ++i)
        lglhs(0, curel++) = posnn0[0].dot(posnntilde0[i]);
    // Entries (0,numalpha::numalpha+numbeta)
    for (int i = 1; i < nc1; ++i)
        lglhs(0, curel++) = posnn1[0].dot(posnntilde1[i]);

    // Entires (1::numalpha,0)
    curel = 1;
    for (int i = 1; i < nc0; ++i)
        lglhs(curel++, 0) = posnntilde0[0].dot(posnn0[i]);
    // Entries (numalpha::numalpha+numbeta,0)
    for (int i = 1; i < nc1; ++i)
        lglhs(curel++, 0) = posnntilde1[0].dot(posnn1[i]);

    // Entries (1::numalpha,1::numalpha)
    for (int i = 1; i < nc0; ++i)
        for (int j = 1; j < nc0; ++j)
            lglhs(i, j) = posnn0[i].dot(posnntilde0[j]);
    // Entries (numalpha::numalpha+numbeta,numalpha::numalpha+numbeta)
    for (int i = 1; i < nc1; ++i)
        for (int j = 1; j < nc1; ++j)
            lglhs(nc0 + i - 1, nc0 + j - 1) = posnn1[i].dot(posnntilde1[j]);

    assert(approxSymmetric(lglhs, 1.0e-6));

    Eigen::VectorXd lgrhs(numalpha);
#ifndef NDEBUG
    lgrhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    // Entry 0
    assert(posnn0[0].size() == m_vnphalf.segment(rod0base, rod0nvdof).size());
    assert(posnn1[0].size() == m_vnphalf.segment(rod1base, rod1nvdof).size());
    lgrhs(0) = posnn0[0].dot(m_vnphalf.segment(rod0base, rod0nvdof)) + posnn1[0].dot(m_vnphalf.segment(rod1base, rod1nvdof));
    assert(approxEq(lgrhs(0), eecol.computeRelativeVelocity(), 1.0e-6));
    assert(lgrhs(0) < 0.0);

    // Entries 1...numalpha
    for (int i = 1; i < nc0; ++i)
    {
        assert(posnn0[i].size() == m_vnphalf.segment(rod0base, rod0nvdof).size());
        lgrhs(i) = posnn0[i].dot(m_vnphalf.segment(rod0base, rod0nvdof));
    }
    // Entries numalpha...numalpha+numbeta
    for (int i = 1; i < nc1; ++i)
    {
        assert(posnn1[i].size() == m_vnphalf.segment(rod1base, rod1nvdof).size());
        lgrhs(nc0 + i - 1) = posnn1[i].dot(m_vnphalf.segment(rod1base, rod1nvdof));
    }
    //TraceStream(g_log, "") << lgrhs << '\n';
    lgrhs *= -1.0;

    // For now, model an inelastic collision
    //lglhs0(0) += 0.0;
    // Entries 1...numalpha
    lgrhs.segment(1, nc0 - 1) += cval0.segment(1, nc0 - 1);
    // Entries numalpha...numalpha+numbeta
    lgrhs.segment(nc0, nc1 - 1) += cval1.segment(1, nc1 - 1);
    //TraceStream(g_log, "") << lgrhs << '\n';

    Eigen::VectorXd alpha(numalpha);
    assert(lglhs.rows() == lglhs.cols());
    assert(lglhs.rows() == lgrhs.size());
    assert(lglhs.rows() == alpha.size());
    alpha = lglhs.lu().solve(lgrhs);

    assert(alpha(0) >= 0.0);

    if (alpha(0) < 0.0)
        DebugStream(g_log, "") << "WARNING NEGATIVE LAGRANGE MULTIPLIER alpha: " << alpha << '\n';

    // BEGIN TEMP
    //if( rod0 == 23 || rod1 == 23 ) DebugStream(g_log, "") << "Collision lagrange multiplier: " << alpha(0) << '\n';
    // END TEMP

    // Add the impulse corresponding to the collision to rod 0
    assert(m_vnphalf.segment(rod0base, rod0nvdof).size() == posnntilde0[0].size());
    m_vnphalf.segment(rod0base, rod0nvdof) += alpha(0) * posnntilde0[0];
    assert(m_vnphalf.segment(rod1base, rod1nvdof).size() == posnntilde1[0].size());
    m_vnphalf.segment(rod1base, rod1nvdof) += alpha(0) * posnntilde1[0];
    // Add the rod 0 specific impulses to rod 0
    for (int i = 1; i < nc0; ++i)
    {
        assert(m_vnphalf.segment(rod0base, rod0nvdof).size() == posnntilde0[i].size());
        m_vnphalf.segment(rod0base, rod0nvdof) += alpha(i) * posnntilde0[i];
    }
    // Add the rod 1 specific impulses to rod 1
    for (int i = 1; i < nc1; ++i)
    {
        assert(m_vnphalf.segment(rod1base, rod1nvdof).size() == posnntilde1[i].size());
        m_vnphalf.segment(rod1base, rod1nvdof) += alpha(nc0 + i - 1) * posnntilde1[i];
    }

    // Ensure that the impulse eliminated the realtive velocity
#ifndef NDEBUG
    //Vec3d postrelvel = computeRelativeVelocity( m_vnphalf, eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t );
    //double postmagrelvel = postrelvel.dot(eecol.n);
    // Ensure the inelastic impulse decreased the realtive velocity
    //assert( fabs(postmagrelvel) <= fabs(magrelvel) );
    // Should add some 'extra kick' to ensure collision gets killed, but for now just be content with small velocity
    //assert( fabs(postmagrelvel) < 1.0e-9 );
    //TraceStream(g_log, "") << magrelvel << "   " << postmagrelvel << '\n';
    //assert( postmagrelvel < 0.0 );
#endif

    // Ensure the 'scripted vertices' achieved the desired velocity for rod 0
#ifndef NDEBUG
    curdof = 1;
    for (size_t i = 0; i < scriptedverts0.size(); ++i)
    {
        assert(scriptedverts0[i] >= 0);
        assert((int) scriptedverts0[i] < m_rods[rod0]->nv());
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+0) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rod0base + 3 * scriptedverts0[i] + 0), cval0(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+1) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rod0base + 3 * scriptedverts0[i] + 1), cval0(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+2) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rod0base + 3 * scriptedverts0[i] + 2), cval0(curdof), 1.0e-6));
        ++curdof;
    }
#endif

    // Ensure the 'scripted vertices' achieved the desired velocity for rod 1
#ifndef NDEBUG
    curdof = 1;
    for (size_t i = 0; i < scriptedverts1.size(); ++i)
    {
        assert(scriptedverts1[i] >= 0);
        assert((int) scriptedverts1[i] < m_rods[rod1]->nv());
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+0) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rod1base + 3 * scriptedverts1[i] + 0), cval1(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+1) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rod1base + 3 * scriptedverts1[i] + 1), cval1(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+2) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rod1base + 3 * scriptedverts1[i] + 2), cval1(curdof), 1.0e-6));
        ++curdof;
    }
#endif
}

void BAGroomingStepper::exertCompliantInelasticEdgeEdgeImpulseOneFixed(const EdgeEdgeCTCollision& eecol)
{
    // Must have one totally fixed and one totally free edge
    // assert(
    //         (YAEdge(eecol.e0_v0, eecol.e0_v1).IsFree(m_geodata) && YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata))
    //                 || (YAEdge(eecol.e1_v0, eecol.e1_v1).IsFree(m_geodata) && YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(
    //                         m_geodata)));

    // TraceStream(g_log, "") << "BAGroomingStepper::exertCompliantInelasticEdgeEdgeImpulseOneFixed: (x[" << eecol.e0_v0 << "]="
    //        << m_geodata.GetPoint(eecol.e0_v0) << " - x[" << eecol.e0_v1 << "]=" << m_geodata.GetPoint(eecol.e0_v1) << ",   x["
    //        << eecol.e1_v0 << "]=" << m_geodata.GetPoint(eecol.e1_v0) << " - x[" << eecol.e1_v1 << "]=" << m_geodata.GetPoint(
    //       eecol.e1_v1) << ")" << " s=" << eecol.s << " t=" << eecol.t << '\n';

    // Determine if either edge is fixed
    bool rod0fixed = YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(m_geodata);
    bool rod1fixed = YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata);
    assert(rod0fixed != rod1fixed);

    double preRelativeVelocity = eecol.computeRelativeVelocity();

    // Compute the relative velocity, which must be negative for a collision to have happened
    //    Vec3d relvel = m_geodata.computeRelativeVelocity(eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t);
    //    double magrelvel = relvel.dot(eecol.n);
    //   double magrelvel = eecol.computeRelativeVelocity(m_geodata);
    assert(eecol.computeRelativeVelocity() < 0.0);

    // Extract the rod index, barycentric coordinate of the collision, and indices of the edge involved in the collision
    int rodidx = -1;
    double u = -1.0;
    int v0 = -1;
    int v1 = -1;
    if (!rod0fixed)
    {
        rodidx = getContainingRod(eecol.e0_v0);
        u = eecol.s;
        v0 = eecol.e0_v0;
        v1 = eecol.e0_v1;
    }
    if (!rod1fixed)
    {
        rodidx = getContainingRod(eecol.e1_v0);
        u = eecol.t;
        v0 = eecol.e1_v0;
        v1 = eecol.e1_v1;
    }
    assert(u >= 0.0);
    assert(u <= 1.0);
    assert(rodidx >= 0);
    assert(rodidx < (int) m_number_of_rods);

    // If the rod has not solved properly, no need to compute its collision response
    if (!m_steppers[rodidx]->HasSolved())
    {
        DebugStream(g_log, "") << "WARNING: attempt to do edge-edge collision with non-dependable rod\n";
        return;
    }

    TraceStream(g_log, "") << "CTcollision: rod " << rodidx << " vs. mesh\n";

    // Convert the vertices' global indices to rod indices
    assert(m_base_dof_indices[rodidx] % 3 == 0);
    v0 = v0 - m_base_dof_indices[rodidx] / 3;
    v1 = v1 - m_base_dof_indices[rodidx] / 3;
    assert(v0 >= 0);
    assert(v0 < m_rods[rodidx]->nv());
    assert(v1 >= 0);
    assert(v1 < m_rods[rodidx]->nv());

    // Get storage for lhs of linear system, get a solver
    LinearSystemSolver* lss = m_solver_collection.getLinearSystemSolver(m_rods[rodidx]->ndof());
    MatrixBase* lhs = lss->m_lhs;
    assert(lhs != NULL);
    LinearSolverBase* solver = lss->m_solver;
    assert(solver != NULL);

    // Compute M - h^2*dF/dx
    computeCompliantLHS(lhs, rodidx);

    int rodbase = m_base_dof_indices[rodidx];

    // Compute the 'base' index of each vertex
    int base0 = m_rods[rodidx]->vertIdx(v0, 0);
    assert(base0 % 4 == 0);
    int base1 = m_rods[rodidx]->vertIdx(v1, 0);
    assert(base1 % 4 == 0);

    int ndof = m_rods[rodidx]->ndof();
    int nvdof = 3 * m_rods[rodidx]->nv();

    // Compute the 'global' normal
    std::vector<VecXd> n;

    n.push_back(VecXd(ndof));
    n.back().setZero();
    n.back().segment<3> (base0) = (1 - u) * eecol.GetNormal();
    n.back().segment<3> (base1) = u * eecol.GetNormal();

    // If rod 0 has the free edge, need to invert sign of normal
    if (!rod0fixed)
    {
        n.back().segment<3> (base0) *= -1.0;
        n.back().segment<3> (base1) *= -1.0;
    }

    // Compute three normals per scripted vertex
    const std::vector<int>& scriptedverts = m_rods[rodidx]->getBoundaryCondition()->scriptedVertices();
    for (size_t i = 0; i < scriptedverts.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            //TraceStream(g_log, "") << m_rods[rodidx]->vertIdx(scriptedverts[i],j) << " ";
            n.push_back(VecXd(ndof));
            n.back().setZero();
            // Zero the velocity of this particular dof of this particular vertex
            n.back()(m_rods[rodidx]->vertIdx(scriptedverts[i], j)) = 1.0;
        }
    }
    //TraceStream(g_log, "") << '\n';

    int numconstraints = (int) (n.size());

    // Determine the desired values for each constraint
    VecXd desired_values(numconstraints);
#ifndef NDEBUG
    desired_values.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    // First constraint, the collision constraint, is 0
    desired_values(0) = 0.0;
    // The 'scripted vertex' constraints are the respective components of the velocity
    int curdof = 1;
    for (size_t i = 0; i < scriptedverts.size(); ++i)
    {
        assert(scriptedverts[i] >= 0);
        assert((int) scriptedverts[i] < m_rods[rodidx]->nv());
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 0);
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 1);
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 2);
    }
    assert(curdof == numconstraints);
    assert(desired_values.size() == numconstraints);

    //TraceStream(g_log, "") << desired_values << '\n';

    //for( size_t i = 0; i < vertexConstraints.size(); ++i ) TraceStream(g_log, "") << vertexConstraints[i] << '\n';

    // Ensure collision constraint adds up to actual normal
#ifndef NDEBUG
    Vec3d testn = n[0].segment<3> (base0) + n[0].segment<3> (base1);
    if (!rod0fixed)
        testn *= -1.0;
    assert(approxEq(testn, eecol.GetNormal(), 1.0e-6));
#endif

    // Currently, all fixed vertex constraints normalized
#ifndef NDEBUG
    for (int i = 1; i < numconstraints; ++i)
        assert(approxEq(n[i].norm(), 1.0, 1.0e-9));
#endif

    std::vector<VecXd> ntilde;
    for (int i = 0; i < numconstraints; ++i)
    {
        ntilde.push_back(VecXd(ndof));
        ntilde.back().setZero();
        int status = solver->solve(ntilde.back(), n[i]);
        if (status < 0)
            DebugStream(g_log, "")
                    << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. exertCompliantInelasticEdgeEdgeImpulseOneFixedThree. Time: "
                    << m_t << '\n';
    }

    //for( size_t i = 0; i < n.size(); ++i ) TraceStream(g_log, "") << ntilde[i].transpose() << '\n';


    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (ElasticRod::edge_iter eit = m_rods[rodidx]->edges_begin(); eit != m_rods[rodidx]->edges_end(); ++eit)
    {
        for (int i = 0; i < numconstraints; ++i)
            assert(ntilde[i](m_rods[rodidx]->edgeIdx(*eit)) == 0.0);
    }
#endif

    // Vectors restriced to just verted DoFs. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn;
    for (int i = 0; i < numconstraints; ++i)
        posnn.push_back(VecXd(nvdof));
#ifndef NDEBUG
    for (int i = 0; i < numconstraints; ++i)
        posnn[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    std::vector<VecXd> posnntilde;
    for (int i = 0; i < numconstraints; ++i)
        posnntilde.push_back(VecXd(nvdof));
#ifndef NDEBUG
    for (int i = 0; i < numconstraints; ++i)
        posnntilde[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (int i = 0; i < numconstraints; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rodidx]->vertices_begin(); vit != m_rods[rodidx]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rodidx]->vertIdx(*vit, 0);
            posnn[i].segment<3> (3 * theidx) = n[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    for (int i = 0; i < numconstraints; ++i)
    {
        int theidx = 0;
        for (ElasticRod::vertex_iter vit = m_rods[rodidx]->vertices_begin(); vit != m_rods[rodidx]->vertices_end(); ++vit)
        {
            int vert_dof = m_rods[rodidx]->vertIdx(*vit, 0);
            posnntilde[i].segment<3> (3 * theidx) = ntilde[i].segment<3> (vert_dof);
            ++theidx;
        }
    }

    // Compute the Lagrange multipliers
    // TODO: This matrix is symmetric, exploit that fact to avoid computations
    MatXd lglhs(numconstraints, numconstraints);
#ifndef NDEBUG
    lglhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    for (int i = 0; i < numconstraints; ++i)
        for (int j = 0; j < numconstraints; ++j)
            lglhs(i, j) = posnn[i].dot(posnntilde[j]);
    assert(approxSymmetric(lglhs, 1.0e-6));

    Eigen::VectorXd lgrhs(numconstraints);
#ifndef NDEBUG
    lgrhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    lgrhs(0) = eecol.computeRelativeVelocity(); // CRAZY SIGN PROBLEM HERE?
    assert(lgrhs(0) < 0.0);
    for (int i = 1; i < numconstraints; ++i)
        lgrhs(i) = posnn[i].dot(m_vnphalf.segment(rodbase, nvdof));
    lgrhs *= -1.0;

    //TraceStream(g_log, "") << lgrhs << '\n';

    lgrhs += desired_values;

    //TraceStream(g_log, "") << lgrhs << '\n';
    //TraceStream(g_log, "") << desired_values << '\n';

    Eigen::VectorXd alpha(numconstraints); // = lglhs.inverse()*lgrhs;
    alpha = lglhs.lu().solve(lgrhs);

    // Contact constraint should 'push not pull'
    assert(alpha(0) >= 0.0);

    // TEMP DEBUG STUFF
    //  if( rodidx == 20 && m_t >= 0.85 )
    //  {
    //    TraceStream(g_log, "") << "---------- ROD 20 ROD 20 ROD 20 ROD 20 ROD 20 ----------" << '\n';
    //    TraceStream(g_log, "") << "alpha: " << alpha << '\n';
    //    TraceStream(g_log, "") << "time: " << m_t << '\n';
    //
    //    if( !rod0fixed ) TraceStream(g_log, "")  << "Free: " << 0 << '\n';
    //    if( !rod1fixed ) TraceStream(g_log, "")  << "Free: " << 1 << '\n';
    //
    //    TraceStream(g_log, "") << "temprelvel: " << posnn[0].dot(m_vnphalf.segment(rodbase,nvdof)) << '\n';
    //  }
    // END TEMP DEBUG STUFF

    // TraceStream(g_log, "") << "Pre-impulse velocity = " << m_vnphalf.segment(rodbase, nvdof) << '\n' << '\n';

    // for (int i = 0; i < nvdof; i += 3)
    // {

    //    double vx = m_vnphalf.segment(rodbase, nvdof)[i];
    //    double vy = m_vnphalf.segment(rodbase, nvdof)[i + 1];
    //   double vz = m_vnphalf.segment(rodbase, nvdof)[i + 2];

    //    double px = m_xn.segment(rodbase, nvdof)[i];
    //    double py = m_xn.segment(rodbase, nvdof)[i + 1];
    //   double pz = m_xn.segment(rodbase, nvdof)[i + 2];

    // TraceStream(g_log, "") << "curve -p " << px << " " << py << " " << pz << " -p " << px + vx << " " << py + vy << " " << pz + vz
    //         << ";" << '\n';
    // }
    // TraceStream(g_log, "") << '\n';

    for (int i = 0; i < numconstraints; ++i)
    {
        //  TraceStream(g_log, "") << "alpha[" << i << "]=" << alpha(i) << " and " << " posnntilde[" << i << "] = " << posnntilde[i]
        //          << '\n' << '\n';

        m_vnphalf.segment(rodbase, nvdof) += alpha(i) * posnntilde[i];
    }

    // TraceStream(g_log, "") << "Post-impulse velocity = " << m_vnphalf.segment(rodbase, nvdof) << '\n' << '\n';

    // for (int i = 0; i < nvdof; i += 3)
    // {
    // double vx = m_vnphalf.segment(rodbase, nvdof)[i];
    // double vy = m_vnphalf.segment(rodbase, nvdof)[i + 1];
    // double vz = m_vnphalf.segment(rodbase, nvdof)[i + 2];

    // double px = m_xn.segment(rodbase, nvdof)[i];
    // double py = m_xn.segment(rodbase, nvdof)[i + 1];
    // double pz = m_xn.segment(rodbase, nvdof)[i + 2];

    // TraceStream(g_log, "") << "curve -p " << px << " " << py << " " << pz << " -p " << px + vx << " " << py + vy << " " << pz + vz
    //         << ";" << '\n';
    //}
    // TraceStream(g_log, "") << '\n';

    applyInextensibilityVelocityFilter(rodidx);

    // ////////////////////////////////////////////////////////
    // // Inextensibility Filter
    // // for my reference: m_xnp1 = m_xn + m_dt * m_vnphalf;
    // //

    // {
    //   Vec3d v0  = m_vnphalf.segment<3> (rodbase+3); // velocity of vertex 1
    //   Vec3d x0  = m_xn.segment<3>      (rodbase+3); // start-of-step position of vertex 1
    //   Vec3d x0N = x0 + m_dt * v0;                   // end-of-step position of vertex 1
    //   for (int i=2; i<m_rods[rodidx]->nv(); ++i)
    //  {
    //    Vec3d x1 = m_xn.segment<3>      (rodbase + 3*i);
    //    Vec3d v1 = m_vnphalf.segment<3> (rodbase + 3*i);

    //    // filter the velocity of vertex i

    //    Vec3d x1N = x1 + m_dt * v1;

    //    double l  = (x1 -x0 ).norm();
    //    double lN = (x1N-x0N).norm();

    //    Vec3d x1Nrevised = x0N + (x1N-x0N)*l/lN;

    //    double lNrevised = (x1Nrevised-x0N).norm();

    //    Vec3d v1revised = (x1Nrevised - x1) / m_dt;

    //    m_vnphalf.segment<3> (rodbase + 3*i) = v1revised;

    //    TraceStream(g_log, "") << "inextensibility: x0N = " << x0N << " x1Nrevised = " << x1Nrevised << " l = " << l << " lNrevised-l = " << (lNrevised-l) << " lNrevised/l = " << (lNrevised/l) << '\n';//x1Nrevised = " << x1Nrevised << " strain = " << (lN/l) << " revised: " << (lNrevised/l) << " l = " << l << " lN = " << lN << " lNrevised = " << lNrevised << " x0N = " << x0N << " x1Nrevised = " << x1Nrevised << " m_vnphalf revised = " << m_vnphalf.segment<3> (rodbase + 3*i) << '\n';

    //    // Vec3d t = (x1-x0).normalized(); // unit tangent

    //    // double relvel = t.dot(v1-v0);

    //    // Vec3d impulse = -relvel*t;

    //    // double relvelAfter = t.dot(v1 + impulse - v0);

    //    // m_vnphalf.segment<3> (rodbase + 3*i) = v1 + impulse;

    //    // TraceStream(g_log, "") << "inextensibility: i=" << i
    //    //        << " x0 = " << x0
    //    //        << " x1 = " << x1
    //    //        << " v0 = " << v0
    //    //        << " v1 = " << v1
    //    //           << " t = " << t
    //    //           << " relvel before = " << relvel
    //    //           << " after = " << relvelAfter << '\n';

    //    x0N = x1Nrevised;
    //    x0  = x1;
    //    v0  = v1;
    //  }
    // }


    ////////////////////////////////////////////////////////


    // TEMP DEBUG STUFF
    //  Vec3d postrelvel = computeRelativeVelocity( m_vnphalf, eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t );
    //  double postmagrelvel = postrelvel.dot(eecol.n);
    //  // Ensure the inelastic impulse decreased the realtive velocity
    //  if( fabs(postmagrelvel) > fabs(magrelvel) )
    //  {
    //    TraceStream(g_log, "") << "Rod: " << rodidx << '\n';
    //    TraceStream(g_log, "") << "Time is: " << m_t << '\n';
    //    TraceStream(g_log, "") << "Incoming relative velocity: " << magrelvel << "      Post-impulse relative velocity: " << postmagrelvel << '\n';
    //
    //    exit(1);
    //  }
    // END TEMP DEBUG STUFF


    // Ensure that the impulse eliminated the realtive velocity
#ifndef NDEBUG
    //Vec3d postrelvel = computeRelativeVelocity( m_vnphalf, eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t );
    //double postmagrelvel = postrelvel.dot(eecol.n);
    // Ensure the inelastic impulse decreased the realtive velocity
    // if( fabs(postmagrelvel) > fabs(magrelvel) )
    // {
    //     TraceStream(g_log, "") << "Rod: " << rodidx << '\n';
    //     TraceStream(g_log, "") << "Time is: " << m_t << '\n';
    //     TraceStream(g_log, "") << "Incoming relative velocity: " << magrelvel << "      Post-impulse relative velocity: " << postmagrelvel << '\n';
    // }
    // assert( fabs(postmagrelvel) <= fabs(magrelvel) );
    // // Should add some 'extra kick' to ensure collision gets killed, but for now just be content with small velocity
    // assert( fabs(postmagrelvel) < 1.0e-9 );
    //assert( postmagrelvel < 0.0 );
#endif

    // Ensure the 'scripted vertices' achieved the desired velocity
#ifndef NDEBUG
    curdof = 1;
    for (size_t i = 0; i < scriptedverts.size(); ++i)
    {
        assert(scriptedverts[i] >= 0);
        assert((int) scriptedverts[i] < m_rods[rodidx]->nv());
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+0) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rodbase + 3 * scriptedverts[i] + 0), desired_values(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+1) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rodbase + 3 * scriptedverts[i] + 1), desired_values(curdof), 1.0e-6));
        ++curdof;
        //TraceStream(g_log, "") << m_vnphalf(rodbase+3*scriptedverts[i]+2) << "   " << desired_values(curdof) << '\n';
        assert(approxEq(m_vnphalf(rodbase + 3 * scriptedverts[i] + 2), desired_values(curdof), 1.0e-6));
        ++curdof;
    }
#endif

    double postRelativeVelocity = eecol.computeRelativeVelocity();

    //  TraceStream(g_log, "") << "BAGroomingStepper::exertCompliantInelasticEdgeEdgeImpulseOneFixed: Relative velocity before = "
    //         << preRelativeVelocity << " after = " << postRelativeVelocity << '\n';
}

bool BAGroomingStepper::checkExplosions(std::vector<bool>& exploding_rods, const std::vector<bool>& failed_collisions_rods,
        const RodSelectionType& selected_rods)
{
    bool explosions_detected = false;
    double maxRate = 0;
    double maxStart = 0;
    double minStart = std::numeric_limits<double>::max();
    int worstViolator = std::numeric_limits<int>::signaling_NaN();
    TraceStream(g_log, "") << "Checking for explosions\n";

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
    {
        if (true || m_steppers[*rod]->HasSolved() && !failed_collisions_rods[*rod])
        {
            TraceStream(g_log, "") << "Rod " << *rod << " Force norms: initial: " << (*(m_startForces[*rod])).norm()
                    << " pre-dynamic: " << (*(m_preDynamicForces[*rod])).norm() << " pre-collisions: "
                    << (*(m_preCollisionForces[*rod])).norm() << " post-collisions: " << (*(m_endForces[*rod])).norm() << "\n";
            for (int j = 0; j < m_rods[*rod]->ndof(); ++j)
            {
                const double s = (*(m_startForces[*rod]))[j];
                const double p = (*(m_preCollisionForces[*rod]))[j];
                const double e = (*(m_endForces[*rod]))[j];
                const double rate = fabs(s - e) / (fabs(s) + m_perf_param.m_explosion_damping);
                maxRate = std::max(maxRate, rate);
                minStart = std::min(fabs(s), minStart);
                maxStart = std::max(fabs(s), maxStart);
                if (maxRate == rate)
                    worstViolator = j;
                if (isnan(rate) || rate > m_perf_param.m_explosion_threshold)
                {
                    explosions_detected = true;
                    exploding_rods[*rod] = true;
                    TraceStream(g_log, "") << "Rod number " << *rod << " had an explosion during collision response: s = " << s
                            << " p = " << p << " e = " << e << " rate = " << rate << " \n";
                }
            }
        }
    }
    if (explosions_detected)
        DebugStream(g_log, "") << "Some rods had explosions\n";

    return explosions_detected;
}

bool BAGroomingStepper::checkLengths(std::vector<bool>& stretching_rods)
{
    bool stretching_detected = false;
    TraceStream(g_log, "") << "Checking for lengths\n";

    for (RodSelectionType::iterator rod = m_simulated_rods.begin(); rod != m_simulated_rods.end(); rod++)
    {
        stretching_rods[*rod] = !checkLength(*rod);

        if (stretching_rods[*rod])
        {
            // Declare the rod collision-immune for the rest of the time step.
            for (int j = 0; j < m_rods[*rod]->nv(); ++j)
                m_collision_immune[m_base_vtx_indices[*rod] + j] = false;
        }

        stretching_detected = stretching_detected || stretching_rods[*rod];
    }

    if (stretching_detected)
        DebugStream(g_log, "") << "Some rods were stretched by a factor > " << m_perf_param.m_stretching_threshold << '\n';

    return stretching_detected;
}

bool BAGroomingStepper::checkLength(int rodIdx)
{
    double length = 0.0;

    for (int j = 1; j < m_rods[rodIdx]->nv(); j++)
        length += (m_rods[rodIdx]->getVertex(j) - m_rods[rodIdx]->getVertex(j - 1)).norm();

    if (length > m_initialLengths[rodIdx] * m_perf_param.m_stretching_threshold)
    {
        TraceStream(g_log, "") << "Rod number " << rodIdx << " was stretched by a factor " << length / m_initialLengths[rodIdx];
        return false;
    }

    return true;
}

void BAGroomingStepper::applyInextensibilityVelocityFilter(int rodidx)
{

    // TraceStream(g_log, "") << "WARNING! SKIPPING INEXTENSIBILITY FILTER!" << '\n';

    // return;


    if (m_level < m_perf_param.m_inextensibility_threshold)
        return;

    int rodbase = m_base_dof_indices[rodidx];

    // if (boundary->isVertexScripted(j))
    //   {
    //  //TraceStream(g_log, "") << "BridsonTimeStepper is calling RodBoundaryCondition at m_t = " << m_t << '\n';
    //  Vec3d desiredposition = boundary->getDesiredVertexPosition(j, m_t);
    //  Vec3d actualvalue = m_xnp1.segment<3> (rodbase + 3 * j);
    //  assert(approxEq(desiredposition, actualvalue, 1.0e-6));
    //   }

    // DebugStream(g_log, "") << "Velocity Filter begin: rod " << rodidx << " vertex 0 = " << m_xn.segment<3>(rodbase) << '\n';

    // m_vnphalf.segment<3> (rodbase + 0) = Vec3d(0, 0, 0);
    // m_vnphalf.segment<3> (rodbase + 3) = Vec3d(0, 0, 0);

    {
        Vec3d v0 = m_vnphalf.segment<3> (rodbase + 3); // velocity of vertex 1
        Vec3d x0 = m_xn.segment<3> (rodbase + 3); // start-of-step position of vertex 1
        Vec3d x0N = x0 + m_dt * v0; // end-of-step position of vertex 1
        for (int i = 2; i < m_rods[rodidx]->nv(); ++i)
        {
            Vec3d x1 = m_xn.segment<3> (rodbase + 3 * i);
            Vec3d v1 = m_vnphalf.segment<3> (rodbase + 3 * i);

            // filter the velocity of vertex i

            Vec3d x1N = x1 + m_dt * v1;

            double l = (x1 - x0).norm();
            double lN = (x1N - x0N).norm();

            Vec3d x1Nrevised = x0N + (x1N - x0N) * l / lN;

            double lNrevised = (x1Nrevised - x0N).norm();

            Vec3d v1revised = (x1Nrevised - x1) / m_dt;

            m_vnphalf.segment<3> (rodbase + 3 * i) = v1revised;

            //TraceStream(g_log, "") << "inextensibility: x0N = " << x0N << " x1N = " << x1N << " x1Nrevised = " << x1Nrevised << " l = " << l << " lNrevised = "
            //            << lNrevised << " l = " << l << '\n';

            //x1Nrevised = " << x1Nrevised << " strain = " << (lN/l) << " revised: " << (lNrevised/l) << " l = " << l << " lN = " << lN << " lNrevised = " << lNrevised << " x0N = " << x0N << " x1Nrevised = " << x1Nrevised << " m_vnphalf revised = " << m_vnphalf.segment<3> (rodbase + 3*i) << '\n';

            // Vec3d t = (x1-x0).normalized(); // unit tangent

            // double relvel = t.dot(v1-v0);

            // Vec3d impulse = -relvel*t;

            // double relvelAfter = t.dot(v1 + impulse - v0);

            // m_vnphalf.segment<3> (rodbase + 3*i) = v1 + impulse;

            // TraceStream(g_log, "") << "inextensibility: i=" << i
            //      << " x0 = " << x0
            //      << " x1 = " << x1
            //      << " v0 = " << v0
            //      << " v1 = " << v1
            //           << " t = " << t
            //           << " relvel before = " << relvel
            //           << " after = " << relvelAfter << '\n';

            x0N = x1Nrevised;
            x0 = x1;
            v0 = v1;
        }
    }

    // TraceStream(g_log, "") << "Edge lengths after inextensibility: ";
    // for (int i = 0; i < m_rods[rodidx]->nv() - 1; ++i)
    // {
    //     Vec3d x0 = m_xn.segment<3> (rodbase + 3 * i);
    //     Vec3d v0 = m_vnphalf.segment<3> (rodbase + 3 * i);
    //     Vec3d x1 = m_xn.segment<3> (rodbase + 3 * i + 3);
    //     Vec3d v1 = m_vnphalf.segment<3> (rodbase + 3 * i + 3);
    //  Vec3d x0N = x0 + m_dt * v0;
    //  Vec3d x1N = x1 + m_dt * v1;
    //  Vec3d eN = (x1N - x0N);
    //  TraceStream(g_log, "") << " " << eN.norm();
    // }
    // TraceStream(g_log, "") << '\n';

    // DebugStream(g_log, "") << "Velocity Filter end: rod " << rodidx << " vertex 0 = " << m_xn.segment<3>(rodbase) << '\n';


}

/**
 * Implicit penalty response
 */
void BAGroomingStepper::setupPenaltyForces(std::list<Collision*>& collisions, const RodSelectionType& selected_rods)
{
    // // Detect proximity collisions
    // m_collision_detector->getCollisions(collisions, Proximity);

    // for (int rod_id = 0; rod_id < m_number_of_rods; rod_id++)
    //     assert(m_implicit_pnlty_forces[rod_id]->cleared());

    // // Store the proximity collision in the RodPenaltyForce
    // for (std::list<Collision*>::const_iterator col = collisions.begin(); col != collisions.end(); col++)
    // {
    //     VertexFaceProximityCollision* vfpcol = dynamic_cast<VertexFaceProximityCollision*> (*col);
    //     if (vfpcol)
    //     {
    //         assert(vfpcol->isAnalysed());
    //         int rod_id = getContainingRod(vfpcol->v0);
    //         TraceStream(g_log, "") << "Creating penalty force for rod " << rod_id << " address " << &m_rods[rod_id] << '\n';
    //         int v_id = vfpcol->v0 - m_base_dof_indices[rod_id] / 3;
    //         m_implicit_pnlty_forces[rod_id]->registerProximityCollision(v_id, vfpcol);
    //     }
    //     // TODO: delete vfpcol once used. Not here though, as vfpcol is used each time RodPenaltyForce::computeForce is called
    // }

}

void BAGroomingStepper::setImplicitPenaltyExtraThickness(const double& h)
{
    m_perf_param.m_implicit_thickness = h;
}

void BAGroomingStepper::setVertexFacePenalty(const double& k)
{
    assert(k >= 0.0);

    m_perf_param.m_implicit_stiffness = k;
}

class RootDistanceComparer
{
    const ElasticRod* m_rod;

public:

    RootDistanceComparer(const ElasticRod* rod) :
        m_rod(rod)
    {
    }

    int operator()(const ElasticRod* rod1, const ElasticRod* rod2)
    {
        return (m_rod->getVertex(0) - rod1->getVertex(0)).norm() < (m_rod->getVertex(0) - rod2->getVertex(0)).norm();
    }
};

void BAGroomingStepper::activateClumpingForce()
{
    selectClumps();

    m_clumpingForce = new RodClumpingForce();

    for (std::vector<GroomingTimeStepper*>::iterator stepper = m_steppers.begin(); stepper != m_steppers.end(); ++stepper)
    {
        (*stepper)->addExternalForce(m_clumpingForce);
    }
}

Scalar SquareDistance(const BoundingBox<Scalar>& bbox, const BoundingBox<Scalar>::PointType point)
{
    Scalar distx = std::max(bbox.min[0] - point[0], 0.0) + std::max(point[0] - bbox.max[0], 0.0);
    Scalar disty = std::max(bbox.min[1] - point[1], 0.0) + std::max(point[1] - bbox.max[1], 0.0);
    Scalar distz = std::max(bbox.min[2] - point[2], 0.0) + std::max(point[2] - bbox.max[2], 0.0);

    return distx * distx + disty * disty + distz * distz;
}

template<typename TreeT>
class RodPointDistance
{
public:
    RodPointDistance(const VecXd& xn, const TreeT& bvh, const Vec3d& point, const std::vector<int>& centerLinePointIndices) :
        m_xn(xn), m_bvh(bvh), m_point(point), m_pointIndices(centerLinePointIndices)
    {
    }

    float NodeDistance(const uint32_t node_index) const
    {
        return SquareDistance(m_bvh.GetNode(node_index).BBox(), m_point);
    }

    float ObjectDistance(const uint32_t obj_index) const
    {
        const Vec3d v = m_xn.segment<3> (3 * m_pointIndices[obj_index]) - m_point;
        return v.dot(v);
    }

private:
    const VecXd& m_xn;
    const TreeT& m_bvh;
    const Vec3d& m_point;
    const std::vector<int>& m_pointIndices;
};

void BAGroomingStepper::findCenterLines(RodSelectionType& centerLineRods)
{
    int decimationTick = m_rods_per_clump;
    for (RodSelectionType::const_iterator rod = m_simulated_rods.begin(); rod != m_simulated_rods.end(); rod++, decimationTick++)
        if (decimationTick == m_rods_per_clump)
        {
            decimationTick = 0;
            centerLineRods.push_back(*rod);
            m_rods[*rod]->setIsClumpCenterLine(true);
        }
        else
        {
            m_rods[*rod]->setIsClumpCenterLine(false);
        }
}

/**
 * Select clumps based on closest rods to a selected few (for now simple decimation of the whole rods set)
 */
void BAGroomingStepper::selectClumps()
{
    extractPositions(m_xn, m_simulated_rods, m_t - m_dt); // Probably superfluous but since we are going to access positions through m_xn let's make double sure it's up-to-date

    // Those will be the clump "center lines".
    RodSelectionType centerLineRods;
    findCenterLines(centerLineRods);

    // Build a vector of (initially sorted) indices (in m_xn) of centerLineRods points
    std::vector<int> centerLinePointIndices;
    for (RodSelectionType::const_iterator rod = centerLineRods.begin(); rod != centerLineRods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
            centerLinePointIndices.push_back(m_base_vtx_indices[*rod] + j);

    // Build a bvh (this reorders rootPointIndices) and the kNN finder
    pantaray::kNN<BVH> knn;
    BVH bvh;
    SimplePointBBoxFunctor bboxes(centerLinePointIndices, m_xn);
    BVHBuilder<SimplePointBBoxFunctor> bvhbuilder;
    bvhbuilder.build(bboxes, &bvh);
    knn.Setup(bvh);

    // For each rod tip find the closest clump center line
    const Scalar searchDistance = std::numeric_limits<Scalar>::max();
    for (RodSelectionType::const_iterator rod = m_simulated_rods.begin(); rod != m_simulated_rods.end(); rod++)
    {
        if (find(centerLineRods.begin(), centerLineRods.end(), *rod) == centerLineRods.end()) // skip if the rod is a centerline
        {
            // const Vec3d& bindingPoint = m_rods[*rod]->getVertex(0); // at the root
            const Vec3d& bindingPoint = m_rods[*rod]->getVertex(m_rods[*rod]->nv() - 1); // at the tip

            const RodPointDistance<BVH> distance(m_xn, bvh, bindingPoint, centerLinePointIndices);
            const size_t numberFound = knn.Run(distance, -1.0, searchDistance, 1);
            assert(numberFound == 1);

            const int foundGlobalIdx = centerLinePointIndices[knn.TopResult().m_node];
            //  DebugStream(g_log, "") << "Rod " << *rod << " found close point " << foundGlobalIdx << " at distance " << sqrt(
            //          knn.TopResult().m_dist) << '\n';

            //  knn.PopResult(); // Probably superfluous since we are only going to use the top result, just need to check that knn.Run correctly resets the result list.
            const int foundRodIdx = getContainingRod(foundGlobalIdx);
            //    int foundLocalIdx = foundGlobalIdx - m_base_vtx_indices[foundRodIdx];
            //  DebugStream(g_log, "") << "That is vertex " << foundLocalIdx << " on rod " << foundRodIdx << " at distance "
            //           << (m_rods[foundRodIdx]->getVertex(foundLocalIdx) - bindingPoint).norm() << '\n';
            assert(find(centerLineRods.begin(), centerLineRods.end(), foundRodIdx) != centerLineRods.end()); // found rod should be in the center lines list

            std::vector<ElasticRod*> neighbours;
            neighbours.push_back(m_rods[foundRodIdx]);
            m_rods[*rod]->setNearestRootNeighbours(neighbours);

            //  DebugStream(g_log, "") << "Rod " << *rod << " is attracted to rod " << foundRodIdx << '\n';
        }
        //  else
        //      DebugStream(g_log, "") << "Rod " << *rod << " is a centerline\n";
    }

}

/*
 void BAGroomingStepper::updateRodsNeighbours() // OBSOLETE
 {
 extractPositions(m_xn, m_simulated_rods, m_t - m_dt); // Probably superfluous

 // Build a vector of (initially sorted) indices of points in m_xn
 std::vector<int> pointIndices;
 for (RodSelectionType::const_iterator rod = m_simulated_rods.begin(); rod != m_simulated_rods.end(); rod++)
 for (int j = 0; j < m_rods[*rod]->nv(); ++j)
 pointIndices.push_back(m_base_dof_indices[*rod] + 3 * j);

 // Build a bvh (this reorders pointIndices)
 BVH bvh;
 SimplePointBBoxFunctor bboxes(pointIndices, m_xn);
 BVHBuilder<SimplePointBBoxFunctor> bvhbuilder;
 bvhbuilder.build(bboxes, &bvh);
 // Prepare the k nearest neighbours finder
 pantaray::kNN<BVH> knn;
 knn.Setup(bvh);

 const Scalar searchDistance = 10 * m_clumpingForce->getDistance();

 // For each point on each rod, find nearest neighbours on other rods
 for (RodSelectionType::const_iterator rod = m_simulated_rods.begin(); rod != m_simulated_rods.end(); rod++)
 {
 std::set<int> nearRodIndices;
 for (int j = 0; j < m_rods[*rod]->nv(); ++j)
 {
 const RodPointDistance<BVH> distance(m_xn, bvh, m_rods[*rod]->getVertex(j));
 const size_t numberFound = knn.Run(distance, -1.0f, searchDistance, m_rods_per_clump);
 for (size_t i = 0; i < numberFound; i++)
 {
 int foundGlobalIdx = knn.TopResult().m_node;
 knn.PopResult();

 int foundRodIdx = getContainingRod(foundGlobalIdx);
 if (foundRodIdx != *rod)
 nearRodIndices.insert(foundRodIdx);
 }
 }

 //       std::ostringstream debugOstr;

 // Turn the set of nearest rod indices into a vector of ElasticRod*
 std::vector<ElasticRod*> neighbours;
 for (std::set<int>::const_iterator idx = nearRodIndices.begin(); idx != nearRodIndices.end(); ++idx)
 {
 neighbours.push_back(m_rods[*idx]);
 //            debugOstr << *idx << ' ';
 }
 m_rods[*rod]->setNearestRootNeighbours(neighbours);

 //       DebugStream(g_log, "") << "Rod " << *rod << " has neighbours " << debugOstr.str() << '\n';
 }
 }
 */

void BAGroomingStepper::setClumpingParameters(const double charge, const double power, const double dist)
{
    assert(m_clumpingForce != NULL);

    DebugStream(g_log, "") << "Changing clumping parameters\n";

    m_clumpingForce->setCharge(charge);
    m_clumpingForce->setPower(power);
    m_clumpingForce->setDistance(dist);

    selectClumps();
}

}
