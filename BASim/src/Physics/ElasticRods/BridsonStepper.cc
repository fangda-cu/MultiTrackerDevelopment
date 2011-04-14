/**
 * \file BridsonStepper.cc
 *
 * \author smith@cs.columbia.edu
 * \date 02/16/2010
 */

#include <typeinfo>
#include "BridsonStepper.hh"
#include "../../Threads/MultithreadedStepper.hh"
#include "../../Core/Timer.hh"

namespace BASim
{

BridsonStepper::BridsonStepper(std::vector<ElasticRod*>& rods, std::vector<TriangleMesh*>& trimeshes,
        std::vector<ScriptingController*>& scripting_controllers, std::vector<RodTimeStepper*>& steppers, const double& dt,
        const double time, int num_threads) :
    m_num_dof(0), m_rods(rods), m_triangle_meshes(trimeshes), m_scripting_controllers(scripting_controllers),
            m_steppers(steppers), m_dt(dt), m_base_indices(), m_base_triangle_indices(), m_edges(), m_faces(),
            m_vertex_radii(), m_edge_radii(), m_face_radii(), m_masses(),
            m_geodata(m_xn, m_vnphalf, m_vertex_radii, m_masses, m_obj_start), m_respns_enbld(true), m_pnlty_enbld(true),
            m_itrv_inlstc_enbld(true), m_num_inlstc_itrns(10), m_vrt_fc_pnlty(200.0), m_nan_enc(false), m_inf_enc(false),
            m_lt0_enc(false), m_gt0_enc(false), m_collision_detector(NULL), // m_bvh(NULL),
            m_obj_start(-1), m_t(time), m_rod_labels(), m_implicit_pnlty_enbld(false), m_implicit_thickness(1.0)
{
#ifdef DEBUG
    for( int i = 0; i < (int) m_rods.size(); ++i ) assert( m_rods[i] != NULL );
    for( int i = 0; i < (int) m_triangle_meshes.size(); ++i ) assert( m_triangle_meshes[i] != NULL );
    for( int i = 0; i < (int) m_steppers.size(); ++i ) assert( m_steppers[i] != NULL );
#endif
    assert(m_dt > 0.0);

    if (num_threads > 0)
        m_num_threads = num_threads;
    else
        m_num_threads = sysconf(_SC_NPROCESSORS_ONLN);

    // Update internal state, prepare for execution
    prepareForExecution();

#ifdef DEBUG
    // Number of degrees of freedom is non-negative multiple of 3 (3 coords per vertex)
    assert( m_num_dof >= 0 );
    assert( m_num_dof%3 == 0 );

    // Base indices refers to rods
    assert( m_base_indices.size() == m_rods.size() );
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
    IntStatTracker::getIntTracker("CONVERGENCE_FAILURES_PROPAGATED_TO_BRIDSONSTEPPER");
#endif
}

BridsonStepper::~BridsonStepper()
{
    for (std::vector<RodPenaltyForce*>::iterator i = m_implicit_pnlty_forces.begin(); i != m_implicit_pnlty_forces.end(); i++)
        delete *i;
    delete m_collision_detector;
}

// TODO: Check the indices here
void BridsonStepper::prepareForExecution()
{
    delete m_collision_detector;
    m_collision_detector = NULL;

    //std::cout << "About to extract rod information" << std::endl;
    for (int i = 0; i < (int) m_rods.size(); ++i)
    {
        assert(m_rods[i] != NULL);
#ifdef DEBUG
        // Sanity checks for collision detection purpouses
        ensureCircularCrossSection( *m_rods[i] );
        ensureNoCollisionsByDefault( *m_rods[i] );
#endif

        // Extract edges from the new rod
        for (int j = 0; j < (int) m_rods[i]->nv() - 1; ++j)
        {
            m_edges.push_back(std::pair<int, int>(getNumVerts() + j, getNumVerts() + j + 1));
            assert(m_rods[i]->getRadiusScale() * m_rods[i]->radiusA(j) > 0.0);
            m_edge_radii.push_back(m_rods[i]->getRadiusScale() * m_rods[i]->radiusA(j));
        }
        assert(m_edges.size() == m_edge_radii.size());

        for (int j = 0; j < (int) m_rods[i]->nv() - 1; ++j)
        {
            assert(m_rods[i]->getRadiusScale() * m_rods[i]->radiusA(j) > 0.0);
            // Radii associated with edges ... what to do if at a vertex with two radii? Average?
            m_vertex_radii.push_back(m_rods[i]->getRadiusScale() * m_rods[i]->radiusA(j));
        }
        m_vertex_radii.push_back(m_vertex_radii.back()); // < TODO: What the $^#! is this call?

        // Extract masses from the new rod
        for (ElasticRod::vertex_iter itr = m_rods[i]->vertices_begin(); itr != m_rods[i]->vertices_end(); ++itr)
        {
            if (m_rods[i]->getBoundaryCondition()->isVertexScripted((*itr).idx()))
            {
                m_masses.push_back(std::numeric_limits<double>::infinity());
            }
            else
            {
                assert(m_rods[i]->getVertexMass(*itr, -1) > 0.0);
                m_masses.push_back(m_rods[i]->getVertexMass(*itr, -1));
            }
        }

        // Update vector that tracks the rod DOF in the system
        m_base_indices.push_back(getNumDof());

        // Update total number of DOF in the system
        m_num_dof += 3 * m_rods[i]->nv();

        assert((int) m_masses.size() == getNumVerts());
        assert((int) m_vertex_radii.size() == getNumVerts());
    }
    assert(m_rods.size() == m_base_indices.size());
    //std::cout << "Extracted rod information" << std::endl;

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
    m_vnphalf.resize(getNumDof());
    m_vnphalf.setConstant(std::numeric_limits<double>::signaling_NaN());

    //std::cout << "About to extract positions" << std::endl;
    // Load positions for initial construction of the BVH
    extractPositions(m_rods, m_base_indices, m_xn);
    extractVelocities(m_rods, m_base_indices, m_vnphalf);
    //std::cout << "Extracted positions" << std::endl;

    std::cout << "About to create CollisionDetector" << std::endl;
    bool skip_rod_rod = true;
    m_collision_detector = new CollisionDetector(m_geodata, m_edges, m_faces, m_dt, skip_rod_rod, m_num_threads);
    std::cout << "Created CollisionDetector" << std::endl;

}

void BridsonStepper::setRodLabels(const std::vector<std::string>& rod_labels)
{
    assert(rod_labels.size() == m_rods.size());
    m_rod_labels = rod_labels;
}

void BridsonStepper::exertInelasticImpluse(EdgeEdgeCTCollision& eecol)
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

    if (eecol.GetRelativeVelocity() >= 0.0)
    {
        std::cerr
                << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Relative velocity computed \
                  incorrectly before applying edge-edge inelastic impulse (bug in normal \
                  computation?). Magnitude of relative velocity: "
                << eecol.GetRelativeVelocity() << ". Rod: " << getContainingRod(eecol.e0_v0) << std::endl;
    }
    assert(eecol.GetRelativeVelocity() < 0.0);

    // Add some extra "kick" to relative velocity to account for FPA errors
    eecol.ApplyRelativeVelocityKick();
    Vec3d I = eecol.computeInelasticImpulse();
    //computeEdgeEdgeInelasticImpulse(m_masses[eecol.e0_v0], m_masses[eecol.e0_v1], m_masses[eecol.e1_v0],
    //m_masses[eecol.e1_v1], eecol.s, eecol.t, eecol.GetRelativeVelocity(), eecol.n);

    exertEdgeImpulse(-I, m_masses[eecol.e0_v0], m_masses[eecol.e0_v1], eecol.s, eecol.e0_v0, eecol.e0_v1, m_vnphalf);
    exertEdgeImpulse(I, m_masses[eecol.e1_v0], m_masses[eecol.e1_v1], eecol.t, eecol.e1_v0, eecol.e1_v1, m_vnphalf);

    assert(eecol.GetRelativeVelocity() >= 0);
}

void BridsonStepper::exertInelasticImpluse(VertexFaceCTCollision& vfcol)
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
    assert(vfcol.GetRelativeVelocity() < 0.0);

    // Add some extra "kick" to relative velocity to account for FPA errors
    vfcol.ApplyRelativeVelocityKick();
    Vec3d I = vfcol.computeInelasticImpulse();
    // computeVertexFaceInelasticImpulse(m_masses[vfcol.v0], m_masses[vfcol.f0], m_masses[vfcol.f1], m_masses[vfcol.f2],
    // vfcol.u, vfcol.v, vfcol.w, eecol.GetRelativeVelocity(), vfcol.n);

    exertFaceImpulse(-I, m_masses[vfcol.f0], m_masses[vfcol.f1], m_masses[vfcol.f2], vfcol.u, vfcol.v, vfcol.w, vfcol.f0,
            vfcol.f1, vfcol.f2, m_vnphalf);
    exertVertexImpulse(I, m_masses[vfcol.v0], vfcol.v0, m_vnphalf);

    assert(vfcol.GetRelativeVelocity() >= 0);
}

void BridsonStepper::executeIterativeInelasticImpulseResponse()
{
    // Detect possible continuous time collisions
    std::list<CTCollision*> collisions;
    m_collision_detector->getContinuousTimeCollisions(collisions);

    // Iterativly apply inelastic impulses
    int itr;
    for (itr = 0; !collisions.empty() && itr < m_num_inlstc_itrns; ++itr)
    {
        // TODO: Add debug checks for repeat collisions.

        // Just sort the collision times to maintain some rough sense of causality
        collisions.sort(CompareTimes);
        while (!collisions.empty())
        {
            CTCollision* collision = collisions.front();
            // collision->Print(std::cerr);
            collisions.pop_front();
            exertCompliantInelasticImpulse(collision);
            delete collision;
            m_collision_detector->updateContinuousTimeCollisions();
        }
        m_collision_detector->getContinuousTimeCollisions(collisions);
    }
    // Just in case we haven't emptied the collisions but exited when itr == m_num_inlstc_itrns
    for (std::list<CTCollision*>::iterator i = collisions.begin(); i != collisions.end(); i++)
        delete *i;

    if (itr > 1)
        std::cerr << "\033[33mIterated collision response " << itr << " times\033[0m" << std::endl;

    if (itr == m_num_inlstc_itrns)
        std::cerr << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Exceeded maximum " << "number of inelastic iterations "
                << m_num_inlstc_itrns << ". Time of warning " << m_t << "." << std::endl;

#ifdef TIMING_ON
    if( itr >= 2 ) IntStatTracker::getIntTracker("STEPS_WITH_MULTIPLE_IMPULSE_ITERATIONS") += 1;
#endif
}

int BridsonStepper::getContainingRod(int vert_idx) const
{
    assert(vert_idx >= 0);
    assert(vert_idx < getNumVerts());
    int rod = -1;

    int accm = 0;
    for (int i = 0; i < (int) m_rods.size(); ++i)
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

bool BridsonStepper::execute()
{
    std::cerr << "Executing time step " << m_t << std::endl;
    Timer::getTimer("BridsonStepper::execute").start();
    bool do_adaptive = false;
    bool result;

    if (do_adaptive)
        result = adaptiveExecute(m_dt);
    else
        result = nonAdaptiveExecute(m_dt);

    Timer::getTimer("BridsonStepper::execute").stop();
    // Timer::report();

    return result;
}

double BridsonStepper::computeTotalForceNorm()
{
    double totalforcenormsqr = 0.0;
    for (int i = 0; i < (int) m_rods.size(); ++i)
    {
        VecXd force(m_rods[i]->ndof());
        force.setZero();
        m_rods[i]->computeForces(force);
        totalforcenormsqr += force.squaredNorm();
    }
    return sqrt(totalforcenormsqr);
}

void BridsonStepper::setDt(double dt)
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

bool BridsonStepper::nonAdaptiveExecute(double dt)
{
    setDt(dt);
    m_t += dt;
    for (int i = 0; i < m_scripting_controllers.size(); ++i)
        m_scripting_controllers[i]->setTime(m_t);
    return step(false);
}

bool BridsonStepper::adaptiveExecute(double dt)
{
    // Backup all rods
    std::vector<MinimalRodStateBackup> rodbackups(m_rods.size());
    for (int i = 0; i < (int) m_rods.size(); ++i)
    {
        rodbackups[i].resize(*m_rods[i]);
        rodbackups[i].backupRod(*m_rods[i]);
    }
    // Backup all objects
    std::vector<MinimalTriangleMeshBackup> objbackups(m_triangle_meshes.size());
    for (int i = 0; i < (int) m_triangle_meshes.size(); ++i)
    {
        objbackups[i].resize(*m_triangle_meshes[i]);
        objbackups[i].backupMesh(*m_triangle_meshes[i]);
    }
    // Backup the current simulation time
    double time = m_t;

    // Set the desired timestep
    std::cout << "Setting dt to: " << dt << std::endl;
    setDt(dt);
    // Advance the current time
    m_t += dt;
    for (int i = 0; i < (int) m_scripting_controllers.size(); ++i)
        m_scripting_controllers[i]->setTime(m_t);

    // Attempt a full time step
    if (step(true))
    {
        // Success!
        return true;
    }

    std::cout << "Adaptive stepping in Bridson stepper" << std::endl;

    // Restore all rods
    for (int i = 0; i < (int) m_rods.size(); ++i)
    {
        rodbackups[i].restoreRod(*m_rods[i]);
        rodbackups[i].clear();
    }
    // Restore all objects
    for (int i = 0; i < (int) m_triangle_meshes.size(); ++i)
    {
        objbackups[i].restoreMesh(*m_triangle_meshes[i]);
        objbackups[i].clear();
    }
    // Restore the time
    m_t = time;
    for (int i = 0; i < (int) m_scripting_controllers.size(); ++i)
        m_scripting_controllers[i]->setTime(m_t);

    // Otherwise attempt two steps of half length
    //setDt(0.5*dt);
    bool first_success = adaptiveExecute(0.5 * dt);
    if (!first_success)
    {
        setDt(dt);
        return false;
    }
    bool second_success = adaptiveExecute(0.5 * dt);
    if (!second_success)
    {
        setDt(dt);
        return false;
    }

    std::cout << "Finished two adaptive steps" << std::endl;
    std::cout << "Restoring dt to: " << dt << std::endl;
    setDt(dt);

    return first_success && second_success;
}

bool BridsonStepper::step(bool check_explosion)
{
    // BEGIN TEMP
    //for( size_t i = 0; i < m_rods.size(); ++i ) if( m_t >= 3.79 ) std::cout << i << ": " << computeMaxEdgeAngle( *m_rods[i] ) << std::endl;
    // END TEMP

    assert(m_edges.size() == m_edge_radii.size());
    assert((int) m_masses.size() == m_xn.size() / 3);
    assert(m_xn.size() == m_xnp1.size());
    assert(m_xn.size() == m_vnphalf.size());
    if (m_rod_labels.size() != 0)
        assert(m_rod_labels.size() == m_rods.size());
    // Sanity check to ensure rods are not "internally colliding" because radius is bigger than edge length
#ifdef DEBUG
    for( int i = 0; i < (int) m_rods.size(); ++i ) ensureNoCollisionsByDefault( *m_rods[i] );
#endif
    // Sanity check to ensure different parts of sim have same time/timetep
#ifdef DEBUG
    for( int i = 0; i < (int) m_scripting_controllers.size(); ++i ) assert( m_scripting_controllers[i]->getTime() == m_t );
    for( int i = 0; i < (int) m_scripting_controllers.size(); ++i ) assert( m_scripting_controllers[i]->getDt() == m_dt );
    for( int i = 0; i < (int) m_steppers.size(); ++i ) assert( m_steppers[i]->getTimeStep() == m_dt );
    for( int i = 0; i < (int) m_rods.size(); ++i ) assert( m_rods[i]->getTimeStep() == m_dt );
#endif

    // Compute the total start of timestep force
    double startforce = std::numeric_limits<double>::signaling_NaN();
    if (check_explosion)
        startforce = computeTotalForceNorm();

    // Save the pre-timestep positions
    extractPositions(m_rods, m_base_indices, m_xn);

    //std::cout << "m_xn: " << m_xn << std::endl;

    // Track whether or not this solve succeeds
    bool dependable_solve = true;

    // Step rods forward ignoring collisions
    START_TIMER("BridsonStepperDynamics");

    // Step scripted objects forward, set boundary conditions
    for (int i = 0; i < (int) m_scripting_controllers.size(); ++i)
        m_scripting_controllers[i]->execute();

    // Jungseock's implicit penalty
    if (m_implicit_pnlty_enbld)
    {
        std::cout << "Disabled this temporarily during testing" << std::endl;
        exit(1);
        // Clear exisiting penalties
        for (int i = 0; i < (int) m_implicit_pnlty_forces.size(); i++)
        {
            m_implicit_pnlty_forces[i]->clearPenaltyForces();
        }
        executeImplicitPenaltyResponse();
    }

    // Launch num_threads threads which will execute all elements of m_steppers.
    MultithreadedStepper<std::vector<RodTimeStepper*> > multithreaded_stepper(m_steppers, m_num_threads);
    dependable_solve = multithreaded_stepper.Execute();

    STOP_TIMER("BridsonStepperDynamics");

    // Post time step position
    extractPositions(m_rods, m_base_indices, m_xnp1);

    // Average velocity over the timestep just completed
    m_vnphalf = (m_xnp1 - m_xn) / m_dt;

    //if( m_pnlty_enbld ) executePenaltyResponse();
    START_TIMER("BridsonStepperResponse");
    if (m_itrv_inlstc_enbld && m_num_inlstc_itrns > 0)
        executeIterativeInelasticImpulseResponse();
    STOP_TIMER("BridsonStepperResponse");

    // Compute final positions from corrected velocities
    m_xnp1 = m_xn + m_dt * m_vnphalf;

    // Ensure boundary conditions respected by corrected positions
#ifdef DEBUG
    // For each rod
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
        RodBoundaryCondition* boundary = m_rods[i]->getBoundaryCondition();
        int rodbase = m_base_indices[i];

        // For each vertex of the current rod
        for( int j = 0; j < m_rods[i]->nv(); ++j )
        {
            // If that vertex has a prescribed position
            if( boundary->isVertexScripted(j) )
            {
                Vec3d desiredposition = boundary->getDesiredVertexPosition(j);
                Vec3d actualvalue = m_xnp1.segment<3>(rodbase+3*j);
                assert( approxEq(desiredposition, actualvalue, 1.0e-6) );
            }
        }
    }
#endif

    // Copy new positions and velocities back to rods
    restorePositions(m_rods, m_xnp1);
    restoreVelocities(m_rods, m_vnphalf);

    // Update frames and such in the rod (Is this correct? Will this do some extra stuff?)
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < (int) m_rods.size(); ++i)
        m_rods[i]->updateProperties();

    // Sanity check to ensure rod's internal state is consistent
#ifdef DEBUG
    for( int i = 0; i < (int) m_rods.size(); ++i ) m_rods[i]->verifyProperties();
#endif

    if (check_explosion)
    {
        // Compute the total end of timestep force
        double endforce = computeTotalForceNorm();

        assert(startforce >= 0.0);
        assert(endforce >= 0.0);
        if (endforce > startforce)
        {
            std::cout << endforce << " " << startforce << "    -    " << (endforce - startforce) / startforce << std::endl;
            // 2.0 and above still get instabilities
            double perrcent_change = 1.5;
            double big_force = 10.0;

            if (!approxEq(startforce, 0.0, 1.0e-6))
            {
                if (endforce - startforce >= perrcent_change * startforce)
                    dependable_solve = false;
            }
            else
            {
                if (endforce > big_force)
                    dependable_solve = false;
            }
        }
    }

    // BEGIN TEMP
    //  if( computeMaxEdgeAngle( *m_rods[23] ) > 1.0 ) std::cout << "Explosion detected after response" << std::endl;
    // END TEMP


    //  #ifdef TIMING_ON
    //    for( int i = 0; i < (int) m_base_indices.size(); ++i )
    //    {
    //      VecXd totalforce(m_rods[i]->ndof());
    //      totalforce.setZero();
    //      m_rods[i]->computeForces(totalforce);
    //      collision_induced_force_change[i].second = totalforce.norm();
    //    }
    //
    //    for( int i = 0; i < (int) m_base_indices.size(); ++i )
    //    {
    //      double fchng = fabs(collision_induced_force_change[i].first-collision_induced_force_change[i].second);
    //      ObjPropHandle<std::pair<double,int> > ophndl;
    //      m_rods[i]->add_property(ophndl,"collision_induced_force_change");
    //      m_rods[i]->property(ophndl).first = fchng;
    //      m_rods[i]->property(ophndl).second = 0;
    //    }
    //  #endif

    //    Timer::getTimer("BridsonStepperDynamics").report();
    //    Timer::getTimer("BridsonStepperResponse").report();
    //    Timer::getTimer("Collision detector").report();

    return dependable_solve;
}

void BridsonStepper::enableReseponse()
{
    m_respns_enbld = true;
}

void BridsonStepper::disableResponse()
{
    m_respns_enbld = false;
}

void BridsonStepper::enableIterativeInelasticImpulses()
{
    m_itrv_inlstc_enbld = true;
}

void BridsonStepper::disableIterativeInelasticImpulses()
{
    m_itrv_inlstc_enbld = false;
}

void BridsonStepper::setNumInelasticIterations(const int& num_itr)
{
    assert(num_itr >= 0);
    m_num_inlstc_itrns = num_itr;
}

// Ensures each rod edge has circular cross section. 
void BridsonStepper::ensureCircularCrossSection(const ElasticRod& rod) const
{
    // Ensure circular cross section
    for (int i = 0; i < (int) rod.ne(); ++i)
    {
        if (rod.getRadiusScale() * rod.radiusA(i) != rod.getRadiusScale() * rod.radiusB(i))
        {
            std::cerr
                    << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m Contact currently not supported for non-circular cross sections. Assuming circular cross section."
                    << m_num_inlstc_itrns << std::endl;
        }
    }
}

// Ensures each internal rod edge has length less than sum of neighbors' radii. 
void BridsonStepper::ensureNoCollisionsByDefault(const ElasticRod& rod) const
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
                    << m_num_inlstc_itrns << std::endl;
        }
    }
}

int BridsonStepper::getNumDof() const
{
    assert(m_num_dof >= 0);
    return m_num_dof;
}

int BridsonStepper::getNumVerts() const
{
    assert(m_num_dof % 3 == 0);
    return m_num_dof / 3;
}

// TODO: pass this triangle vector too
void BridsonStepper::extractPositions(const std::vector<ElasticRod*>& rods, const std::vector<int>& base_indices,
        VecXd& positions)
{
    assert(rods.size() == base_indices.size());
    assert(getNumDof() == positions.size());

    if (getNumDof() == 0)
        return;

#ifdef DEBUG
    positions.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (int i = 0; i < (int) rods.size(); ++i)
        for (int j = 0; j < rods[i]->nv(); ++j)
        {
            assert(base_indices[i] + 3 * j + 2 < positions.size());
            positions.segment<3> (base_indices[i] + 3 * j) = rods[i]->getVertex(j);
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
    for( int i = 0; i < (int) m_rods.size(); ++i )
    {
        RodBoundaryCondition* boundary = m_rods[i]->getBoundaryCondition();
        int rodbase = m_base_indices[i];

        // For each vertex of the current rod, if that vertex has a prescribed position
        for( int j = 0; j < m_rods[i]->nv(); ++j ) if( boundary->isVertexScripted(j) )
        {
            Vec3d desiredposition = boundary->getDesiredVertexPosition(j);
            Vec3d actualvalue = positions.segment<3>(rodbase+3*j);
            assert( approxEq(desiredposition, actualvalue, 1.0e-6) );
        }
    }
#endif
}

void BridsonStepper::extractVelocities(const std::vector<ElasticRod*>& rods, const std::vector<int>& base_indices,
        VecXd& velocities)
{
    assert(rods.size() == base_indices.size());
    assert(getNumDof() == velocities.size());

    if (getNumDof() == 0)
        return;

#ifdef DEBUG
    velocities.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    for (int i = 0; i < (int) rods.size(); ++i)
        for (int j = 0; j < rods[i]->nv(); ++j)
        {
            assert(base_indices[i] + 3 * j + 2 < velocities.size());
            velocities.segment<3> (base_indices[i] + 3 * j) = rods[i]->getVelocity(j);
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

void BridsonStepper::restorePositions(std::vector<ElasticRod*>& rods, const VecXd& positions)
{
    assert(rods.size() == m_base_indices.size());

    for (int i = 0; i < (int) m_base_indices.size(); ++i)
        for (int j = 0; j < rods[i]->nv(); ++j)
            if (!m_rods[i]->getBoundaryCondition()->isVertexScripted(j))
                rods[i]->setVertex(j, positions.segment<3> (m_base_indices[i] + 3 * j));
}

void BridsonStepper::restoreVelocities(std::vector<ElasticRod*>& rods, const VecXd& velocities)
{
    assert(rods.size() == m_base_indices.size());

    for (int i = 0; i < (int) m_base_indices.size(); ++i)
        for (int j = 0; j < rods[i]->nv(); ++j)
            if (!m_rods[i]->getBoundaryCondition()->isVertexScripted(j))
                rods[i]->setVelocity(j, velocities.segment<3> (m_base_indices[i] + 3 * j));
}

bool BridsonStepper::isRodVertex(int vert) const
{
    assert(vert >= 0);
    assert(vert < getNumVerts());

    // Is a vertex if index is less than start of object vertices in global array
    return vert < m_obj_start;
}

bool BridsonStepper::vertexAndFaceShareVertex(const int& v, const int& f0, const int& f1, const int& f2) const
{
    return v == f0 || v == f1 || v == f2;
}

bool BridsonStepper::vertexAndFaceShareVertex(const int& vertex, const int& face) const
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

bool BridsonStepper::isProperCollisionTime(double time)
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

/*
 Vec3d BridsonStepper::computeEdgeEdgeInelasticImpulse(const double& ma0, const double& ma1, const double& mb0,
 const double& mb1, const double& s, const double& t, const double& relvel, const Vec3d& n)
 {
 assert(ma0 > 0.0);
 assert(ma1 > 0.0);
 assert(mb0 > 0.0);
 assert(mb1 > 0.0);
 assert(s >= 0.0);
 assert(s <= 1.0);
 assert(t >= 0.0);
 assert(t <= 1.0);
 assert(relvel <= 0.0);

 // Assumes negative relative velocity
 Vec3d numerator = -relvel * n;
 double denominator = (1 - s) * (1 - s) / ma0 + s * s / ma1 + (1 - t) * (1 - t) / mb0 + t * t / mb1;
 assert(denominator != 0.0);
 return numerator / denominator;
 }
 */
/*
 Vec3d BridsonStepper::computeVertexFaceInelasticImpulse(const double& mvrt, const double& mfc0, const double& mfc1,
 const double& mfc2, const double& u, const double& v, const double& w, const double& relvel, const Vec3d& n)
 {
 assert(mvrt >= 0.0);
 assert(mfc0 >= 0.0);
 assert(mfc1 >= 0.0);
 assert(mfc2 >= 0.0);
 assert(approxEq(u + v + w, 1.0));
 //assert( u >= 0.0 ); assert( u <= 1.0 ); assert( v >= 0.0 ); assert( v <= 1.0 ); assert( w >= 0.0 ); assert( w <= 1.0 );
 assert(relvel <= 0.0);

 Vec3d numerator = -relvel * n;
 double denominator = 1 / mvrt + u * u / mfc0 + v * v / mfc1 + w * w / mfc2;
 assert(denominator != 0.0);
 return numerator / denominator;
 }
 */

void BridsonStepper::exertVertexImpulse(const Vec3d& I, const double& m, const int& idx, VecXd& v)
{
    assert(m > 0.0);
    assert(idx >= 0);
    assert(idx < getNumVerts());
    assert(v.size() == getNumDof());

    v.segment<3> (3 * idx) += I / m;
}

void BridsonStepper::exertEdgeImpulse(const Vec3d& I, const double& m0, const double& m1, const double& alpha, const int& idx0,
        const int& idx1, VecXd& v)
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

void BridsonStepper::exertFaceImpulse(const Vec3d& I, const double& m0, const double& m1, const double& m2, const double& u,
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

void BridsonStepper::computeCompliantLHS(MatrixBase* lhs, int rodidx)
{
    assert(lhs != NULL);
    assert(rodidx >= 0);
    assert(rodidx < (int) m_rods.size());

    // lhs = -h^2*dF/dx
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
    //    std::cout << "Force: " << (*fIt)->getName() << std::endl;
    // END TEMP
}

void BridsonStepper::exertCompliantInelasticImpulse(const CTCollision* cllsn)
{
    const EdgeEdgeCTCollision* eecol = dynamic_cast<const EdgeEdgeCTCollision*> (cllsn);
    const VertexFaceCTCollision* vfcol = dynamic_cast<const VertexFaceCTCollision*> (cllsn);
    if (eecol)
    {
        exertCompliantInelasticEdgeEdgeImpulse(*eecol);
        //exertInelasticImpluse(cllsn.getEdgeEdge());
    }
    else if (vfcol)
    {
        exertCompliantInelasticVertexFaceImpulse(*vfcol);
        //exertCompliantInelasticVertexFaceImpulse(cllsn.getVertexFace());
        //exertInelasticImpluse(cllsn.getVertexFace());
    }
}

void BridsonStepper::exertCompliantInelasticVertexFaceImpulse(const VertexFaceCTCollision& vfcol)
{
    //   std::cerr << "Vertex-face compliant inelastic impulse" << std::endl;
    //   std::cerr << vfcol << std::endl;

    // For now, assume vertex is free and entire face is fixed
    assert(!m_geodata.isVertexFixed(vfcol.v0));
    assert(YATriangle(vfcol.f0, vfcol.f1, vfcol.f2).IsFixed(m_geodata));

    // Determine which rod the free vertex belongs to
    int rodidx = getContainingRod(vfcol.v0);

    // Ensure the free vertex belongs to a rod
    assert(rodidx >= 0);
    assert(rodidx < (int) m_rods.size());

    // Determine which vertex of the rod the free vertex is
    assert(m_base_indices[rodidx] % 3 == 0);
    int v0 = vfcol.v0 - m_base_indices[rodidx] / 3;
    assert(v0 >= 0);
    assert(v0 < m_rods[rodidx]->nv());

    // Compute the relative velocity of the collision
    //    double magrelvel = vfcol.computeRelativeVelocity(m_geodata);
    assert(vfcol.GetRelativeVelocity() < 0.0);

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

    int rodbase = m_base_indices[rodidx];

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

    //std::cout << "nc: " << numconstraints << std::endl;
    //for( size_t i = 0; i < n.size(); ++i ) std::cout << n[i] << std::endl;

    // Determine the desired values for each constraint
    VecXd desired_values(numconstraints);
#ifdef DEBUG
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

    //std::cout << "Constraint values: " << desired_values << std::endl;

    // Currently, all fixed vertex constraints normalized
#ifdef DEBUG
    for( int i = 0; i < numconstraints; ++i ) assert( approxEq(normal[i].norm(),1.0,1.0e-9) );
#endif

    std::vector<VecXd> ntilde;
    for (int i = 0; i < numconstraints; ++i)
    {
        ntilde.push_back(VecXd(ndof));
        ntilde.back().setZero();
        int status = solver->solve(ntilde.back(), normal[i]);
        if (status < 0)
            std::cerr
                    << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. exertCompliantInelasticVertexFaceImpulse. Time: "
                    << m_t << std::endl;
    }

    // Ensure the edge degrees of freedom experience no impulse
#ifdef DEBUG
    for( ElasticRod::edge_iter eit = m_rods[rodidx]->edges_begin(); eit != m_rods[rodidx]->edges_end(); ++eit )
    {
        for( int i = 0; i < numconstraints; ++i ) assert( ntilde[i](m_rods[rodidx]->edgeIdx(*eit)) == 0.0 );
    }
#endif

    // Vectors restriced to just verted DoFs. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn;
    for (int i = 0; i < numconstraints; ++i)
        posnn.push_back(VecXd(nvdof));
#ifdef DEBUG
    for( int i = 0; i < numconstraints; ++i ) posnn[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    std::vector<VecXd> posnntilde;
    for (int i = 0; i < numconstraints; ++i)
        posnntilde.push_back(VecXd(nvdof));
#ifdef DEBUG
    for( int i = 0; i < numconstraints; ++i ) posnntilde[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

#ifdef DEBUG
    for( int i = 0; i < numconstraints; ++i ) assert( (normal[i].cwise()==normal[i]).all() );
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
#ifdef DEBUG
    for( int i = 0; i < numconstraints; ++i ) assert( (posnn[i].cwise()==posnn[i]).all() );
#endif

#ifdef DEBUG
    for( int i = 0; i < numconstraints; ++i ) assert( (ntilde[i].cwise()==ntilde[i]).all() );
#endif
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
#ifdef DEBUG
    for( int i = 0; i < numconstraints; ++i ) assert( (posnntilde[i].cwise()==posnntilde[i]).all() );
#endif

    // Compute the Lagrange multipliers
    // TODO: This matrix is symmetric, exploit that fact to avoid computations
    MatXd lglhs(numconstraints, numconstraints);
#ifdef DEBUG
    lglhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    for (int i = 0; i < numconstraints; ++i)
        for (int j = 0; j < numconstraints; ++j)
            lglhs(i, j) = posnn[i].dot(posnntilde[j]);
    assert(approxSymmetric(lglhs, 1.0e-6));

    Eigen::VectorXd lgrhs(numconstraints);
#ifdef DEBUG
    lgrhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    lgrhs(0) = vfcol.GetRelativeVelocity(); //posnN0.dot(m_vnphalf.segment(m_base_indices[rodidx],posnN0.size()));
    assert(lgrhs(0) < 0.0);
    for (int i = 1; i < numconstraints; ++i)
        lgrhs(i) = posnn[i].dot(m_vnphalf.segment(rodbase, nvdof));
    lgrhs *= -1.0;

    //std::cout << lgrhs << std::endl;

    lgrhs += desired_values;

    //std::cout << lgrhs << std::endl;
    //std::cout << desired_values << std::endl;

    Eigen::VectorXd alpha(numconstraints); // = lglhs.inverse()*lgrhs;
    alpha = lglhs.lu().solve(lgrhs);

    assert(alpha(0) >= 0.0);

    for (int i = 0; i < numconstraints; ++i)
        m_vnphalf.segment(rodbase, nvdof) += alpha(i) * posnntilde[i];

    // Ensure that the impulse eliminated the realtive velocity
#ifdef DEBUG
    Vec3d postrelvel = computeRelativeVelocity(m_vnphalf,vfcol.v0,vfcol.f0,vfcol.f1,vfcol.f2,vfcol.u,vfcol.v,vfcol.w);
    double postmagrelvel = postrelvel.dot(vfcol.normal);
    // Ensure the inelastic impulse decreased the realtive velocity
    assert( fabs(postmagrelvel) <= fabs(magrelvel) );
    // Should add some 'extra kick' to ensure collision gets killed, but for now just be content with small velocity
    assert( fabs(postmagrelvel) < 1.0e-9 );
    //std::cout << magrelvel << "   " << postmagrelvel << std::endl;
    //assert( postmagrelvel < 0.0 );
#endif

    // Ensure the 'scripted vertices' achieved the desired velocity
#ifdef DEBUG
    curdof = 1;
    for( size_t i = 0; i < scriptedverts.size(); ++i )
    {
        assert( scriptedverts[i] >= 0 );
        assert( (int) scriptedverts[i] < m_rods[rodidx]->nv() );
        //std::cout << m_vnphalf(rodbase+3*scriptedverts[i]+0) << "   " << desired_values(curdof) << std::endl;
        assert( approxEq(m_vnphalf(rodbase+3*scriptedverts[i]+0),desired_values(curdof),1.0e-6) );
        ++curdof;
        //std::cout << m_vnphalf(rodbase+3*scriptedverts[i]+1) << "   " << desired_values(curdof) << std::endl;
        assert( approxEq(m_vnphalf(rodbase+3*scriptedverts[i]+1),desired_values(curdof),1.0e-6) );
        ++curdof;
        //std::cout << m_vnphalf(rodbase+3*scriptedverts[i]+2) << "   " << desired_values(curdof) << std::endl;
        assert( approxEq(m_vnphalf(rodbase+3*scriptedverts[i]+2),desired_values(curdof),1.0e-6) );
        ++curdof;
    }
#endif
}

void BridsonStepper::exertCompliantInelasticEdgeEdgeImpulse(const EdgeEdgeCTCollision& eecol)
{
    //   std::cerr << "Edge-edge compliant inelastic impulse" << std::endl;
    //   std::cerr << eecol << std::endl;

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

void BridsonStepper::exertCompliantInelasticEdgeEdgeImpulseBothFree(const EdgeEdgeCTCollision& eecol)
{
    // Ensure both edges have two free vertices
    //assert( isEntireEdgeFree(eecol.e0_v0,eecol.e0_v1) );
    //assert( isEntireEdgeFree(eecol.e1_v0,eecol.e1_v1) );

    assert(!YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(m_geodata) && !YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata));

    // Determine which rod each edge belongs to
    int rod0 = getContainingRod(eecol.e0_v0);
    assert(rod0 == getContainingRod(eecol.e0_v1));
    int rod1 = getContainingRod(eecol.e1_v0);
    assert(rod1 == getContainingRod(eecol.e1_v1));

    // Don't do self-collisions, for now
    assert(rod0 != rod1);

    // Compute the relative velocity, which must be negative
    //    Vec3d relvel = m_geodata.computeRelativeVelocity(eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t);
    //    double magrelvel = relvel.dot(eecol.n);
    //   double magrelvel = eecol.computeRelativeVelocity(m_geodata);
    assert(eecol.GetRelativeVelocity() < 0.0);

    // Determine the total number of degrees of freedom in each rod
    int rod0ndof = m_rods[rod0]->ndof();
    int rod1ndof = m_rods[rod1]->ndof();

    // Determine the number of vertices in each rod
    int rod0nv = m_rods[rod0]->nv();
    int rod1nv = m_rods[rod1]->nv();

    // Determine the number of vertices in each rod
    int rod0nvdof = 3 * rod0nv;
    int rod1nvdof = 3 * rod1nv;

    //std::cout << "Ndof 0: " << rod0ndof << std::endl;
    //std::cout << "Ndof 1: " << rod1ndof << std::endl;
    //std::cout << "NV 0: " << rod0nv << std::endl;
    //std::cout << "NV 1: " << rod1nv << std::endl;
    //std::cout << "NVdof 0: " << rod0nvdof << std::endl;
    //std::cout << "NVdof 1: " << rod1nvdof << std::endl;

    // Determine where in the 'global' vertex dof pool each rod begins
    int rod0base = m_base_indices[rod0];
    assert(rod0base % 3 == 0);
    int rod1base = m_base_indices[rod1];
    assert(rod1base % 3 == 0);

    //std::cout << "rod0base: " << rod0base << std::endl;
    //std::cout << "rod1base: " << rod1base << std::endl;

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
#ifdef DEBUG
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
#ifdef DEBUG
    Vec3d testn0 = n0[0].segment<3>(base_i0)+n0[0].segment<3>(base_i1);
    Vec3d actln0 = -eecol.n;
    assert( approxEq(testn0, actln0, 1.0e-6) );
#endif
    // Currently, all fixed vertex constraints normalized
#ifdef DEBUG
    for( int i = 1; i < nc0; ++i ) assert( approxEq(n0[i].norm(),1.0,1.0e-9) );
#endif
    // Ensure the edge degrees of freedom experience no impulse
#ifdef DEBUG
    for( int i = 0; i < nc0; ++i )
    for( ElasticRod::edge_iter eit = m_rods[rod0]->edges_begin(); eit != m_rods[rod0]->edges_end(); ++eit )
    assert( n0[i](m_rods[rod0]->edgeIdx(*eit)) == 0.0 );
#endif

    //for( int i = 0; i < nc0; ++i ) std::cout << "n0:    " << n0[i] << std::endl;
    //for( int i = 0; i < nc0; ++i ) std::cout << "cval0: " << cval0[i] << std::endl;

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
#ifdef DEBUG
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
#ifdef DEBUG
    Vec3d testn1 = n1[0].segment<3>(base_j0)+n1[0].segment<3>(base_j1);
    Vec3d actln1 = eecol.n;
    assert( approxEq(testn1, actln1, 1.0e-6) );
#endif
    // Currently, all fixed vertex constraints normalized
#ifdef DEBUG
    for( int i = 1; i < nc1; ++i ) assert( approxEq(n1[i].norm(),1.0,1.0e-9) );
#endif
    // Ensure the edge degrees of freedom experience no impulse
#ifdef DEBUG
    for( int i = 0; i < nc1; ++i )
    for( ElasticRod::edge_iter eit = m_rods[rod1]->edges_begin(); eit != m_rods[rod1]->edges_end(); ++eit )
    assert( n1[i](m_rods[rod1]->edgeIdx(*eit)) == 0.0 );
#endif

    //for( int i = 0; i < nc1; ++i ) std::cout << "n1:    " << n1[i] << std::endl;
    //for( int i = 0; i < nc1; ++i ) std::cout << "cval1: " << cval1[i] << std::endl;

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
            std::cerr << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. " << std::endl;
    }
    // Ensure the edge degrees of freedom experience no impulse
#ifdef DEBUG
    for( int i = 0; i < nc0; ++i )
    for( ElasticRod::edge_iter eit = m_rods[rod0]->edges_begin(); eit != m_rods[rod0]->edges_end(); ++eit )
    assert( ntilde0[i](m_rods[rod0]->edgeIdx(*eit)) == 0.0 );
#endif

    //for( int i = 0; i < nc0; ++i ) std::cout << "ntilde0: " << ntilde0[i] << std::endl;

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
            std::cerr << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. " << std::endl;
    }
    // Ensure the edge degrees of freedom experience no impulse
#ifdef DEBUG
    for( int i = 0; i < nc1; ++i )
    for( ElasticRod::edge_iter eit = m_rods[rod1]->edges_begin(); eit != m_rods[rod1]->edges_end(); ++eit )
    assert( ntilde1[i](m_rods[rod1]->edgeIdx(*eit)) == 0.0 );
#endif

    //for( int i = 0; i < nc1; ++i ) std::cout << "ntilde1: " << ntilde1[i] << std::endl;

    // Vectors restriced to just vertex DoFs for rod 0. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn0;
    for (int i = 0; i < nc0; ++i)
        posnn0.push_back(VecXd(rod0nvdof));
#ifdef DEBUG
    for( int i = 0; i < nc0; ++i ) posnn0[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    std::vector<VecXd> posnntilde0;
    for (int i = 0; i < nc0; ++i)
        posnntilde0.push_back(VecXd(rod0nvdof));
#ifdef DEBUG
    for( int i = 0; i < nc0; ++i ) posnntilde0[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

#ifdef DEBUG
    for( int i = 0; i < nc0; ++i ) assert( (n0[i].cwise()==n0[i]).all() );
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
#ifdef DEBUG
    for( int i = 0; i < nc0; ++i ) assert( (posnn0[i].cwise()==posnn0[i]).all() );
#endif

    //for( int i = 0; i < nc0; ++i ) std::cout << "posnn0: " << posnn0[i] << std::endl;

#ifdef DEBUG
    for( int i = 0; i < nc0; ++i ) assert( (ntilde0[i].cwise()==ntilde0[i]).all() );
#endif
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
#ifdef DEBUG
    for( int i = 0; i < nc0; ++i ) assert( (posnntilde0[i].cwise()==posnntilde0[i]).all() );
#endif

    //for( int i = 0; i < nc0; ++i ) std::cout << "posnntilde0: " << posnntilde0[i] << std::endl;


    // Vectors restriced to just vertex DoFs for rod 1. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn1;
    for (int i = 0; i < nc1; ++i)
        posnn1.push_back(VecXd(rod1nvdof));
#ifdef DEBUG
    for( int i = 0; i < nc1; ++i ) posnn1[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    std::vector<VecXd> posnntilde1;
    for (int i = 0; i < nc1; ++i)
        posnntilde1.push_back(VecXd(rod1nvdof));
#ifdef DEBUG
    for( int i = 0; i < nc1; ++i ) posnntilde1[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

#ifdef DEBUG
    for( int i = 0; i < nc1; ++i ) assert( (n1[i].cwise()==n1[i]).all() );
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
#ifdef DEBUG
    for( int i = 0; i < nc1; ++i ) assert( (posnn1[i].cwise()==posnn1[i]).all() );
#endif

    //for( int i = 0; i < nc1; ++i ) std::cout << "posnn1: " << posnn1[i] << std::endl;

#ifdef DEBUG
    for( int i = 0; i < nc1; ++i ) assert( (ntilde1[i].cwise()==ntilde1[i]).all() );
#endif
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
#ifdef DEBUG
    for( int i = 0; i < nc1; ++i ) assert( (posnntilde1[i].cwise()==posnntilde1[i]).all() );
#endif

    //for( int i = 0; i < nc1; ++i ) std::cout << "posnntilde1: " << posnntilde1[i] << std::endl;


    // Compute the Lagrange multipliers for all constraints
    // TODO: This matrix is symmetric, exploit that fact to avoid computations
    int numalpha = nc0 + nc1 - 1;
    MatXd lglhs(numalpha, numalpha);
    //#ifdef DEBUG
    //  lglhs.setConstant(std::numeric_limits<double>::signaling_NaN());
    //#endif
    lglhs.setZero();

    //std::cout << "System size: " << numalpha << std::endl;

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
#ifdef DEBUG
    lgrhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    // Entry 0
    assert(posnn0[0].size() == m_vnphalf.segment(rod0base, rod0nvdof).size());
    assert(posnn1[0].size() == m_vnphalf.segment(rod1base, rod1nvdof).size());
    lgrhs(0) = posnn0[0].dot(m_vnphalf.segment(rod0base, rod0nvdof)) + posnn1[0].dot(m_vnphalf.segment(rod1base, rod1nvdof));
    assert(approxEq(lgrhs(0), eecol.GetRelativeVelocity(), 1.0e-6));
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
    //std::cout << lgrhs << std::endl;
    lgrhs *= -1.0;

    // For now, model an inelastic collision
    //lglhs0(0) += 0.0;
    // Entries 1...numalpha
    lgrhs.segment(1, nc0 - 1) += cval0.segment(1, nc0 - 1);
    // Entries numalpha...numalpha+numbeta
    lgrhs.segment(nc0, nc1 - 1) += cval1.segment(1, nc1 - 1);
    //std::cout << lgrhs << std::endl;

    Eigen::VectorXd alpha(numalpha);
    assert(lglhs.rows() == lglhs.cols());
    assert(lglhs.rows() == lgrhs.size());
    assert(lglhs.rows() == alpha.size());
    alpha = lglhs.lu().solve(lgrhs);

    assert(alpha(0) >= 0.0);

    if (alpha(0) < 0.0)
        std::cout << "WARNING NEGATIVE LAGRANGE MULTIPLIER alpha: " << alpha << std::endl;

    // BEGIN TEMP
    //if( rod0 == 23 || rod1 == 23 ) std::cout << "Collision lagrange multiplier: " << alpha(0) << std::endl;
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
#ifdef DEBUG
    Vec3d postrelvel = computeRelativeVelocity( m_vnphalf, eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t );
    double postmagrelvel = postrelvel.dot(eecol.n);
    // Ensure the inelastic impulse decreased the realtive velocity
    assert( fabs(postmagrelvel) <= fabs(magrelvel) );
    // Should add some 'extra kick' to ensure collision gets killed, but for now just be content with small velocity
    assert( fabs(postmagrelvel) < 1.0e-9 );
    //std::cout << magrelvel << "   " << postmagrelvel << std::endl;
    //assert( postmagrelvel < 0.0 );
#endif

    // Ensure the 'scripted vertices' achieved the desired velocity for rod 0
#ifdef DEBUG
    curdof = 1;
    for( size_t i = 0; i < scriptedverts0.size(); ++i )
    {
        assert( scriptedverts0[i] >= 0 );
        assert( (int) scriptedverts0[i] < m_rods[rod0]->nv() );
        //std::cout << m_vnphalf(rodbase+3*scriptedverts[i]+0) << "   " << desired_values(curdof) << std::endl;
        assert( approxEq(m_vnphalf(rod0base+3*scriptedverts0[i]+0),cval0(curdof),1.0e-6) );
        ++curdof;
        //std::cout << m_vnphalf(rodbase+3*scriptedverts[i]+1) << "   " << desired_values(curdof) << std::endl;
        assert( approxEq(m_vnphalf(rod0base+3*scriptedverts0[i]+1),cval0(curdof),1.0e-6) );
        ++curdof;
        //std::cout << m_vnphalf(rodbase+3*scriptedverts[i]+2) << "   " << desired_values(curdof) << std::endl;
        assert( approxEq(m_vnphalf(rod0base+3*scriptedverts0[i]+2),cval0(curdof),1.0e-6) );
        ++curdof;
    }
#endif

    // Ensure the 'scripted vertices' achieved the desired velocity for rod 1
#ifdef DEBUG
    curdof = 1;
    for( size_t i = 0; i < scriptedverts1.size(); ++i )
    {
        assert( scriptedverts1[i] >= 0 );
        assert( (int) scriptedverts1[i] < m_rods[rod1]->nv() );
        //std::cout << m_vnphalf(rodbase+3*scriptedverts[i]+0) << "   " << desired_values(curdof) << std::endl;
        assert( approxEq(m_vnphalf(rod1base+3*scriptedverts1[i]+0),cval1(curdof),1.0e-6) );
        ++curdof;
        //std::cout << m_vnphalf(rodbase+3*scriptedverts[i]+1) << "   " << desired_values(curdof) << std::endl;
        assert( approxEq(m_vnphalf(rod1base+3*scriptedverts1[i]+1),cval1(curdof),1.0e-6) );
        ++curdof;
        //std::cout << m_vnphalf(rodbase+3*scriptedverts[i]+2) << "   " << desired_values(curdof) << std::endl;
        assert( approxEq(m_vnphalf(rod1base+3*scriptedverts1[i]+2),cval1(curdof),1.0e-6) );
        ++curdof;
    }
#endif
}

void BridsonStepper::exertCompliantInelasticEdgeEdgeImpulseOneFixed(const EdgeEdgeCTCollision& eecol)
{
    // Must have one totally fixed and one totally free edge
    assert(
            (YAEdge(eecol.e0_v0, eecol.e0_v1).IsFree(m_geodata) && YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata))
                    || (YAEdge(eecol.e1_v0, eecol.e1_v1).IsFree(m_geodata) && YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(
                            m_geodata)));

    // Determine if either edge is fixed
    bool rod0fixed = YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(m_geodata);
    bool rod1fixed = YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata);
    assert(rod0fixed != rod1fixed);

    // Compute the relative velocity, which must be negative for a collision to have happened
    //    Vec3d relvel = m_geodata.computeRelativeVelocity(eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t);
    //    double magrelvel = relvel.dot(eecol.n);
    //   double magrelvel = eecol.computeRelativeVelocity(m_geodata);
    assert(eecol.GetRelativeVelocity() < 0.0);

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
    assert(rodidx < (int) m_rods.size());

    // Convert the vertices' global indices to rod indices
    assert(m_base_indices[rodidx] % 3 == 0);
    v0 = v0 - m_base_indices[rodidx] / 3;
    v1 = v1 - m_base_indices[rodidx] / 3;
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

    int rodbase = m_base_indices[rodidx];

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
            //std::cout << m_rods[rodidx]->vertIdx(scriptedverts[i],j) << " ";
            n.push_back(VecXd(ndof));
            n.back().setZero();
            // Zero the velocity of this particular dof of this particular vertex
            n.back()(m_rods[rodidx]->vertIdx(scriptedverts[i], j)) = 1.0;
        }
    }
    //std::cout << std::endl;

    int numconstraints = (int) (n.size());

    // Determine the desired values for each constraint
    VecXd desired_values(numconstraints);
#ifdef DEBUG
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

    //std::cout << desired_values << std::endl;

    //for( size_t i = 0; i < vertexConstraints.size(); ++i ) std::cout << vertexConstraints[i] << std::endl;

    // Ensure collision constraint adds up to actual normal
#ifdef DEBUG
    Vec3d testn = n[0].segment<3>(base0)+n[0].segment<3>(base1);
    if( !rod0fixed ) testn *= -1.0;
    assert( approxEq(testn, eecol.n, 1.0e-6) );
#endif

    // Currently, all fixed vertex constraints normalized
#ifdef DEBUG
    for( int i = 1; i < numconstraints; ++i ) assert( approxEq(n[i].norm(),1.0,1.0e-9) );
#endif

    std::vector<VecXd> ntilde;
    for (int i = 0; i < numconstraints; ++i)
    {
        ntilde.push_back(VecXd(ndof));
        ntilde.back().setZero();
        int status = solver->solve(ntilde.back(), n[i]);
        if (status < 0)
            std::cerr
                    << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. exertCompliantInelasticEdgeEdgeImpulseOneFixedThree. Time: "
                    << m_t << std::endl;
    }

    //for( size_t i = 0; i < n.size(); ++i ) std::cout << ntilde[i].transpose() << std::endl;


    // Ensure the edge degrees of freedom experience no impulse
#ifdef DEBUG
    for( ElasticRod::edge_iter eit = m_rods[rodidx]->edges_begin(); eit != m_rods[rodidx]->edges_end(); ++eit )
    {
        for( int i = 0; i < numconstraints; ++i ) assert( ntilde[i](m_rods[rodidx]->edgeIdx(*eit)) == 0.0 );
    }
#endif

    // Vectors restriced to just verted DoFs. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn;
    for (int i = 0; i < numconstraints; ++i)
        posnn.push_back(VecXd(nvdof));
#ifdef DEBUG
    for( int i = 0; i < numconstraints; ++i ) posnn[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

    std::vector<VecXd> posnntilde;
    for (int i = 0; i < numconstraints; ++i)
        posnntilde.push_back(VecXd(nvdof));
#ifdef DEBUG
    for( int i = 0; i < numconstraints; ++i ) posnntilde[i].setConstant(std::numeric_limits<double>::signaling_NaN());
#endif

#ifdef DEBUG
    for( int i = 0; i < numconstraints; ++i ) assert( (n[i].cwise()==n[i]).all() );
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
#ifdef DEBUG
    for( int i = 0; i < numconstraints; ++i ) assert( (posnn[i].cwise()==posnn[i]).all() );
#endif

#ifdef DEBUG
    for( int i = 0; i < numconstraints; ++i ) assert( (ntilde[i].cwise()==ntilde[i]).all() );
#endif
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
#ifdef DEBUG
    for( int i = 0; i < numconstraints; ++i ) assert( (posnntilde[i].cwise()==posnntilde[i]).all() );
#endif

    // Compute the Lagrange multipliers
    // TODO: This matrix is symmetric, exploit that fact to avoid computations
    MatXd lglhs(numconstraints, numconstraints);
#ifdef DEBUG
    lglhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    for (int i = 0; i < numconstraints; ++i)
        for (int j = 0; j < numconstraints; ++j)
            lglhs(i, j) = posnn[i].dot(posnntilde[j]);
    //  assert(approxSymmetric(lglhs, 1.0e-6)); // FIXME: why is this failing?

    Eigen::VectorXd lgrhs(numconstraints);
#ifdef DEBUG
    lgrhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    lgrhs(0) = eecol.GetRelativeVelocity(); // CRAZY SIGN PROBLEM HERE?
    assert(lgrhs(0) < 0.0);
    for (int i = 1; i < numconstraints; ++i)
        lgrhs(i) = posnn[i].dot(m_vnphalf.segment(rodbase, nvdof));
    lgrhs *= -1.0;

    //std::cout << lgrhs << std::endl;

    lgrhs += desired_values;

    //std::cout << lgrhs << std::endl;
    //std::cout << desired_values << std::endl;

    Eigen::VectorXd alpha(numconstraints); // = lglhs.inverse()*lgrhs;
    alpha = lglhs.lu().solve(lgrhs);

    // Contact constraint should 'push not pull'
    assert(alpha(0) >= 0.0);

    // TEMP DEBUG STUFF
    //  if( rodidx == 20 && m_t >= 0.85 )
    //  {
    //    std::cout << "---------- ROD 20 ROD 20 ROD 20 ROD 20 ROD 20 ----------" << std::endl;
    //    std::cout << "alpha: " << alpha << std::endl;
    //    std::cout << "time: " << m_t << std::endl;
    //
    //    if( !rod0fixed ) std::cout  << "Free: " << 0 << std::endl;
    //    if( !rod1fixed ) std::cout  << "Free: " << 1 << std::endl;
    //
    //    std::cout << "temprelvel: " << posnn[0].dot(m_vnphalf.segment(rodbase,nvdof)) << std::endl;
    //  }
    // END TEMP DEBUG STUFF

    for (int i = 0; i < numconstraints; ++i)
        m_vnphalf.segment(rodbase, nvdof) += alpha(i) * posnntilde[i];

    // TEMP DEBUG STUFF
    //  Vec3d postrelvel = computeRelativeVelocity( m_vnphalf, eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t );
    //  double postmagrelvel = postrelvel.dot(eecol.n);
    //  // Ensure the inelastic impulse decreased the realtive velocity
    //  if( fabs(postmagrelvel) > fabs(magrelvel) )
    //  {
    //    std::cout << "Rod: " << rodidx << std::endl;
    //    std::cout << "Time is: " << m_t << std::endl;
    //    std::cout << "Incoming relative velocity: " << magrelvel << "      Post-impulse relative velocity: " << postmagrelvel << std::endl;
    //
    //    exit(1);
    //  }
    // END TEMP DEBUG STUFF


    // Ensure that the impulse eliminated the realtive velocity
#ifdef DEBUG
    Vec3d postrelvel = computeRelativeVelocity( m_vnphalf, eecol.e0_v0, eecol.e0_v1, eecol.e1_v0, eecol.e1_v1, eecol.s, eecol.t );
    double postmagrelvel = postrelvel.dot(eecol.n);
    // Ensure the inelastic impulse decreased the realtive velocity
    if( fabs(postmagrelvel) > fabs(magrelvel) )
    {
        std::cout << "Rod: " << rodidx << std::endl;
        std::cout << "Time is: " << m_t << std::endl;
        std::cout << "Incoming relative velocity: " << magrelvel << "      Post-impulse relative velocity: " << postmagrelvel << std::endl;
    }
    assert( fabs(postmagrelvel) <= fabs(magrelvel) );
    // Should add some 'extra kick' to ensure collision gets killed, but for now just be content with small velocity
    assert( fabs(postmagrelvel) < 1.0e-9 );
    //assert( postmagrelvel < 0.0 );
#endif

    // Ensure the 'scripted vertices' achieved the desired velocity
#ifdef DEBUG
    curdof = 1;
    for( size_t i = 0; i < scriptedverts.size(); ++i )
    {
        assert( scriptedverts[i] >= 0 );
        assert( (int) scriptedverts[i] < m_rods[rodidx]->nv() );
        //std::cout << m_vnphalf(rodbase+3*scriptedverts[i]+0) << "   " << desired_values(curdof) << std::endl;
        assert( approxEq(m_vnphalf(rodbase+3*scriptedverts[i]+0),desired_values(curdof),1.0e-6) );
        ++curdof;
        //std::cout << m_vnphalf(rodbase+3*scriptedverts[i]+1) << "   " << desired_values(curdof) << std::endl;
        assert( approxEq(m_vnphalf(rodbase+3*scriptedverts[i]+1),desired_values(curdof),1.0e-6) );
        ++curdof;
        //std::cout << m_vnphalf(rodbase+3*scriptedverts[i]+2) << "   " << desired_values(curdof) << std::endl;
        assert( approxEq(m_vnphalf(rodbase+3*scriptedverts[i]+2),desired_values(curdof),1.0e-6) );
        ++curdof;
    }
#endif
}

void BridsonStepper::detectVertexFaceImplicitPenaltyCollisions(const VecXd& x,
        std::vector<VertexFaceProximityCollision>& pssbl_cllsns,
        std::vector<VertexFaceImplicitPenaltyCollision>& vetex_face_collisions) const
{
    assert((int) m_masses.size() == getNumVerts());
    assert(x.size() == getNumDof());
    assert(m_faces.size() == m_face_radii.size());

    //std::vector<VertexFaceProximityCollision> pssbl_cllsns;
    //generateAllVertexFaceProximityCollisionPairs( pssbl_cllsns );

    for (int i = 0; i < (int) pssbl_cllsns.size(); ++i)
    {
        assert(pssbl_cllsns[i].v0 < getNumVerts());
        assert(pssbl_cllsns[i].f0 < getNumVerts());
        assert(pssbl_cllsns[i].f1 < getNumVerts());
        assert(pssbl_cllsns[i].f2 < getNumVerts());
        assert(pssbl_cllsns[i].r0 >= 0.0);
        assert(pssbl_cllsns[i].r1 >= 0.0);

        int vrtidx = pssbl_cllsns[i].v0;
        int fcidx0 = pssbl_cllsns[i].f0;
        int fcidx1 = pssbl_cllsns[i].f1;
        int fcidx2 = pssbl_cllsns[i].f2;

        if (vertexAndFaceShareVertex(vrtidx, fcidx0, fcidx1, fcidx2))
            continue;

        // TODO: Add check for both having fixed vertices

        const Vec3d p1 = x.segment<3> (3 * vrtidx);
        const Vec3d t0 = x.segment<3> (3 * fcidx0);
        const Vec3d t1 = x.segment<3> (3 * fcidx1);
        const Vec3d t2 = x.segment<3> (3 * fcidx2);

        Vec3d cp = ClosestPtPointTriangle(p1, t0, t1, t2);
        double sqrdist = (p1 - cp).norm();

        if (sqrdist < (m_implicit_thickness + pssbl_cllsns[i].r0 + pssbl_cllsns[i].r1) * (m_implicit_thickness
                + pssbl_cllsns[i].r0 + pssbl_cllsns[i].r1))
        {
            VertexFaceImplicitPenaltyCollision real_cllsn;

            real_cllsn.v0 = vrtidx;
            real_cllsn.f0 = pssbl_cllsns[i].f0;
            real_cllsn.f1 = pssbl_cllsns[i].f1;
            real_cllsn.f2 = pssbl_cllsns[i].f2;

            real_cllsn.r0 = pssbl_cllsns[i].r0;
            real_cllsn.r1 = pssbl_cllsns[i].r1;
            real_cllsn.h = m_implicit_thickness;

            real_cllsn.k = m_vrt_fc_pnlty;

            //      Barycentric( t0, t1, t2, p1, pssbl_cllsns[i].u, pssbl_cllsns[i].v, pssbl_cllsns[i].w );

            // u,v,w CAN be outside of 0,1. We're really looking at minkowski sum of triangle with sphere here.
            //      assert( approxEq(pssbl_cllsns[i].u+pssbl_cllsns[i].v+pssbl_cllsns[i].w,1.0) );

            //      pssbl_cllsns[i].pen = pssbl_cllsns[i].r0+pssbl_cllsns[i].r1-sqrt(sqrdist);
            //      assert( pssbl_cllsns[i].pen > 0.0 );

            //      pssbl_cllsns[i].n = p1-cp;
            real_cllsn.n = (t1 - t0).cross(t2 - t0);
            assert(real_cllsn.n.norm() > 0.0);

            real_cllsn.n.normalize();
            assert(fabs(real_cllsn.n.norm() - 1.0) < 1.0e-6);

            real_cllsn.cp = cp;

            // TODO: Add some checks that if the closest point is inside the triangle, normal is normal to each edge

            vetex_face_collisions.push_back(real_cllsn);
        }
    }
}

void BridsonStepper::executeImplicitPenaltyResponse()
{
    // Detect possible "proximity" collisions based on pre-step positions
    std::vector<EdgeEdgeProximityCollision> edge_edge_collisions;
    std::vector<VertexFaceProximityCollision> vertex_face_collisions;
    m_collision_detector->getImplicitPenaltyCollisions(edge_edge_collisions, vertex_face_collisions);

    std::vector<VertexFaceImplicitPenaltyCollision> vertex_face_for_real;
    detectVertexFaceImplicitPenaltyCollisions(m_xn, vertex_face_collisions, vertex_face_for_real);

    // Add penalty forces
    //  exertImplicitPenalty( vertex_face_for_real, m_vnphalf );
    for (int i = 0; i < (int) vertex_face_for_real.size(); i++)
    {
        int rod_id = getContainingRod(vertex_face_for_real[i].v0);

        if (rod_id >= 0)
        {
            int v_id = vertex_face_for_real[i].v0 - m_base_indices[rod_id] / 3;

            m_implicit_pnlty_forces[rod_id]->addRodPenaltyForce(v_id, vertex_face_for_real[i]);
        }
    }
}

void BridsonStepper::enableImplicitPenaltyImpulses()
{
    m_implicit_pnlty_enbld = true;

    for (int i = 0; i < (int) m_rods.size(); i++)
    {
        RodPenaltyForce *pnlty = new RodPenaltyForce();
        m_implicit_pnlty_forces.push_back(pnlty);
        m_steppers[i]->addExternalForce(pnlty);
    }
}

void BridsonStepper::disableImplicitPenaltyImpulses()
{
    m_implicit_pnlty_enbld = false;

    for (int i = 0; i < (int) m_rods.size(); i++)
    {
        std::vector<RodExternalForce*>& forces = m_steppers[i]->getExternalForces();
        for (int j = 0; j < (int) forces.size(); j++)
        {
            if (typeid(*(forces[j])) == typeid(RodPenaltyForce))
                    {
                        forces.erase(forces.begin() + j);
                        break;
                    }
                }
            }

            m_implicit_pnlty_forces.clear();
        }

        void BridsonStepper::setImplicitPenaltyExtraThickness(const double& h)
        {
            m_implicit_thickness = h;
        }

        void BridsonStepper::setVertexFacePenalty(const double& k)
        {
            assert(k >= 0.0);
            m_vrt_fc_pnlty = k;
        }

        }

