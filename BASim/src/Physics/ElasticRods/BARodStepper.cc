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
#include "../../Collisions/CTCollision.hh"

#include <iostream>
#include <fstream>
#include "../../Eigen/Eigenvalues"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

namespace BASim
{

BARodStepper::BARodStepper(std::vector<ElasticRod*>& rods, std::vector<TriangleMesh*>& trimeshes,
        std::vector<ScriptingController*>& scripting_controllers, std::vector<RodTimeStepper*>& steppers, const double& dt,
        const double time, const int num_threads, const PerformanceTuningParameters perf_param,
        std::vector<LevelSet*>* levelSets) :
            m_num_dof(0),
            m_rods(rods),
            m_number_of_rods(rods.size()),
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
            m_level_sets(levelSets ? *levelSets : std::vector<LevelSet*>(trimeshes.size(), NULL)),
            m_geodata(m_xn, m_vnphalf, m_vertex_radii, m_masses, m_collision_immune, m_obj_start,
                    m_perf_param.m_implicit_thickness, m_perf_param.m_implicit_stiffness)
{
    g_log = new TextLog(std::cerr, MsgInfo::kDebug, true);
    InfoStream(g_log, "") << "Started logging BARodStepper\n";

#ifndef NDEBUG
    checkDataConsistency();
#endif

    setNumThreads(num_threads);

    prepareForExecution();

    allocateBackups();

    buildCollisionDetector();

    initializeSimulationList();

#ifndef NDEBUG
    checkInternalConsistency();
#endif

    CopiousStream(g_log, "") << "Finished BARodStepper construction\n";
}

BARodStepper::~BARodStepper()
{
    delete m_collision_detector;

    for (int i = 0; i < m_number_of_rods; i++)
    {
        delete m_startForces[i];
        delete m_preDynamicForces[i];
        delete m_preCollisionForces[i];
        delete m_endForces[i];
    }

    delete g_log;
}

/**
 * Preparation
 */
void BARodStepper::checkDataConsistency()
{
    assert(m_rods.size() == m_number_of_rods);
    assert(m_steppers.size() == m_number_of_rods);

    for (int i = 0; i < (int) m_steppers.size(); ++i)
        assert(m_steppers[i] != NULL);
    for (int i = 0; i < (int) m_number_of_rods; ++i)
    {
        assert(m_rods[i] != NULL);
        // Sanity checks for collision detection purposes
        ensureCircularCrossSection(*m_rods[i]);
        ensureNoCollisionsByDefault(*m_rods[i]);
    }
    for (int i = 0; i < (int) m_triangle_meshes.size(); ++i)
        assert(m_triangle_meshes[i] != NULL);

    // Do not check if any level sets are null as that may be valid as not all triangle meshes
    // may have a level set.
    assert(m_level_sets.size() == m_triangle_meshes.size());

    assert(m_dt > 0.0);
}

void BARodStepper::prepareForExecution()
{
    for (int i = 0; i < m_number_of_rods; ++i)
        m_rods[i]->globalRodIndex = i;

    CopiousStream(g_log, "") << "About to extract rod information\n";

    for (int i = 0; i < m_number_of_rods; ++i)
    {
        // Extract edges from the new rod
        for (int j = 0; j < m_rods[i]->nv() - 1; ++j)
        {
            m_edges.push_back(std::pair<int, int>(getNumVerts() + j, getNumVerts() + j + 1));
            assert(m_rods[i]->getRadiusScale() * m_rods[i]->radiusA(j) > 0.0);
            m_edge_radii.push_back(m_rods[i]->getRadiusScale() * m_rods[i]->radiusA(j));
        }

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
        CopiousStream(g_log, "") << "Finished extracting face stuff: " << i + 1 << " mesh" << (i > 0 ? "es\n" : "\n");

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

    for (std::vector<RodTimeStepper*>::iterator stepper = m_steppers.begin(); stepper != m_steppers.end(); ++stepper)
    {
        (*stepper)->setMaxIterations(m_perf_param.m_solver.m_max_iterations);
    }
}

void BARodStepper::allocateBackups()
{
    // Resize the internal storage
    m_xn.resize(getNumDof());
    m_xn.setConstant(std::numeric_limits<double>::signaling_NaN());
    m_xnp1.resize(getNumDof());
    m_xnp1.setConstant(std::numeric_limits<double>::signaling_NaN());
    m_xdebug.resize(getNumDof());
    m_xdebug.setConstant(std::numeric_limits<double>::signaling_NaN());
    m_vnphalf.resize(getNumDof());
    m_vnphalf.setConstant(std::numeric_limits<double>::signaling_NaN());

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
}

void BARodStepper::checkInternalConsistency()
{
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

    // TODO: Further, check that they all belong to same rod or triangle obj
}

void BARodStepper::setNumThreads(int num_threads)
{
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
}

void BARodStepper::buildCollisionDetector()
{
    // Load positions for initial construction of the BVH
    RodSelectionType selected_rods;
    for (int i = 0; i < m_number_of_rods; i++)
        selected_rods.push_back(i);

    assert(m_xn.size() == m_num_dof);
    extractPositions(m_xn, selected_rods, 0.0);
    assert(m_vnphalf.size() == m_num_dof);
    extractVelocities(m_vnphalf, selected_rods);

    m_collision_detector = new CollisionDetectorType(m_geodata, m_edges, m_faces, m_dt, m_perf_param.m_skipRodRodCollisions,
            m_num_threads);

    m_collision_immune.resize(getNumVerts());

    if (m_perf_param.m_enable_penalty_response)
        enableImplicitPenaltyImpulses();

    m_initialLengths.resize(m_number_of_rods);
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 1; j < m_rods[*rod]->nv(); j++)
            m_initialLengths[*rod] += (m_rods[*rod]->getVertex(j) - m_rods[*rod]->getVertex(j - 1)).norm();
}

void BARodStepper::initializeSimulationList()
{
    assert(m_simulated_rods.empty());

    // For debugging purposes
#ifdef KEEP_ONLY_SOME_RODS
    WarningStream(g_log, "", MsgInfo::kOncePerMessage)
    << "WARNING: KEEP_ONLY_SOME_RODS: Simulating only a specified subset of rods!\n***********************************************************\n";
    std::set<int> keep_only;

    keep_only.insert(11);

    // Only the rods in the keep_only set are kept, the others are killed.
    for (int i = 0; i < m_number_of_rods; i++)
    if (keep_only.find(i) == keep_only.end())
    for (int j = 0; j < m_rods[i]->nv(); j++)
    m_rods[i]->setVertex(j, 0 * m_rods[i]->getVertex(j));
    else
    m_simulated_rods.push_back(i);
#else
    // Initially all rods passed from Maya will be simulated
    for (int i = 0; i < m_number_of_rods; i++)
        m_simulated_rods.push_back(i);
#endif
    DebugStream(g_log, "") << "This BARodStepper (" << this << ") will simulate " << m_simulated_rods.size() << " rods\n";

    m_killed_rods.clear();
    failed_collisions_rods.resize(m_number_of_rods);
    stretching_rods.resize(m_number_of_rods);
}

/**
 * Execution
 */
bool BARodStepper::execute()
{
    START_TIMER("BARodStepper::execute")

    m_perf_param.m_solver.resetNum();
    m_perf_param.m_collision.resetNum();
    m_perf_param.m_explosion.resetNum();
    m_perf_param.m_stretching.resetNum();

    DebugStream(g_log, "") << "Executing time step " << m_t << '\n';

    m_collision_detector->buildBVH();
    DebugStream(g_log, "") << "BVH has been rebuilt\n";

    if (!m_collision_disabled_rods.empty())
    {
        TraceStream(g_log, "") << "The following rods were collision-disabled in the previous time step: ";
        for (std::set<int>::const_iterator rod = m_collision_disabled_rods.begin(); rod != m_collision_disabled_rods.end(); rod++)
            TraceStream(g_log, "") << *rod << " ";
        TraceStream(g_log, "") << '\n';
    }
    m_collision_disabled_rods.clear();

    bool do_adaptive = true;
    bool result;

    int k = 0;
    for (int i = 0; i < m_number_of_rods; ++i)
    {
        // Extract masses from the new rod
        for (ElasticRod::vertex_iter itr = m_rods[i]->vertices_begin(); itr != m_rods[i]->vertices_end(); ++itr)
        {
            if (m_rods[i]->getBoundaryCondition()->isVertexScripted((*itr).idx()))
                m_masses[k++] = std::numeric_limits<double>::infinity();
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

    DebugStream(g_log, "") << "Time step finished, " << selected_rods.size() << " rods remaining out of " << m_rods.size()
            << '\n';
    DebugStream(g_log, "") << m_perf_param.m_solver.sumMessage() << '\n';
    DebugStream(g_log, "") << m_perf_param.m_collision.sumMessage() << '\n';
    DebugStream(g_log, "") << m_perf_param.m_explosion.sumMessage() << '\n';
    DebugStream(g_log, "") << m_perf_param.m_stretching.sumMessage() << '\n';

    STOP_TIMER("BARodStepper::execute")

#ifdef TIMING_ON
    // This is not using TextLog because std::setw is not supported. TODO: you know.
    std::cout << "Cumulative timing results (entire run up to this point)\n";
    std::cout << "========================================\n";
    Timer::report();
    std::cout << "========================================\n";
#endif

    return result;
}

bool BARodStepper::nonAdaptiveExecute(double dt, RodSelectionType& selected_rods)
{
    setDt(dt);
    setTime(m_t + dt);
    //for (int i = 0; i < m_scripting_controllers.size(); ++i)
    //  m_scripting_controllers[i]->setTime(m_t);
    step(selected_rods);
    return (selected_rods.size() == 0);
}

bool BARodStepper::adaptiveExecute(double dt, RodSelectionType& selected_rods)
{
    START_TIMER("BARodStepper::adaptiveExecute")

    DebugStream(g_log, "") << "adaptiveExecute at level " << m_level << " with " << selected_rods.size() << " rod(s), m_t = "
            << m_t << ", dt = " << dt << '\n';

    // Backup all selected rods
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        m_rodbackups[*rod].backupRod(*m_rods[*rod]);

    // Backup all objects
    for (int i = 0; i < (int) m_triangle_meshes.size(); ++i)
        m_objbackups[i].backupMesh(*m_triangle_meshes[i]);

    // Backup the current simulation time
    double time = m_t;

    // Set the desired timestep
    setDt(dt);
    // Advance the current time
    setTime(m_t + dt);

    // Attempt a full time step
    step(selected_rods);

    if (m_simulationFailed)
    {
        WarningStream(g_log, "", MsgInfo::kOncePerMessage) << "t = " << m_t
                << ": **** SIMULATION FAILED AND IS NOW STOPPED! ****\n***********************************************************\n";
        STOP_TIMER("BARodStepper::adaptiveExecute")
        return true;
    }
    if (selected_rods.empty()) // Success!
    {
        TraceStream(g_log, "") << "t = " << m_t << ": adaptiveExecute has simulated (or killed) all rods\n";
        STOP_TIMER("BARodStepper::adaptiveExecute")
        return true;
    }
    // Otherwise do two half time steps
    DebugStream(g_log, "") << "t = " << m_t << ": adaptiveExecute left " << selected_rods.size() << " rods for substepping\n";
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

    DebugStream(g_log, "") << "t = " << m_t << " selected_rods: adaptiveExecute substepping (part 1) " << selected_rods.size()
            << " rods\n";

    bool first_success = adaptiveExecute(0.5 * dt, selected_rods);
    if (!first_success)
    {
        setDt(dt);
        STOP_TIMER("BARodStepper::adaptiveExecute")
        return false;
    }
    DebugStream(g_log, "") << "t = " << m_t << " selected_rods: adaptiveExecute substepping (part 2) "
            << selected_rods_2.size() << " rods\n";

    // Remove from the rod selection any one that might have been killed during the first time step
    for (RodSelectionType::iterator rod = selected_rods_2.begin(); rod != selected_rods_2.end(); rod++)
        if (find(m_simulated_rods.begin(), m_simulated_rods.end(), *rod) == m_simulated_rods.end())
        {
            DebugStream(g_log, "") << "Erasing from second time step rod number " << *rod << '\n';
            selected_rods_2.erase(rod--);
        }
    bool second_success = adaptiveExecute(0.5 * dt, selected_rods_2);
    if (!second_success)
    {
        setDt(dt);
        STOP_TIMER("BARodStepper::adaptiveExecute")
        return false;
    }
    TraceStream(g_log, "") << "Finished two adaptive steps\n";
    setDt(dt);
    m_level--;

    STOP_TIMER("BARodStepper::adaptiveExecute")
    return first_success && second_success;
}

void BARodStepper::step(RodSelectionType& selected_rods)
{
    START_TIMER("BARodStepper::step")

    if (m_simulationFailed) // We stopped simulating already
    {
        STOP_TIMER("BARodStepper::step")
        return;
    }

    DebugStream(g_log, "") << "t = " << m_t << ": BARodStepper::step() begins with " << selected_rods.size() << " rods\n";

    step_setup(selected_rods);

    step_dynamic(selected_rods);

    step_collision(selected_rods);

    step_failure(selected_rods);

    START_TIMER("BARodStepper::step/penalty");
    // Clean up penalty collisions list
    for (std::list<Collision*>::iterator i = m_penalty_collisions.begin(); i != m_penalty_collisions.end(); i++)
        delete *i;
    m_penalty_collisions.clear();
    // Clear existing penalties
    for (std::vector<RodPenaltyForce*>::const_iterator penalty_force = m_implicit_pnlty_forces.begin(); penalty_force
            != m_implicit_pnlty_forces.end(); penalty_force++)
        (*penalty_force)->clearProximityCollisions();
    STOP_TIMER("BARodStepper::step/penalty");

    STOP_TIMER("BARodStepper::step");
}

void BARodStepper::step_setup(const RodSelectionType& selected_rods)
{
    assert(m_edges.size() == m_edge_radii.size());
    assert((int) m_masses.size() == m_xn.size() / 3);
    assert(m_xn.size() == m_xnp1.size());
    assert(m_xn.size() == m_vnphalf.size());
    assert(m_rod_labels.size() == 0 || m_rod_labels.size() == m_number_of_rods);

    // Sanity check to ensure rods are not "internally colliding" because radius is bigger than edge length
#ifndef NDEBUG
    for (int i = 0; i < (int) m_number_of_rods; ++i)
        ensureNoCollisionsByDefault(*m_rods[i]);

    // Sanity check to ensure different parts of sim have same time/timetep
    for (int i = 0; i < (int) m_scripting_controllers.size(); ++i)
        assert(m_scripting_controllers[i]->getTime() == m_t);
    for (int i = 0; i < (int) m_scripting_controllers.size(); ++i)
        assert(m_scripting_controllers[i]->getDt() == m_dt);
    for (int i = 0; i < (int) m_steppers.size(); ++i)
        assert(m_steppers[i]->getTimeStep() == m_dt);
    for (int i = 0; i < (int) m_number_of_rods; ++i)
        assert(m_rods[i]->getTimeStep() == m_dt);
#endif

    START_TIMER("BARodStepper::step/explo");
    if (m_perf_param.m_enable_explosion_detection)
    {
        TraceStream(g_log, "") << "BARodStepper::step: computing start forces" << '\n';
        computeForces(m_startForces, selected_rods);
    }STOP_TIMER("BARodStepper::step/explo");

    START_TIMER("BARodStepper::step/setup");
    // Save the pre-timestep positions
    extractPositions(m_xn, selected_rods, m_t - m_dt);
    STOP_TIMER("BARodStepper::step/setup");

    START_TIMER("BARodStepper::step/immune");
    // Determine which vertex are to be considered collision-immune for this step
    computeImmunity(selected_rods);
    STOP_TIMER("BARodStepper::step/immune");

    // Step rods forward ignoring collisions
    START_TIMER("BARodStepper::step/scripting");
    // Step scripted objects forward, set boundary conditions
    for (std::vector<ScriptingController*>::const_iterator scripting_controller = m_scripting_controllers.begin(); scripting_controller
            != m_scripting_controllers.end(); scripting_controller++)
        (*scripting_controller)->execute();
    STOP_TIMER("BARodStepper::step/scripting");

    START_TIMER("BARodStepper::step/penalty");
    // Jungseock's implicit penalty
    assert(m_penalty_collisions.empty());
    // The penalty collisions list is used to create penalty forces. All that is deleted and cleared at the end of this step.
    if (m_perf_param.m_enable_penalty_response)
        setupPenaltyForces(m_penalty_collisions, selected_rods);
    STOP_TIMER("BARodStepper::step/penalty");

    if (m_perf_param.m_enable_explosion_detection)
    {
        TraceStream(g_log, "") << "BARodStepper::step: computing pre-dynamic forces" << '\n';
        computeForces(m_preDynamicForces, selected_rods);
    }
}

void BARodStepper::step_dynamic(const RodSelectionType& selected_rods)
{
    START_TIMER("BARodStepper::step/steppers");

    bool dependable_solve = true;

    START_TIMER("BARodStepper::step/setup");
    // Prepare list of steppers to be executed.
    std::vector<RodTimeStepper*> selected_steppers;
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        selected_steppers.push_back(m_steppers[*rod]);
    STOP_TIMER("BARodStepper::step/setup");

#ifdef HAVE_OPENMP
    //#pragma omp parallel for
#endif
    for (int i = 0; i < selected_steppers.size(); i++)
    {
        RodTimeStepper* const stepper = selected_steppers[i];

        bool result = stepper->execute();
        if (!result)
            TraceStream(g_log, "") << stepper->getDiffEqSolver().getName() << " solver for rod "
                    << stepper->getRod()->globalRodIndex << " failed to converge after " << stepper->getMaxIterations()
                    << " iterations\n";
        dependable_solve = dependable_solve && result;
    }

    STOP_TIMER("BARodStepper::step/steppers");

    TraceStream(g_log, "") << "Dynamic step is " << (dependable_solve ? "" : "not ") << "entirely dependable"
            << (dependable_solve ? " :-)\n" : "!\n");

    // If we do rod-rod collisions (meaning no selective adaptivity) and global dependability failed, we might as well stop here.
    if (!m_perf_param.m_skipRodRodCollisions && !dependable_solve)
    {
        WarningStream(g_log, "", MsgInfo::kOncePerMessage) << "t = " << m_t
                << " selected_rods: step() failed (due to rod-rod) for " << selected_rods.size() << " rods\n";
        STOP_TIMER("BARodStepper::step")
        return;
    }

    START_TIMER("BARodStepper::step/immune");
    // Mark invalid rods as entirely collision-immune, so we don't waste time on colliding them.
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        if (!m_steppers[*rod]->HasSolved())
            setRodImmunity(*rod, true);
    STOP_TIMER("BARodStepper::step/immune");

    START_TIMER("BARodStepper::step/setup");
    // Post time step position
    extractPositions(m_xnp1, selected_rods, m_t);
    // Average velocity over the timestep just completed
    m_vnphalf = (m_xnp1 - m_xn) / m_dt;
    STOP_TIMER("BARodStepper::step/setup");

    START_TIMER("BARodStepper::step/inexten");
    for (int i = 0; i < selected_steppers.size(); i++)
        applyInextensibilityVelocityFilter(selected_steppers[i]->getRod()->globalRodIndex);
    STOP_TIMER("BARodStepper::step/inexten");
}

void BARodStepper::step_collision(const RodSelectionType& selected_rods)
{
    for (int i = 0; i < m_number_of_rods; i++)
        failed_collisions_rods[i] = stretching_rods[i] = false;

    START_TIMER("BARodStepper::step/response");
    DebugStream(g_log, "") << "Starting collision response\n";

    if (m_perf_param.m_enable_explosion_detection)
    {
        TraceStream(g_log, "") << "Computing pre-collision forces\n";
        computeForces(m_preCollisionForces, selected_rods);
    }

    if (m_perf_param.m_collision.m_max_iterations > 0)
        if (!executeIterativeInelasticImpulseResponse(failed_collisions_rods, stretching_rods))
            TraceStream(g_log, "") << "Some collision responses failed!\n";

    TraceStream(g_log, "") << "Finished collision response\n";
    STOP_TIMER("BARodStepper::step/response");

    START_TIMER("BARodStepper::step/setup");

    // Store the response part for visualization
    m_vnresp = m_vnphalf - (m_xnp1 - m_xn) / m_dt;

    // Compute final positions from corrected velocities
    m_xnp1 = m_xn + m_dt * m_vnphalf;

#ifndef NDEBUG
    // Ensure boundary conditions respected by corrected positions
    // For each selected rod
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
    {
        RodBoundaryCondition* boundary = m_rods[*rod]->getBoundaryCondition();
        int rodbase = m_base_dof_indices[*rod];

        // For each vertex of the current rod
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
        {
            // If that vertex has a prescribed position
            if (boundary->isVertexScripted(j))
            {
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

    TraceStream(g_log, "") << "About to update properties" << '\n';

    // Update frames and such in the rod (Is this correct? Will this do some extra stuff?)
    for (RodSelectionType::const_iterator selected_rod = selected_rods.begin(); selected_rod != selected_rods.end(); selected_rod++)
        m_rods[*selected_rod]->updateProperties();

    // Sanity check to ensure rod's internal state is consistent
#ifndef NDEBUG
    for (RodSelectionType::const_iterator selected_rod = selected_rods.begin(); selected_rod != selected_rods.end(); selected_rod++)
        if (m_steppers[*selected_rod]->HasSolved())
            m_rods[*selected_rod]->verifyProperties();
#endif

    STOP_TIMER("BARodStepper::step/setup");
}

void BARodStepper::step_failure(RodSelectionType& selected_rods)
{
    START_TIMER("BARodStepper::step/explo");

    // Explosion detection
    std::vector<bool> exploding_rods(m_number_of_rods);
    if (m_perf_param.m_enable_explosion_detection)
    {
        DebugStream(g_log, "") << "About to detect explosions" << '\n';
        computeForces(m_endForces, selected_rods);
        checkExplosions(exploding_rods, failed_collisions_rods, selected_rods);
    }
    // Decide whether to substep or kill some rods

    STOP_TIMER("BARodStepper::step/explo");

    // Check lengths again... Why is this necessary?
    checkLengths(stretching_rods);

    START_TIMER("BARodStepper::step/exception");

    int rod_kill = 0;
    for (RodSelectionType::iterator rodit = selected_rods.begin(); rodit != selected_rods.end(); rodit++)
    {
        int rodidx = *rodit;

        bool solveFailure = !m_steppers[rodidx]->HasSolved();
        bool explosion = exploding_rods[rodidx];
        bool collisionFailure = failed_collisions_rods[rodidx];
        bool stretching = stretching_rods[rodidx];

        bool substep = (solveFailure && m_level < m_perf_param.m_solver.m_max_substeps) //
                || (explosion && m_level < m_perf_param.m_explosion.m_max_substeps) //
                || (collisionFailure && m_level < m_perf_param.m_collision.m_max_substeps) //
                || (stretching && m_level < m_perf_param.m_stretching.m_max_substeps);

        bool killRod = (solveFailure && m_perf_param.m_solver.m_in_case_of == FailureMode::KillTheRod) //
                || (explosion && m_perf_param.m_explosion.m_in_case_of == FailureMode::KillTheRod)//
                || (collisionFailure && m_perf_param.m_collision.m_in_case_of == FailureMode::KillTheRod) //
                || (stretching && m_perf_param.m_stretching.m_in_case_of == FailureMode::KillTheRod);

        bool haltSim = (solveFailure && m_perf_param.m_solver.m_in_case_of == FailureMode::HaltSimulation) //
                || (explosion && m_perf_param.m_explosion.m_in_case_of == FailureMode::HaltSimulation) //
                || (collisionFailure && m_perf_param.m_collision.m_in_case_of == FailureMode::HaltSimulation) //
                || (stretching && m_perf_param.m_stretching.m_in_case_of == FailureMode::HaltSimulation);

        if (substep) // Only in that case keep the rod in the selected list
            continue;
        else if (killRod)
        {
            // DEBUG
            if (solveFailure && m_perf_param.m_solver.m_in_case_of == FailureMode::KillTheRod)
                ++m_perf_param.m_solver;
            if (collisionFailure && m_perf_param.m_collision.m_in_case_of == FailureMode::KillTheRod)
                ++m_perf_param.m_collision;
            if (explosion && m_perf_param.m_explosion.m_in_case_of == FailureMode::KillTheRod)
                ++m_perf_param.m_explosion;
            if (stretching && m_perf_param.m_stretching.m_in_case_of == FailureMode::KillTheRod)
                ++m_perf_param.m_stretching;

            rod_kill++;
            killTheRod(*rodit);
        }
        else if (haltSim)
            m_simulationFailed = true;
        else if (m_useKineticDamping)
        {
            // TODO (sainsley) : add flag check here
            //     std::cout << "treatment: accept this step as-is" << '\n';
            // at this point, the step is either successful, or includes only ignorable errors

            // Accept this step

            ElasticRod* rod = m_rods[rodidx];

            // std::cout << "KE[" << rodidx << "] = " << rod->computeKineticEnergy() << '\n';

            // Apply kinetic damping
            rod->recordKineticEnergy();
            if (rod->isKineticEnergyPeaked())
            {
                // std::cout << "Zeroing energy for rod " << rodidx << '\n';
                for (int i = 0; i < rod->nv(); ++i)
                {
                    rod->setVelocity(i, Vec3d(0, 0, 0));
                }
            }
        }
        selected_rods.erase(rodit--);
        // the -- compensates for the erased rod; this is dangerous since it assumes array (rather than linked-list) semantics for selected_rods
        // Yes I know, I'm surprised this even works.
    }

    STOP_TIMER("BARodStepper::step/exception");

    if (rod_kill)
        NoticeStream(g_log, "") << "This step killed " << rod_kill << " rods\n";

    if (selected_rods.size() > 0)
        TraceStream(g_log, "") << "Step finished, " << selected_rods.size() << " rods must be substepped\n";
    else
        TraceStream(g_log, "") << "Step finished, all rods treated (either successful step, removed, or errors ignored)\n";

    if (m_killed_rods.size() > 0)
    {
        std::ostringstream ost;
        for (RodSelectionType::const_iterator rod = m_killed_rods.begin(); rod != m_killed_rods.end(); rod++)
            ost << *rod << ' ';
        DebugStream(g_log, "") << "List of rods killed: " << ost.str() << '\n';
    }
}

/*
 * Extracting/Restoring
 */
void BARodStepper::extractPositions(VecXd& positions, const RodSelectionType& selected_rods, const double time) const
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
                const Vec3d desiredposition = boundary->getDesiredVertexPosition(j, time);
                const Vec3d actualvalue = positions.segment<3> (rodbase + 3 * j);
                assert(approxEq(desiredposition, actualvalue, 1.0e-12));
            }
    }
#endif
}

void BARodStepper::extractVelocities(VecXd& velocities, const RodSelectionType& selected_rods) const
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

void BARodStepper::restorePositions(const VecXd& positions, const RodSelectionType& selected_rods)
{
    assert(m_number_of_rods == m_base_dof_indices.size());

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
            if (!m_rods[*rod]->getBoundaryCondition()->isVertexScripted(j))
                m_rods[*rod]->setVertex(j, positions.segment<3> (m_base_dof_indices[*rod] + 3 * j));
}

void BARodStepper::restoreVelocities(const VecXd& velocities, const RodSelectionType& selected_rods)
{
    assert(m_number_of_rods == m_base_dof_indices.size());

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
            if (!m_rods[*rod]->getBoundaryCondition()->isVertexScripted(j))
                m_rods[*rod]->setVelocity(j, velocities.segment<3> (m_base_dof_indices[*rod] + 3 * j));
}

void BARodStepper::restoreResponses(const VecXd& responses, const RodSelectionType& selected_rods)
{
    assert(m_number_of_rods == m_base_dof_indices.size());

    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        for (int j = 0; j < m_rods[*rod]->nv(); ++j)
            if (!m_rods[*rod]->getBoundaryCondition()->isVertexScripted(j))
                m_rods[*rod]->setResponse(j, responses.segment<3> (m_base_dof_indices[*rod] + 3 * j));
}

/**
 * Enabling/Disabling
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

    TraceStream(g_log, "") << "Implicit penalty response is now enabled\n";

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

void BARodStepper::setNumInelasticIterations(const int& num_itr)
{
    assert(num_itr >= 0);
    m_perf_param.m_collision.m_max_iterations = num_itr;
}

void BARodStepper::computeImmunity(const RodSelectionType& selected_rods)
{
    // Initially, the rods are all collision-immune
    for (std::vector<bool>::iterator i = m_collision_immune.begin(); i != m_collision_immune.begin() + m_obj_start; i++)
        *i = true;
    // Except for the ones in the selected rods list
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
        // Unless they are in the collision disabled set
        if (m_collision_disabled_rods.find(*rod) == m_collision_disabled_rods.end())
            setRodImmunity(*rod, false);

    // On the top of that, find the initial vertex-face intersections and mark them collision-immune
    std::list<Collision*> collisions;
    m_collision_detector->getCollisions(collisions, EdgeFace);
    while (!collisions.empty())
    {
        EdgeFaceIntersection* intersection = dynamic_cast<EdgeFaceIntersection*> (collisions.front());
        assert(intersection != NULL);

        m_collision_immune[intersection->v0] = true;

        collisions.pop_front();
        delete intersection;
    }
}

/**
 * Utilities
 */
const double BARodStepper::computeTotalForceNorm() const
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
    // std::cout << "settingTime in BARodStepper to be " << m_t << '\n';

    // Set the time for the rod controllers
    for (int i = 0; i < (int) m_steppers.size(); ++i)
        m_steppers[i]->setTime(m_t);

    // Set the time for the scripted object controllers
    for (int i = 0; i < (int) m_scripting_controllers.size(); ++i)
        m_scripting_controllers[i]->setTime(m_t);
}

const double BARodStepper::getDt() const
{
    return m_dt;
}

const double BARodStepper::getTime() const
{
    //std::cout << "BARodStepper::getTime() = " << m_t << '\n';
    return m_t;
}

void BARodStepper::skipRodRodCollisions(bool skipRodRodCollisions)
{
    TraceStream(g_log, "") << "Switching rod-rod collisions " << (skipRodRodCollisions ? "OFF" : "ON") << '\n';
    m_perf_param.m_skipRodRodCollisions = skipRodRodCollisions;

    if (m_collision_detector)
        m_collision_detector->setSkipRodRodCollisions(skipRodRodCollisions);
}

void BARodStepper::setRodLabels(const std::vector<std::string>& rod_labels)
{
    assert(rod_labels.size() == m_number_of_rods);
    m_rod_labels = rod_labels;
}

const int BARodStepper::getContainingRod(int vert_idx) const
{
    assert(vert_idx >= 0);
    assert(vert_idx < getNumVerts());

    return upper_bound(m_base_vtx_indices.begin(), m_base_vtx_indices.end(), vert_idx) - m_base_vtx_indices.begin() - 1;
}

bool BARodStepper::isRodVertex(int vert) const
{
    assert(vert >= 0);
    assert(vert < getNumVerts());

    // Is a vertex if index is less than start of object vertices in global array
    return vert < m_obj_start;
}

const int BARodStepper::getNumDof() const
{
    assert(m_num_dof >= 0);
    return m_num_dof;
}

const int BARodStepper::getNumVerts() const
{
    assert(m_num_dof % 3 == 0);
    return m_num_dof / 3;
}

void BARodStepper::ensureCircularCrossSection(const ElasticRod& rod) const
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

void BARodStepper::ensureNoCollisionsByDefault(const ElasticRod& rod) const
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

void BARodStepper::setRodImmunity(const int rodIdx, const bool immune)
{
    for (int j = 0; j < m_rods[rodIdx]->nv(); ++j)
        m_collision_immune[m_base_vtx_indices[rodIdx] + j] = immune;
}

void BARodStepper::killTheRod(int rod) // TODO: remove the rod properly in Maya
{
    m_killed_rods.push_back(rod);
    for (int j = 0; j < m_rods[rod]->nv(); j++)
        m_rods[rod]->setVertex(j, 0 * m_rods[rod]->getVertex(j));
    m_simulated_rods.erase(find(m_simulated_rods.begin(), m_simulated_rods.end(), rod));

#ifndef KEEP_ONLY_SOME_RODS
    assert(m_simulated_rods.size() + m_killed_rods.size() == m_number_of_rods);
#endif
}

void BARodStepper::computeForces(std::vector<VecXd*> Forces, const RodSelectionType& selected_rods) const
{
    for (RodSelectionType::const_iterator rod = selected_rods.begin(); rod != selected_rods.end(); rod++)
    {
#ifndef NDEBUG
        if (m_steppers[*rod]->HasSolved())
            m_rods[*rod]->verifyProperties(); // Sanity check to ensure rod's internal state is consistent
#endif
        Forces[*rod]->setZero();
        m_rods[*rod]->computeForces(*Forces[*rod]); // If we want the internal forces only
        //  m_steppers[*rod]->evaluatePDot(*Forces[*rod]); // If we want the external forces too
    }
}

void BARodStepper::addRod(ElasticRod* rod, RodTimeStepper* stepper)
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

    failed_collisions_rods.resize(m_number_of_rods);
    stretching_rods.resize(m_number_of_rods);

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

void BARodStepper::removeRod(int rodIdx)
{
    killTheRod(rodIdx);
}

/**
 * Iterative collision resolution
 */
bool BARodStepper::executeIterativeInelasticImpulseResponse(std::vector<bool>& failed_collisions_rods,
        std::vector<bool>& stretching_rods)
{
    bool all_rods_collisions_ok = true;

    // As a safeguard, check whether the solver left some rods stretched, in which case they will be declared collision-immune
    checkLengths(stretching_rods);

    // Detect continuous time collisions
    std::list<Collision*> collisions_list;
    TraceStream(g_log, "") << "Detecting collisions...\n";
    m_collision_detector->getCollisions(collisions_list, ContinuousTime);
    TraceStream(g_log, "") << "Initial potential collisions: " << m_collision_detector->m_potential_collisions << "\n";

    // Iteratively apply inelastic impulses
    for (int itr = 0; !collisions_list.empty() && itr < m_perf_param.m_collision.m_max_iterations; ++itr)
    {
        TraceStream(g_log, "") << "CTcollision response iteration " << itr << '\n';
        TraceStream(g_log, "") << "Detected " << collisions_list.size() << " continuous time collisions (potential: "
                << m_collision_detector->m_potential_collisions << ")\n";

        // Just sort the collision times to maintain some rough sense of causality
        collisions_list.sort(CompareTimes);

        for (; !collisions_list.empty();)
        {
            // Find the first collision on a unique rod, put it in first_collisions_list and remove it from collisions_list.
            // Subsequent collisions on the same rod are kept in collisions_list.
            std::set<int> already_collided_rods;
            std::vector<Collision*> first_collisions_list;
            for (std::list<Collision*>::iterator col_it = collisions_list.begin(); col_it != collisions_list.end(); col_it++)
            {
                CTCollision* collision = dynamic_cast<CTCollision*> (*col_it);
                if (collision)
                {
                    // Only keep collisions involving a rod we see for the first time
                    const int collidingRodIdx = getContainingRod(collision->GetRodVertex());
                    if (already_collided_rods.find(collidingRodIdx) == already_collided_rods.end())
                    {
                        already_collided_rods.insert(collidingRodIdx);
                        first_collisions_list.push_back(collision);
                        collisions_list.erase(col_it--);
                    }
                }
            }
            TraceStream(g_log, "") << "of which " << collisions_list.size() << " are on different rods\n";

            // Now apply response to the first collisions.
#ifdef HAVE_OPENMP
#pragma omp parallel for
#endif
            for (std::vector<Collision*>::iterator col_it = first_collisions_list.begin(); col_it
                    != first_collisions_list.end(); col_it++)
            {
                const CTCollision* collision = dynamic_cast<CTCollision*> (*col_it);
                assert(collision != NULL);

                // So, which rod was it? Here we are assuming that this is not a rod-rod collision
                const int collidingRodIdx = getContainingRod(collision->GetRodVertex());
                // If the rod has not yet failed collision (this time step)
                if (!failed_collisions_rods[collidingRodIdx])
                {
                    ElasticRod* const collidingRod = m_rods[collidingRodIdx];
                    RodSelectionType oneRodList; // So we can call existing routines that take a rod list. Not very efficient.
                    oneRodList.push_back(collidingRodIdx);

                    // Save pre-impulse velocities
                    VecXd velBackup(3 * collidingRod->nv());
                    for (int i = 0; i < 3 * collidingRod->nv(); i++)
                        velBackup[i] = m_vnphalf[m_base_dof_indices[collidingRodIdx] + i];

                    // Position the rod at end of time step so the compliance works on the "right" configuration
                    restorePositions(m_xn + m_dt * m_vnphalf, oneRodList);
                    collidingRod->updateProperties();

                    if (!exertCompliantInelasticImpulse(collision))
                        failed_collisions_rods[collidingRodIdx] = true;

                    if (m_perf_param.m_enable_explosion_detection)
                    {
                        // Prepare the rod for explosion checking
                        restorePositions(m_xn + m_dt * m_vnphalf, oneRodList);
                        collidingRod->updateProperties();
                        m_endForces[collidingRodIdx]->setZero();
                        collidingRod->computeForces(*m_endForces[collidingRodIdx]);

                        // Line search on the impulse size: if full impulse causes explosion, halve it and so on.
                        int splitCounter = 5;
                        double splitFactor = 1.0;
                        while (hadExplosion(collidingRodIdx))
                        {
                            splitFactor *= 0.5;
                            splitCounter--;
                            TraceStream(g_log, "") << "Downsizing impulse for rod " << collidingRodIdx << " to " << splitFactor
                                    * 100.0 << "%\n";
                            // Collision is marked as failed if the original impulse was exploding NOT
                            // failed_collisions_rods[collidingRodIdx] = true;
                            // Interpolate velocity between pre-impulse (velBackup) and (resized) post-impulse (m_vnphalf)
                            for (int v = 0; v < collidingRod->nv(); ++v)
                            {
                                m_vnphalf.segment<3> (m_base_dof_indices[collidingRodIdx] + 3 * v) = 0.5
                                        * m_vnphalf.segment<3> (m_base_dof_indices[collidingRodIdx] + 3 * v) + 0.5
                                        * velBackup.segment<3> (3 * v);
                                m_collision_immune[m_base_vtx_indices[collidingRodIdx] + v] = true;
                            }
                            if (splitCounter == 0)
                                break;

                            // Prepare the rod again for explosion checking
                            restorePositions(m_xn + m_dt * m_vnphalf, oneRodList);
                            collidingRod->updateProperties();
                            m_endForces[collidingRodIdx]->setZero();
                            collidingRod->computeForces(*m_endForces[collidingRodIdx]);
                        }
                    }

                    // Test for stretching and if so, mark the rod immune
                    stretching_rods[collidingRodIdx] = !checkLength(collidingRodIdx);
                    if (stretching_rods[collidingRodIdx])
                    {
                        // Declare the rod collision-immune for the rest of the time step
                        for (int j = 0; j < collidingRod->nv(); ++j)
                            m_collision_immune[m_base_vtx_indices[collidingRodIdx] + j] = true;
                    }
                }
                delete collision;
            }
            // All first collisions have been treated, see what remains
            first_collisions_list.clear();
            already_collided_rods.clear();
            m_collision_detector->updateCollisions(collisions_list);
        }

        // Detect remaining collisions (including at the end of the last iteration, so we know what failed)
        TraceStream(g_log, "") << "Detecting collisions...\n";
        m_collision_detector->getCollisions(collisions_list, ContinuousTime, false); // No need to update the mesh bvh bounding boxes.
    }

    // If collisions were still detected before exiting the loop
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

    return all_rods_collisions_ok;
}

/**
 * Non-compliant Inelastic response
 */
void BARodStepper::exertVertexImpulse(const Vec3d& I, const double& m, const int& idx, VecXd& v)
{
    assert(m > 0.0);
    assert(idx >= 0);
    assert(idx < getNumVerts());
    assert(v.size() == getNumDof());

    v.segment<3> (3 * idx) += I / m;
}

void BARodStepper::exertEdgeImpulse(const Vec3d& I, const double& m0, const double& m1, const double& alpha, const int& idx0,
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

void BARodStepper::exertFaceImpulse(const Vec3d& I, const double& m0, const double& m1, const double& m2, const double& u,
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

void BARodStepper::exertInelasticImpulse(EdgeEdgeCTCollision& eecol)
{
    assert(eecol.isAnalysed());
    assert(eecol.e0_v0 >= 0);
    assert(eecol.e0_v0 < getNumVerts());
    assert(eecol.e0_v1 >= 0);
    assert(eecol.e0_v1 < getNumVerts());
    assert(eecol.e1_v0 >= 0);
    assert(eecol.e1_v0 < getNumVerts());
    assert(eecol.e1_v1 >= 0);
    assert(eecol.e1_v1 < getNumVerts());

    assert(eecol.GetCachedRelativeVelocity() < 0.0);

    // Add some extra "kick" to relative velocity to account for FPA errors
    //eecol.ApplyRelativeVelocityKick();

    Vec3d I = eecol.computeInelasticImpulse();

    exertEdgeImpulse(-I, m_masses[eecol.e0_v0], m_masses[eecol.e0_v1], eecol.s, eecol.e0_v0, eecol.e0_v1, m_vnphalf);
    exertEdgeImpulse(I, m_masses[eecol.e1_v0], m_masses[eecol.e1_v1], eecol.t, eecol.e1_v0, eecol.e1_v1, m_vnphalf);

    assert(eecol.computeRelativeVelocity() >= 0);
}

void BARodStepper::exertInelasticImpulse(VertexFaceCTCollision& vfcol)
{
    assert(vfcol.isAnalysed());
    assert(vfcol.v0 >= 0);
    assert(vfcol.v0 < getNumVerts());
    assert(vfcol.f0 >= 0);
    assert(vfcol.f0 < getNumVerts());
    assert(vfcol.f1 >= 0);
    assert(vfcol.f1 < getNumVerts());
    assert(vfcol.f2 >= 0);
    assert(vfcol.f2 < getNumVerts());

    assert(vfcol.GetCachedRelativeVelocity() < 0.0);

    Vec3d I = vfcol.computeInelasticImpulse();

    exertFaceImpulse(-I, m_masses[vfcol.f0], m_masses[vfcol.f1], m_masses[vfcol.f2], vfcol.u, vfcol.v, vfcol.w, vfcol.f0,
            vfcol.f1, vfcol.f2, m_vnphalf);
    exertVertexImpulse(I, m_masses[vfcol.v0], vfcol.v0, m_vnphalf);

    assert(vfcol.computeRelativeVelocity() >= 0);
}

/**
 * Compliant inelastic response
 */
static const double NORMAL_SCALING = 10.0;
static const double MATRIX_ASYMMETRY = 1.0e-6 * NORMAL_SCALING;

bool BARodStepper::exertCompliantInelasticImpulse(const CTCollision* cllsn)
{
    assert(cllsn->isAnalysed());
    const EdgeEdgeCTCollision* eecol = dynamic_cast<const EdgeEdgeCTCollision*> (cllsn);
    const VertexFaceCTCollision* vfcol = dynamic_cast<const VertexFaceCTCollision*> (cllsn);
    //TraceStream(g_log, "") << "BARodStepper:exertCompliantInelasticImpulse: pre-impulse e-e relative velocity = " << cllsn->GetRelativeVelocity() << '\n';

    if (eecol)
    {
        return exertCompliantInelasticEdgeEdgeImpulse(*eecol);
        //exertInelasticImpulse(cllsn.getEdgeEdge());
    }
    else if (vfcol)
    {
        return exertCompliantInelasticVertexFaceImpulse(*vfcol);
        //exertCompliantInelasticVertexFaceImpulse(cllsn.getVertexFace());
        //exertInelasticImpulse(cllsn.getVertexFace());
    }
    //TraceStream(g_log, "") << "BARodStepper:exertCompliantInelasticImpulse: post-impulse e-e relative velocity = " << cllsn->GetRelativeVelocity() << '\n';
}

bool BARodStepper::exertCompliantInelasticEdgeEdgeImpulse(const EdgeEdgeCTCollision& eecol)
{
    //   TraceStream(g_log, "") << "Edge-edge compliant inelastic impulse" << '\n';
    //   TraceStream(g_log, "") << eecol << '\n';

    // Determine if either edge is totally fixed
    bool rod0fixed = YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(m_geodata);
    bool rod1fixed = YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata);

    if (rod0fixed && rod1fixed) // Skip collisions in which both edges are fixed
        return true;
    else if (!rod0fixed && !rod1fixed)
        return exertCompliantInelasticEdgeEdgeImpulseBothFree(eecol);
    else
        return exertCompliantInelasticEdgeEdgeImpulseOneFree(eecol);
}

bool BARodStepper::exertCompliantInelasticVertexFaceImpulse(const VertexFaceCTCollision& vfcol)
{
    assert(vfcol.isAnalysed());

    // For now, assume vertex is free and entire face is fixed
    assert(!m_geodata.isVertexFixed(vfcol.v0));
    assert(YATriangle(vfcol.f0, vfcol.f1, vfcol.f2).IsFixed(m_geodata));

    // Determine which rod the free vertex belongs to
    int rodidx = getContainingRod(vfcol.v0);

    // Ensure the free vertex belongs to a rod
    assert(rodidx >= 0);
    assert(rodidx < (int) m_number_of_rods);

    ElasticRod* const colliding_rod = m_rods[rodidx];

    // If the rod has not solved properly, no need to compute its collision response
    if (!m_steppers[rodidx]->HasSolved())
    {
        DebugStream(g_log, "") << "WARNING: attempt to do vertex-face collision with non-dependable rod\n";
        return false;
    }

    // Determine which vertex of the rod the free vertex is
    int rodbase = m_base_dof_indices[rodidx];
    assert(rodbase % 3 == 0);
    int v0 = vfcol.v0 - rodbase / 3;
    assert(v0 >= 0);
    assert(v0 < colliding_rod->nv());

    TraceStream(g_log, "") << "CTcollision: vertex " << v0 << " of rod " << rodidx << '\n';

    // Compute the relative velocity of the collision
    assert(vfcol.computeRelativeVelocity() < 0.0);

    // Get storage for lhs of linear system, get a solver
    LinearSystemSolver* lss = m_solver_collection.getLinearSystemSolver(colliding_rod->ndof());
    MatrixBase* lhs = lss->m_lhs;
    assert(lhs != NULL);
    LinearSolverBase* solver = lss->m_solver;
    assert(solver != NULL);

    // Compute M - h^2*dF/dx
    computeCompliantLHS(lhs, colliding_rod);

    // Compute the 'base' index of the vertex
    int base0 = colliding_rod->vertIdx(v0, 0);
    assert(base0 % 4 == 0);

    int ndof = colliding_rod->ndof();

    // Compute the 'global' normal
    std::vector<VecXd> normal;

    normal.push_back(VecXd(ndof));
    normal.back().setZero();
    normal.back().segment<3> (base0) = vfcol.GetNormal() * NORMAL_SCALING;

    // Compute three normals per scripted vertex
    const std::vector<int>& scriptedverts = colliding_rod->getBoundaryCondition()->scriptedVertices();
    for (std::vector<int>::const_iterator vertex = scriptedverts.begin(); vertex != scriptedverts.end(); ++vertex)
    {
        for (int j = 0; j < 3; ++j)
        {
            normal.push_back(VecXd(ndof));
            normal.back().setZero();
            // Zero the velocity of this particular dof of this particular vertex
            normal.back()(colliding_rod->vertIdx(*vertex, j)) = 1.0;
        }
    }
    int numconstraints = normal.size();
    int nvdof = 3 * colliding_rod->nv();

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
        assert((int) scriptedverts[i] < colliding_rod->nv());
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 0);
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 1);
        desired_values(curdof++) = m_vnphalf(rodbase + 3 * scriptedverts[i] + 2);
    }
    assert(curdof == numconstraints);
    assert(desired_values.size() == numconstraints);

    // Currently, all fixed vertex constraints normalized
#ifndef NDEBUG
    assert(approxEq(normal[0].norm(), NORMAL_SCALING));
    for (int i = 1; i < numconstraints; ++i)
        assert(approxEq(normal[i].norm(), 1.0, 1.0e-9));
#endif

    std::vector<VecXd> ntilde(numconstraints);
    for (std::vector<VecXd>::iterator constraint = ntilde.begin(); constraint != ntilde.end(); ++constraint)
    {
        if (solver->solve(*constraint, normal[constraint - ntilde.begin()]) < 0)
        {
            WarningStream(g_log, "")
                    << "BARodStepper::exertCompliantInelasticVertexFaceImpulse: linear solve for ntilde failed\n";
            return false;
        }
    }

    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (ElasticRod::edge_iter edge = colliding_rod->edges_begin(); edge != colliding_rod->edges_end(); ++edge)
        for (int i = 0; i < numconstraints; ++i)
            assert(ntilde[i](colliding_rod->edgeIdx(*edge)) == 0.0);
#endif

    // Vectors restricted to just vertex DoFs. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnnormal(numconstraints);
    extractVertexDOF(posnnormal, normal, colliding_rod);
    std::vector<VecXd> posnntilde(numconstraints);
    extractVertexDOF(posnntilde, ntilde, colliding_rod);

    bool success = changeVelocityOneFree(posnnormal, posnntilde, scriptedverts, desired_values, numconstraints, rodbase, nvdof,
            vfcol);

    // applyInextensibilityVelocityFilter(rodidx);

    return success;
}

bool BARodStepper::exertCompliantInelasticEdgeEdgeImpulseOneFree(const EdgeEdgeCTCollision& eecol)
{
    assert(eecol.isAnalysed());

    // Determine if either edge is fixed. Exactly one of them must be.
    bool rod0fixed = YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(m_geodata);
    bool rod1fixed = YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata);
    assert(rod0fixed != rod1fixed);

    // Compute the relative velocity, which must be negative for a collision to have happened
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

    ElasticRod* const colliding_rod = m_rods[rodidx];

    // If the rod has not solved properly, no need to compute its collision response
    if (!m_steppers[rodidx]->HasSolved())
    {
        DebugStream(g_log, "") << "WARNING: attempt to do edge-edge collision with non-dependable rod\n";
        return;
    }

    TraceStream(g_log, "") << "CTcollision: rod " << rodidx << " vs. mesh\n";

    // Convert the vertices' global indices to rod indices
    int rodbase = m_base_dof_indices[rodidx];

    assert(rodbase % 3 == 0);
    v0 -= rodbase / 3;
    v1 -= rodbase / 3;
    assert(v0 >= 0);
    assert(v0 < colliding_rod->nv());
    assert(v1 >= 0);
    assert(v1 < colliding_rod->nv());

    // Get storage for lhs of linear system, get a solver
    LinearSystemSolver* lss = m_solver_collection.getLinearSystemSolver(colliding_rod->ndof());
    MatrixBase* lhs = lss->m_lhs;
    assert(lhs != NULL);
    LinearSolverBase* solver = lss->m_solver;
    assert(solver != NULL);

    // Compute M - h^2*dF/dx
    computeCompliantLHS(lhs, colliding_rod);

    // Compute the 'base' index of each vertex
    int base0 = colliding_rod->vertIdx(v0, 0);
    assert(base0 % 4 == 0);
    int base1 = colliding_rod->vertIdx(v1, 0);
    assert(base1 % 4 == 0);

    int ndof = colliding_rod->ndof();
    int nvdof = 3 * colliding_rod->nv();

    // Compute the 'global' normal
    std::vector<VecXd> normal;

    normal.push_back(VecXd(ndof));
    normal.back().setZero();
    normal.back().segment<3> (base0) = (1 - u) * eecol.GetNormal() * NORMAL_SCALING;
    normal.back().segment<3> (base1) = u * eecol.GetNormal() * NORMAL_SCALING;

    // If rod 0 has the free edge, need to invert sign of normal
    if (!rod0fixed)
    {
        normal.back().segment<3> (base0) *= -1.0;
        normal.back().segment<3> (base1) *= -1.0;
    }

    // Compute three normals per scripted vertex
    const std::vector<int>& scriptedverts = colliding_rod->getBoundaryCondition()->scriptedVertices();
    for (size_t i = 0; i < scriptedverts.size(); ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            normal.push_back(VecXd(ndof));
            normal.back().setZero();
            // Zero the velocity of this particular dof of this particular vertex
            normal.back()(colliding_rod->vertIdx(scriptedverts[i], j)) = 1.0;
        }
    }

    int numconstraints = (int) (normal.size());

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
        assert((int) scriptedverts[i] < colliding_rod->nv());
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
    Vec3d testn = normal[0].segment<3> (base0) + normal[0].segment<3> (base1);
    if (!rod0fixed)
        testn *= -1.0;
    Vec3d actln = NORMAL_SCALING * eecol.GetNormal();
    assert(approxEq(testn, actln, 1.0e-6));
#endif

    // Currently, all fixed vertex constraints normalized
#ifndef NDEBUG
    for (int i = 1; i < numconstraints; ++i)
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
                    << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. exertCompliantInelasticEdgeEdgeImpulseOneFixedThree. Time: "
                    << m_t << '\n';
    }

    // Ensure the edge degrees of freedom experience no impulse
#ifndef NDEBUG
    for (ElasticRod::edge_iter eit = colliding_rod->edges_begin(); eit != colliding_rod->edges_end(); ++eit)
    {
        for (int i = 0; i < numconstraints; ++i)
            assert(ntilde[i](colliding_rod->edgeIdx(*eit)) == 0.0);
    }
#endif

    // Vectors restriced to just verted DoFs. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnnormal(numconstraints);
    extractVertexDOF(posnnormal, normal, colliding_rod);
    std::vector<VecXd> posnntilde(numconstraints);
    extractVertexDOF(posnntilde, ntilde, colliding_rod);

    bool success = changeVelocityOneFree(posnnormal, posnntilde, scriptedverts, desired_values, numconstraints, rodbase, nvdof,
            eecol);

    //  applyInextensibilityVelocityFilter(rodidx);

    return success;
}

bool BARodStepper::changeVelocityOneFree(const std::vector<VecXd>& posnnormal, const std::vector<VecXd>& posnntilde,
        const std::vector<int>& scriptedverts, const VecXd& desired_values, int numconstraints, int rodbase, int nvdof,
        const CTCollision& ctcol)
{
    // Compute the Lagrange multiplier alpha for an asymmetric (e.g. rod vs. mesh) collision.
    MatXd lglhs(numconstraints, numconstraints);
#ifndef NDEBUG
    lglhs.setConstant(std::numeric_limits<double>::signaling_NaN()); // TODO: This matrix is symmetric, exploit that fact to avoid computations
#endif
    for (int i = 0; i < numconstraints; ++i)
        for (int j = 0; j < numconstraints; ++j)
            lglhs(i, j) = posnnormal[i].dot(posnntilde[j]);
    assert(approxSymmetric(lglhs, MATRIX_ASYMMETRY));

#ifndef NDEBUG
    Eigen::EigenSolver<MatXd> es;
    es.compute(lglhs, false); // compute the eigenvalues, not the eigenvectors
    TraceStream(g_log, "") << "Compliance matrix eigenvalues: " << es.eigenvalues() << '\n';

    // TODO: Test the conditioning number here
#endif

    VecXd lgrhs(desired_values);
    assert(lgrhs(0) == 0.0);
    assert(ctcol.computeRelativeVelocity() < 0);
    lgrhs(0) -= ctcol.computeRelativeVelocity() * NORMAL_SCALING;
    for (int i = 1; i < numconstraints; ++i)
        lgrhs(i) -= posnnormal[i].dot(m_vnphalf.segment(rodbase, nvdof));

    VecXd alpha(lglhs.lu().solve(lgrhs));

    // Contact constraint should 'push not pull'
    assert(alpha(0) >= 0.0);

    for (int i = 0; i < numconstraints; ++i)
        m_vnphalf.segment(rodbase, nvdof) += alpha(i) * posnntilde[i];

#ifndef NDEBUG
    // Ensure that the impulse eliminated the relative velocity
    assert(fabs(ctcol.computeRelativeVelocity()) < 1.0e-8);

    // Ensure the 'scripted vertices' achieved the desired velocity
    int curdof = 1;
    for (std::vector<int>::const_iterator vertex = scriptedverts.begin(); vertex != scriptedverts.end(); vertex++, curdof += 3)
        assert((m_vnphalf.segment<3> (rodbase + 3 * (*vertex)) - desired_values.segment<3> (curdof)).norm() < 1.0e-8);
#endif

    return true;
}

bool BARodStepper::exertCompliantInelasticEdgeEdgeImpulseBothFree(const EdgeEdgeCTCollision& eecol)
{
    assert(eecol.isAnalysed());
    // Ensure both edges have two free vertices
    //assert( isEntireEdgeFree(eecol.e0_v0,eecol.e0_v1) );
    //assert( isEntireEdgeFree(eecol.e1_v0,eecol.e1_v1) );

    assert(!YAEdge(eecol.e0_v0, eecol.e0_v1).IsFixed(m_geodata) && !YAEdge(eecol.e1_v0, eecol.e1_v1).IsFixed(m_geodata));

    // TraceStream(g_log, "") << "BARodStepper::exertCompliantInelasticEdgeEdgeImpulseBothFree: (x[" << eecol.e0_v0 << "]="
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
    n0.back().segment<3> (base_i0) = -(1.0 - eecol.s) * eecol.GetNormal() * NORMAL_SCALING;
    n0.back().segment<3> (base_i1) = -eecol.s * eecol.GetNormal() * NORMAL_SCALING;
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
    Vec3d actln0 = -NORMAL_SCALING * eecol.GetNormal();
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
    n1.back().segment<3> (base_j0) = (1.0 - eecol.t) * eecol.GetNormal() * NORMAL_SCALING;
    n1.back().segment<3> (base_j1) = eecol.t * eecol.GetNormal() * NORMAL_SCALING;
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
    Vec3d actln1 = NORMAL_SCALING * eecol.GetNormal();
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
    computeCompliantLHS(lhs0, m_rods[rod0]);
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
    computeCompliantLHS(lhs1, m_rods[rod1]);
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

    // Vectors restricted to just vertex DoFs for rod 0. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn0(nc0);
    extractVertexDOF(posnn0, n0, m_rods[rod0]);
    std::vector<VecXd> posnntilde0(nc0);
    extractVertexDOF(posnntilde0, ntilde0, m_rods[rod0]);

    // Vectors restricted to just vertex DoFs for rod 1. TODO: do all of this in place in the vectors I already have
    std::vector<VecXd> posnn1(nc1);
    extractVertexDOF(posnn1, n1, m_rods[rod1]);
    std::vector<VecXd> posnntilde1(nc1);
    extractVertexDOF(posnntilde1, ntilde1, m_rods[rod1]);

    bool success = changeVelocityBothFree(posnn0, posnn1, posnntilde0, posnntilde1, cval0, cval1, nc0, nc1, rod0base, rod1base,
            rod0nvdof, rod1nvdof, eecol);

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

    return success;
}

bool BARodStepper::changeVelocityBothFree(const std::vector<VecXd>& posnnormal0, const std::vector<VecXd>& posnnormal1,
        const std::vector<VecXd>& posnntilde0, const std::vector<VecXd>& posnntilde1, const VecXd& cval0, const VecXd& cval1,
        int nc0, int nc1, int rod0base, int rod1base, int rod0nvdof, int rod1nvdof, const EdgeEdgeCTCollision& eecol)
{
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
    lglhs(0, 0) = posnnormal0[0].dot(posnntilde0[0]) + posnnormal1[0].dot(posnntilde1[0]);
    // Entires (0,1::numalpha)
    int curel = 1;
    for (int i = 1; i < nc0; ++i)
        lglhs(0, curel++) = posnnormal0[0].dot(posnntilde0[i]);
    // Entries (0,numalpha::numalpha+numbeta)
    for (int i = 1; i < nc1; ++i)
        lglhs(0, curel++) = posnnormal1[0].dot(posnntilde1[i]);

    // Entires (1::numalpha,0)
    curel = 1;
    for (int i = 1; i < nc0; ++i)
        lglhs(curel++, 0) = posnntilde0[0].dot(posnnormal0[i]);
    // Entries (numalpha::numalpha+numbeta,0)
    for (int i = 1; i < nc1; ++i)
        lglhs(curel++, 0) = posnntilde1[0].dot(posnnormal1[i]);

    // Entries (1::numalpha,1::numalpha)
    for (int i = 1; i < nc0; ++i)
        for (int j = 1; j < nc0; ++j)
            lglhs(i, j) = posnnormal0[i].dot(posnntilde0[j]);
    // Entries (numalpha::numalpha+numbeta,numalpha::numalpha+numbeta)
    for (int i = 1; i < nc1; ++i)
        for (int j = 1; j < nc1; ++j)
            lglhs(nc0 + i - 1, nc0 + j - 1) = posnnormal1[i].dot(posnntilde1[j]);

    assert(approxSymmetric(lglhs, MATRIX_ASYMMETRY));

#ifndef NDEBUG
    Eigen::EigenSolver<MatXd> es;
    es.compute(lglhs, false); // compute the eigenvalues, not the eigenvectors
    DebugStream(g_log, "") << "Compliance matrix eigenvalues: " << es.eigenvalues() << '\n';

    // TODO: Test the conditioning number here
#endif

    Eigen::VectorXd lgrhs(numalpha);
#ifndef NDEBUG
    lgrhs.setConstant(std::numeric_limits<double>::signaling_NaN());
#endif
    // Entry 0
    assert(posnnormal0[0].size() == m_vnphalf.segment(rod0base, rod0nvdof).size());
    assert(posnnormal1[0].size() == m_vnphalf.segment(rod1base, rod1nvdof).size());
    lgrhs(0) = (posnnormal0[0].dot(m_vnphalf.segment(rod0base, rod0nvdof)) + posnnormal1[0].dot(
            m_vnphalf.segment(rod1base, rod1nvdof)));
    assert(approxEq(lgrhs(0), eecol.computeRelativeVelocity() * NORMAL_SCALING, 1.0e-6));
    assert(lgrhs(0) < 0.0);

    // Entries 1...numalpha
    for (int i = 1; i < nc0; ++i)
    {
        assert(posnnormal0[i].size() == m_vnphalf.segment(rod0base, rod0nvdof).size());
        lgrhs(i) = posnnormal0[i].dot(m_vnphalf.segment(rod0base, rod0nvdof));
    }
    // Entries numalpha...numalpha+numbeta
    for (int i = 1; i < nc1; ++i)
    {
        assert(posnnormal1[i].size() == m_vnphalf.segment(rod1base, rod1nvdof).size());
        lgrhs(nc0 + i - 1) = posnnormal1[i].dot(m_vnphalf.segment(rod1base, rod1nvdof));
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

}

void BARodStepper::extractVertexDOF(std::vector<VecXd>& posnnormal, const std::vector<VecXd>& normal,
        const ElasticRod* const rod)
{
    assert(rod != NULL);

    for (std::vector<VecXd>::iterator constraint = posnnormal.begin(); constraint != posnnormal.end(); constraint++)
    {
        *constraint = VecXd(3 * rod->nv());
        int theidx = 0;
        for (ElasticRod::vertex_iter vertex = rod->vertices_begin(); vertex != rod->vertices_end(); ++vertex, theidx += 3)
            constraint->segment<3> (theidx) = normal[constraint - posnnormal.begin()].segment<3> (rod->vertIdx(*vertex, 0));
    }
}

void BARodStepper::computeCompliantLHS(MatrixBase* lhs, ElasticRod* const rod)
{
    assert(lhs != NULL);
    assert(rod != NULL);

    //double emphasizeStaticEquilibium = 0;
    //InfoStream(g_log,"") << "BARodStepper::computeCompliantLHS: emphasizing static equilibrium = " << emphasizeStaticEquilibium;

    // lhs = -h^2*dF/dx
    lhs->setZero();
    //TraceStream(g_log, "") << "WARNING: COMPLIANCE is disabled!" << '\n';
    //colliding_rod->computeJacobian(0, -m_dt * m_dt - emphasizeStaticEquilibium, *lhs);
    rod->computeJacobian(0, -m_dt * m_dt, *lhs);
    lhs->finalize();

    // lhs = M - h^2*dF/dx
    for (int i = 0; i < rod->ndof(); ++i)
        lhs->add(i, i, rod->getMass(i));
    lhs->finalize();

    // Set rows/cols of twist DOFs to identity
    // TODO: Handle this better
    EPropHandle<int> edgeIdx;
    rod->property_handle(edgeIdx, "edge index");
    std::vector<int> dofs_to_clear;
    for (ElasticRod::edge_iter eit = rod->edges_begin(); eit != rod->edges_end(); ++eit)
    {
        int twist_dof = rod->property(edgeIdx)[*eit];
        dofs_to_clear.push_back(twist_dof);
    }

    lhs->zeroRows(dofs_to_clear, 1.0);
    lhs->finalize();
    lhs->zeroCols(dofs_to_clear, 1.0);
    lhs->finalize();

    assert(isSymmetric(*lhs));
}

/**
 * Post-collision checks and hacks
 */
bool BARodStepper::checkExplosions(std::vector<bool>& exploding_rods, const std::vector<bool>& failed_collisions_rods,
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
            if (hadExplosion(*rod))
            {
                explosions_detected = true;
                exploding_rods[*rod] = true;
            }
        }
    }
    if (explosions_detected)
        DebugStream(g_log, "") << "Some rods had explosions\n";

    return explosions_detected;
}

bool BARodStepper::hadExplosion(int rodIdx) const // All forces must have been computed
{
    TraceStream(g_log, "") << "Rod " << rodIdx << " Force norms: initial: " << (*(m_startForces[rodIdx])).norm()
            << " pre-dynamic: " << (*(m_preDynamicForces[rodIdx])).norm() << " pre-collisions: "
            << (*(m_preCollisionForces[rodIdx])).norm() << " post-collisions: " << (*(m_endForces[rodIdx])).norm() << "\n";

    for (int j = 0; j < m_rods[rodIdx]->ndof(); ++j)
    {
        const double s = (*(m_startForces[rodIdx]))[j];
        const double p = (*(m_preCollisionForces[rodIdx]))[j];
        const double e = (*(m_endForces[rodIdx]))[j];
        const double rate = fabs(e - p) / (fabs(p) + m_perf_param.m_explosion_damping);
        if (isnan(rate) || rate > m_perf_param.m_explosion_threshold)
        {
            TraceStream(g_log, "") << "Rod " << rodIdx << " had an explosion on DOF " << j << ": s = " << s << " p = " << p
                    << " e = " << e << " rate = " << rate << " \n";
            return true;
        }
    }
    return false;
}

bool BARodStepper::checkLengths(std::vector<bool>& stretching_rods)
{
    bool stretching_detected = false;
    DebugStream(g_log, "") << "Checking for lengths: " << m_simulated_rods.size() << " simulated rods\n";

    for (RodSelectionType::iterator rod = m_simulated_rods.begin(); rod != m_simulated_rods.end(); rod++)
    {
        if (stretching_rods[*rod] = !checkLength(*rod)) // Not a mistake
        {
            stretching_detected = true;
            setRodImmunity(*rod, true);
        }
    }

    if (stretching_detected)
        DebugStream(g_log, "") << "Some rods were stretched by a factor > " << m_perf_param.m_stretching_threshold << '\n';

    return stretching_detected;
}

bool BARodStepper::checkLength(const int rodIdx) const
{
    double length = 0.0;
    for (int j = 1; j < m_rods[rodIdx]->nv(); j++)
        length += (m_rods[rodIdx]->getVertex(j) - m_rods[rodIdx]->getVertex(j - 1)).norm();

    if (length > m_initialLengths[rodIdx] * m_perf_param.m_stretching_threshold)
    {
        TraceStream(g_log, "") << "Rod number " << rodIdx << " was stretched by a factor " << length / m_initialLengths[rodIdx]
                << '\n';
        return false;
    }

    return true;
}

void BARodStepper::applyInextensibilityVelocityFilter(int rodidx)
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
void BARodStepper::setupPenaltyForces(std::list<Collision*>& collisions, const RodSelectionType& selected_rods)
{
    // Detect proximity collisions
    m_collision_detector->getCollisions(collisions, Proximity);

    for (int rod_id = 0; rod_id < m_number_of_rods; rod_id++)
        assert(m_implicit_pnlty_forces[rod_id]->cleared());

    // Store the proximity collision in the RodPenaltyForce
    for (std::list<Collision*>::const_iterator col = collisions.begin(); col != collisions.end(); col++)
    {
        VertexFaceProximityCollision* vfpcol = dynamic_cast<VertexFaceProximityCollision*> (*col);
        if (vfpcol)
        {
            assert(vfpcol->isAnalysed());
            int rod_id = getContainingRod(vfpcol->v0);
            TraceStream(g_log, "") << "Creating penalty force for rod " << rod_id << " address " << &m_rods[rod_id] << '\n';
            int v_id = vfpcol->v0 - m_base_dof_indices[rod_id] / 3;
            m_implicit_pnlty_forces[rod_id]->registerProximityCollision(v_id, vfpcol);
        }
        // TODO: delete vfpcol once used. Not here though, as vfpcol is used each time RodPenaltyForce::computeForce is called
    }

}

void BARodStepper::setImplicitPenaltyThickness(const double& h)
{
    assert(h >= 0);
    m_perf_param.m_implicit_thickness = h;
}

void BARodStepper::setImplicitPenaltyStiffness(const double& k)
{
    assert(k >= 0.0);
    m_perf_param.m_implicit_stiffness = k;
}

}

