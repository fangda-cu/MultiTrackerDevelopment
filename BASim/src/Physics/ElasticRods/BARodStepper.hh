/**
 * \file BARodStepper.hh
 *
 * \author smith@cs.columbia.edu
 * \date 02/16/2010
 */

#ifndef BARODSTEPPER_HH
#define BARODSTEPPER_HH

#ifdef WETA
#include "../../Core/ScriptingController.hh"
#include "ElasticRod.hh"
#include "RodTimeStepper.hh"
#include "RodMassDamping.hh"
#include "RodGravity.hh"
#include "../../Math/Math.hh"
#include "../../Math/LinearSystemSolver.hh"
#include "../../Core/TriangleMesh.hh"
//#include "../../Collisions/BVHAABB.hh"
#include "../../Collisions/CollisionDetector.hh"
#include "../../Collisions/CollisionUtils.hh"
//#include "../../../../Apps/BASimulator/Problems/ProblemBase.hh"
#include "../../Core/StatTracker.hh"
#include "MinimalTriangleMeshBackup.hh"
#include "RodPenaltyForce.hh"
#include "PerformanceTuningParameters.hh"
#else
#include "BASim/src/Core/ScriptingController.hh"
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"
#include "BASim/src/Physics/ElasticRods/RodTimeStepper.hh"
#include "BASim/src/Physics/ElasticRods/RodMassDamping.hh"
#include "BASim/src/Physics/ElasticRods/RodGravity.hh"
#include "BASim/src/Math/Math.hh"

#include "BASim/src/Core/TriangleMesh.hh"

#include "BASim/src/Collisions/CollisionUtils.hh"
#include "BASim/src/Collisions/CollisionDetector.hh"

#include "Apps/BASimulator/Problems/ProblemBase.hh"

#include "BASim/src/Core/StatTracker.hh"

#include "BASim/src/Physics/ElasticRods/MinimalTriangleMeshBackup.hh"
#endif

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <limits>
#include <list>
#include <queue>
#include <typeinfo>

namespace BASim
{

/**
 * Class to evolve a collection of rods forward in time, resolving collisions using
 * a "velocity filter" in the spirit of Bridson's 2002 paper "Robust Treatment of
 * Collisions, Contact, and Friction for Cloth Animation."
 */
class BARodStepper: public ObjectControllerBase
{
    typedef std::list<int> RodSelectionType;

public:
    /**
     * Creates a BARodStepper with user-supplied options.
     *
     * \param[in] intgrtr Integrator (class RodTimeStepper) to use. Assumes implicit euler, for now.
     * \param[in] max_implct_itrtns If an implicit integrator is selected, the maximum iterations allowed.
     * \param[in] dt Timestep to use.
     * \param[in] mass_dmpng Amount of damping that acts in opposition to vertex velocities (I think? Miklos has been mucking with the damping :)).
     * \param[in] grav Three dimensional vector that specifies gravity.
     */
    // Parameter num_threads = -1 will cause the number of threads to be set equal to the number of available processors.
    BARodStepper(std::vector<ElasticRod*>& rods, std::vector<TriangleMesh*>& trimeshes,
            std::vector<ScriptingController*>& scripting_controllers, std::vector<RodTimeStepper*>& steppers, const double& dt,
            const double time = 0.0, const int num_threads = -1,
            const PerformanceTuningParameters perf_param = PerformanceTuningParameters());

    /**
     * Destructor.
     */
    virtual ~BARodStepper();

    /**
     * Evolves all inserted rods forward in time.
     */
    bool execute();

    /**
     * Returns the timestep.
     **/
    double getDt() const;

    /**
     * Returns the simulation time
     **/
    double getTime() const;

    /**
     *  Enable or disable self collisions between all rods
     */
    void skipRodRodCollisions(bool skipRodRodCollisions)
    {
        std::cerr << "Switching rod-rod collisions " << (skipRodRodCollisions ? "OFF" : "ON") << std::endl;
        m_perf_param.m_skipRodRodCollisions = skipRodRodCollisions;

        if (m_collision_detector)
            m_collision_detector->setSkipRodRodCollisions(skipRodRodCollisions);
    }

    void setStopOnRodError(bool stopOnRodError)
    {
        if (!m_stopOnRodError && stopOnRodError)
        {
            std::cerr << "BARodStepper::m_stopOnError changed to \033[33mtrue\033[0m" << std::endl;
            // If we change from non-stopping to stopping, reset m_simulationFailed so we take only future errors into account.
            m_simulationFailed = false;
        }
        if (m_stopOnRodError && !stopOnRodError)
            std::cerr << "BARodStepper::m_stopOnError changed to \033[33mfalse\033[0m" << std::endl;

        m_stopOnRodError = stopOnRodError;
    }

private:
    /**
     * Modifies the timestep.
     */
    void setDt(double dt);

    /**
     * After adding new rods or objects, this method must be called.
     */
    void prepareForExecution();

    double computeMaxEdgeAngle(const ElasticRod& rod) const
    {
        double maxangle = -std::numeric_limits<double>::infinity();
        for (int i = 0; i < rod.ne() - 1; ++i)
        {
            Vec3d edge0 = rod.getEdge(i);
            Vec3d edge1 = rod.getEdge(i + 1);
            double numer = edge0.cross(edge1).norm();
            double denom = edge0.dot(edge1);
            double angle = atan2(numer, denom);
            if (angle < 0.0)
                std::cout << "NEGATIVE ANGLE AGHHHHHH" << std::endl;
            if (angle > maxangle)
                maxangle = angle;
        }
        return maxangle;
    }

    /**
     * Adds a rod that will be evolved in time using this BARodStepper.
     */
    //void addRod( ElasticRod* rod );

    /**
     * Adds a non-simulated triangle mesh for objects to collide with.
     */
    //void addTriangleMesh( TriangleMesh* tri_mesh );

    /**
     * Disables all response.
     */
    void disableResponse();

    /**
     * Enables response, subject to state of enableIterativeInelasticImpulses().
     */
    void enableResponse();

    /**
     * Enables penalty response.
     */
    void enablePenaltyImpulses();

    /**
     * Disables penalty response.
     */
    void disablePenaltyImpulses();

    /**
     * Enables iterative impulse response.
     */
    void enableIterativeInelasticImpulses();

    /**
     * Disables iterative impulse response.
     */
    void disableIterativeInelasticImpulses();

    void computeImmunity(const RodSelectionType& selected_rods);

    /**
     * Sets the maximum number of inelastic impulses to apply iterativly.
     */
    void setNumInelasticIterations(const int& num_itr);

    /**
     * Number of rods this controller is responsible for.
     */
    int getNumRods() const
    {
        return (int) (m_rods.size());
    }

    /**
     * Number of triangle meshes this controller is responsible for.
     */
    int getNumTriangleMeshes() const
    {
        return (int) (m_triangle_meshes.size());
    }
    ;

    void setTime(double time);

    // TODO: Move these to some kind of automated test suite
    //void testCoplanarityTime();

    /**
     * Options names for rods used in output.
     */
    void setRodLabels(const std::vector<std::string>& rod_labels);

    double computeTotalForceNorm() const;
    void step(RodSelectionType& selected_rods);
    bool nonAdaptiveExecute(double dt, RodSelectionType& selected_rods);
    bool adaptiveExecute(double dt, RodSelectionType& selected_rods);

    /////////////////////////////////////////////////////
    // Methods for checking the sanity of input rods

    // Currently we do not support collisions for anisotropic cross sections
    void ensureCircularCrossSection(const ElasticRod& rod) const;

    // If the cross-sectional radius is too large and edge lengths too small,
    // non-adjacent portions of the rod will be in contact by default. We can
    // probably add some special case code to handle this later, but just
    // disallow the situation for now.
    void ensureNoCollisionsByDefault(const ElasticRod& rod) const;

    /////////////////////////////////////////////////////
    // Helper methods

    // Returns the total number of degrees of freedom in the system
    int getNumDof() const;

    // Returns the total number of vertices in the system
    int getNumVerts() const;

    void extractPositions(VecXd& positions, const RodSelectionType& selected_rods) const;
    void extractVelocities(VecXd& velocities, const RodSelectionType& selected_rods) const;

    void restorePositions(const VecXd& positions, const RodSelectionType& selected_rods);
    void restoreVelocities(const VecXd& velocities, const RodSelectionType& selected_rods);
    void restoreResponses(const VecXd& responses, const RodSelectionType& selected_rods);

    bool isRodVertex(int vert) const;
    bool isRodRodCollision(const EdgeEdgeCTCollision& collision) const;

    int getContainingRod(int vert_idx) const;

    // Determines if a vertex and a face share a vertex
    bool vertexAndFaceShareVertex(const int& vertex, const int& face) const;
    bool vertexAndFaceShareVertex(const int& v, const int& f0, const int& f1, const int& f2) const;
    bool isProperCollisionTime(double time);

    void applyInextensibilityVelocityFilter(int rodidx);

    /////////////////////////////////////////////////////
    // Collision response routines

    void executePenaltyResponse();
    bool executeIterativeInelasticImpulseResponse(std::vector<bool>& rods_failed_because_of_iterated_collisions);
    //	void filterCollisions(std::list<ContinuousTimeCollision>& cllsns);

    void exertPenaltyImpulses(std::vector<EdgeEdgeProximityCollision>& edg_edg_cllsns,
            std::vector<VertexFaceProximityCollision>& vrtx_fce_cllsns, VecXd& v);

    void exertInelasticImpluse(EdgeEdgeCTCollision& clssn);
    void exertInelasticImpluse(VertexFaceCTCollision& clssn);
    void exertInelasticImpulses(std::vector<EdgeEdgeCTCollision>& edg_edg_cllsns,
            std::vector<VertexFaceCTCollision>& vrtx_fce_cllsns, VecXd& v);

    void exertVertexImpulse(const Vec3d& I, const double& m, const int& idx, VecXd& v);
    void exertEdgeImpulse(const Vec3d& I, const double& m0, const double& m1, const double& alpha, const int& idx0,
            const int& idx1, VecXd& v);
    void exertFaceImpulse(const Vec3d& I, const double& m0, const double& m1, const double& m2, const double& u,
            const double& v, const double& w, const int& idx0, const int& idx1, const int& idx2, VecXd& vel);

    void computeCompliantLHS(MatrixBase* lhs, int rodidx);

    void exertCompliantInelasticImpulse(const CTCollision* cllsn);
    void exertCompliantInelasticVertexFaceImpulse(const VertexFaceCTCollision& vfcol);
    void exertCompliantInelasticEdgeEdgeImpulse(const EdgeEdgeCTCollision& eecol);
    void exertCompliantInelasticEdgeEdgeImpulseOneFixed(const EdgeEdgeCTCollision& eecol);
    void exertCompliantInelasticEdgeEdgeImpulseBothFree(const EdgeEdgeCTCollision& eecol);
    bool checkExplosions(std::vector<bool>& exploding_rods, const std::vector<bool>& failed_collisions_rods,
            const RodSelectionType& selected_rods);

    //////////////////////////////////
    // Jungseock's penalty response
    //  void detectVertexFaceImplicitPenaltyCollisions(const VecXd& x, std::vector<VertexFaceProximityCollision>& pssbl_cllsns,
    //          std::vector<VertexFaceImplicitPenaltyCollision>& vetex_face_collisions) const;
    void executeImplicitPenaltyResponse(std::list<Collision*>& collisions, const RodSelectionType& selected_rods);

    void enableImplicitPenaltyImpulses();
    void disableImplicitPenaltyImpulses();
    void setImplicitPenaltyExtraThickness(const double& h);
    void setVertexFacePenalty(const double& k);

    void killTheRod(int rod);
    void computeForces(VecXd ** Forces, const RodSelectionType& selected_rods);

    /*
     * Member variables
     */

    // Total number of degrees of freedom in the system
    int m_num_dof;
    // Time steppers to evolve rods forward (ignoring collisions)
    std::vector<RodTimeStepper*>& m_steppers;
#ifdef KEEP_ONLY_SOME_RODS
    // Vector of rods this BARodStepper evolves in time
    std::vector<ElasticRod*>& m_rods;
    size_t m_number_of_rods; // set to m_rods.size()
#else
    // Vector of rods this BARodStepper evolves in time
    const std::vector<ElasticRod*>& m_rods;
    const size_t m_number_of_rods; // set to m_rods.size()
#endif
    // Vector of ScriptedTriangleObjects in the system
    const std::vector<TriangleMesh*>& m_triangle_meshes;
    // Controllers to move scripted geometry/rods forward in time and to set boundary conditions
    const std::vector<ScriptingController*>& m_scripting_controllers;
    // Integrator selected by user
    // RodTimeStepper::Method m_method;
    // Timestep selected by user
    double m_dt;
    // Current level of substepping
    int m_level;

    // Entry i is base index of ith rod in global array of position dofs
    std::vector<int> m_base_indices;

    // Entry i is base index of ith ScriptedTriangleMesh in global array of position dofs
    std::vector<int> m_base_triangle_indices;

    // Vector of edges in the system. FREE (not part of a face) edges
    std::vector<std::pair<int, int> > m_edges;

    // Vector of triangular faces in the system
    std::vector<TriangularFace> m_faces;

    std::vector<double> m_vertex_radii;
    // TODO: Possibly get rid of these, just pull radii from m_vertex_radii.
    std::vector<double> m_edge_radii;
    std::vector<double> m_face_radii;
    // Vector of masses per vertex
    std::vector<double> m_masses;
    // Vector of booleans indicating whether a vertex should be considered for collisions
    std::vector<bool> m_collision_immune;
    // Structure of references to the geometry
    const GeometricData m_geodata;

    // Enable/Disable portions of the collision response
    bool m_respns_enbld;
    //    bool m_pnlty_enbld;
    //    bool m_itrv_inlstc_enbld;
    //    int m_num_inlstc_itrns;

    // Some debug stuff in for now.
    bool m_nan_enc;
    bool m_inf_enc;
    bool m_lt0_enc;
    bool m_gt0_enc;

    // Collision detector
    CollisionDetector* m_collision_detector;

    // Assuming rods are stored first (which they now are), the index that points one past the last rods data
    int m_obj_start;

    // Current time
    double m_t;

    std::vector<std::string> m_rod_labels;

    LinearSystemSolverCollection m_solver_collection;

    //////////////////////////////////
    // Jungseock's penalty response
    std::vector<RodPenaltyForce*> m_implicit_pnlty_forces;

    // Number of threads to be used for dynamics and collisions
    int m_num_threads;

    // Toggle selective adaptivity
    bool m_simulationFailed;
    // Flag indicating whether we should freeze the simulation on first failure.
    bool m_stopOnRodError;
    // Set of rods that will no longer be collided during a time step
    std::set<int> m_collision_disabled_rods;

    /**
     * Backup structures
     */
    // Positions
    VecXd m_xn;
    VecXd m_xnp1;
    VecXd m_xdebug;
    // Velocities
    VecXd m_vnphalf;
    VecXd m_vnresp;
    // Forces
    VecXd** m_startForces;
    VecXd** m_preCollisionForces;
    VecXd** m_endForces;
    // Global back-ups for adaptive substepping
    std::vector<MinimalRodStateBackup> m_rodbackups;
    std::vector<MinimalTriangleMeshBackup> m_objbackups;

    RodSelectionType m_simulated_rods;

    PerformanceTuningParameters m_perf_param;
};

} // namespace BASim

#endif // RODTIMESTEPPER_HH
