/**
 * \file BridsonStepper.hh
 *
 * \author smith@cs.columbia.edu
 * \date 02/16/2010
 */

// This class is currently a little hacked together, and
// does some really inefficient stuff (copies all DOFs 
// into flat arrays and such) just because it was
// easier than using BASim iterators, for now :).

// Note also that the detection acceleartion is as basic
// as possible (AABB Hierarchy built with no knowledge of
// connectivity). Lots of low hanging fruit for acceleration.

#ifndef BRIDSONSTEPPER_HH
#define BRIDSONSTEPPER_HH

#ifdef WETA
#include "../../Core/ScriptingController.hh"
#include "ElasticRod.hh"
#include "RodTimeStepper.hh"
#include "RodMassDamping.hh"
#include "RodGravity.hh"
#include "../../Math/Math.hh"
#include "../../Core/TriangleMesh.hh"
//#include "../../Collisions/BVHAABB.hh"
#include "../../Collisions/CollisionDetector.hh"
#include "../../Collisions/CollisionUtils.hh"
//#include "../../../../Apps/BASimulator/Problems/ProblemBase.hh"
#include "../../Core/StatTracker.hh"
#include "MinimalTriangleMeshBackup.hh"
#include "RodPenaltyForce.hh"
#else
#include "BASim/src/Core/ScriptingController.hh"
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"
#include "BASim/src/Physics/ElasticRods/RodTimeStepper.hh"
#include "BASim/src/Physics/ElasticRods/RodMassDamping.hh"
#include "BASim/src/Physics/ElasticRods/RodGravity.hh"
#include "BASim/src/Math/Math.hh"

#include "BASim/src/Core/TriangleMesh.hh"

//#include "BASim/src/Collisions/BVHAABB.hh"
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

// Linear solver and matrix that will work together
class LinearSystemSolver
{
public:

    LinearSystemSolver(MatrixBase* lhs, LinearSolverBase* solver) :
        m_sys_size(lhs->rows()), m_lhs(lhs), m_solver(solver)
    {
        assert(lhs->rows() == lhs->cols());
        assert(lhs->rows() == m_sys_size);
    }

    ~LinearSystemSolver()
    {
        if (m_lhs != NULL)
        {
            delete m_lhs;
            m_lhs = NULL;
        }
        if (m_solver != NULL)
        {
            delete m_solver;
            m_solver = NULL;
        }
    }

    int m_sys_size;
    MatrixBase* m_lhs;
    LinearSolverBase* m_solver;
};

// Collection of linear system solvers of different sizes
//  NOT THREAD SAFE
//  ASSUMES YOU WANT A BAND MATRIX OF SIZE 10 :)
class LinearSystemSolverCollection
{
public:

    ~LinearSystemSolverCollection()
    {
        std::map<int, LinearSystemSolver*>::iterator it = m_solver_map.begin();
        for (; it != m_solver_map.end(); ++it)
        {
            assert(it->second != NULL);
            delete it->second;
            it->second = NULL;
        }
    }

    LinearSystemSolver* getLinearSystemSolver(int size)
    {
        assert(size > 0);

        // Attempt to locate a solver of the requested size
        std::map<int, LinearSystemSolver*>::iterator it = m_solver_map.find(size);
        // If a solver of the size exists, return it
        if (it != m_solver_map.end())
        {
            return it->second;
        }
        // Otherwise create a new solver
        int band = 10;
        MatrixBase* lhs = SolverUtils::instance()->createBandMatrix(size, size, band, band);
        LinearSolverBase* solver = SolverUtils::instance()->createLinearSolver(lhs);
        m_solver_map.insert(std::pair<int, LinearSystemSolver*>(size, new LinearSystemSolver(lhs, solver)));

        it = m_solver_map.find(size);
        return it->second;
    }

private:
    std::map<int, LinearSystemSolver*> m_solver_map;
};

/**
 * Class to evolve a collection of rods forward in time, resolving collisions using
 * a "velocity filter" in the spirit of Bridson's 2002 paper "Robust Treatment of
 * Collisions, Contact, and Friction for Cloth Animation."
 */
class BridsonStepper: public ObjectControllerBase
{

public:
    /**
     * Default constructor.
     */
    //BridsonStepper();

    /**
     * Creates a BridsonStepper with user-supplied options.
     *
     * \param[in] intgrtr Integrator (class RodTimeStepper) to use. Assumes implicit euler, for now.
     * \param[in] max_implct_itrtns If an implicit integrator is selected, the maximum iterations allowed.
     * \param[in] dt Timestep to use.
     * \param[in] mass_dmpng Amount of damping that acts in opposition to vertex velocities (I think? Miklos has been mucking with the damping :)).
     * \param[in] grav Three dimensional vector that specifies gravity.
     */
    // Parameter num_threads = -1 will cause the number of threads to be set equal to the number of available processors.
    BridsonStepper(std::vector<ElasticRod*>& rods, std::vector<TriangleMesh*>& trimeshes,
            std::vector<ScriptingController*>& scripting_controllers, std::vector<RodTimeStepper*>& steppers, const double& dt,
            const double time = 0.0, int num_threads = -1);

    /**
     * Destructor.
     */
    virtual ~BridsonStepper();

private:
    /**
     * After adding new rods or objects, this method must be called.
     */
    void prepareForExecution();

public:

    double computeMaxEdgeAngle(const ElasticRod& rod)
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
     * Adds a rod that will be evolved in time using this BridsonStepper.
     */
    //void addRod( ElasticRod* rod );

    /**
     * Adds a non-simulated triangle mesh for objects to collide with.
     */
    //void addTriangleMesh( TriangleMesh* tri_mesh );

    /**
     * Evolves all inserted rods forward in time.
     */
    bool execute();

    /**
     * Modifies the timestep.
     */
    void setDt(double dt);

    /**
     * Disables all response.
     */
    void disableResponse();

    /**
     * Enables response, subject to state of enableIterativeInelasticImpulses().
     */
    void enableReseponse();

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

    /**
     *  Enable or disable self collisions between all rods
     */
    void skipRodRodCollisions( bool skipRodRodCollisions )
    {
        m_skipRodRodCollisions = true;
        
        if ( m_collision_detector )
            m_collision_detector->skipRodRodCollisions( skipRodRodCollisions );
    }
        
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

    // TODO: Move these to some kind of automated test suite
    //void testCoplanarityTime();

    /**
     * Options names for rods used in output.
     */
    void setRodLabels(const std::vector<std::string>& rod_labels);

    double computeTotalForceNorm();
    bool step(bool check_explosion);
    bool nonAdaptiveExecute(double dt);
    bool adaptiveExecute(double dt);

private:
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

    void extractPositions(const std::vector<ElasticRod*>& rods, const std::vector<int>& base_indices, VecXd& positions);
    void extractVelocities(const std::vector<ElasticRod*>& rods, const std::vector<int>& base_indices, VecXd& velocities);
    //void extractMasses( const std::vector<ElasticRod*>& rods, std::vector<double>& masses ) const;
    //void extractFixedVertices( const std::vector<ElasticRod*>& rods, std::vector<bool>& fixed ) const;
    //void extractEdges( const std::vector<ElasticRod*>& rods, const std::vector<int>& base_indices, std::vector<std::pair<int,int> >& edges, std::vector<double>& radii );

    void restorePositions(std::vector<ElasticRod*>& rods, const VecXd& positions);
    void restoreVelocities(std::vector<ElasticRod*>& rods, const VecXd& velocities);
    void restoreResponses(std::vector<ElasticRod*>& rods, const VecXd& responses);

    bool isRodVertex(int vert) const;
    bool isRodRodCollision(const EdgeEdgeCTCollision& collision) const;

    int getContainingRod(int vert_idx) const;

    /////////////////////////////////////////////////////
    // Collision detection routines

    /*
    bool isVertexFixed(int vert_idx) const;
    bool isEntireFaceFixed(int v0, int v1, int v2) const;
    bool isEntireEdgeFree(int v0, int v1) const;
    bool isEntireEdgeFixed(int v0, int v1) const;
    bool isOneVertexFixed(int v0, int v1) const;

    // Determines if two edges share a vertex
    bool edgesShareVertex(const std::pair<int, int>& edgei, const std::pair<int, int>& edgej) const;
    bool edgesSharevertex(const int& e0v0, const int& e0v1, const int& e1v0, const int& e1v1) const;
*/

    // Determines if a vertex and a face share a vertex
    bool vertexAndFaceShareVertex(const int& vertex, const int& face) const;
    bool vertexAndFaceShareVertex(const int& v, const int& f0, const int& f1, const int& f2) const;
    // Generates a list of ALL possible edge-edge collisions
    void generateAllEdgeEdgeProximityCollisionPairs(std::vector<EdgeEdgeProximityCollision>& edge_edge_collisions) const;
    void
    generateAllEdgeEdgeContinuousTimeCollisionPairs(std::vector<EdgeEdgeCTCollision>& edge_edge_collisions) const;

    // Generates a list of ALL possible vertex-face collisions
    void generateAllVertexFaceProximityCollisionPairs(std::vector<VertexFaceProximityCollision>& vertex_face_collisions) const;
    void generateAllVertexFaceContinuousTimeCollisionPairs(std::vector<VertexFaceCTCollision>& vertex_face_collisions) const;

    // Returns a list of all edges that are in "proximity" for the positions specified in the vector x.
    //void detectEdgeEdgeProximityCollisions( const VecXd& x, std::vector<EdgeEdgeProximityCollision>& edge_edge_collisions ) const;
    void detectEdgeEdgeProximityCollisions(const VecXd& x, std::vector<EdgeEdgeProximityCollision>& pssbl_cllsns,
            std::vector<EdgeEdgeProximityCollision>& cllsns) const;

    // Returns a list of all vertex-face pairs that are in "proximity" for the positions specified in the vector x.
    void detectVertexFaceProximityCollisions(const VecXd& x, std::vector<VertexFaceProximityCollision>& pssbl_cllsns,
            std::vector<VertexFaceProximityCollision>& vetex_face_collisions) const;

    bool isProperCollisionTime(double time);

    // Computes the impulse necessary to eliminate all relative velocity at given points on two edges
    Vec3d computeEdgeEdgeInelasticImpulse(const double& ma0, const double& ma1, const double& mb0, const double& mb1,
            const double& s, const double& t, const double& relvel, const Vec3d& n);

    // Computes the impulse necessary to eliminate all relative velocity at a given point on a face and a vertex
    Vec3d computeVertexFaceInelasticImpulse(const double& mvrt, const double& mfc0, const double& mfc1, const double& mfc2,
            const double& u, const double& v, const double& w, const double& relvel, const Vec3d& n);

    /////////////////////////////////////////////////////
    // Collision response routines

    void executePenaltyResponse();
    bool executeIterativeInelasticImpulseResponse();
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

    //////////////////////////////////
    // Jungseock's penalty response
    void detectVertexFaceImplicitPenaltyCollisions(const VecXd& x, std::vector<VertexFaceProximityCollision>& pssbl_cllsns,
            std::vector<VertexFaceImplicitPenaltyCollision>& vetex_face_collisions) const;
    void executeImplicitPenaltyResponse();

public:
    void enableImplicitPenaltyImpulses();
    void disableImplicitPenaltyImpulses();
    void setImplicitPenaltyExtraThickness(const double& h);
    void setVertexFacePenalty(const double& k);

private:

    // Total number of degrees of freedom in the system
    int m_num_dof;
    // Vector of rods this BridsonStepper evolves in time
    std::vector<ElasticRod*>& m_rods;
    // Vector of ScriptedTriangleObjects in the system
    const std::vector<TriangleMesh*>& m_triangle_meshes;
    // Controllers to move scripted geometry/rods forward in time and to set boundary conditions
    const std::vector<ScriptingController*>& m_scripting_controllers;
    // Time steppers to evolve rods forward (ignoring collisions)
    const std::vector<RodTimeStepper*>& m_steppers;
    // Integrator selected by user
    RodTimeStepper::Method m_method;
    // Timestep selected by user
    double m_dt;

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

    std::vector<double> m_masses;

    const GeometricData m_geodata;

    VecXd m_xn;
    VecXd m_xnp1;
    VecXd m_vnphalf;
    VecXd m_vnresp;

    // Enable/Disable portions of the collision response
    bool m_respns_enbld;
    bool m_pnlty_enbld;
    bool m_itrv_inlstc_enbld;
    int m_num_inlstc_itrns;
    double m_vrt_fc_pnlty;

    // Some debug stuff in for now.
    bool m_nan_enc;
    bool m_inf_enc;
    bool m_lt0_enc;
    bool m_gt0_enc;

    //	BVHAABB* m_bvh;
    CollisionDetector* m_collision_detector;

    // Assuming rods are stored first (which they now are), the index that points one past the last rods data
    int m_obj_start;

    //Problem* m_problem;
    double m_t;

    std::vector<std::string> m_rod_labels;

    LinearSystemSolverCollection m_solver_collection;

    //////////////////////////////////
    // Jungseock's penalty response
    bool m_implicit_pnlty_enbld;
    double m_implicit_thickness;
    std::vector<RodPenaltyForce*> m_implicit_pnlty_forces;

    int m_num_threads;
    
    // Toggle self collisions on or off
    bool m_skipRodRodCollisions;

};

} // namespace BASim

#endif // RODTIMESTEPPER_HH
