#ifndef WMSWEENEYRODMANAGER_HH_
#define WMSWEENEYRODMANAGER_HH_

#ifdef WETA
#include <weta/Wfigaro/Core/EigenIncludes.hh>
#include <weta/Wfigaro/Physics/World.hh>
#include <weta/Wfigaro/Core/ObjectControllerBase.hh>
#include <weta/Wfigaro/Render/RodRenderer.hh>
#include <weta/Wfigaro/Core/TriangleMesh.hh>
#include <weta/Wfigaro/Core/ScriptingController.hh>
#include <weta/Wfigaro/Collisions/LevelSet.hh>
#include <weta/Wfigaro/Physics/ElasticRods/BAGroomingStepper.hh>
#include <weta/Wfigaro/Physics/ElasticRods/ElasticRod.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodUtils.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodMayaForces.hh>
#include <weta/Wfigaro/Physics/World.hh>
#include <weta/Wfigaro/Physics/ElasticRods/PerformanceTuningParameters.hh>
#else
#include <BASim/src/Core/EigenIncludes.hh>
#include <BASim/src/Core/ObjectControllerBase.hh>
#include <BASim/src/Physics/World.hh>
#include <BASim/src/Render/RodRenderer.hh>
#include <BASim/src/Physics/ElasticRods/BAGroomingStepper.hh>
#endif

#include "../WmFigMeshController.hh"

#include <tr1/unordered_map>

/** \class WmSweeneyRodManager
 * \brief A class to create and control a simulation involving one of more rods.
 * 
 * This class has methods to allow the creation of a simulation world, populate it with rods
 * and then control the rods with various forces.
 */
class WmSweeneyRodManager
{
public:
    WmSweeneyRodManager();
    ~WmSweeneyRodManager();
    
    bool addRod( const std::vector< BASim::Vec3d >& i_vertices, 
                 const double i_time, const double i_youngsModulus =  1e7 * 1000.0,
                 const double i_shearModulus = 1e7 * 340.0, const double i_viscosity = 10.0, 
                 const double i_density = 1.3, const double i_radiusA = 1e-1 * 0.05, 
                 const double i_radiusB = 1e-1 * 0.05, 
                 const BASim::ElasticRod::RefFrameType i_referenceFrame = BASim::ElasticRod::TimeParallel,
                 const double i_massDamping = 10.0, 
                 const BASim::Vec3d i_gravity = BASim::Vec3d( 0.0, -980.0, 0.0 ),                 
                 const BASim::GroomingTimeStepper::Method i_solverType = BASim::GroomingTimeStepper::IMPL_EULER );

    void setUseKineticDamping ( bool i_useKinecticDamping );

    void initialiseSimulation( const double i_timeStep, const double i_startTime, BASim::PerformanceTuningParameters perfParams,
                               double i_atol, double i_stol, double i_rtol, double i_inftol, int i_numLineSearchIters );

    void updateSolverSettings( double i_atol, double i_stol, double i_rtol, double i_inftol, int i_numLineSearchIters, double i_penaltyStiffness);


    void addCollisionMesh( BASim::TriangleMesh* i_triangleMesh, 
                           BASim::LevelSet* i_levelSet,
                           WmFigMeshController* i_figMeshController );

    void takeStep();
    
    void drawAllRods();
    
    size_t numberOfRods()
    {
        return m_rods.size();
    }
    
    BASim::ElasticRod* rod( const size_t i_rodIndex )
    {
        return m_rods[ i_rodIndex ];
    }

    std::vector< BASim::ElasticRod* > m_rods;
  
private:
    BASim::BAGroomingStepper* m_bridsonStepper;
    std::vector< BASim::GroomingTimeStepper* > m_rodTimeSteppers;
    std::vector< BASim::RodRenderer* > m_rodRenderers;
    std::vector< BASim::TriangleMesh* > m_triangleMeshes;
    std::vector< BASim::LevelSet* > m_levelSets;
    std::vector< BASim::ScriptingController* > m_scriptingControllers;    
};

#endif
