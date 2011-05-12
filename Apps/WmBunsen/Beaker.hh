/**
 * \file Beaker.hh
 *
 * \author acoull@wetafx.co.nz
 * \date 10/26/2009
 */

#ifndef BEAKER_HH_
#define BEAKER_HH_

#ifdef WETA
#include <weta/Wfigaro/Core/EigenIncludes.hh>
#include <weta/Wfigaro/Physics/World.hh>
//#include <weta/Wfigaro/Collisions/CollisionMeshData.hh>
#include <weta/Wfigaro/Core/ObjectControllerBase.hh>
//#include <weta/Wfigaro/Physics/ElasticRods/RodCollisionTimeStepper.hh>
#include <weta/Wfigaro/Render/RodRenderer.hh>
#include <weta/Wfigaro/Core/TriangleMesh.hh>
#include <weta/Wfigaro/Core/ScriptingController.hh>
#include <weta/Wfigaro/Physics/ElasticRods/BridsonStepper.hh>
#else
#include <BASim/src/Core/EigenIncludes.hh>
//#include <BASim/src/Collisions/CollisionMeshData.hh>
#include <BASim/src/Core/ObjectControllerBase.hh>
#include <BASim/src/Physics/World.hh>
//#include <BASim/src/Physics/ElasticRods/RodCollisionTimeStepper.hh>
#include <BASim/src/Render/RodRenderer.hh>
#include <BASim/src/Physics/ElasticRods/BridsonStepper.hh>
#endif

#include "RodData.hh"
#include <tr1/unordered_map>
#include <iostream>
#include <ext/hash_map>
#include <iostream>
#include <fstream>
#include "WmFigRodGroup.hh"
//#include <weta/Wfigaro/SceneXML/SceneXML.hh>

#undef USE_MPI

/*#include <Geometry/SEGMENTED_CURVE.h>
#include <Grids/GRID_3D.h>
#include <Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Matrices_And_Vectors/VECTOR.h>*/
//#include "VolumetricCollisions/VolumetricCollisionsCPU.hh"

using namespace BASim;
using namespace std;

typedef std::tr1::unordered_map< int, BASim::Vec3d > FixedVertexMap ;
typedef std::tr1::unordered_map< int, FixedVertexMap > FixedRodVertexMap ;

typedef std::tr1::unordered_map< int, BASim::Vec3d > LockedVertexMap ;
typedef std::tr1::unordered_map< int, LockedVertexMap > LockedRodVertexMap ;

//typedef std::tr1::unordered_map< int, TriangleMesh* > CollisionMeshHashMap;
//typedef std::tr1::unordered_map< int, ScriptedController* > ScriptedControllerHashMap;

class CollisionMeshData
{
public:
    CollisionMeshData( TriangleMesh* i_triangleMesh, ScriptingController* i_scriptingController )
    {
        m_triangleMesh = i_triangleMesh;
        m_scriptingController = i_scriptingController;
    }

    TriangleMesh* triangleMesh()
    {
        return m_triangleMesh;
    }

    ScriptingController* scriptingController()
    {
        return m_scriptingController;
    }

private:
    TriangleMesh* m_triangleMesh;
    ScriptingController* m_scriptingController;
};

typedef std::tr1::unordered_map< int, CollisionMeshData* > CollisionMeshDataHashMap;


class Beaker
{
public:
    Beaker();
    ~Beaker();

    void BaseAtEachTimestep();
    void BaseAfterLoad();
    
    // Accessor functions
    
    // TODO - remove world as it's not needed any more I don't think!
    
    const World& getWorld() const { return *m_world; }
    World& getWorld() { return *m_world; }
    
    Scalar getTime() const
    {
        //return m_world->property( m_timeHandle );
        // 
        return m_bridsonStepper->getTime();
    }
    
    Scalar getDt() const
    {
        //return m_world->property( m_dtHandle );
        return m_bridsonStepper->getDt();
    }

    void setDt( const Scalar& i_dt )
    {
        // Does setting this property on the world mean anything?
        m_world->property( m_dtHandle ) = i_dt;
        
        m_bridsonStepper->setDt( i_dt );
    }

    const BASim::Vec3d& getGravity() const
    {
      return m_world->property( m_gravityHandle );
    }
    
    void setGravity( const BASim::Vec3d& gravity )
    {
        m_gravity = gravity;
        m_world->property( m_gravityHandle ) = m_gravity;
    }

    const int& getMaxIter() const
    {
        return m_world->property( m_maxIterHandle );
    }

    void setMaxIter(const int& maxIter)
    {
        m_world->property( m_maxIterHandle ) = std::max(maxIter,1);
    }

    /*vector<RodData*>* rodData( int i_rodGroup )
    {
        return &( m_rodDataMap[ i_rodGroup ] );
    }*/

    void setPlasticDeformations( bool i_plasticDeformations )
    {
        m_plasticDeformations = i_plasticDeformations;
    }
    
    void shouldDrawSubsteppedVertices( bool i_shouldDrawSubsteppedVertices )
    {
        m_shouldDrawSubsteppedVertices = i_shouldDrawSubsteppedVertices;
    }

    void clumpingEnabled( bool i_isClumpingEnabled )
    {
        m_isClumpingEnabled = i_isClumpingEnabled;
    }

    void clumpingCoefficient( double i_clumpingCoefficient )
    {
        m_clumpingCoefficient = i_clumpingCoefficient;
    }

    void set_stol( double tols )
    {
        m_stol = tols;
    }

    void set_atol( double tola )
    {
        m_atol = tola;
    }

    void set_rtol( double tolr )
    {
        m_rtol = tolr;
    }

    void set_inftol( double tolinf )
    {
        m_inftol = tolinf;
    }

    ////////////////////////////////////////////////////////////////////////////////////////
    //
    // Volumetric collision stuff
    //
    ////////////////////////////////////////////////////////////////////////////////////////

    void setFlip( double i_flip )
    {
        m_flip = i_flip;
    }
    
    void setSlip( double i_slip )
    {
        m_slip = i_slip;
    }
    
    void setDoVolumetricCollisions( bool i_doVolumetricCollisions )
    {
        m_doVolumetricCollisions = i_doVolumetricCollisions;
    }

    void setTargetEdgeDensity( double i_targetEdgeDensity )
    {
        m_targetEdgeDensity = i_targetEdgeDensity;
    }

    void setVolumetricRadius( double i_volumetricRadius )
    {
        m_volumetricRadius = i_volumetricRadius;
    }

    void setGridDX( double i_gridDX )
    {
        m_gridDX = i_gridDX;
    }

    void setSeparationCondition( double i_separationConditionX, double i_separationConditionY, double i_separationConditionZ )
    {
        m_separationCondition[ 0 ] = i_separationConditionX;
        m_separationCondition[ 1 ] = i_separationConditionY;
        m_separationCondition[ 2 ] = i_separationConditionZ;
    }

    void setDisplayGrid( bool i_displayGrid )
    {
        m_displayGrid = i_displayGrid;
    }    

    void setDisplayGridVelocitiesMultiplier( double i_displayGridVelocitiesMultiplier )
    {
        m_displayGridVelocitiesMultiplier = i_displayGridVelocitiesMultiplier;
    }
    
    void setMaxDisplayDensity( double i_maxDisplayDensity )
    {
        m_maxDisplayDensity = i_maxDisplayDensity;
    }
    
    void setDisplayCollisionBoundary( bool i_displayCollisionBoundary )
    {
        m_displayCollisionBoundary = i_displayCollisionBoundary;
    }

    void setDisplayAirBoundary( bool i_displayAirBoundary )
    {
        m_displayAirBoundary = i_displayAirBoundary;
    }
    
   // RodCollisionTimeStepper* setupRodTimeStepper( RodData* i_rodData );
    
    void draw(void);

    int calculateNumSubSteps(int i_subSteps,Scalar i_stepSize, double i_subDistanceMax);

    void takeTimeStep( int i_numberOfThreadsToUse, Scalar i_stepSize, 
                       int i_subSteps, double i_subDistanceMax, bool i_collisionsEnabled,
                       bool i_selfCollisionPenaltyForcesEnabled,
                       bool i_fullSelfCollisionsEnabled,
                       int i_fullSelfCollisionIters,
                       double i_selfCollisionCOR,
                       FixedRodVertexMap* i_fixedVertices = NULL,
                       bool i_zeroAllTwist = false, double i_constraintSrength = 10.0,
                       LockedRodVertexMap* i_lockedRodVertexMap = NULL );
    
    /*(void addRod( int i_rodGroup,
                 vector<BASim::Vec3d>& i_initialVertexPositions, 
                 vector<BASim::Vec3d>& i_undeformedVertexPositions,
                 RodOptions& i_options );*/
    
    void initialiseWorld( const double i_time, const double i_dt );
    void resetEverything(  const double i_time, const double i_dt );
    //void createSpaceForRods( int i_rodGroup, int i_numRods );
    void addRodsToWorld( int i_rodGroupIndex, WmFigRodGroup* i_rodGroup, double startTime );
    bool collisionMeshInitialised( const int i_collisionMeshIndex );
    void initialiseCollisionMesh( TriangleMesh* i_collisionMesh, 
                        ScriptingController* i_scriptingController, const int i_collisionMeshIndex );
    void removeCollisionMesh( const int i_collisionMeshIndex );
    
  // void checkAllRodForces(); 
    void startTimer( timeval& i_startTimer );
    double stopTimer( timeval& i_startTimer );
    void resetTimers();
    void printTimingInfo();
    void setTimingsFile( std::string i_fileName );
    void setTimingEnabled( bool i_timingsEnabled );
    void setStopOnRodError( bool i_stopOnRodError );

    void startXMLLogging( std::string& i_xmlFilePath, std::string& i_mayaSceneFilename );
    void writeXMLFileToDisk();


private:
    //! Returns true if any rods are active, false otherwise.
    /*!
      All rods can be too short to be simulated or the rod node may just be being used as a
      way to store simulated nurbs curve data in a handy file on disk. In which case there will
      be no rods active.
    */
    bool anyRodsActive();

    //! Calculates the interpolated position of collision meshes for a specific substep
    /*!
      The simulation steps at a different speed from Maya animations. Since we only get the 
      positions of the meshes at the Maya frame rate we need to work out where the meshes are
      at each simulation substep. This function does that interpolation and stores it in
      m_collisionMeshMap.
      @param i_interpolateFactor (0.0 - 1.0) The time to evaluate at between the start mesh and end mesh positions.   
      @param o_timeTaken Optional parameter which is filled in with the time the function takes to evaluate.
    */
   // void interpolateCollisionMeshesForSubstep( const float i_interpolateFactor, float* o_timeTaken = NULL );

    //! Sets a few parameters on the RodCollisionTimeStepper ready for the sim
    /*!
      FIXME: This function is very nearly pointless. It could all be done once before the substeps
      start.
      @param i_pRodGroup Pointer to a WmFigRodGroup that contains the rods to be initialised
      @param i_subStep The currently being evalued sumStep
      @param i_collisionsEnabled Whether collisions are turned on or not
    */    
    void setupRodTimeStepperForSubStep( WmFigRodGroup* i_pRodGroup, const int i_subStep,
                                        const bool i_collisionsEnabled );

    
    //void storeMaterialFrames();
    
    World* m_world;
    RodDataMap m_rodDataMap;
    //CollisionMeshDataHashMap m_collisionMeshMap;
    CollisionMeshDataHashMap m_collisionMeshDataHashMap;    
    
    ObjPropHandle<Scalar> m_timeHandle;
    ObjPropHandle<Scalar> m_dtHandle;
    ObjPropHandle<BASim::Vec3d> m_gravityHandle;
    ObjPropHandle<int> m_maxIterHandle;
    bool m_plasticDeformations;    
    BASim::Vec3d m_gravity;

    // FIXME:
    // Pointless vector with pointers to the rods. Get rid of it. It 
    // is here simply because the self collision code needed it and I have
    // to get some numbers out of this before I leave Columbia.
    vector<ElasticRod*> m_rods;

    bool m_stopOnRodError;

    bool m_timingEnabled;
    std::string m_timingsFile;
    timeval m_timerStart;
    
    double m_meshInterpolationTime;
    double m_vertexInterpolationTime;
    double m_objectCollisionForces;
    double m_objectCollisionResponse;
    double m_collisionStructuresTidyup;
    double m_selfCollisionPenaltyForceTime;
    double m_selfCollisionsResponseTime;
    double m_fastestSelfCollisionPenaltyForceTime;
    double m_slowestSelfCollisionPenaltyForceTime;
    double m_fastestSelfCollisionsResponseTime;
    double m_slowestSelfCollisionsResponseTime;
    double m_integrationStepTime;
    double m_slowestIntegrationTime;
    double m_fastestIntegrationTime;
    double m_slowestCollisionForcesTime;
    double m_fastestCollisionForcesTime;
    double m_slowestCollisionResponseTime;
    double m_fastestCollisionResponseTime;
    double m_fastestFrameTime;
    double m_slowestFrameTime;
    double m_totalSimTime;
    int m_numberOfFramesSimulated;
    
    ofstream m_timingsFP;
    
    int m_numberofThreadsUsed;
    int m_numRods;
    
    vector<MaterialFrame> m_rodRootMaterialFrame;
    vector<MaterialFrame> m_strandRootMaterialFrame;
    vector<MaterialFrame> m_rodRefMaterialFrame;

    vector< vector < vector < BASim::Vec3d > > > m_subSteppedVertexPositions;
    bool m_shouldDrawSubsteppedVertices;

    bool m_isClumpingEnabled;
    double m_clumpingCoefficient;

    double m_stol;
    double m_atol;
    double m_rtol;
    double m_inftol;

    bool m_isXMLLoggingEnabled;
    //SceneXML* m_sceneXML;

    //std::vector< InitialRodConfiguration > m_initialRodConfigurations;

    //VolumetricCollisions* m_volumetricCollisions;
    
    double m_flip;
    double m_slip;
    bool m_doVolumetricCollisions;
    double m_targetEdgeDensity;
    double m_volumetricRadius;
    double m_gridDX;    
    bool m_displayGrid;
    double m_displayGridVelocitiesMultiplier;
    double m_maxDisplayDensity;
    bool m_displayCollisionBoundary;
    bool m_displayAirBoundary;
    BASim::Vec3d m_separationCondition;

    BridsonStepper* m_bridsonStepper;

    // TODO: This is dumb, the meshes and controllers are stored in a hash map of collisionmeshdata
    // now we have this representation just to send to bridsonStepper. Fix one or the other
    vector< TriangleMesh* > m_triangleMeshes;
    vector< ScriptingController* > m_scriptingControllers;
    vector< RodTimeStepper* > m_rodTimeSteppers;
};

#endif // BEAKER_HH_
