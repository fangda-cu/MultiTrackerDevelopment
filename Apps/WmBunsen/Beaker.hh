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
#include <weta/Wfigaro/Collisions/CollisionMeshData.hh>
#include <weta/Wfigaro/Core/ObjectControllerBase.hh>
#include <weta/Wfigaro/Physics/World.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodCollisionTimeStepper.hh>
#include <weta/Wfigaro/Render/RodRenderer.hh>
#else
#include <BASim/src/Core/EigenIncludes.hh>
#include <BASim/src/Collisions/CollisionMeshData.hh>
#include <BASim/src/Core/ObjectControllerBase.hh>
#include <BASim/src/Physics/World.hh>
#include <BASim/src/Physics/ElasticRods/RodCollisionTimeStepper.hh>
#include <BASim/src/Render/RodRenderer.hh>
#endif

#include "RodData.hh"
#include <tr1/unordered_map>
#include <iostream>
#include <ext/hash_map>
#include <iostream>
#include <fstream>

using namespace BASim;
using namespace std;

typedef tr1::unordered_map<size_t, vector<RodData*> > RodDataMap;
typedef RodDataMap::iterator RodDataMapIterator;

class Beaker
{
public:
    Beaker();
    ~Beaker();

    void BaseAtEachTimestep();
    void BaseAfterLoad();
    
    // Accessor functions
    
    const World& getWorld() const { return *m_world; }
    World& getWorld() { return *m_world; }
    
    const Scalar& getTime() const
    {
        return m_world->property( m_timeHandle );
    }

    void setTime( const Scalar& time )
    {
        m_world->property( m_timeHandle ) = time;
    }

    const Scalar& getDt() const
    {
        return m_world->property( m_dtHandle );
    }

    void setDt( const Scalar& dt )
    {
        m_world->property( m_dtHandle ) = dt;
    }

    const Vec3d& getGravity() const
    {
      return m_world->property( m_gravityHandle );
    }
    
    void setGravity( const Vec3d& gravity )
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

    vector<RodData*>* rodData( size_t i_rodGroup )
    {
        return &( m_rodDataMap[ i_rodGroup ] );
    }

    void setPlasticDeformations( bool i_plasticDeformations )
    {
        m_plasticDeformations = i_plasticDeformations;
    }
    
    RodCollisionTimeStepper* setupRodTimeStepper( RodData* i_rodData );
    
    void draw(void);

    void takeTimeStep( int i_numberOfThreadsToUse, Scalar i_stepSize, 
                       int i_subSteps, bool i_collisionsEnabled,
                       bool i_selfCollisionPenaltyForcesEnabled,
                       bool i_fullSelfCollisionsEnabled,
                       int i_fullSelfCollisionIters,
                       double i_selfCollisionCOR );
    
    /*(void addRod( size_t i_rodGroup,
                 vector<Vec3d>& i_initialVertexPositions, 
                 vector<Vec3d>& i_undeformedVertexPositions,
                 RodOptions& i_options );*/
    
    void initialiseWorld();
    void resetEverything();
    void createSpaceForRods( size_t i_rodGroup, size_t i_numRods );
    void createRods( size_t i_rodGroup, ObjectControllerBase::SolverLibrary solverLibrary );
    bool collisionMeshInitialised( const size_t id );
    void initialiseCollisionMesh( BASim::CollisionMeshData *collisionMeshData, size_t id );
    void removeCollisionMesh( const size_t id );
    void checkAllRodForces(); 
    void startTimer( timeval& i_startTimer );
    double stopTimer( timeval& i_startTimer );
    void resetTimers();
    void printTimingInfo();
    void setTimingsFile( std::string i_fileName );
    void setTimingEnabled( bool i_timingsEnabled );
    std::string makeString( double i_val );

private:
    //void storeMaterialFrames();
    
    World* m_world;
    RodDataMap m_rodDataMap;
    CollisionMeshDataHashMap m_collisionMeshMap;
    
    ObjPropHandle<Scalar> m_timeHandle;
    ObjPropHandle<Scalar> m_dtHandle;
    ObjPropHandle<Vec3d> m_gravityHandle;
    ObjPropHandle<int> m_maxIterHandle;
    bool m_plasticDeformations;    
    Vec3d m_gravity;

    // FIXME:
    // Pointless vector with pointers to the rods. Get rid of it. It 
    // is here simply because the self collision code needed it and I have
    // to get some numbers out of this before I leave Columbia.
    vector<ElasticRod*> m_rods;
    
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
    size_t m_numberOfFramesSimulated;
    
    ofstream m_timingsFP;
    
    int m_numberofThreadsUsed;
    size_t m_numRods;
    
    vector<MaterialFrame> m_rodRootMaterialFrame;
    vector<MaterialFrame> m_strandRootMaterialFrame;
    vector<MaterialFrame> m_rodRefMaterialFrame;
};

#endif // BEAKER_HH_
