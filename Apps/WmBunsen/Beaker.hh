/**
 * \file Beaker.hh
 *
 * \author acoull@wetafx.co.nz
 * \date 10/26/2009
 */

#ifndef BEAKER_HH_
#define BEAKER_HH_

#include <tr1/unordered_map>
#include <BASim/BASim>
#include <iostream>
#include <ext/hash_map>
#include "RodData.hh"
#include <BASim/src/Collisions/CollisionMeshData.hh>

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
    
    RodCollisionTimeStepper* setupRodTimeStepper( BASim::ElasticRod& rod, 
        ObjectControllerBase::SolverLibrary solverLibrary );
    
    void draw(void);

    void takeTimeStep( int i_numberOfThreadsToUse, Scalar i_stepSize, 
                       int i_subSteps, bool i_collisionsEnabled );
    
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
    
private:
    World* m_world;
    RodDataMap m_rodDataMap;
    CollisionMeshDataHashMap m_collisionMeshMap;
    
    ObjPropHandle<Scalar> m_timeHandle;
    ObjPropHandle<Scalar> m_dtHandle;
    ObjPropHandle<Vec3d> m_gravityHandle;
    ObjPropHandle<int> m_maxIterHandle;
    
    Vec3d m_gravity;
};

#endif // BEAKER_HH_
