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

using namespace BASim;
using namespace std;
//using namespace tr1;


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
        return m_world->property( m_time );
    }

    void setTime( const Scalar& time )
    {
        m_world->property( m_time ) = time;
    }

    const Scalar& getDt() const
    {
        return m_world->property( m_dt );
    }

    void setDt( const Scalar& dt )
    {
        m_world->property( m_dt ) = dt;
    }

    const Vec3d& getGravity() const
    {
      return m_world->property( m_gravity );
    }
    
    // FIXME:
    // Changing this does nothing!
    void setGravity( const Vec3d& gravity )
    {
      m_world->property( m_gravity ) = gravity;
    }

    vector<RodData*>* rodData( size_t i_rodGroup )
    {
        return &( m_rodDataMap[ i_rodGroup ] );
    }
    
    RodTimeStepper* setupRodTimeStepper( ElasticRod& rod );
    
    void draw(void);
    void takeTimeStep();
    
    void addRod( size_t i_rodGroup,
                 vector<Vec3d>& i_initialVertexPositions, 
                 vector<Vec3d>& i_undeformedVertexPositions,
                 RodOptions& i_options );
    
    void initialiseWorld();
    void resetEverything();
    void createSpaceForRods( size_t i_rodGroup, size_t i_numRods );
    void createRods( size_t i_rodGroup );

private:
    World* m_world;
    RodDataMap m_rodDataMap;
    
    ObjPropHandle<Scalar> m_time;
    ObjPropHandle<Scalar> m_dt;
    ObjPropHandle<Vec3d> m_gravity;

};

#endif // BEAKER_HH_
