/**
 * \file Beaker.hh
 *
 * \author acoull@wetafx.co.nz
 * \date 10/26/2009
 */

#ifndef BEAKER_HH_
#define BEAKER_HH

#include <BASim/BASim>
#include <iostream>

using namespace BASim;
using namespace std;

// This class holds all the info needed to simulate and render a single rod.
// It feels like stepper and RodRenderer should perhaps be members of the
// ElasticRod class.
class RodData
{
public:
    RodData() { rod = NULL; stepper = NULL; rodRenderer = NULL; }
    RodData( ElasticRod* i_rod, RodTimeStepper* i_stepper, RodRenderer* i_rodRenderer ) 
    { 
        rod = i_rod; stepper = i_stepper; rodRenderer = i_rodRenderer; 
    }
    ~RodData()
    {
        if ( rod != NULL )
            delete rod;
        if ( stepper != NULL )
            delete stepper;
        if ( rodRenderer != NULL )
            delete rodRenderer;
    }
    
    ElasticRod* rod;
    RodTimeStepper* stepper;
    RodRenderer* rodRenderer;
};

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
    
    void setGravity( const Vec3d& gravity )
    {
      m_world->property( m_gravity ) = gravity;
    }

    RodTimeStepper* setupRodTimeStepper(ElasticRod& rod);
    
    void draw(void);
    void takeTimeStep();
    
    void addRod( vector<Vec3d>& i_initialVertexPositions, 
                 vector<Vec3d>& i_undeformedVertexPositions,
                 RodOptions& i_options );

private:
    World* m_world;
    vector<RodData*> m_rods;
    
    ObjPropHandle<Scalar> m_time;
    ObjPropHandle<Scalar> m_dt;
    ObjPropHandle<Vec3d> m_gravity;

};

#endif // BASIMULATOR_HH
