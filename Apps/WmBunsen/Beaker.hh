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

    
    RodTimeStepper* getRodTimeStepper(ElasticRod& rod);
    
    void display(void);
    void idle();
    
    void addRod( vector<Vec3d>& i_initialVertexPositions, 
                 vector<Vec3d>& i_undeformedVertexPositions,
                 RodOptions& i_options );

private:
    World* m_world;
    
    bool m_dynamicsProps;
    ObjPropHandle<Scalar> m_time;
    ObjPropHandle<Scalar> m_dt;
    ObjPropHandle<Vec3d> m_gravity;
    
    void AtEachTimestep();
    
    ElasticRod* rod;
    RodTimeStepper* stepper;
    
    RodRenderer* m_rod_renderer;
    
    Scalar m_maxTwist;
    Scalar m_twistRate;
    Scalar m_currentTwist;
};

#endif // BASIMULATOR_HH
