/**
 * \file RodModelForce.h
 *
 * \author fang@cs.columbia.edu
 * \date May 12, 2011
 */

#ifndef RODMODELFORCE_HH
#define RODMODELFORCE_HH

#include "BASim/src/Physics/DeformableObjects/Rods/ElasticRodModel.hh"

namespace BASim
{
  // interface following ElasticShellForce
  class RodModelForce
  {
  public:
    RodModelForce(ElasticRodModel & rod, Scalar timestep, const std::string & name = "RodModelForce") : 
      m_rod(rod), 
      m_name(name),
      m_timestep(timestep)
    { }
    virtual ~RodModelForce() { }

    std::string getName() const { return m_name; }
    
  public:
    virtual Scalar globalEnergy() = 0;
    virtual void globalForce(VecXd & force) = 0;
    virtual void globalJacobian(Scalar scale, MatrixBase & Jacobian) = 0;
    
  public:
    virtual void updateStiffness() { }               // called whenever rod radii change, or time step changes (for viscous stiffness)
    virtual void updateViscousReferenceStrain() { }  // called at the beginning of every time step
    virtual void updateProperties() { }              // called at every solver iteration (rod updateProperties()), updating cached properties
    virtual void updateAfterRemesh() { }             // called after remeshing to correct stencils and reference configs

  public:
    ElasticRodModel & rod() { return m_rod; }
    const ElasticRodModel & rod() const { return m_rod; }
    
    Scalar & timeStep() { return m_timestep; }
    const Scalar & timeStep() const { return m_timestep; }
    
  protected:
    ElasticRodModel & m_rod;
    std::string m_name;

    Scalar m_timestep;  // this is needed for computing viscous forces' stiffnesses
    
  };
  
}


#endif // RODMODELFORCE_HH
