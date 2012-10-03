/**
 * \file DefoObjForce.h
 *
 * \author fang@cs.columbia.edu
 * \date 09/07/2012
 */

#ifndef DEFOOBJFORCE_H
#define DEFOOBJFORCE_H

#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"

namespace BASim 
{

// Force that acts on one or multiple models
class DefoObjForce
{
public:
  DefoObjForce(DeformableObject & obj, Scalar timestep, const std::string & name = "DefoObjForce") :
    m_obj(obj),
    m_name(name),
    m_timestep(timestep)
  { }
  virtual ~DefoObjForce() { }

  virtual std::string getName() const { return m_name; }
  
public:
  virtual Scalar globalEnergy() = 0;
  virtual void globalForce(VecXd & force) = 0;
  virtual void globalJacobian(Scalar scale, MatrixBase & Jacobian) = 0;

public:
  virtual void updateStiffness() { }               // called whenever rod radii change, or time step changes (for viscous stiffness)
  virtual void updateViscousReferenceStrain() { }  // called at the beginning of every time step
  virtual void updateProperties() { }              // called at every solver iteration (rod updateProperties()), updating cached properties
  
public:
  virtual void startStep(Scalar time, Scalar timestep) { }
  virtual void endStep(Scalar time, Scalar timestep) { }
  virtual void startIteration(Scalar time, Scalar timestep) { }
  virtual void endIteration(Scalar time, Scalar timestep) { }

public:
  DeformableObject & deformableObject() { return m_obj; }
  const DeformableObject & deformableObject() const { return m_obj; }
  
  Scalar & timeStep() { return m_timestep; }
  const Scalar & timeStep() const { return m_timestep; }
  
protected:
  DeformableObject & m_obj;
  std::string m_name;

  Scalar m_timestep;  // this is needed for computing viscous forces' stiffnesses

};



} // namespace BASim

#endif // RODFORCE_HH
