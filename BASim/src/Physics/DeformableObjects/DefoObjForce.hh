/**
 * \file DefoObjForce.h
 *
 * \author batty@cs.columbia.edu
 * \date 12/04/2011
 */

#ifndef DEFOOBJFORCE_H
#define DEFOOBJFORCE_H

#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Physics/DeformableObjects/PositionDofsModel.hh"

namespace BASim {

/** Base class for a force that acts on topological objects. */
class DefoObjForce
{
public:

  DefoObjForce(PositionDofsModel& model, const std::string& name = "DefoObjForce"): m_name(name), m_model(model) {}
  virtual ~DefoObjForce() {}

  virtual std::string getName() const { return m_name; }

  virtual Scalar globalEnergy() const = 0;
  virtual void globalForce(VecXd& force) const = 0;
  virtual void globalJacobian(Scalar scale, MatrixBase& Jacobian) const = 0;
 
protected:

  std::string m_name;

  PositionDofsModel& m_model;
};



} // namespace BASim

#endif // RODFORCE_HH
