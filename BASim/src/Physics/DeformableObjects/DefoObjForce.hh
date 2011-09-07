/**
 * \file DefoObjForce.h
 *
 * \author batty@cs.columbia.edu
 * \date 12/04/2011
 */

#ifndef DEFOOBJFORCE_H
#define DEFOOBJFORCE_H

#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Physics/SimplexObjects/DeformableObject.h"
namespace BASim {

/** Base class for a force that acts on topological objects. */
class DefoObjForce
{
public:

  DefoObjForce(DefoObjForce& object, const std::string& name = "DefoObjForce");
  virtual ~DefoObjForce() {}

  virtual std::string getName() const;

  virtual Scalar globalEnergy() = 0;
  virtual void globalForce(VecXd& force) = 0;
  virtual void globalJacobian(int baseidx, Scalar scale, MatrixBase& Jacobian) = 0;
 

protected:


  std::string m_name;

  DeformableObject& m_obj;

};



} // namespace BASim

#endif // RODFORCE_HH
