/**
 * \file DSBendingForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date April 18, 2011
 */

#ifndef GRAVITYFORCE_H
#define GRAVITYFORCE_H

#include "BASim/src/Physics/DeformableObjects/DefoObjForce.hh"

namespace BASim {

class GravityForce : public DefoObjForce {

public:

  GravityForce (PositionDofsModel& shell, const std::string& name = "GravityForce", const Vec3d& gravityVector = Vec3d(0,-9.81,0));
  virtual ~GravityForce () {}

  

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;

protected:

  Vec3d m_gravity; //acceleration due to gravity

};




}


#endif //GRAVITYFORCE_H
