/**
 * \file GravityForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date Oct 3, 2012
 */

#ifndef GRAVITYFORCE_H
#define GRAVITYFORCE_H

#include "BASim/src/Physics/DeformableObjects/DefoObjForce.hh"

namespace BASim {

class GravityForce : public DefoObjForce {

public:

  GravityForce (DeformableObject& obj, Scalar timestep, const std::string& name = "GravityForce", const Vec3d& gravityVector = Vec3d(0,-9.81,0));
  virtual ~GravityForce () {}

  Scalar globalEnergy() ;
  void globalForce(VecXd& force) ;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) ;

protected:

  Vec3d m_gravity; //acceleration due to gravity

};




}


#endif //GRAVITYFORCE_H
