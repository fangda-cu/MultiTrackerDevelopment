/**
 * \file DSBendingForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date April 18, 2011
 */

#ifndef SHELLGRAVITYFORCE_H
#define SHELLGRAVITYFORCE_H

#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

namespace BASim {

class ShellGravityForce : public ElasticShellForce {

public:

  ShellGravityForce (ElasticShell& shell, const std::string& name = "ShellGravityForce", const Vec3d& gravityVector = Vec3d(0,-9.81,0));
  virtual ~ShellGravityForce () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;

protected:

  Vec3d m_gravity; //acceleration due to gravity

};




}


#endif //DSBENDINGFORCE_H
