/**
 * \file ShellBathForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date Jan 6, 2011
 */

#ifndef SHELLBATHFORCE_H
#define SHELLBATHFORCE_H

//A force representing the effect of the fluid bath

#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

namespace BASim {

class ShellBathForce : public ElasticShellForce {

public:

  ShellBathForce  (ElasticShell& shell, const std::string& name = "ShellBathForce", const Vec3d& gravityVector = Vec3d(0,-9.81,0), Scalar density = 1, Scalar bathHeight = 0);
  virtual ~ShellBathForce  () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;

protected:

  Vec3d m_gravity; //acceleration due to gravity
  Scalar m_density; //density
  Scalar m_bath_height;
};




}


#endif //SHELLBATHFORCE_H
