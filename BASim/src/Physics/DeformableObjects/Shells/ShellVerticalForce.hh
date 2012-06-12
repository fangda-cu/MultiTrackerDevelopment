/**
 * \file DSBendingForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date April 18, 2011
 */

#ifndef SHELLVERTICALFORCE_H
#define SHELLVERTICALFORCE_H

#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

namespace BASim {

class ShellVerticalForce : public ElasticShellForce {

public:

  ShellVerticalForce  (ElasticShell& shell, const std::string& name = "ShellVerticalForce ", const Vec3d& gravityVector = Vec3d(0,-9.81,0), Scalar forcePerUnitArea = 0);
  virtual ~ShellVerticalForce  () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;

protected:

  Vec3d m_gravity; //acceleration direction for gravity
  Scalar m_strength; //vertical force strength per unit area
};




}


#endif //SHELLVERTICALFORCE_H
