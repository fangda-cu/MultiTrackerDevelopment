/**
 * \file ShellRadialForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date Sept 20, 2011
 */

#ifndef SHELLRADIALFORCE_H
#define SHELLRADIALFORCE_H

#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

namespace BASim {

class ShellRadialForce : public ElasticShellForce {

public:

  ShellRadialForce (ElasticShell& shell, const std::string& name = "ShellRadialForce", const Vec3d& centrePos=Vec3d(0,0,0), Scalar strength=0, bool constPressure=false);
  virtual ~ShellRadialForce () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;

protected:

  Vec3d m_centre; //centre point of force
  Scalar m_strength; //outward strength of the force
  bool m_constant_pressure; //whether to divide the force by rad^2
};




}


#endif //SHELLRADIALFORCE_H
