/**
 * \file ShellPointForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date June something, 2012
 */

#ifndef SHELLPOINTFORCE_H
#define SHELLPOINTFORCE_H

#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

namespace BASim {

class ShellPointForce : public ElasticShellForce {

public:

  ShellPointForce   (ElasticShell& shell, const std::string& name, const std::vector<VertexHandle>& verts, const std::vector<Vec3d>& forces);
  virtual ~ShellPointForce   () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;

protected:

  std::vector<VertexHandle> m_vertices;
  std::vector<Vec3d> m_forces;
};




}


#endif //SHELLVERTICALFORCE_H
