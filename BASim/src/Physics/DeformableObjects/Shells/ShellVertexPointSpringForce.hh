/**
 * \file ShellVertexPointSpringForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date Nov 22, 2011
 */

#ifndef SHELLVERTEXPOINTSPRINGFORCE_HH
#define SHELLVERTEXPOINTSPRINGFORCE_HH


#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

#include "BASim/src/Math/ADT/adreal.h"
#include "BASim/src/Math/ADT/advec.h"


//A set of springs connecting vertices to fixed positions

namespace BASim {


typedef Scalar Real;
typedef CVec3T<Real> Vector3d;

const int NumSpringDofs = 3;
const int NumSpringVerts = 1;

class ShellVertexPointSpringForce : public ElasticShellForce {

public:

  ShellVertexPointSpringForce (ElasticShell& shell, const std::string& name = "ShellVertexPointSpringForce", Scalar timestep = 1.0);
  virtual ~ShellVertexPointSpringForce () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;
  
  void addSpring(const VertexHandle& vh, const Vec3d& position, Scalar stiffness, Scalar damping, Scalar restlen);

  bool hasSpring(const VertexHandle& vh);
protected:

  bool gatherDOFs(const VertexHandle& vh, std::vector<Vec3d>& deformed, std::vector<Vec3d>& undeformed_damp, std::vector<int>& indices) const;

  Scalar elementEnergy(const std::vector<Vec3d>& deformed, const Vec3d& position, Scalar strength, Scalar restlen) const;
  void elementForce(const std::vector<Vec3d>& deformed, const Vec3d& position, Scalar strength, Scalar restlen,
                    Eigen::Matrix<Scalar, NumSpringDofs, 1>& force) const;
  void elementJacobian(const std::vector<Vec3d>& deformed, const Vec3d& position, Scalar strength, Scalar restlen,
                       Eigen::Matrix<Scalar, NumSpringDofs, NumSpringDofs>& J) const;
  
  //List of springs
  std::vector<VertexHandle> m_vertices;
  std::vector<Vec3d> m_positions;
  std::vector<Scalar> m_stiffnesses;
  std::vector<Scalar> m_damping;
  std::vector<Scalar> m_restlen;

  Scalar m_timestep; //for damping/viscosity
};




}


#endif //SHELLVERTEXPOINTSPRINGFORCE_HH
