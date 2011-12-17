/**
 * \file ShellVertexTriSpringForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date Nov 22, 2011
 */

#ifndef SHELLVERTEXTRISPRINGFORCE_HH
#define SHELLVERTEXTRISPRINGFORCE_HH


#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

#include "BASim/src/Math/ADT/adreal.h"
#include "BASim/src/Math/ADT/advec.h"


//A set of zero-rest-length springs, each between a tri and a vert.

namespace BASim {


typedef Scalar Real;
typedef CVec3T<Real> Vector3d;

const int NumVTDof = 12;
const int NumVTVerts = 4;

class ShellVertexTriSpringForce : public ElasticShellForce {

public:

  ShellVertexTriSpringForce (ElasticShell& shell, const std::string& name = "ShellVertexTriSpringForce", Scalar strength = 0);
  virtual ~ShellVertexTriSpringForce () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;
  
  void addSpring(const FaceHandle& fh, const VertexHandle& vh, const Vec3d& baryCoords);

protected:

  bool gatherDOFs(const FaceHandle& fh, const VertexHandle& vh, std::vector<Vec3d>& deformed, std::vector<int>& indices) const;

  Scalar elementEnergy(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords) const;
  void elementForce(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords, 
                    Eigen::Matrix<Scalar, NumVTDof, 1>& force) const;
  void elementJacobian(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords, 
                       Eigen::Matrix<Scalar, NumVTDof, NumVTDof>& J) const;
  
  //List of springs
  std::vector<FaceHandle> m_faces;
  std::vector<VertexHandle> m_vertices;
  std::vector<Vec3d> m_barycoords;

  Scalar m_strength;

};




}


#endif //SHELLVERTEXTRISPRINGFORCE_HH
