/**
 * \file ShellLinearSurfaceTensionForce2.hh
 *
 * \author batty@cs.columbia.edu
 * \date Aug 24, 2012
 */

#ifndef SHELLLINEARSURFACETENSIONFORCE2_HH
#define SHELLLINEARSURFACETENSIONFORCE2_HH


#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

#include "BASim/src/Math/ADT/adreal.h"
#include "BASim/src/Math/ADT/advec.h"


//A surface tension force based on minimizing surface energy, assuming piecewise linear
//across triangles, based on interpolated thicknesses at vertices

namespace BASim {

typedef Scalar Real;
typedef CVec3T<Real> Vector3d;

const int NumLST2Dof = 24*3;
const int NumLST2VertsMax = 24;

class ShellLinearSurfaceTensionForce2  : public ElasticShellForce {

public:

  ShellLinearSurfaceTensionForce2 (ElasticShell& shell, const std::string& name = "ShellLinearSurfaceTensionForce2", Scalar surfCoeff = 0);
  virtual ~ShellLinearSurfaceTensionForce2  () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;
  
protected:

  bool gatherDOFs(const FaceHandle& fh, 
    std::vector<Vec3d>& vertices, 
    std::vector<Vec3i>& faces,
    std::vector<Scalar>& volumes,
    Vec3i& incidentFacesPerVertex,
    Vec3i& isBdryVertex, 
    std::vector<int>& indices) const;
  
  Scalar elementEnergy(const std::vector<Vec3d>& vertices, const std::vector<Vec3i>& faces, const std::vector<Scalar>& volumes, const Vec3i& incidentFacesPerMainVertex, const Vec3i& isBdryVertex) const;
  void elementForce(const std::vector<Vec3d>& vertices, const std::vector<Vec3i>& faces, const std::vector<Scalar>& volumes, const Vec3i& incidentFacesPerMainVertex, const Vec3i& isBdryVertex, Eigen::Matrix<Scalar, NumLST2Dof, 1>& force) const;
  void elementJacobian(const std::vector<Vec3d>& vertices, const std::vector<Vec3i>& faces, const std::vector<Scalar>& volumes, const Vec3i& incidentFacesPerMainVertex, const Vec3i& isBdryVertex, Eigen::Matrix<Scalar,NumLST2Dof,NumLST2Dof>& jac) const;
  
  Scalar m_surface_tension_coeff;

};




}


#endif //SHELLLINEARSURFACETENSIONFORCE2_HH
