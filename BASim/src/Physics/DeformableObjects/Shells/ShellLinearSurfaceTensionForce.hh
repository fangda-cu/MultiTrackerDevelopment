/**
 * \file ShellLinearSurfaceTensionForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date Aug 21, 2012
 */

#ifndef SHELLLINEARSURFACETENSIONFORCE_HH
#define SHELLLINEARSURFACETENSIONFORCE_HH


#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

#include "BASim/src/Math/ADT/adreal.h"
#include "BASim/src/Math/ADT/advec.h"


//A surface tension force based on minimizing surface energy, assuming piecewise linear
//liquid thickness between triangle barycentres. See the writeup by Christopher.

namespace BASim {


typedef Scalar Real;
typedef CVec3T<Real> Vector3d;

const int NumLSTDof = 12;

class ShellLinearSurfaceTensionForce  : public ElasticShellForce {

public:

  ShellLinearSurfaceTensionForce (ElasticShell& shell, const std::string& name = "ShellLinearSurfaceTensionForce", Scalar surfCoeff = 0);
  virtual ~ShellLinearSurfaceTensionForce  () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;
  
protected:

  bool gatherDOFs(const EdgeHandle& eh, std::vector<Vec3d>& deformed, FaceHandle& f1, FaceHandle& f2, std::vector<int>& indices) const;

  Scalar elementEnergy(const std::vector<Vec3d>& deformed, Scalar vol1, Scalar vol2) const;
  void elementForce(const std::vector<Vec3d>& deformed, Scalar vol1, Scalar vol2, 
                    Eigen::Matrix<Scalar, NumLSTDof, 1>& force) const;
  void elementJacobian(const std::vector<Vec3d>& deformed, Scalar vol1, Scalar vol2, 
                       Eigen::Matrix<Scalar, NumLSTDof, NumLSTDof>& J) const;
  
  Scalar m_surface_tension_coeff;

};




}


#endif //SHELLLINEARSURFACETENSIONFORCE_HH
