/**
 * \file ShellSurfaceTensionForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date Nov 22, 2011
 */

#ifndef SHELLSURFACETENSIONFORCE_HH
#define SHELLSURFACETENSIONFORCE_HH


#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

#include "BASim/src/Math/ADT/adreal.h"
#include "BASim/src/Math/ADT/advec.h"


//A surface tension force based on minimizing surface energy, assuming piecewise constant
//liquid thickness when computing the surface area (as a shallow triangular prism)

namespace BASim {


typedef Scalar Real;
typedef CVec3T<Real> Vector3d;

const int NumSTDof = 9;

class ShellSurfaceTensionForce : public ElasticShellForce {

public:

  ShellSurfaceTensionForce (ElasticShell& shell, const std::string& name = "ShellSurfaceTensionForce", Scalar surfCoeff = 0);
  virtual ~ShellSurfaceTensionForce () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;
  
protected:

  bool gatherDOFs(const FaceHandle& fh, std::vector<Vec3d>& deformed, std::vector<int>& indices) const;

  Scalar elementEnergy(const std::vector<Vec3d>& deformed) const;
  void elementForce(const std::vector<Vec3d>& deformed, 
                    Eigen::Matrix<Scalar, 9, 1>& force) const;
  void elementJacobian(const std::vector<Vec3d>& deformed, 
                       Eigen::Matrix<Scalar, 9, 9>& J) const;
  
  Scalar m_surface_tension_coeff;

};




}


#endif //SHELLSURFACETENSIONFORCE_HH
