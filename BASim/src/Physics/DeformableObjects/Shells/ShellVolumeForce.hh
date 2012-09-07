/**
 * \file ShellSurfaceTensionForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date Nov 22, 2011
 */

#ifndef SHELLVOLUMEFORCE_HH
#define SHELLVOLUMEFORCE_HH


#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

#include "BASim/src/Math/ADT/adreal.h"
#include "BASim/src/Math/ADT/advec.h"


//A force that maintains the volume of a closed region of shell.

namespace BASim {


typedef Scalar Real;
typedef CVec3T<Real> Vector3d;

const int NumVolDof = 9;

class ShellVolumeForce : public ElasticShellForce {

public:

  ShellVolumeForce (ElasticShell& shell, const std::string& name = "ShellVolumeForce", Scalar strength = 0);
  virtual ~ShellVolumeForce () {}

  std::string getName() const;

  void computeReferenceVolume();
  void computeRefPoint();

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;
  
  void update();

protected:

  bool gatherDOFs(const FaceHandle& fh, std::vector<Vec3d>& deformed, std::vector<int>& indices) const;

  Scalar elementEnergy(const std::vector<Vec3d>& deformed) const;
  void elementForce(const std::vector<Vec3d>& deformed, 
                    Eigen::Matrix<Scalar, 9, 1>& force) const;
  void elementJacobian(const std::vector<Vec3d>& deformed, 
                       Eigen::Matrix<Scalar, 9, 9>& J) const;
  
  std::vector<Scalar> m_target_volumes;
  Vec3d m_ref_point;
  Scalar m_strength;

};




}


#endif //SHELLVOLUMEFORCE_HH
