/**
 * \file DSBendingForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date April 18, 2011
 */

#ifndef DSBENDINGFORCE_H
#define DSBENDINGFORCE_H

#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/TransferFunction.hh"
#include <memory>
#include "BASim/src/Physics/DeformableObjects/Shells/BendingDerivatives.hh"

namespace BASim {

class DSBendingForce : public ElasticShellForce {

public:

  DSBendingForce (ElasticShell& shell, const std::string& name = "DSBendingForce", Scalar stiffness = 0, Scalar damping = 0, Scalar timestep = 1.0);
  virtual ~DSBendingForce () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;

private:

  Scalar elementEnergy(const std::vector<Vec3d>& undeformed,
                    const std::vector<Vec3d>& deformed) const;

  void elementForce(const std::vector<Vec3d>& undeformed,
                    const std::vector<Vec3d>& deformed,
                    Eigen::Matrix<Scalar,12,1>& force) const;

 void elementJacobian(const std::vector<Vec3d>& undeformed, 
                      const std::vector<Vec3d>& deformed, 
                      Eigen::Matrix<Scalar,12,12>& jac) const;

  Scalar ComputeAngle(Vec3d BA, Vec3d  BC) const;

  void ComputeDihedralAngleDerivatives(Vec3d& del_p1_theta,
    Vec3d& del_p2_theta,
    Vec3d& del_q1_theta,
    Vec3d& del_q2_theta,
    const Vec3d& p1,
    const Vec3d& p2,
    const Vec3d& q1,
    const Vec3d& q2) const;

  void Symmetrize(Mat3d& m) const;

  bool gatherDOFs(const EdgeHandle& edge, std::vector<Vec3d>& undeformed, std::vector<Vec3d>& undeformed_damp, std::vector<Vec3d>& deformed, std::vector<int>& indices ) const;

  void ComputeDihedralAngleSecondDerivatives(EnergyHessian& J,
    Scalar Kb,
    Vec3d& q1,
    Vec3d& p2,
    Vec3d& q2,
    Vec3d& p1,
    int q1Index,
    int p2Index,
    int q2Index,
    int p1Index) const;

  Scalar getEdgeThickness(const EdgeHandle& edge) const;

  Scalar m_stiffness, m_damping, m_timestep;

  std::tr1::shared_ptr<TransferFunction> m_func;
};




}


#endif //DSBENDINGFORCE_H
