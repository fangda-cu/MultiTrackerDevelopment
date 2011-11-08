/**
 * \file MNBendingForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date April 18, 2011
 */

#ifndef MNBENDINGFORCE_HH
#define MNBENDINGFORCE_HH

#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"
#include "BASim/src/Core/TopologicalObject/TopologicalObject.hh"
#include <memory>

namespace BASim {

struct MNPrecomputed;

const int NumMNBendDof = 21;

class MNBendingForce : public ElasticShellForce {

public:

  MNBendingForce (ElasticShell& shell, const std::string& name = "MNBendingForce", Scalar Youngs = 0, Scalar Poisson = 0, Scalar Youngs_damping = 0, Scalar Poisson_damping = 0, Scalar timestep = 1);
  virtual ~MNBendingForce () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;

private:

  Scalar elementEnergy(const std::vector<Scalar>& undeformed,
                    const std::vector<Scalar>& deformed, MNPrecomputed* pre) const;

  void elementForce(const std::vector<Scalar>& undeformed,
                    const std::vector<Scalar>& deformed,
                    Eigen::Matrix<Scalar,NumMNBendDof,1>& force, MNPrecomputed* pre) const;

 void elementJacobian(const std::vector<Scalar>& undeformed, 
                      const std::vector<Scalar>& deformed, 
                      Eigen::Matrix<Scalar,NumMNBendDof,NumMNBendDof>& jac, MNPrecomputed* pre) const;
 
 void initializePrecomp( const FaceHandle& face, const std::vector<Scalar>& undeformed, const std::vector<Scalar>& deformed, MNPrecomputed* pre) const;
 void updatePrecomp(const FaceHandle& face, const std::vector<Scalar>& undeformed, const std::vector<Scalar>& deformed, MNPrecomputed* pre) const;

 bool gatherDOFs(const FaceHandle& edge, 
                std::vector<Scalar>& undeformed, 
                std::vector<Scalar>& undeformed_damp, 
                std::vector<Scalar>& deformed, 
                std::vector<int>& indices ) const;
 
  Scalar m_timestep;
  Scalar m_Youngs, m_Poisson;
  Scalar m_Youngs_damp, m_Poisson_damp;

};




}


#endif //MNBENDINGFORCE_HH
