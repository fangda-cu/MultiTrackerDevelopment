/**
 * \file CSTMembraneForce.h
 *
 * \author batty@cs.columbia.edu
 * \date April 18, 2011
 */

#ifndef CSTMEMBRANEFORCE_HH
#define CSTMEMBRANEFORCE_HH

#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

//A Constant Strain Triangle discretization of membrane forces in an elastic shell
//See the (unpublished) document by Denis Zorin for details for this specific derivation.
//Also see the update document by Christopher Batty that corrects some scaling errors.

namespace BASim {

class CSTMembraneForce : public ElasticShellForce {

public:

  CSTMembraneForce(ElasticShell& shell, const std::string& name = "CSTMembraneForceNew", Scalar Youngs = 0, Scalar Poisson = 0, Scalar Youngs_damping = 0, Scalar Poisson_damping = 0, Scalar timestep = 1);
  virtual ~CSTMembraneForce() {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;
  
protected:

  bool gatherDOFs(const FaceHandle& fh, std::vector<Vec3d>& undeformed, std::vector<Vec3d>& undeformed_damp, std::vector<Vec3d>& deformed, std::vector<int>& indices) const;

  Scalar elementEnergy(const std::vector<Vec3d>& undeformed, 
                      const std::vector<Vec3d>& deformed,
                      Scalar Young, Scalar Poisson, Scalar thickness) const;
  
  void elementForce(const std::vector<Vec3d>& undeformed, 
                    const std::vector<Vec3d>& deformed, 
                    Eigen::Matrix<Scalar, 9, 1>& force,
                    Scalar Young, Scalar Poisson, Scalar thickness) const;
  
  void elementJacobian(const std::vector<Vec3d>& undeformed, 
                       const std::vector<Vec3d>& deformed, 
                       Eigen::Matrix<Scalar, 9, 9>& J,
                       Scalar Young, Scalar Poisson, Scalar thickness) const;

  void computeHash(const std::vector<Vec3d>& undeformed, 
                   Eigen::Matrix<Scalar, 3, 3>& Tm, 
                   Scalar* lenSqv,
                   Scalar Youngs, Scalar Poisson, Scalar thickness) const;
  
  Scalar m_timestep;
  Scalar m_Youngs, m_Poisson;
  Scalar m_Youngs_damp, m_Poisson_damp;

};




}


#endif //CSTMEMBRANEFORCE_HH
