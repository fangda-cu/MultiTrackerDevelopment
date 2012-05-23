/**
 * \file SolidElasticityForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date May 23, 2012
 */

#ifndef SOLIDELASTICITYFORCE_HH
#define SOLIDELASTICITYFORCE_HH


#include "BASim/src/Physics/DeformableObjects/Solids/ElasticSolidForce.hh"

#include "BASim/src/Math/ADT/adreal.h"
#include "BASim/src/Math/ADT/advec.h"


//A linear elastic solid deformation force

namespace BASim {


typedef Scalar Real;
typedef CVec3T<Real> Vector3d;

const int NumElasticityDof = 12;

class SolidElasticityForce : public ElasticSolidForce {

public:

  SolidElasticityForce (ElasticSolid& shell, const std::string& name = "SolidElasticityForce", 
      Scalar Youngs = 0, Scalar Poisson = 0,
      Scalar Youngs_damp = 0, Scalar Poisson_damp = 0, Scalar timestep = 0);
  virtual ~SolidElasticityForce () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;
  
protected:

  bool gatherDOFs(const TetHandle& th, std::vector<Vec3d>& deformed, 
                                       std::vector<Vec3d>& undeformed_damp, 
                                       std::vector<Vec3d>& undeformed, 
                                       std::vector<int>& indices) const;

  Scalar elementEnergy(const std::vector<Vec3d>& deformed, const std::vector<Vec3d>& undeformed, Scalar Youngs, Scalar Poisson) const;
  void elementForce(const std::vector<Vec3d>& deformed, const std::vector<Vec3d>& undeformed, Scalar Youngs, Scalar Poisson,
                    Eigen::Matrix<Scalar, 12, 1>& force) const;
  void elementJacobian(const std::vector<Vec3d>& deformed, const std::vector<Vec3d>& undeformed, Scalar Youngs, Scalar Poisson,
                       Eigen::Matrix<Scalar, 12, 12>& J) const;
  
  Scalar m_Youngs, m_Poisson;
  Scalar m_Youngs_damp, m_Poisson_damp;
  Scalar m_timestep;

};




}


#endif //CSTMEMBRANEFORCE_HH
