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
#include "BASim/src/Math/ADT/adreal.h"
#include "BASim/src/Math/ADT/advec.h"

//TODO This is not fully debugged or verified. Free boundaries behave oddly, and the xi values are not updated
//to reflect changes in the reference coordinates required for the updated Lagrangian approach. - Christopher

namespace BASim {

const int NumTriPoints = 3;

//typedefs and definitions from Meshopt
typedef Scalar Real;
typedef CVec3T<Real> Vector3d;

struct MNPrecomputed { 

  // 2 parameters in the energy formula related to ref config. ->Per Face!
  Vector3d tau0[NumTriPoints];  // perp to the average edge normal for ref config
  Real     c0[NumTriPoints];    // 1/dot(tunit0[0],tau0[0])    

  // 3 parameters related to undeformed config that need to be computed
  // only once ->Per Face!
  Real   T[NumTriPoints*NumTriPoints]; // matrix T used in the energy computation
  int    s[NumTriPoints];              // +/- 1 indicating ownership
  Real   w_undef[NumTriPoints];        // coefficients of the shape operator 
                                       // in the undeformed state in t_i X t_i basis
};

const int NumMNBendDof = 12;  //3 vertices x dofs each, + 3 edges x 1 dof each
const int MNBendStencilSize = 21; //same as above + 3 flap vertices x 3 dofs each

class MNBendingForce : public ElasticShellForce {

public:

  MNBendingForce (ElasticShell& shell, const std::string& name = "MNBendingForce", Scalar Youngs = 0, Scalar Poisson = 0, Scalar Youngs_damping = 0, Scalar Poisson_damping = 0, Scalar timestep = 1);
  virtual ~MNBendingForce () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;

  void update();

private:

  Scalar elementEnergy(const std::vector<Scalar>& undeformed,
                    const std::vector<Scalar>& deformed, MNPrecomputed* pre) const;

  void elementForce(const std::vector<Scalar>& undeformed,
                    const std::vector<Scalar>& deformed,
                    Eigen::Matrix<Scalar,NumMNBendDof,1>& force, MNPrecomputed* pre) const;

 void elementJacobian(const std::vector<Scalar>& undeformed, 
                      const std::vector<Scalar>& deformed, 
                      Eigen::Matrix<Scalar,NumMNBendDof,NumMNBendDof>& jac, MNPrecomputed* pre) const;
 
 void doPrecomputation() const;
 void initializePrecomp( const FaceHandle& face, const std::vector<Scalar>& undeformed, MNPrecomputed* pre) const;
 void updatePrecomp(const FaceHandle& face, const std::vector<Scalar>& deformed, MNPrecomputed* pre) const;

 bool gatherDOFs(const FaceHandle& edge, 
                std::vector<Scalar>& undeformed, 
                std::vector<Scalar>& undeformed_damp, 
                std::vector<Scalar>& deformed, 
                std::vector<int>& indices ) const;
 
  Scalar m_timestep;
  Scalar m_Youngs, m_Poisson;
  Scalar m_Youngs_damp, m_Poisson_damp;
  mutable FaceProperty<MNPrecomputed> m_precomputed;
  mutable bool m_initialized;
  mutable bool m_updated;
};


}

#endif //MNBENDINGFORCE_HH
