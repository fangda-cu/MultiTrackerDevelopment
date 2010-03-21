/**
 * \file RodBendingForceSym.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 02/18/2010
 */

#ifndef RODBENDINGFORCESYM_HH
#define RODBENDINGFORCESYM_HH

#ifdef WETA
#include "RodForce.hh"
#include "VertexStencil.hh"
#endif

namespace BASim {

class RodBendingForceSym : public RodForceT<VertexStencil>
{
public:
  RodBendingForceSym(ElasticRod& rod);

  virtual Scalar globalEnergy();
  virtual void globalForce(VecXd& force);
  virtual void globalJacobian(MatrixBase& J);

  Scalar localEnergy(const vertex_handle& vh);
  void localForce(VecXd& F, const vertex_handle& vh);
  void localJacobian(MatXd& J, const vertex_handle& vh);

  const Mat2d& getB(const vertex_handle& vh) const;
  void setB(const vertex_handle& vh, const Mat2d& B);

  const Vec2d& getKappa(const vertex_handle& vh) const;
  void setKappa(const vertex_handle& vh, const Vec2d& kappa);

  const Vec2d& getKappaBar(const vertex_handle& vh) const;
  void setKappaBar(const vertex_handle& vh, const Vec2d& kappaBar);

  Scalar getRefVertexLength(const vertex_handle& vh) const;
  void setRefVertexLength(const vertex_handle& vh, const Scalar& length);

  virtual void updateProperties();
  virtual void updateStiffness();
  virtual void updateUndeformedStrain();
  virtual void updateReferenceDomain();

protected:

  const MatXd& getGradKappa(const vertex_handle& vh) const;
  const std::pair<MatXd, MatXd>& getHessKappa(const vertex_handle& vh) const;

  void computeGradKappa();
  void computeHessKappa();

  bool m_gradKappaValid;
  bool m_hessKappaValid;

  VPropHandle<Vec2d> m_kappa;
  VPropHandle<Vec2d> m_kappaBar;
  VPropHandle<MatXd> m_gradKappa; ///< each entry is a 11x2 matrix
  VPropHandle< std::pair<MatXd, MatXd> > m_hessKappa;

  VPropHandle<Mat2d> m_B;
  VPropHandle<Scalar> m_refVertexLength;
};

} // namespace BASim

#endif // RODBENDINGFORCESYM_HH
