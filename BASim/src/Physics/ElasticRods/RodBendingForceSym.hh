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
#else
#include "BASim/src/Physics/ElasticRods/RodForce.hh"
#include "BASim/src/Physics/ElasticRods/VertexStencil.hh"
#endif

namespace BASim
{

class RodBendingForceSym: public RodForceT<VertexStencil>
{
public:
    typedef Eigen::Matrix<Scalar, 11, 1> ElementForce;
    typedef Eigen::Matrix<Scalar, 11, 2> ElementBiForce;
    typedef Eigen::Matrix<Scalar, 11, 11> ElementJacobian;
    typedef std::pair<ElementJacobian, ElementJacobian> ElementBiJacobian;

		 explicit RodBendingForceSym(ElasticRod& rod, bool vscs = false, bool runinit = true);

  virtual Scalar globalEnergy();
  virtual void globalForce(VecXd& force);
  virtual void globalJacobian(int baseidx, Scalar scale, MatrixBase& J);
  virtual void globalForceEnergy(VecXd& force, Scalar& energy);
  virtual void globalJacobianForceEnergy(int baseidx, Scalar scale, MatrixBase& Jacobian, 
					 VecXd& force, Scalar& energy);

  Scalar localEnergy(const vertex_handle& vh);
  void localForce(ElementForce& force, const vertex_handle& vh);
  void localJacobian(ElementJacobian& localJ, const vertex_handle& vh);
  void localForceEnergy(VecXd& F, Scalar& energy, const vertex_handle& vh);
  void localJacobianForceEnergy(MatXd& J, VecXd& F, Scalar& energy, const vertex_handle& vh);

    virtual void globalReverseJacobian(MatrixBase& Jacobian);
    virtual void updateReverseUndeformedStrain(const VecXd& e);

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

    virtual void updateUndeformedConfiguration(std::vector<Scalar>& vals);
    virtual void setReferenceLengths(std::vector<Scalar>& vals);
    virtual void reattatchProperties();

protected:
    const ElementBiForce& getGradKappa(const vertex_handle& vh) const;
    const ElementBiJacobian& getHessKappa(const vertex_handle& vh) const;

    void computeGradKappa();
    void computeHessKappa();

    VPropHandle<Vec2d> m_kappa;
    VPropHandle<Vec2d> m_kappaBar;

    ObjPropHandle<bool> m_gradKappaValid;
    ObjPropHandle<bool> m_hessKappaValid;
    VPropHandle<ElementBiForce> m_gradKappa; ///< each entry is a 11x2 matrix
    VPropHandle<ElementBiJacobian> m_hessKappa; ///< each entry is a pair of 11x11 matrix

    VPropHandle<Mat2d> m_B;
    VPropHandle<Scalar> m_refVertexLength;
};

} // namespace BASim

#endif // RODBENDINGFORCESYM_HH
