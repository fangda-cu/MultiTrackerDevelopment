/**
 * \file RodModelBendingForce.h
 *
 * \author fang@cs.columbia.edu
 * \date May 12, 2011
 */

#ifndef RODMODELBENDINGFORCE_HH
#define RODMODELBENDINGFORCE_HH

#include "BASim/src/Physics/DeformableObjects/Rods/RodModelForce.hh"

namespace BASim 
{
  class RodModelBendingForce : public RodModelForce
  {
  public:
    typedef Eigen::Matrix<Scalar, 11, 1> ElementForce;
    typedef Eigen::Matrix<Scalar, 11, 2> ElementBiForce;
    typedef Eigen::Matrix<Scalar, 11, 11> ElementJacobian;
    typedef std::pair<ElementJacobian, ElementJacobian> ElementBiJacobian;

    struct Stencil : public ElasticRodModel::JointStencil
    {
      Stencil(const ElasticRodModel::JointStencil & s) : ElasticRodModel::JointStencil(s) { }
      
      // cached stiffness
      Mat2d stiffness;
      Mat2d viscous_stiffness;
      
      // reference strain
      Vec2d undeformed_kappa;
      Vec2d damping_undeformed_kappa;
      Scalar reference_voronoi_length;
      
      // cached properties
      Vec2d kappa;
    };
    
  public:
    RodModelBendingForce(ElasticRodModel & rod, const std::vector<ElasticRodModel::JointStencil> & stencils, Scalar youngs_modulus, Scalar youngs_modulus_damping, Scalar timestep);
    virtual ~RodModelBendingForce();
    
  public:
    void addStencil(Stencil & s) { m_stencils.push_back(s); }
    std::vector<Stencil> & stencils() { return m_stencils; }
    const std::vector<Stencil> & stencils() const { return m_stencils; }
    
  public:
    void updateStiffness();
    void updateViscousReferenceStrain();
    void updateProperties();
    void updateAfterRemesh();

  public:
    Scalar globalEnergy();
    void globalForce(VecXd & force);
    void globalJacobian(Scalar scale, MatrixBase & Jacobian);
    
  protected:
    Scalar localEnergy(Stencil & s, bool viscous);
    void localForce(ElementForce & force, Stencil & s, bool viscous);
    void localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous);
    
    void computeReferenceStrain();
    
    ElementBiForce computeGradKappa(Stencil & s);
    ElementBiJacobian computeHessKappa(Stencil & s);
    
  protected:
    std::vector<Stencil> m_stencils;
    
    Scalar m_youngs_modulus;
    Scalar m_youngs_modulus_damping;
    
  };
  
}


#endif // RODMODELBENDINGFORCE_HH
