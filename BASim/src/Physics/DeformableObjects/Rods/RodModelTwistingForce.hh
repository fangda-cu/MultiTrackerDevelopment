/**
 * \file RodModelTwistingForce.h
 *
 * \author fang@cs.columbia.edu
 * \date May 12, 2011
 */

#ifndef RODMODELTWISTINGFORCE_HH
#define RODMODELTWISTINGFORCE_HH

#include "BASim/src/Physics/DeformableObjects/Rods/RodModelForce.hh"

namespace BASim 
{
  class RodModelTwistingForce : public RodModelForce
  {
  public:
    typedef Eigen::Matrix<Scalar, 11, 1> ElementForce;
    typedef Eigen::Matrix<Scalar, 11, 11> ElementJacobian;
    
    struct Stencil : public ElasticRodModel::JointStencil
    {
      Stencil(const ElasticRodModel::JointStencil & s) : ElasticRodModel::JointStencil(s) { }
      
      // cached stiffness
      Scalar stiffness;
      Scalar viscous_stiffness;
      
      // reference strain
      Scalar undeformed_twist;
      Scalar damping_undeformed_twist;
      Scalar reference_voronoi_length;
      
      // cached properties
      Scalar twist;
    };
  
  public:
    RodModelTwistingForce(ElasticRodModel & rod, const std::vector<ElasticRodModel::JointStencil> & stencils, Scalar shear_modulus, Scalar shear_modulus_damping, Scalar timestep);
    virtual ~RodModelTwistingForce();
    
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

    ElementForce computeGradTwist(Stencil & s);
    ElementJacobian computeHessTwist(Stencil & s);
    
  protected:
    std::vector<Stencil> m_stencils;
    
    Scalar m_shear_modulus;
    Scalar m_shear_modulus_damping;
    
  };
  
}


#endif // RODMODELTWISTINGFORCE_HH
