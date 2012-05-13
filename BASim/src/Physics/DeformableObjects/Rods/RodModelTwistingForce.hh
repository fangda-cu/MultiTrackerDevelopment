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
    
    typedef ElasticRodModel::JointStencil Stencil;
    
  public:
    RodModelTwistingForce(ElasticRodModel & rod, Scalar shear_modulus, Scalar shear_modulus_damping, Scalar timestep);
    virtual ~RodModelTwistingForce();
    
  public:
    void addStencil(Stencil & s) { m_stencils.push_back(s); }
    std::vector<Stencil> & stencils() { return m_stencils; }
    const std::vector<Stencil> & stencils() const { return m_stencils; }
    
  public:
    void updateStiffness();
    void updateViscousReferenceStrain();
    void updateProperties();
    
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
    
    Scalar m_timestep;
    
    // cached stiffnesses
    VertexProperty<Scalar> m_stiffness;
    VertexProperty<Scalar> m_viscous_stiffness;
    
    // reference strains
    VertexProperty<Scalar> m_undeformed_twist;
    VertexProperty<Scalar> m_damping_undeformed_twist;
    VertexProperty<Scalar> m_reference_voronoi_length;
    
    // cached properties
    VertexProperty<Scalar> m_twist;
    
  };
  
}


#endif // RODMODELTWISTINGFORCE_HH
