/**
 * \file RodModelStretchingForce.h
 *
 * \author fang@cs.columbia.edu
 * \date May 12, 2011
 */

#ifndef RODMODELSTRETCHINGFORCE_HH
#define RODMODELSTRETCHINGFORCE_HH

#include "BASim/src/Physics/DeformableObjects/Rods/RodModelForce.hh"

namespace BASim 
{
  class RodModelStretchingForce : public RodModelForce
  {
  public:
    typedef Eigen::Matrix<Scalar, 6, 1> ElementForce;
    typedef Eigen::Matrix<Scalar, 6, 6> ElementJacobian;

    struct Stencil : public ElasticRodModel::EdgeStencil
    {
      Stencil(const ElasticRodModel::EdgeStencil & s) : ElasticRodModel::EdgeStencil(s) { }
      
      // cached stiffness
      Scalar stiffness;
      Scalar viscous_stiffness;
      
      // reference strain
      Scalar undeformed_length;
      Scalar damping_undeformed_length;
      
      // cached properties (none)
    };

  public:
    RodModelStretchingForce(ElasticRodModel & rod, const std::vector<ElasticRodModel::EdgeStencil> & stencils, Scalar youngs_modulus, Scalar youngs_modulus_damping, Scalar timestep);
    virtual ~RodModelStretchingForce();

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
    
  protected:
    std::vector<Stencil> m_stencils;
    
    Scalar m_youngs_modulus;
    Scalar m_youngs_modulus_damping;
  };
  
}


#endif // RODMODELSTRETCHINGFORCE_HH
