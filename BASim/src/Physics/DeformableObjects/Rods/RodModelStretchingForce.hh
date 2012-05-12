/**
 * \file RodModelStretchingForce.h
 *
 * \author fang@cs.columbia.edu
 * \date May 12, 2011
 */

#ifndef RODMODELSTRETCHINGFORCE_HH
#define RODMODELSTRETCHINGFORCE_HH

#include "BASim/src/Physics/DeformableObjects/Rods/RodModelForce.hh"
//#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"

namespace BASim 
{
  class RodModelStretchingForce : public RodModelForce
  {
  public:
    typedef Eigen::Matrix<Scalar, 6, 1> ElementForce;
    typedef Eigen::Matrix<Scalar, 6, 6> ElementJacobian;

    struct Stencil
    {
      EdgeHandle e;
      VertexHandle v1;
      VertexHandle v2;
      IntArray dofindices;
    };

  public:
    RodModelStretchingForce(ElasticRodModel & rod, Scalar youngs_modulus, Scalar youngs_modulus_damping, Scalar timestep);
    virtual ~RodModelStretchingForce();

  public:
    void addStencil(Stencil & s) { m_stencils.push_back(s); }
    std::vector<Stencil> & stencils() { return m_stencils; }
    const std::vector<Stencil> & stencils() const { return m_stencils; }
    
  public:
    Scalar globalEnergy();
    void globalForce(VecXd & force);
    void globalJacobian(Scalar scale, MatrixBase & Jacobian);
  
  protected:
    Scalar computeStiffness(Stencil & s, bool viscous);
    
    Scalar localEnergy(Stencil & s, bool viscous);
    void localForce(ElementForce & force, Stencil & s, bool viscous);
    void localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous);

  protected:
    std::vector<Stencil> m_stencils;
    
    Scalar m_youngs_modulus;
    Scalar m_youngs_modulus_damping;
    
    Scalar m_timestep;
    
  };
  
}


#endif // RODMODELSTRETCHINGFORCE_HH
