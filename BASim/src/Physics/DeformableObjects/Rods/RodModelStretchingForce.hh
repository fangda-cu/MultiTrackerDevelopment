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
    RodModelStretchingForce(ElasticRodModel & rod, Scalar youngs_modulus, Scalar youngs_modulus_damping);
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
    Scalar localEnergy(Stencil & s);
    void localForce(ElementForce & f, Stencil & s);
    void localJacobian(ElementJacobian & f, Stencil & s);

  protected:
    std::vector<Stencil> m_stencils;
    
    Scalar m_youngs_modulus;
    Scalar m_youngs_modulus_damping;
    
  };
  
}


#endif // RODMODELSTRETCHINGFORCE_HH
