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
    
    struct Stencil
    {
      VertexHandle v;
      EdgeHandle e1;
      EdgeHandle e2;
    };
    
  public:
    RodModelTwistingForce(ElasticRodModel & rod, Scalar shear_modulus, Scalar shear_modulus_damping);
    virtual ~RodModelTwistingForce();
    
  public:
    Scalar globalEnergy();
    void globalForce(VecXd & force);
    void globalJacobian(Scalar scale, MatrixBase & Jacobian);
    
  protected:
    Scalar localEnergy(Stencil & s);
    void localForce(ElementForce & f, Stencil & s);
    void localJacobian(ElementJacobian & f, Stencil & s);
    
  protected:
    Scalar m_shear_modulus;
    Scalar m_shear_modulus_damping;
    
  };
  
}


#endif // RODMODELTWISTINGFORCE_HH
