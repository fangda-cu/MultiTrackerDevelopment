#include "RodModelTwistingForce.hh"

using namespace BASim;

RodModelTwistingForce::RodModelTwistingForce(ElasticRodModel & rod, Scalar shear_modulus, Scalar shear_modulus_damping) :
  RodModelForce(rod, "RodModelStretchingForce"),
  m_shear_modulus(shear_modulus),
  m_shear_modulus_damping(shear_modulus_damping)
{
  
}

RodModelTwistingForce::~RodModelTwistingForce()
{
  
}

Scalar RodModelTwistingForce::globalEnergy()
{
  
}

void RodModelTwistingForce::globalForce(VecXd & force)
{
  
}

void RodModelTwistingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
{
  
}

Scalar RodModelTwistingForce::localEnergy(Stencil & s)
{
  
}

void RodModelTwistingForce::localForce(ElementForce & f, Stencil & s)
{
  
}

void RodModelTwistingForce::localJacobian(ElementJacobian & f, Stencil & s)
{
  
}
