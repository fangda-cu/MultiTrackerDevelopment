#include "RodModelStretchingForce.hh"

using namespace BASim;

RodModelStretchingForce::RodModelStretchingForce(ElasticRodModel & rod, Scalar youngs_modulus, Scalar youngs_modulus_damping) :
  RodModelForce(rod, "RodModelStretchingForce"),
  m_youngs_modulus(youngs_modulus),
  m_youngs_modulus_damping(youngs_modulus_damping)
{
  
}

RodModelStretchingForce::~RodModelStretchingForce()
{
  
}

Scalar RodModelStretchingForce::globalEnergy()
{
  
}

void RodModelStretchingForce::globalForce(VecXd & force)
{
  
}

void RodModelStretchingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
{
  
}

Scalar RodModelStretchingForce::localEnergy(Stencil & s)
{
  
}

void RodModelStretchingForce::localForce(ElementForce & f, Stencil & s)
{
  
}

void RodModelStretchingForce::localJacobian(ElementJacobian & f, Stencil & s)
{
  
}
