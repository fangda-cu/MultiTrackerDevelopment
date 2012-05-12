#include "RodModelBendingForce.hh"

using namespace BASim;

RodModelBendingForce::RodModelBendingForce(ElasticRodModel & rod, Scalar youngs_modulus, Scalar youngs_modulus_damping) :
  RodModelForce(rod, "RodModelStretchingForce"),
  m_youngs_modulus(youngs_modulus),
  m_youngs_modulus_damping(youngs_modulus_damping)
{
  
}

RodModelBendingForce::~RodModelBendingForce()
{
  
}

Scalar RodModelBendingForce::globalEnergy()
{
  
}

void RodModelBendingForce::globalForce(VecXd & force)
{
  
}

void RodModelBendingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
{
  
}

Scalar RodModelBendingForce::localEnergy(Stencil & s)
{
  
}

void RodModelBendingForce::localForce(ElementForce & f, Stencil & s)
{
  
}

void RodModelBendingForce::localJacobian(ElementJacobian & f, Stencil & s)
{
  
}
