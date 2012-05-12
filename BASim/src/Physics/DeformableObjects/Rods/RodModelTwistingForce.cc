#include "RodModelTwistingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"

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
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    energy += localEnergy(m_stencils[i]);
  }
  return energy;
}

void RodModelTwistingForce::globalForce(VecXd & force)
{
  ElementForce localforce;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    localForce(localforce, m_stencils[i]);    
    for (size_t j = 0; j < m_stencils[i].dofindices.size(); j++)
      force(m_stencils[i].dofindices[j]) += localforce(j);
  }
}

void RodModelTwistingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
{
  ElementJacobian localjacobian;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    localJacobian(localjacobian, m_stencils[i]);    
    Jacobian.add(m_stencils[i].dofindices, m_stencils[i].dofindices, scale * localjacobian);
  }
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
