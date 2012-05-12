#include "RodModelBendingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"

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
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    energy += localEnergy(m_stencils[i]);
  }
  return energy;
}

void RodModelBendingForce::globalForce(VecXd & force)
{
  ElementForce localforce;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    localForce(localforce, m_stencils[i]);    
    for (size_t j = 0; j < m_stencils[i].dofindices.size(); j++)
      force(m_stencils[i].dofindices[j]) += localforce(j);
  }
}

void RodModelBendingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
{
  ElementJacobian localjacobian;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    localJacobian(localjacobian, m_stencils[i]);    
    Jacobian.add(m_stencils[i].dofindices, m_stencils[i].dofindices, scale * localjacobian);
  }
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
