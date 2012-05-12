#include "RodModelStretchingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"

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
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    energy += localEnergy(m_stencils[i]);
  }
  return energy;
}

void RodModelStretchingForce::globalForce(VecXd & force)
{
  ElementForce localforce;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    localForce(localforce, m_stencils[i]);    
    for (size_t j = 0; j < m_stencils[i].dofindices.size(); j++)
      force(m_stencils[i].dofindices[j]) += localforce(j);
  }
}

void RodModelStretchingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
{
  ElementJacobian localjacobian;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    localJacobian(localjacobian, m_stencils[i]);    
    Jacobian.add(m_stencils[i].dofindices, m_stencils[i].dofindices, scale * localjacobian);
  }
}

Scalar RodModelStretchingForce::localEnergy(Stencil & s)
{
//  Scalar ks = getKs(eh);
//  if (ks == 0.0)
//    return 0;
//  
//  Scalar refLength = getRefLength(eh);
//  Scalar len = m_rod.getEdgeLength(eh);
//  
//  return ks / 2.0 * square(len / refLength - 1.0) * refLength;
}

void RodModelStretchingForce::localForce(ElementForce & f, Stencil & s)
{
  
}

void RodModelStretchingForce::localJacobian(ElementJacobian & f, Stencil & s)
{
  
}
