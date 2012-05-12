#include "RodModelStretchingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"

using namespace BASim;

RodModelStretchingForce::RodModelStretchingForce(ElasticRodModel & rod, Scalar youngs_modulus, Scalar youngs_modulus_damping, Scalar timestep) :
  RodModelForce(rod, "RodModelStretchingForce"),
  m_youngs_modulus(youngs_modulus),
  m_youngs_modulus_damping(youngs_modulus_damping),
  m_timestep(timestep)
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
    // energy only include non-viscous forces
    energy += localEnergy(m_stencils[i], false);
  }
  return energy;
}

void RodModelStretchingForce::globalForce(VecXd & force)
{
  ElementForce localforce;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    // non-viscous force
    localForce(localforce, m_stencils[i], false);
    for (size_t j = 0; j < m_stencils[i].dofindices.size(); j++)
      force(m_stencils[i].dofindices[j]) += localforce(j);
    
    // viscous force
    localForce(localforce, m_stencils[i], true);
    for (size_t j = 0; j < m_stencils[i].dofindices.size(); j++)
      force(m_stencils[i].dofindices[j]) += localforce(j) / m_timestep;
  }
}

void RodModelStretchingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
{
  ElementJacobian localjacobian;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    // non-viscous force
    localJacobian(localjacobian, m_stencils[i], false);
    Jacobian.add(m_stencils[i].dofindices, m_stencils[i].dofindices, scale * localjacobian);
    
    // viscous force
    localJacobian(localjacobian, m_stencils[i], true);
    Jacobian.add(m_stencils[i].dofindices, m_stencils[i].dofindices, scale / m_timestep * localjacobian);
  }
}

Scalar RodModelStretchingForce::localEnergy(Stencil & s, bool viscous)
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

void RodModelStretchingForce::localForce(ElementForce & f, Stencil & s, bool viscous)
{
  
}

void RodModelStretchingForce::localJacobian(ElementJacobian & f, Stencil & s, bool viscous)
{
  
}
