#include "RodModelStretchingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"

using namespace BASim;

RodModelStretchingForce::RodModelStretchingForce(ElasticRodModel & rod, Scalar youngs_modulus, Scalar youngs_modulus_damping, Scalar timestep) :
  RodModelForce(rod, timestep, "RodModelStretchingForce"),
  m_youngs_modulus(youngs_modulus),
  m_youngs_modulus_damping(youngs_modulus_damping),
  m_stiffness(&rod.getDefoObj()),
  m_viscous_stiffness(&rod.getDefoObj()),
  m_undeformed_length(&rod.getDefoObj()),
  m_damping_undeformed_length(&rod.getDefoObj())
{
  updateProperties();
  updateStiffness();
  updateViscousReferenceStrain();
  computeReferenceStrain();
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
      force(m_stencils[i].dofindices[j]) += localforce(j) / timeStep();
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
    Jacobian.add(m_stencils[i].dofindices, m_stencils[i].dofindices, scale / timeStep() * localjacobian);
  }
}

Scalar RodModelStretchingForce::localEnergy(Stencil & s, bool viscous)
{
  Scalar ks = (viscous ? m_viscous_stiffness[s.e] : m_stiffness[s.e]);
  
  Scalar reflen = (viscous ? m_damping_undeformed_length[s.e] : m_undeformed_length[s.e]);
  Scalar len = rod().getEdgeLength(s.e);
  
  return ks / 2.0 * square(len / reflen - 1.0) * reflen;
}

void RodModelStretchingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
  Scalar ks = (viscous ? m_viscous_stiffness[s.e] : m_stiffness[s.e]);
  
  Scalar reflen = (viscous ? m_damping_undeformed_length[s.e] : m_undeformed_length[s.e]);
  Scalar len = rod().getEdgeLength(s.e);
  Vec3d tangent = rod().getEdgeTangent(s.e);
  
  Vec3d f = ks * (len / reflen - 1.0) * tangent;
  force.segment<3> (0) = f;
  force.segment<3> (3) = -f;
}

void RodModelStretchingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
  Scalar ks = (viscous ? m_viscous_stiffness[s.e] : m_stiffness[s.e]);
  
  Vec3d e = rod().getEdge(s.e);
  Scalar reflen = (viscous ? m_damping_undeformed_length[s.e] : m_undeformed_length[s.e]);
  Scalar len = rod().getEdgeLength(s.e);
  Mat3d M = ks * ((1.0 / reflen - 1.0 / len) * Mat3d::Identity() + 1.0 / len * outerProd(e, e) / square(len));
  
  //TODO: this is copied from RodStretchingForce::elemengJacobian(). can't this be implemented using blocks?
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      jacobian(i, j) = jacobian(3 + i, 3 + j) = -M(i, j);
      jacobian(3 + i, j) = jacobian(i, 3 + j) = M(i, j);
    }
  }
  
  assert(isSymmetric(j));
}

void RodModelStretchingForce::updateStiffness()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    Vec2d r = rod().getRadii(s.e);
    Scalar cross_section = M_PI * r(0) * r(1);
    m_stiffness[s.e] = m_youngs_modulus * cross_section;
    m_viscous_stiffness[s.e] = m_youngs_modulus_damping * cross_section;
  }
}

void RodModelStretchingForce::updateViscousReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    m_damping_undeformed_length[s.e] = rod().getEdgeLength(s.e);
  }
}

void RodModelStretchingForce::updateProperties()
{
  
}

void RodModelStretchingForce::computeReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    m_undeformed_length[s.e] = rod().getEdgeLength(s.e);
  }
}

