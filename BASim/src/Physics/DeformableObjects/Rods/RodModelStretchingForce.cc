#include "RodModelStretchingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"

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
  //TODO: this can use optimization (caching quantities like edge length, like BASim does with updateProperties)
  Scalar E = (viscous ? m_youngs_modulus_damping : m_youngs_modulus);
  if (E == 0) return 0;

  Vec2d r = rod().getRadii(s.e);
  Scalar ks = E * M_PI * r(0) * r(1);
  
  DeformableObject & obj = rod().getDefoObj();
  Scalar reflen = (viscous ? 
                      (obj.getVertexDampingUndeformedPosition(s.v2) - obj.getVertexDampingUndeformedPosition(s.v1)).norm() : 
                      (obj.getVertexUndeformedPosition(s.v2) - obj.getVertexUndeformedPosition(s.v1)).norm());
  Scalar len = (obj.getVertexPosition(s.v2) - obj.getVertexPosition(s.v1)).norm();
  
  return ks / 2.0 * square(len / reflen - 1.0) * reflen;
}

void RodModelStretchingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
  //TODO: this can use optimization (caching quantities like edge length, like BASim does with updateProperties)
  Scalar E = (viscous ? m_youngs_modulus_damping : m_youngs_modulus);
  if (E == 0) return;
  
  Vec2d r = rod().getRadii(s.e);
  Scalar ks = E * M_PI * r(0) * r(1);
  
  DeformableObject & obj = rod().getDefoObj();
  Scalar reflen = (viscous ? 
                   (obj.getVertexDampingUndeformedPosition(s.v2) - obj.getVertexDampingUndeformedPosition(s.v1)).norm() : 
                   (obj.getVertexUndeformedPosition(s.v2) - obj.getVertexUndeformedPosition(s.v1)).norm());
  Scalar len = (obj.getVertexPosition(s.v2) - obj.getVertexPosition(s.v1)).norm();
  Vec3d tangent = (obj.getVertexPosition(s.v2) - obj.getVertexPosition(s.v1)).normalized();
  
  Vec3d f = ks * (len / reflen - 1.0) * tangent;
  force.segment<3> (0) = f;
  force.segment<3> (3) = -f;
}

void RodModelStretchingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
  //TODO: this can use optimization (caching quantities like edge length, like BASim does with updateProperties)
  Scalar E = (viscous ? m_youngs_modulus_damping : m_youngs_modulus);
  if (E == 0) return;
  
  Vec2d r = rod().getRadii(s.e);
  Scalar ks = E * M_PI * r(0) * r(1);
  
  DeformableObject & obj = rod().getDefoObj();
  Vec3d e = obj.getVertexPosition(s.v2) - obj.getVertexPosition(s.v1);
  Scalar reflen = (viscous ? 
                   (obj.getVertexDampingUndeformedPosition(s.v2) - obj.getVertexDampingUndeformedPosition(s.v1)).norm() : 
                   (obj.getVertexUndeformedPosition(s.v2) - obj.getVertexUndeformedPosition(s.v1)).norm());
  Scalar len = (obj.getVertexPosition(s.v2) - obj.getVertexPosition(s.v1)).norm();
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
