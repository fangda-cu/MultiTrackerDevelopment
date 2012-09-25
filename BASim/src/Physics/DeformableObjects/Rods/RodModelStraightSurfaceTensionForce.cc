#include "RodModelStraightSurfaceTensionForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"

using namespace BASim;

RodModelStraightSurfaceTensionForce::RodModelStraightSurfaceTensionForce  (ElasticRodModel & rod, const std::vector<ElasticRodModel::EdgeStencil> & stencils, Scalar surface_tension_coeff) :
  RodModelForce(rod, 0, "RodModelStraightSurfaceTensionForce"),
  m_stencils(),
  m_surface_tension_coeff(surface_tension_coeff)
{
  for (size_t i = 0; i < stencils.size(); i++)
  {
    Stencil s(stencils[i]);
    m_stencils.push_back(s);
  }

}

RodModelStraightSurfaceTensionForce::~RodModelStraightSurfaceTensionForce  ()
{
  
}

Scalar RodModelStraightSurfaceTensionForce::globalEnergy()
{
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    energy += localEnergy(m_stencils[i]);
  }
  return energy;
}

void RodModelStraightSurfaceTensionForce::globalForce(VecXd & force)
{
  ElementForce localforce;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    localForce(localforce, m_stencils[i]);
    for (size_t j = 0; j < m_stencils[i].dofindices.size(); j++)
      force(m_stencils[i].dofindices[j]) += localforce(j);
  }
}

void RodModelStraightSurfaceTensionForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
{
  ElementJacobian localjacobian;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    localJacobian(localjacobian, m_stencils[i]);
    Jacobian.add(m_stencils[i].dofindices, m_stencils[i].dofindices, scale * localjacobian);
  }
}

Scalar RodModelStraightSurfaceTensionForce::localEnergy(Stencil & s)
{

  Scalar len = rod().getEdgeLength(s.e);
  Vec2d radii = rod().getRadii(s.e);
  Scalar vol = rod().getVolume(s.e);
  
  //assume radii are the same on each axis (i.e. an isotropic thread)
  //use the cylinder surface area formula: 2*pi*r*h, and simplify
  return m_surface_tension_coeff * 2 * sqrt(M_PI * vol * len);
}

void RodModelStraightSurfaceTensionForce::localForce(ElementForce & force, Stencil & s)
{
  Vec2d radii = rod().getRadii(s.e);
  Scalar len = rod().getEdgeLength(s.e);
  Vec3d edgeVec = rod().getEdge(s.e);
  Scalar vol = rod().getVolume(s.e);

  Vec3d f = m_surface_tension_coeff * sqrt(M_PI * vol) * edgeVec / pow(len, 1.5);
  force.segment<3> (0) = f;
  force.segment<3> (3) = -f;
}

void RodModelStraightSurfaceTensionForce::localJacobian(ElementJacobian & jacobian, Stencil & s)
{
  Vec2d radii = rod().getRadii(s.e);
  Scalar len = rod().getEdgeLength(s.e);
  Vec3d edgeVec = rod().getEdge(s.e);
  Scalar vol = rod().getVolume(s.e);

  Mat3d M = m_surface_tension_coeff * sqrt(M_PI*vol) * ( 1.0 / pow(len, 1.5) * Mat3d::Identity() - 1.5 * outerProd(edgeVec, edgeVec) / pow(len, 3.5));

  //TODO: this is copied from RodStretchingForce::elementJacobian(). can't this be implemented using blocks?
  for (int i = 0; i < 3; i++)
  {
    for (int j = 0; j < 3; j++)
    {
      jacobian(i, j) = jacobian(3 + i, 3 + j) = -M(i, j);
      jacobian(3 + i, j) = jacobian(i, 3 + j) = M(i, j);
    }
  }
  
  assert(isSymmetric(jacobian));
}


