#include "RodModelBendingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"

using namespace BASim;

RodModelBendingForce::RodModelBendingForce(ElasticRodModel & rod, Scalar youngs_modulus, Scalar youngs_modulus_damping, Scalar timestep) :
  RodModelForce(rod, "RodModelStretchingForce"),
  m_youngs_modulus(youngs_modulus),
  m_youngs_modulus_damping(youngs_modulus_damping),
  m_timestep(timestep),
  m_stiffness(&rod.getDefoObj()),
  m_viscous_stiffness(&rod.getDefoObj()),
  m_undeformed_kappa(&rod.getDefoObj()),
  m_damping_undeformed_kappa(&rod.getDefoObj()),
  m_reference_voronoi_length(&rod.getDefoObj()),
  m_kappa(&rod.getDefoObj())
{
  updateProperties();
  updateStiffness();
  updateViscousReferenceStrain();
  computeReferenceStrain();
}

RodModelBendingForce::~RodModelBendingForce()
{
  
}

Scalar RodModelBendingForce::globalEnergy()
{
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    // energy only include non-viscous forces
    energy += localEnergy(m_stencils[i], false);
  }
  return energy;
}

void RodModelBendingForce::globalForce(VecXd & force)
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

void RodModelBendingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
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

Scalar RodModelBendingForce::localEnergy(Stencil & s, bool viscous)
{
//  //TODO: this can use optimization (caching quantities like edge length, like BASim does with updateProperties)
//  Mat2d B = computeStiffness(s, viscous);
//  Scalar len = getRefVertexLength(s.v);
//  
//  const Vec2d& kappa = getKappa(vh);
//  const Vec2d& kappaBar = getKappaBar(vh);
//  
//  return 0.5 / len * (kappa - kappaBar).dot(B * (kappa - kappaBar));
}

void RodModelBendingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
  
}

void RodModelBendingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
  
}

void RodModelBendingForce::updateStiffness()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    
    Vec2d ra = rod().getRadii(s.e1);
    Vec2d rb = rod().getRadii(s.e2);
    
    Scalar a = (ra(0) + rb(0)) / 2.0;
    Scalar b = (ra(1) + rb(1)) / 2.0;
    Mat2d B(Mat2d::Zero());
    B(0, 0) = M_PI * cube(a) * b / 4.0;
    B(1, 1) = M_PI * a * cube(b) / 4.0;
    
// base rotation is not supported
//    // rotate cross section
//    Mat2d rot(Mat2d::Zero());
//    rot(0, 0) = cos(m_rod.baseRotation());
//    rot(1, 0) = sin(m_rod.baseRotation());
//    rot(0, 1) = -1 * rot(1, 0);
//    rot(1, 1) = rot(0, 0);
//
//    B = rot * B * rot.transpose();
    
    m_stiffness[s.v2] = B * m_youngs_modulus;
    m_viscous_stiffness[s.v2] = B * m_youngs_modulus_damping;
  }
}

void RodModelBendingForce::updateViscousReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    m_damping_undeformed_kappa[s.v2] = m_kappa[s.v2];
  }
}

void RodModelBendingForce::updateProperties()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    const Vec3d & kb = rod().getCurvatureBinormal(s.v2);
    const Vec3d & m1e = rod().getMaterialDirector1(s.e1);
    const Vec3d & m2e = rod().getMaterialDirector2(s.e1);
    const Vec3d & m1f = rod().getMaterialDirector1(s.e2);
    const Vec3d & m2f = rod().getMaterialDirector2(s.e2);
    m_kappa[s.v2] = Vec2d(0.5 * kb.dot(m2e + m2f), -0.5 * kb.dot(m1e + m1f));
  }
}

void RodModelBendingForce::computeReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    m_undeformed_kappa[s.v2] = m_kappa[s.v2];
    m_reference_voronoi_length[s.v2] = rod().getVoronoiLength(s.v2);
  }
}


