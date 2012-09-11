#include "RodTwistShellFaceCouplingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"

using namespace BASim;

RodTwistShellFaceCouplingForce::RodTwistShellFaceCouplingForce(ElasticRodModel & rod, ElasticShell & shell, const std::vector<Stencil> & stencils, Scalar stiffness, Scalar stiffness_damp, Scalar timestep) :
  DefoObjForce(rod.getDefoObj(), timestep, "RodTwistShellFaceCouplingForce"),
  m_stencils(),
  m_stiffness(stiffness),
  m_stiffness_damp(stiffness_damp)
{
  for (size_t i = 0; i < stencils.size(); i++)
  {
    Stencil s(stencils[i]);
    s.stiffness = 0;
    s.viscous_stiffness = 0;
    s.undeformed_delta = 0;
    s.damping_undeformed_delta = 0;
    
    m_stencils.push_back(s);
  }

  updateProperties();
  updateStiffness();
  updateViscousReferenceStrain();
  computeReferenceStrain();
}

RodTwistShellFaceCouplingForce::~RodTwistShellFaceCouplingForce()
{
  
}

RodTwistShellFaceCouplingForce::Vector3d RodTwistShellFaceCouplingForce::vec2vector(const Vec3d & input)
{
  Vector3d output;
  output.x() = input.x();
  output.y() = input.y();
  output.z() = input.z();
  return output;
}

template <int DO_HESS>
adreal<RodTwistShellFaceCouplingForce::NumDof, DO_HESS, Scalar> 
RodTwistShellFaceCouplingForce::adEnergy(const RodTwistShellFaceCouplingForce & mn, const Vec3d & A, const Vec3d & B, const Vec3d & C, Scalar theta, const Vec3d & ref1, const Vec3d & ref2, Scalar undeformed_delta, Scalar stiffness) 
{  
  // typedefs to simplify code below
  typedef adreal<RodTwistShellFaceCouplingForce::NumDof, DO_HESS, Scalar> adrealElast;
  typedef CVec3T<adrealElast> advecElast;
  Mat3T<adrealElast> temp;
  typedef Mat3T<adrealElast> admatElast;

  Vector3d vA = vec2vector(A);
  Vector3d vB = vec2vector(B);
  Vector3d vC = vec2vector(C);
  Vector3d vRef1 = vec2vector(ref1);
  Vector3d vRef2 = vec2vector(ref2);
  
  advecElast p[3];
  adrealElast t;
  set_independent(p[0], vA, 0);
  set_independent(p[1], vB, 3);
  set_independent(p[2], vC, 6);
  t.set_independent(theta, 9);

  adrealElast e(0);
  
  adrealElast theta_deformed = atan2(dot(p[0] - p[1], vRef2), dot(p[0] - p[1], vRef1));
  adrealElast delta = t - theta_deformed;
  
  e = stiffness * ((delta - undeformed_delta) * (delta - undeformed_delta));
  
  return e;
}

Scalar RodTwistShellFaceCouplingForce::globalEnergy()
{
  Scalar energy = 0;
//  for (size_t i = 0; i < m_stencils.size(); i++)
//  {
//    // energy only include non-viscous forces
//    energy += localEnergy(m_stencils[i], false);
//  }
  return energy;
}

void RodTwistShellFaceCouplingForce::globalForce(VecXd & force)
{
//  ElementForce localforce;
//  for (size_t i = 0; i < m_stencils.size(); i++)
//  {
//    // non-viscous force
//    localForce(localforce, m_stencils[i], false);
//    for (size_t j = 0; j < m_stencils[i].dofindices.size(); j++)
//      force(m_stencils[i].dofindices[j]) += localforce(j);
//    
//    // viscous force
//    localForce(localforce, m_stencils[i], true);
//    for (size_t j = 0; j < m_stencils[i].dofindices.size(); j++)
//      force(m_stencils[i].dofindices[j]) += localforce(j) / timeStep();
//  }
}

void RodTwistShellFaceCouplingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
{
//  ElementJacobian localjacobian;
//  for (size_t i = 0; i < m_stencils.size(); i++)
//  {
//    // non-viscous force
//    localJacobian(localjacobian, m_stencils[i], false);
//    Jacobian.add(m_stencils[i].dofindices, m_stencils[i].dofindices, scale * localjacobian);
//    
//    // viscous force
//    localJacobian(localjacobian, m_stencils[i], true);
//    Jacobian.add(m_stencils[i].dofindices, m_stencils[i].dofindices, scale / timeStep() * localjacobian);
//  }
}

Scalar RodTwistShellFaceCouplingForce::localEnergy(Stencil & s, bool viscous)
{
//  Scalar ks = (viscous ? s.viscous_stiffness : s.stiffness);
//  
//  Scalar reflen = (viscous ? s.damping_undeformed_length : s.undeformed_length);
//  Scalar len = rod().getEdgeLength(s.e);
//  
//  return ks / 2.0 * square(len / reflen - 1.0) * reflen;
  
  Vec3d A = defoObj().getVertexPosition(s.v);
  Vec3d B = defoObj().getVertexPosition(defoObj().fromVertex(s.e));
  Vec3d C = defoObj().getVertexPosition(defoObj().toVertex(s.e));
  Scalar theta = rod().getEdgeTheta(s.e);
  
  Vec3d & ref1 = rod().getReferenceDirector1(s.e);
  Vec3d & ref2 = rod().getReferenceDirector2(s.e);
  
  adreal<NumDof, 0, Scalar> e = adEnergy<0>(*this, A, B, C, theta, ref1, ref2, (viscous ? s.damping_undeformed_delta : s.undeformed_delta), (viscous ? m_stiffness_damp : m_stiffness));
//  Scalar energy = e.value();
  Scalar energy = 0 ;

  return energy;
}

void RodTwistShellFaceCouplingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
//  Scalar ks = (viscous ? s.viscous_stiffness : s.stiffness);
//  
//  Scalar reflen = (viscous ? s.damping_undeformed_length : s.undeformed_length);
//  Scalar len = rod().getEdgeLength(s.e);
//  Vec3d tangent = rod().getEdgeTangent(s.e);
//  
//  Vec3d f = ks * (len / reflen - 1.0) * tangent;
//  force.segment<3> (0) = f;
//  force.segment<3> (3) = -f;
}

void RodTwistShellFaceCouplingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
//  Scalar ks = (viscous ? s.viscous_stiffness : s.stiffness);
//  
//  Vec3d e = rod().getEdge(s.e);
//  Scalar reflen = (viscous ? s.damping_undeformed_length : s.undeformed_length);
//  Scalar len = rod().getEdgeLength(s.e);
//  Mat3d M = ks * ((1.0 / reflen - 1.0 / len) * Mat3d::Identity() + 1.0 / len * outerProd(e, e) / square(len));
//  
//  //TODO: this is copied from RodStretchingForce::elemengJacobian(). can't this be implemented using blocks?
//  for (int i = 0; i < 3; i++)
//  {
//    for (int j = 0; j < 3; j++)
//    {
//      jacobian(i, j) = jacobian(3 + i, 3 + j) = -M(i, j);
//      jacobian(3 + i, j) = jacobian(i, 3 + j) = M(i, j);
//    }
//  }
//  
//  assert(isSymmetric(jacobian));
}

void RodTwistShellFaceCouplingForce::updateStiffness()
{
//  for (size_t i = 0; i < m_stencils.size(); i++)
//  {
//    Stencil & s = m_stencils[i];
//    Vec2d r = rod().getRadii(s.e);
//    Scalar cross_section = M_PI * r(0) * r(1);
//    s.stiffness = m_youngs_modulus * cross_section;
//    s.viscous_stiffness = m_youngs_modulus_damping * cross_section;
//  }
}

void RodTwistShellFaceCouplingForce::updateViscousReferenceStrain()
{
//  for (size_t i = 0; i < m_stencils.size(); i++)
//  {
//    Stencil & s = m_stencils[i];
//    s.damping_undeformed_length = rod().getEdgeLength(s.e);
//  }
}

void RodTwistShellFaceCouplingForce::updateProperties()
{
  
}

void RodTwistShellFaceCouplingForce::computeReferenceStrain()
{
//  for (size_t i = 0; i < m_stencils.size(); i++)
//  {
//    Stencil & s = m_stencils[i];
//    s.undeformed_length = rod().getEdgeLength(s.e);
//  }
}

