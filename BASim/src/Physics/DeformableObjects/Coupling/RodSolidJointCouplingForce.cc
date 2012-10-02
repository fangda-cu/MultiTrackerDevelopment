#include "RodSolidJointCouplingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"

using namespace BASim;

RodSolidJointCouplingForce::RodSolidJointCouplingForce(ElasticRodModel & rod, ElasticSolid & solid, const std::vector<Stencil> & stencils, Scalar stiffness, Scalar stiffness_damp, Scalar timestep) :
  DefoObjForce(rod.getDefoObj(), timestep, "RodSolidJointCouplingForce"),
  m_rod(&rod),
  m_solid(&solid),
  m_stencils(),
  m_stiffness(stiffness),
  m_stiffness_damp(stiffness_damp)
{
  for (size_t i = 0; i < stencils.size(); i++)
  {
    Stencil s(stencils[i]);
    s.stiffness = 0;
    s.viscous_stiffness = 0;
    s.undeformed_AB.setZero();
    s.undeformed_AC.setZero();
    s.undeformed_AD.setZero();
    s.damping_undeformed_AB.setZero();
    s.damping_undeformed_AC.setZero();
    s.damping_undeformed_AD.setZero();
    
    std::vector<VertexHandle> vh = getVertices(s);
    s.dofindices.resize(NumDof);
    int dofbase = defoObj().getPositionDofBase(vh[0]);
    s.dofindices[0] = dofbase;
    s.dofindices[1] = dofbase + 1;
    s.dofindices[2] = dofbase + 2;
    dofbase = defoObj().getPositionDofBase(vh[1]);
    s.dofindices[3] = dofbase;
    s.dofindices[4] = dofbase + 1;
    s.dofindices[5] = dofbase + 2;
    dofbase = defoObj().getPositionDofBase(vh[2]);
    s.dofindices[6] = dofbase;
    s.dofindices[7] = dofbase + 1;
    s.dofindices[8] = dofbase + 2;
    dofbase = defoObj().getPositionDofBase(vh[3]);
    s.dofindices[9] = dofbase;
    s.dofindices[10] = dofbase + 1;
    s.dofindices[11] = dofbase + 2;
    dofbase = defoObj().getPositionDofBase(vh[4]);
    s.dofindices[12] = dofbase;
    s.dofindices[13] = dofbase + 1;
    s.dofindices[14] = dofbase + 2;
    
    dofbase = rod.getEdgeDofBase(s.e);
    s.dofindices[15] = dofbase;
    
    s.AB.setZero();  // will be computed by updateProperties() below
    s.AC.setZero();
    s.AD.setZero();
    
    m_stencils.push_back(s);
  }

  updateProperties();
  updateStiffness();
  updateViscousReferenceStrain();
  computeReferenceStrain();
}

RodSolidJointCouplingForce::~RodSolidJointCouplingForce()
{
  
}

std::vector<VertexHandle> RodSolidJointCouplingForce::getVertices(const Stencil & s)
{
  std::vector<VertexHandle> vh(4);
  EdgeVertexIterator evit = defoObj().ev_iter(s.e);
  VertexHandle ve1 = *evit; ++evit; assert(evit);
  VertexHandle ve2 = *evit; ++evit; assert(!evit);
  TetVertexIterator tvit = defoObj().tv_iter(s.t);
  VertexHandle vt1 = *tvit; ++tvit; assert(tvit);
  VertexHandle vt2 = *tvit; ++tvit; assert(tvit);
  VertexHandle vt3 = *tvit; ++tvit; assert(tvit);
  VertexHandle vt4 = *tvit; ++tvit; assert(!tvit);
  if (ve1 == vt1 || ve1 == vt2 || ve1 == vt3 || ve1 == vt4) 
    vh[0] = ve1, vh[3] = ve2;
  else
    vh[0] = ve2, vh[3] = ve1;
  if (vh[0] == vt1) vh[1] = vt2, vh[2] = vt3;
  if (vh[0] == vt2) vh[1] = vt3, vh[2] = vt4;
  if (vh[0] == vt3) vh[1] = vt4, vh[2] = vt1;
  if (vh[0] == vt4) vh[1] = vt1, vh[2] = vt2;
  return vh;
}

RodSolidJointCouplingForce::Vector3d RodSolidJointCouplingForce::vec2vector(const Vec3d & input)
{
  Vector3d output;
  output.x() = input.x();
  output.y() = input.y();
  output.z() = input.z();
  return output;
}

template <int DO_HESS>
adreal<RodSolidJointCouplingForce::NumDof, DO_HESS, Scalar> 
RodSolidJointCouplingForce::adEnergy(const RodSolidJointCouplingForce & mn, const Vec3d & A, const Vec3d & B, const Vec3d & C, const Vec3d & D, const Vec3d & E, Scalar theta, const Vec3d & ref1, const Vec3d & ref2, const Vec3d & undeformed_AB, const Vec3d & undeformed_AC, const Vec3d & undeformed_AD, Scalar stiffness) 
{  
  // typedefs to simplify code below
  typedef adreal<RodSolidJointCouplingForce::NumDof, DO_HESS, Scalar> adrealElast;
  typedef CVec3T<adrealElast> advecElast;
  Mat3T<adrealElast> temp;
  typedef Mat3T<adrealElast> admatElast;

  Vector3d vA = vec2vector(A);
  Vector3d vB = vec2vector(B);
  Vector3d vC = vec2vector(C);
  Vector3d vD = vec2vector(D);
  Vector3d vE = vec2vector(E);
  Vector3d vRef1 = vec2vector(ref1);
  Vector3d vRef2 = vec2vector(ref2);
  Vector3d vRef3 = vec2vector(ref1.cross(ref2));
  Vector3d vUndeformedAB = vec2vector(undeformed_AB);
  Vector3d vUndeformedAC = vec2vector(undeformed_AC);
  Vector3d vUndeformedAD = vec2vector(undeformed_AD);

  advecElast p[4];
  adrealElast t;
  set_independent(p[0], vA, 0);
  set_independent(p[1], vB, 3);
  set_independent(p[2], vC, 6);
  set_independent(p[3], vD, 9);
  set_independent(p[4], vE, 12);
  t.set_independent(theta, 15);

  adrealElast e(0);
  
  advecElast md1 = vRef1 * cos(t) + vRef2 * sin(t);
  advecElast md2 = -vRef1 * sin(t) + vRef2 * cos(t);
  Vector3d md3 = vRef3;
  
  advecElast vAB = advecElast(dot(p[1] - p[0], md1), dot(p[1] - p[0], md2), dot(p[1] - p[0], md3));
  advecElast vAC = advecElast(dot(p[2] - p[0], md1), dot(p[2] - p[0], md2), dot(p[2] - p[0], md3));
  advecElast vAD = advecElast(dot(p[3] - p[0], md1), dot(p[3] - p[0], md2), dot(p[3] - p[0], md3));

  e = stiffness * (dot(vAB - vUndeformedAB, vAB - vUndeformedAB) + dot(vAC - vUndeformedAC, vAC - vUndeformedAC) + dot(vAD - vUndeformedAD, vAD - vUndeformedAD));

  return e;
}

Scalar RodSolidJointCouplingForce::globalEnergy()
{
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    // energy only include non-viscous forces
    energy += localEnergy(m_stencils[i], false);
  }
  return energy;
}

void RodSolidJointCouplingForce::globalForce(VecXd & force)
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

void RodSolidJointCouplingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
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

Scalar RodSolidJointCouplingForce::localEnergy(Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Vec3d E = defoObj().getVertexPosition(vh[4]);
  Scalar theta = rod().getEdgeTheta(s.e);
  
  Vec3d & ref1 = rod().getReferenceDirector1(s.e);
  Vec3d & ref2 = rod().getReferenceDirector2(s.e);
  
  adreal<NumDof, 1, Scalar> e = adEnergy<1>(*this, A, B, C, D, E, theta, ref1, ref2, (viscous ? s.damping_undeformed_AB : s.undeformed_AB), (viscous ? s.damping_undeformed_AC : s.undeformed_AC), (viscous ? s.damping_undeformed_AD : s.undeformed_AD), (viscous ? m_stiffness_damp : m_stiffness));
  Scalar energy = e.value();

  return energy;
}

void RodSolidJointCouplingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Vec3d E = defoObj().getVertexPosition(vh[4]);
  Scalar theta = rod().getEdgeTheta(s.e);
  
  Vec3d & ref1 = rod().getReferenceDirector1(s.e);
  Vec3d & ref2 = rod().getReferenceDirector2(s.e);
  
  adreal<NumDof, 1, Scalar> e = adEnergy<1>(*this, A, B, C, D, E, theta, ref1, ref2, (viscous ? s.damping_undeformed_AB : s.undeformed_AB), (viscous ? s.damping_undeformed_AC : s.undeformed_AC), (viscous ? s.damping_undeformed_AD : s.undeformed_AD), (viscous ? m_stiffness_damp : m_stiffness));
  for (int i = 0; i < NumDof; i++)
  {
    force[i] = -e.gradient(i);
  }
}

void RodSolidJointCouplingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Vec3d E = defoObj().getVertexPosition(vh[4]);
  Scalar theta = rod().getEdgeTheta(s.e);
  
  Vec3d & ref1 = rod().getReferenceDirector1(s.e);
  Vec3d & ref2 = rod().getReferenceDirector2(s.e);
  
  adreal<NumDof, 1, Scalar> e = adEnergy<1>(*this, A, B, C, D, E, theta, ref1, ref2, (viscous ? s.damping_undeformed_AB : s.undeformed_AB), (viscous ? s.damping_undeformed_AC : s.undeformed_AC), (viscous ? s.damping_undeformed_AD : s.undeformed_AD), (viscous ? m_stiffness_damp : m_stiffness));
  for (int i = 0; i < NumDof; i++)
    for (int j = 0; j < NumDof; j++)
    {
      jacobian(i, j) = -e.hessian(i, j);
    }
}

void RodSolidJointCouplingForce::updateStiffness()
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

void RodSolidJointCouplingForce::updateViscousReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.damping_undeformed_AB = s.AB;
    s.damping_undeformed_AC = s.AC;
    s.damping_undeformed_AD = s.AD;
  }
}

void RodSolidJointCouplingForce::updateProperties()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    
    std::vector<VertexHandle> vh = getVertices(s);
    Vec3d A = defoObj().getVertexPosition(vh[0]);
    Vec3d B = defoObj().getVertexPosition(vh[1]);
    Vec3d C = defoObj().getVertexPosition(vh[2]);
    Vec3d D = defoObj().getVertexPosition(vh[3]);
    Vec3d E = defoObj().getVertexPosition(vh[4]);
    Vec3d md1 = rod().getMaterialDirector1(s.e);
    Vec3d md2 = rod().getMaterialDirector2(s.e);
    Vec3d md3 = rod().getEdgeTangent(s.e);
    Vec3d AB = B - A;
    Vec3d AC = C - A;
    Vec3d AD = D - A;
    s.AB = Vec3d(AB.dot(md1), AB.dot(md2), AB.dot(md3));
    s.AC = Vec3d(AC.dot(md1), AC.dot(md2), AC.dot(md3));
    s.AD = Vec3d(AD.dot(md1), AD.dot(md2), AD.dot(md3));
  }
}

void RodSolidJointCouplingForce::computeReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.undeformed_AB = s.AB;
    s.undeformed_AC = s.AC;
    s.undeformed_AD = s.AD;
  }
}

void RodSolidJointCouplingForce::startStep(Scalar time, Scalar timestep)
{
  updateViscousReferenceStrain();
}

void RodSolidJointCouplingForce::endStep(Scalar time, Scalar timestep)
{
  updateStiffness();
}

void RodSolidJointCouplingForce::startIteration(Scalar time, Scalar timestep)
{
  
}

void RodSolidJointCouplingForce::endIteration(Scalar time, Scalar timestep)
{
  updateProperties();
}

