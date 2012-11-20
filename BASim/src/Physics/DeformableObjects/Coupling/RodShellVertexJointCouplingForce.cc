#include "RodShellVertexJointCouplingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"

using namespace BASim;

RodShellVertexJointCouplingForce::RodShellVertexJointCouplingForce(ElasticRodModel & rod, ElasticShell & shell, const std::vector<Stencil> & stencils, Scalar stiffness, Scalar stiffness_damp, Scalar timestep) :
  DefoObjForce(rod.getDefoObj(), timestep, "RodShellVertexJointCouplingForce"),
  m_rod(&rod),
  m_shell(&shell),
  m_stencils(),
  m_stiffness(stiffness),
  m_stiffness_damp(stiffness_damp)
{
  for (size_t i = 0; i < stencils.size(); i++)
  {
    Stencil s(stencils[i]);
//    s.stiffness = 0;
//    s.viscous_stiffness = 0;
    s.undeformed_AP.setZero();
    s.undeformed_delta = 0;
    s.damping_undeformed_AP.setZero();
    s.damping_undeformed_delta = 0;
    
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
    
    dofbase = rod.getEdgeDofBase(s.e);
    s.dofindices[12] = dofbase;
    
    s.AP.setZero();  // will be computed by updateProperties() below
    s.delta = 0;
    
    m_stencils.push_back(s);
  }

  updateProperties();
  updateStiffness();
  updateViscousReferenceStrain();
  computeReferenceStrain();
}

RodShellVertexJointCouplingForce::~RodShellVertexJointCouplingForce()
{
  
}

std::vector<VertexHandle> RodShellVertexJointCouplingForce::getVertices(const Stencil & s)
{
  std::vector<VertexHandle> vh(4);
  EdgeVertexIterator evit = defoObj().ev_iter(s.e);
  VertexHandle ve1 = *evit; ++evit; assert(evit);
  VertexHandle ve2 = *evit; ++evit; assert(!evit);
  FaceVertexIterator fvit = defoObj().fv_iter(s.f);
  VertexHandle vf1 = *fvit; ++fvit; assert(fvit);
  VertexHandle vf2 = *fvit; ++fvit; assert(fvit);
  VertexHandle vf3 = *fvit; ++fvit; assert(!fvit);
  if (ve1 == vf1 || ve1 == vf2 || ve1 == vf3) vh[0] = ve1, vh[3] = ve2;
  if (ve2 == vf1 || ve2 == vf2 || ve2 == vf3) vh[0] = ve2, vh[3] = ve1;
  if (vh[0] == vf1) vh[1] = vf2, vh[2] = vf3;
  if (vh[0] == vf2) vh[1] = vf3, vh[2] = vf1;
  if (vh[0] == vf3) vh[1] = vf1, vh[2] = vf2;
  return vh;
}

RodShellVertexJointCouplingForce::Vector3d RodShellVertexJointCouplingForce::vec2vector(const Vec3d & input)
{
  Vector3d output;
  output.x() = input.x();
  output.y() = input.y();
  output.z() = input.z();
  return output;
}

template <int DO_HESS>
adreal<RodShellVertexJointCouplingForce::NumDof, DO_HESS, Scalar> 
RodShellVertexJointCouplingForce::adEnergy(const RodShellVertexJointCouplingForce & mn, const Vec3d & A, const Vec3d & B, const Vec3d & C, const Vec3d & D, Scalar delta, Scalar theta, const Vec3d & ref1, const Vec3d & ref2, const Vec3d & undeformed_AP, Scalar undeformed_delta, Scalar stiffness) 
{  
  // typedefs to simplify code below
  typedef adreal<RodShellVertexJointCouplingForce::NumDof, DO_HESS, Scalar> adrealElast;
  typedef CVec3T<adrealElast> advecElast;
  Mat3T<adrealElast> temp;
  typedef Mat3T<adrealElast> admatElast;

  Vector3d vA = vec2vector(A);
  Vector3d vB = vec2vector(B);
  Vector3d vC = vec2vector(C);
  Vector3d vD = vec2vector(D);
  Vector3d vRef1 = vec2vector(ref1);
  Vector3d vRef2 = vec2vector(ref2);
  Vector3d vRef3 = vec2vector(ref1.cross(ref2));
  Vector3d vUndeformedAP = vec2vector(undeformed_AP);

  advecElast p[4];
  adrealElast t;
  set_independent(p[0], vA, 0);
  set_independent(p[1], vB, 3);
  set_independent(p[2], vC, 6);
  set_independent(p[3], vD, 9);
  t.set_independent(theta, 12);

  adrealElast e(0);
  
  advecElast m11 = vRef1 * cos(t) + vRef2 * sin(t);
  advecElast m12 = -vRef1 * sin(t) + vRef2 * cos(t);
  Vector3d m13 = vRef3;
  if (ref1.cross(ref2).dot(A - D) < 0) // align the rod direction with AP
    m11 = -m11, m12 = -m12, m13 = -m13;

  advecElast vAP = ((p[1] + p[2]) - p[0] * 2.0);  vAP /= len(vAP);
  
  advecElast m23 = vAP;
  advecElast m21 = ((p[2] - p[1]) - dot(p[2] - p[1], vAP) * vAP); m21 /= len(m21);
  advecElast m22 = cross(m23, m21);
  advecElast m21r = m21 * cos(delta) + m22 * sin(delta);
  
  // implementation of parallel transport
  advecElast m21r_in_m1;
  advecElast b = cross(m23, m13);
  if (len(b).value() == 0)
  {
    m21r_in_m1 = m21r;
  } else
  {
    b /= len(b);
    advecElast n1 = cross(m23, b);
    advecElast n2 = cross(m13, b);
    m21r_in_m1 = dot(m21r, m23) * m13 + dot(m21r, n1) * n2 + dot(m21r, b) * b;
  }
  
  adrealElast newdelta = delta + atan2(dot(cross(m21r_in_m1, m11), m13), dot(m21r_in_m1, m11));
  advecElast vAPm = advecElast(dot(vAP, m11), dot(vAP, m12), dot(vAP, m13));

  e = stiffness * (dot(vAPm - vUndeformedAP, vAPm - vUndeformedAP) + sqr(newdelta - undeformed_delta));

  return e;
}

Scalar RodShellVertexJointCouplingForce::globalEnergy()
{
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    // energy only include non-viscous forces
    energy += localEnergy(m_stencils[i], false);
  }
  return energy;
}

void RodShellVertexJointCouplingForce::globalForce(VecXd & force)
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

void RodShellVertexJointCouplingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
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

Scalar RodShellVertexJointCouplingForce::localEnergy(Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Scalar theta = rod().getEdgeTheta(s.e);

  Vec3d & ref1 = rod().getReferenceDirector1(s.e);
  Vec3d & ref2 = rod().getReferenceDirector2(s.e);
  
  adreal<NumDof, 0, Scalar> e = adEnergy<0>(*this, A, B, C, D, s.delta, theta, ref1, ref2, (viscous ? s.damping_undeformed_AP : s.undeformed_AP), (viscous ? s.damping_undeformed_delta : s.undeformed_delta), (viscous ? m_stiffness_damp : m_stiffness));
  Scalar energy = e.value();

  return energy;
}

void RodShellVertexJointCouplingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Scalar theta = rod().getEdgeTheta(s.e);
  
  Vec3d & ref1 = rod().getReferenceDirector1(s.e);
  Vec3d & ref2 = rod().getReferenceDirector2(s.e);
  
  adreal<NumDof, 0, Scalar> e = adEnergy<0>(*this, A, B, C, D, s.delta, theta, ref1, ref2, (viscous ? s.damping_undeformed_AP : s.undeformed_AP), (viscous ? s.damping_undeformed_delta : s.undeformed_delta), (viscous ? m_stiffness_damp : m_stiffness));
  for (int i = 0; i < NumDof; i++)
  {
    force[i] = -e.gradient(i);
  }
}

void RodShellVertexJointCouplingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Scalar theta = rod().getEdgeTheta(s.e);
  
  Vec3d & ref1 = rod().getReferenceDirector1(s.e);
  Vec3d & ref2 = rod().getReferenceDirector2(s.e);
  
  adreal<NumDof, 1, Scalar> e = adEnergy<1>(*this, A, B, C, D, s.delta, theta, ref1, ref2, (viscous ? s.damping_undeformed_AP : s.undeformed_AP), (viscous ? s.damping_undeformed_delta : s.undeformed_delta), (viscous ? m_stiffness_damp : m_stiffness));
  for (int i = 0; i < NumDof; i++)
    for (int j = 0; j < NumDof; j++)
    {
      jacobian(i, j) = -e.hessian(i, j);
    }
}

void RodShellVertexJointCouplingForce::updateStiffness()
{

}

void RodShellVertexJointCouplingForce::updateViscousReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.damping_undeformed_AP = s.AP;
    s.damping_undeformed_delta = s.delta;
  }
}

void RodShellVertexJointCouplingForce::updateProperties()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    
    std::vector<VertexHandle> vh = getVertices(s);
    Vec3d A = defoObj().getVertexPosition(vh[0]);
    Vec3d B = defoObj().getVertexPosition(vh[1]);
    Vec3d C = defoObj().getVertexPosition(vh[2]);
    Vec3d D = defoObj().getVertexPosition(vh[3]);
    Vec3d m11 = rod().getMaterialDirector1(s.e);
    Vec3d m12 = rod().getMaterialDirector2(s.e);
    Vec3d m13 = rod().getEdgeTangent(s.e);
    if (m13.dot(A - D) < 0) // align the rod orientation with AP
      m11 = -m11, m12 = -m12, m13 = -m13;
    Vec3d AP = ((B + C) - A * 2).normalized();
    s.AP = Vec3d(AP.dot(m11), AP.dot(m12), AP.dot(m13));
    Vec3d m23 = AP;
    Vec3d m21 = ((C - B) - (C - B).dot(AP) * AP).normalized();
    Vec3d m22 = m23.cross(m21);
    Vec3d m21r = m21 * cos(s.delta) + m22 * sin(s.delta);
    Vec3d m21r_in_m1 = parallel_transport(m21r, m23, m13);
    s.delta += atan2(m21r_in_m1.cross(m11).dot(m13), m21r_in_m1.dot(m11));
  }
}

void RodShellVertexJointCouplingForce::computeReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.undeformed_AP = s.AP;
    s.undeformed_delta = s.delta;
  }
}

void RodShellVertexJointCouplingForce::startStep(Scalar time, Scalar timestep)
{
  updateViscousReferenceStrain();
}

void RodShellVertexJointCouplingForce::endStep(Scalar time, Scalar timestep)
{
  updateStiffness();
}

void RodShellVertexJointCouplingForce::startIteration(Scalar time, Scalar timestep)
{
  
}

void RodShellVertexJointCouplingForce::endIteration(Scalar time, Scalar timestep)
{
  updateProperties();
}

