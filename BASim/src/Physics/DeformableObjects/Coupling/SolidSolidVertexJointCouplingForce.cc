#include "SolidSolidVertexJointCouplingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"

using namespace BASim;

SolidSolidVertexJointCouplingForce::SolidSolidVertexJointCouplingForce(ElasticSolid & solid, const std::vector<Stencil> & stencils, Scalar stiffness, Scalar stiffness_damp, Scalar timestep) :
  DefoObjForce(solid.getDefoObj(), timestep, "SolidSolidVertexJointCouplingForce"),
  m_solid(&solid),
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
    dofbase = defoObj().getPositionDofBase(vh[4]);
    s.dofindices[12] = dofbase;
    s.dofindices[13] = dofbase + 1;
    s.dofindices[14] = dofbase + 2;
    dofbase = defoObj().getPositionDofBase(vh[5]);
    s.dofindices[15] = dofbase;
    s.dofindices[16] = dofbase + 1;
    s.dofindices[17] = dofbase + 2;
    dofbase = defoObj().getPositionDofBase(vh[6]);
    s.dofindices[18] = dofbase;
    s.dofindices[19] = dofbase + 1;
    s.dofindices[20] = dofbase + 2;
        
    s.AP.setZero();  // will be computed by updateProperties() below
    s.delta = 0;
    
    m_stencils.push_back(s);
  }

  updateProperties();
  updateStiffness();
  updateViscousReferenceStrain();
  computeReferenceStrain();
}

SolidSolidVertexJointCouplingForce::~SolidSolidVertexJointCouplingForce()
{
  
}

std::vector<VertexHandle> SolidSolidVertexJointCouplingForce::getVertices(const Stencil & s)
{
  std::vector<VertexHandle> vh(7);
  TetVertexIterator tvit = defoObj().tv_iter(s.t1);
  VertexHandle vt11 = *tvit; ++tvit; assert(tvit);
  VertexHandle vt12 = *tvit; ++tvit; assert(tvit);
  VertexHandle vt13 = *tvit; ++tvit; assert(tvit);
  VertexHandle vt14 = *tvit; ++tvit; assert(!tvit);
  TetVertexIterator tvit2 = defoObj().tv_iter(s.t2);
  VertexHandle vt21 = *tvit2; ++tvit2; assert(tvit2);
  VertexHandle vt22 = *tvit2; ++tvit2; assert(tvit2);
  VertexHandle vt23 = *tvit2; ++tvit2; assert(tvit2);
  VertexHandle vt24 = *tvit2; ++tvit2; assert(!tvit2);
  if (vt11 == vt21 || vt11 == vt22 || vt11 == vt23 || vt11 == vt24) vh[0] = vt11, vh[4] = vt12, vh[5] = vt13, vh[6] = vt14;
  if (vt12 == vt21 || vt12 == vt22 || vt12 == vt23 || vt12 == vt24) vh[0] = vt12, vh[4] = vt13, vh[5] = vt14, vh[6] = vt11;
  if (vt13 == vt21 || vt13 == vt22 || vt13 == vt23 || vt13 == vt24) vh[0] = vt13, vh[4] = vt14, vh[5] = vt11, vh[6] = vt12;
  if (vt14 == vt21 || vt14 == vt22 || vt14 == vt23 || vt14 == vt24) vh[0] = vt14, vh[4] = vt11, vh[5] = vt12, vh[6] = vt13;
  if (vh[0] == vt21) vh[1] = vt22, vh[2] = vt23, vh[3] = vt24;
  if (vh[0] == vt22) vh[1] = vt23, vh[2] = vt24, vh[3] = vt21;
  if (vh[0] == vt23) vh[1] = vt24, vh[2] = vt21, vh[3] = vt22;
  if (vh[0] == vt24) vh[1] = vt21, vh[2] = vt22, vh[3] = vt23;
  return vh;
}

SolidSolidVertexJointCouplingForce::Vector3d SolidSolidVertexJointCouplingForce::vec2vector(const Vec3d & input)
{
  Vector3d output;
  output.x() = input.x();
  output.y() = input.y();
  output.z() = input.z();
  return output;
}

template <int DO_HESS>
adreal<SolidSolidVertexJointCouplingForce::NumDof, DO_HESS, Scalar> 
SolidSolidVertexJointCouplingForce::adEnergy(const SolidSolidVertexJointCouplingForce & mn, const Vec3d & A, const Vec3d & B, const Vec3d & C, const Vec3d & D, const Vec3d & E, const Vec3d & F, const Vec3d & G, Scalar delta, const Vec3d & undeformed_AP, Scalar undeformed_delta, Scalar stiffness)
{  
  // typedefs to simplify code below
  typedef adreal<SolidSolidVertexJointCouplingForce::NumDof, DO_HESS, Scalar> adrealElast;
  typedef CVec3T<adrealElast> advecElast;
  Mat3T<adrealElast> temp;
  typedef Mat3T<adrealElast> admatElast;

  Vector3d vA = vec2vector(A);
  Vector3d vB = vec2vector(B);
  Vector3d vC = vec2vector(C);
  Vector3d vD = vec2vector(D);
  Vector3d vE = vec2vector(E);
  Vector3d vF = vec2vector(F);
  Vector3d vG = vec2vector(G);
  Vector3d vUndeformedAP = vec2vector(undeformed_AP);

  advecElast p[7];
  set_independent(p[0], vA, 0);
  set_independent(p[1], vB, 3);
  set_independent(p[2], vC, 6);
  set_independent(p[3], vD, 9);
  set_independent(p[4], vE, 12);
  set_independent(p[5], vF, 15);
  set_independent(p[6], vG, 18);

  adrealElast e(0);
  
  advecElast vQA = (p[0] * 3.0 - (p[4] + p[5] + p[6]));  vQA /= len(vQA);
  
  advecElast m13 = vQA;
  advecElast m11 = ((p[5] - p[4]) - dot(p[5] - p[4], vQA) * vQA); m11 /= len(m11);
  advecElast m12 = cross(m13, m11);
  
  advecElast vAP = ((p[1] + p[2] + p[3]) - p[0] * 3.0);  vAP /= len(vAP);
  
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

Scalar SolidSolidVertexJointCouplingForce::globalEnergy()
{
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    // energy only include non-viscous forces
    energy += localEnergy(m_stencils[i], false);
  }
  return energy;
}

void SolidSolidVertexJointCouplingForce::globalForce(VecXd & force)
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

void SolidSolidVertexJointCouplingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
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

Scalar SolidSolidVertexJointCouplingForce::localEnergy(Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Vec3d E = defoObj().getVertexPosition(vh[4]);
  Vec3d F = defoObj().getVertexPosition(vh[5]);
  Vec3d G = defoObj().getVertexPosition(vh[6]);
  
  adreal<NumDof, 0, Scalar> e = adEnergy<0>(*this, A, B, C, D, E, F, G, s.delta, (viscous ? s.damping_undeformed_AP : s.undeformed_AP), (viscous ? s.damping_undeformed_delta : s.undeformed_delta), (viscous ? m_stiffness_damp : m_stiffness));
  Scalar energy = e.value();

  return energy;
}

void SolidSolidVertexJointCouplingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Vec3d E = defoObj().getVertexPosition(vh[4]);
  Vec3d F = defoObj().getVertexPosition(vh[5]);
  Vec3d G = defoObj().getVertexPosition(vh[6]);
  
  adreal<NumDof, 0, Scalar> e = adEnergy<0>(*this, A, B, C, D, E, F, G, s.delta, (viscous ? s.damping_undeformed_AP : s.undeformed_AP), (viscous ? s.damping_undeformed_delta : s.undeformed_delta), (viscous ? m_stiffness_damp : m_stiffness));
  for (int i = 0; i < NumDof; i++)
  {
    force[i] = -e.gradient(i);
  }
}

void SolidSolidVertexJointCouplingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Vec3d E = defoObj().getVertexPosition(vh[4]);
  Vec3d F = defoObj().getVertexPosition(vh[5]);
  Vec3d G = defoObj().getVertexPosition(vh[6]);
  
  adreal<NumDof, 1, Scalar> e = adEnergy<1>(*this, A, B, C, D, E, F, G, s.delta, (viscous ? s.damping_undeformed_AP : s.undeformed_AP), (viscous ? s.damping_undeformed_delta : s.undeformed_delta), (viscous ? m_stiffness_damp : m_stiffness));
  for (int i = 0; i < NumDof; i++)
    for (int j = 0; j < NumDof; j++)
    {
      jacobian(i, j) = -e.hessian(i, j);
    }
}

void SolidSolidVertexJointCouplingForce::updateStiffness()
{

}

void SolidSolidVertexJointCouplingForce::updateViscousReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.damping_undeformed_AP = s.AP;
    s.damping_undeformed_delta = s.delta;
  }
}

void SolidSolidVertexJointCouplingForce::updateProperties()
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
    Vec3d F = defoObj().getVertexPosition(vh[5]);
    Vec3d G = defoObj().getVertexPosition(vh[6]);
    Vec3d QA = (A * 3 - (E + F + G)).normalized();
    Vec3d m13 = QA;
    Vec3d m11 = ((F - E) - (F - E).dot(QA) * QA).normalized();
    Vec3d m12 = m13.cross(m11);
    Vec3d AP = ((B + C + D) - A * 3).normalized();
    s.AP = Vec3d(AP.dot(m11), AP.dot(m12), AP.dot(m13));
    Vec3d m23 = AP;
    Vec3d m21 = ((C - B) - (C - B).dot(AP) * AP).normalized();
    Vec3d m22 = m23.cross(m21);
    Vec3d m21r = m21 * cos(s.delta) + m22 * sin(s.delta);
    Vec3d m21r_in_m1 = parallel_transport(m21r, m23, m13);
    s.delta += atan2(m21r_in_m1.cross(m11).dot(m13), m21r_in_m1.dot(m11));
  }
}

void SolidSolidVertexJointCouplingForce::computeReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.undeformed_AP = s.AP;
    s.undeformed_delta = s.delta;
  }
}

void SolidSolidVertexJointCouplingForce::startStep(Scalar time, Scalar timestep)
{
  updateViscousReferenceStrain();
}

void SolidSolidVertexJointCouplingForce::endStep(Scalar time, Scalar timestep)
{
  updateStiffness();
}

void SolidSolidVertexJointCouplingForce::startIteration(Scalar time, Scalar timestep)
{
  
}

void SolidSolidVertexJointCouplingForce::endIteration(Scalar time, Scalar timestep)
{
  updateProperties();
}

