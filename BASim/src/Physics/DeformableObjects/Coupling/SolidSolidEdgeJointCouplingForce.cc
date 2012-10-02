#include "SolidSolidEdgeJointCouplingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"

using namespace BASim;

SolidSolidEdgeJointCouplingForce::SolidSolidEdgeJointCouplingForce(ElasticSolid & solid, const std::vector<Stencil> & stencils, Scalar stiffness, Scalar stiffness_damp, Scalar timestep) :
  DefoObjForce(solid.getDefoObj(), timestep, "SolidSolidEdgeJointCouplingForce"),
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
    s.undeformed_delta = 0;
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
        
    s.delta = 0;  // will be computed by updateProperties() below
    
    m_stencils.push_back(s);
  }

  updateProperties();
  updateStiffness();
  updateViscousReferenceStrain();
  computeReferenceStrain();
}

SolidSolidEdgeJointCouplingForce::~SolidSolidEdgeJointCouplingForce()
{
  
}

std::vector<VertexHandle> SolidSolidEdgeJointCouplingForce::getVertices(const Stencil & s)
{
  std::vector<VertexHandle> vh(6);
  TetVertexIterator tvit = defoObj().tv_iter(s.t1);
  VertexHandle vt11 = *tvit; ++tvit; assert(tvit);
  VertexHandle vt12 = *tvit; ++tvit; assert(tvit);
  VertexHandle vt13 = *tvit; ++tvit; assert(tvit);
  VertexHandle vt14 = *tvit; ++tvit; assert(!tvit);
  tvit = defoObj().tv_iter(s.t2);
  VertexHandle vt21 = *tvit; ++tvit; assert(tvit);
  VertexHandle vt22 = *tvit; ++tvit; assert(tvit);
  VertexHandle vt23 = *tvit; ++tvit; assert(tvit);
  VertexHandle vt24 = *tvit; ++tvit; assert(!tvit);
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

SolidSolidEdgeJointCouplingForce::Vector3d SolidSolidEdgeJointCouplingForce::vec2vector(const Vec3d & input)
{
  Vector3d output;
  output.x() = input.x();
  output.y() = input.y();
  output.z() = input.z();
  return output;
}

template <int DO_HESS>
adreal<SolidSolidEdgeJointCouplingForce::NumDof, DO_HESS, Scalar> 
SolidSolidEdgeJointCouplingForce::adEnergy(const SolidSolidEdgeJointCouplingForce & mn, const Vec3d & A, const Vec3d & B, const Vec3d & C, const Vec3d & D, const Vec3d & E, const Vec3d & F, Scalar delta, Scalar undeformed_delta, Scalar stiffness) 
{  
  // typedefs to simplify code below
  typedef adreal<SolidSolidEdgeJointCouplingForce::NumDof, DO_HESS, Scalar> adrealElast;
  typedef CVec3T<adrealElast> advecElast;
  Mat3T<adrealElast> temp;
  typedef Mat3T<adrealElast> admatElast;

  Vector3d vA = vec2vector(A);
  Vector3d vB = vec2vector(B);
  Vector3d vC = vec2vector(C);
  Vector3d vD = vec2vector(D);
  Vector3d vE = vec2vector(E);
  Vector3d vF = vec2vector(F);

  advecElast p[6];
  set_independent(p[0], vA, 0);
  set_independent(p[1], vB, 3);
  set_independent(p[2], vC, 6);
  set_independent(p[3], vD, 9);
  set_independent(p[4], vE, 12);
  set_independent(p[5], vF, 15);

  adrealElast e(0);
  
  advecElast vAB = (p[1] - p[0]); vAB /= len(vAB);
  advecElast vM = (p[0] + p[1]) * 0.5;
  advecElast vP = (p[2] + p[3]) * 0.5;
  advecElast vQ = (p[4] + p[5]) * 0.5;
  advecElast vMQ = ((vQ - vM) - dot(vQ - vM, vAB) * vAB);  vMQ /= len(vMQ);
  advecElast vMP = ((vP - vM) - dot(vP - vM, vAB) * vAB);  vMP /= len(vMP);
  advecElast vMQrot = vMQ * cos(delta) + cross(vAB, vMQ) * sin(delta);
  adrealElast newdelta = delta + atan2(dot(cross(vMQrot, vMP), vAB), dot(vMQrot, vMP));
  
  e = stiffness * (newdelta - undeformed_delta) * (newdelta - undeformed_delta);

  return e;
}

Scalar SolidSolidEdgeJointCouplingForce::globalEnergy()
{
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    // energy only include non-viscous forces
    energy += localEnergy(m_stencils[i], false);
  }
  return energy;
}

void SolidSolidEdgeJointCouplingForce::globalForce(VecXd & force)
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

void SolidSolidEdgeJointCouplingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
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

Scalar SolidSolidEdgeJointCouplingForce::localEnergy(Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Vec3d E = defoObj().getVertexPosition(vh[4]);
  Vec3d F = defoObj().getVertexPosition(vh[5]);
  Vec3d G = defoObj().getVertexPosition(vh[6]);
  
  adreal<NumDof, 1, Scalar> e = adEnergy<1>(*this, A, B, C, D, E, F, s.delta, (viscous ? s.damping_undeformed_delta : s.undeformed_delta), (viscous ? m_stiffness_damp : m_stiffness));
  Scalar energy = e.value();

  return energy;
}

void SolidSolidEdgeJointCouplingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Vec3d E = defoObj().getVertexPosition(vh[4]);
  Vec3d F = defoObj().getVertexPosition(vh[5]);
  Vec3d G = defoObj().getVertexPosition(vh[6]);
  
  adreal<NumDof, 1, Scalar> e = adEnergy<1>(*this, A, B, C, D, E, F, s.delta, (viscous ? s.damping_undeformed_delta : s.undeformed_delta), (viscous ? m_stiffness_damp : m_stiffness));
  for (int i = 0; i < NumDof; i++)
  {
    force[i] = -e.gradient(i);
  }
}

void SolidSolidEdgeJointCouplingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Vec3d E = defoObj().getVertexPosition(vh[4]);
  Vec3d F = defoObj().getVertexPosition(vh[5]);
  Vec3d G = defoObj().getVertexPosition(vh[6]);
  
  adreal<NumDof, 1, Scalar> e = adEnergy<1>(*this, A, B, C, D, E, F, s.delta, (viscous ? s.damping_undeformed_delta : s.undeformed_delta), (viscous ? m_stiffness_damp : m_stiffness));
  for (int i = 0; i < NumDof; i++)
    for (int j = 0; j < NumDof; j++)
    {
      jacobian(i, j) = -e.hessian(i, j);
    }
}

void SolidSolidEdgeJointCouplingForce::updateStiffness()
{

}

void SolidSolidEdgeJointCouplingForce::updateViscousReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.damping_undeformed_delta = s.delta;
  }
}

void SolidSolidEdgeJointCouplingForce::updateProperties()
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
    Vec3d AB = (B - A).normalized();
    Vec3d M = (A + B) * 0.5;
    Vec3d P = (C + D) * 0.5;
    Vec3d Q = (E + F) * 0.5;
    Vec3d MQ = ((Q - M) - (Q - M).dot(AB) * AB).normalized();
    Vec3d MP = ((P - M) - (P - M).dot(AB) * AB).normalized();
    Vec3d MQrot = MQ * cos(s.delta) + AB.cross(MQ) * sin(s.delta);
    s.delta += atan2(MQrot.cross(MP).dot(AB), MQrot.dot(MP));
  }
}

void SolidSolidEdgeJointCouplingForce::computeReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.undeformed_delta = s.delta;
  }
}

void SolidSolidEdgeJointCouplingForce::startStep(Scalar time, Scalar timestep)
{
  updateViscousReferenceStrain();
}

void SolidSolidEdgeJointCouplingForce::endStep(Scalar time, Scalar timestep)
{
  updateStiffness();
}

void SolidSolidEdgeJointCouplingForce::startIteration(Scalar time, Scalar timestep)
{
  
}

void SolidSolidEdgeJointCouplingForce::endIteration(Scalar time, Scalar timestep)
{
  updateProperties();
}

