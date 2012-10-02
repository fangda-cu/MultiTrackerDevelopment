#include "ShellShellVertexJointCouplingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"

using namespace BASim;

ShellShellVertexJointCouplingForce::ShellShellVertexJointCouplingForce(ElasticShell & shell, const std::vector<Stencil> & stencils, Scalar stiffness, Scalar stiffness_damp, Scalar timestep) :
  DefoObjForce(shell.getDefoObj(), timestep, "ShellShellVertexJointCouplingForce"),
  m_shell(&shell),
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
    s.damping_undeformed_AB.setZero();
    s.damping_undeformed_AC.setZero();
    
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
        
    s.AB.setZero();  // will be computed by updateProperties() below
    s.AC.setZero();
    
    m_stencils.push_back(s);
  }

  updateProperties();
  updateStiffness();
  updateViscousReferenceStrain();
  computeReferenceStrain();
}

ShellShellVertexJointCouplingForce::~ShellShellVertexJointCouplingForce()
{
  
}

std::vector<VertexHandle> ShellShellVertexJointCouplingForce::getVertices(const Stencil & s)
{
  std::vector<VertexHandle> vh(5);
  FaceVertexIterator fvit = defoObj().fv_iter(s.f1);
  VertexHandle vf11 = *fvit; ++fvit; assert(fvit);
  VertexHandle vf12 = *fvit; ++fvit; assert(fvit);
  VertexHandle vf13 = *fvit; ++fvit; assert(!fvit);
  fvit = defoObj().fv_iter(s.f2);
  VertexHandle vf21 = *fvit; ++fvit; assert(fvit);
  VertexHandle vf22 = *fvit; ++fvit; assert(fvit);
  VertexHandle vf23 = *fvit; ++fvit; assert(!fvit);
  if (vf11 == vf21 || vf11 == vf22 || vf11 == vf23) vh[0] = vf11, vh[3] = vf12, vh[4] = vf13;
  if (vf12 == vf21 || vf12 == vf22 || vf12 == vf23) vh[0] = vf12, vh[3] = vf13, vh[4] = vf11;
  if (vf13 == vf21 || vf13 == vf22 || vf13 == vf23) vh[0] = vf13, vh[3] = vf11, vh[4] = vf12;
  if (vh[0] == vf21) vh[1] = vf22, vh[2] = vf23;
  if (vh[0] == vf22) vh[1] = vf23, vh[2] = vf21;
  if (vh[0] == vf23) vh[1] = vf21, vh[2] = vf22;
  return vh;
}

ShellShellVertexJointCouplingForce::Vector3d ShellShellVertexJointCouplingForce::vec2vector(const Vec3d & input)
{
  Vector3d output;
  output.x() = input.x();
  output.y() = input.y();
  output.z() = input.z();
  return output;
}

template <int DO_HESS>
adreal<ShellShellVertexJointCouplingForce::NumDof, DO_HESS, Scalar> 
ShellShellVertexJointCouplingForce::adEnergy(const ShellShellVertexJointCouplingForce & mn, const Vec3d & A, const Vec3d & B, const Vec3d & C, const Vec3d & D, const Vec3d & E, const Vec3d & undeformed_AB, const Vec3d & undeformed_AC, Scalar stiffness) 
{  
  // typedefs to simplify code below
  typedef adreal<ShellShellVertexJointCouplingForce::NumDof, DO_HESS, Scalar> adrealElast;
  typedef CVec3T<adrealElast> advecElast;
  Mat3T<adrealElast> temp;
  typedef Mat3T<adrealElast> admatElast;

  Vector3d vA = vec2vector(A);
  Vector3d vB = vec2vector(B);
  Vector3d vC = vec2vector(C);
  Vector3d vD = vec2vector(D);
  Vector3d vE = vec2vector(E);
  Vector3d vUndeformedAB = vec2vector(undeformed_AB);
  Vector3d vUndeformedAC = vec2vector(undeformed_AC);

  advecElast p[6];
  set_independent(p[0], vA, 0);
  set_independent(p[1], vB, 3);
  set_independent(p[2], vC, 6);
  set_independent(p[3], vD, 9);
  set_independent(p[4], vE, 12);

  adrealElast e(0);
  
  advecElast md1 = normalized(p[3] + p[4] - p[0] * 2);
  advecElast md2 = normalized(p[3] - p[4] - dot(p[3] - p[4], md1) * md1);
  advecElast md3 = cross(md1, md2);
  
  advecElast vAB = advecElast(dot(p[1] - p[0], md1), dot(p[1] - p[0], md2), dot(p[1] - p[0], md3));
  advecElast vAC = advecElast(dot(p[2] - p[0], md1), dot(p[2] - p[0], md2), dot(p[2] - p[0], md3));

  e = stiffness * (dot(vAB - vUndeformedAB, vAB - vUndeformedAB) + dot(vAC - vUndeformedAC, vAC - vUndeformedAC));

  return e;
}

Scalar ShellShellVertexJointCouplingForce::globalEnergy()
{
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    // energy only include non-viscous forces
    energy += localEnergy(m_stencils[i], false);
  }
  return energy;
}

void ShellShellVertexJointCouplingForce::globalForce(VecXd & force)
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

void ShellShellVertexJointCouplingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
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

Scalar ShellShellVertexJointCouplingForce::localEnergy(Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Vec3d E = defoObj().getVertexPosition(vh[4]);
  
  adreal<NumDof, 1, Scalar> e = adEnergy<1>(*this, A, B, C, D, E, (viscous ? s.damping_undeformed_AB : s.undeformed_AB), (viscous ? s.damping_undeformed_AC : s.undeformed_AC), (viscous ? m_stiffness_damp : m_stiffness));
  Scalar energy = e.value();

  return energy;
}

void ShellShellVertexJointCouplingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Vec3d E = defoObj().getVertexPosition(vh[4]);
  
  adreal<NumDof, 1, Scalar> e = adEnergy<1>(*this, A, B, C, D, E, (viscous ? s.damping_undeformed_AB : s.undeformed_AB), (viscous ? s.damping_undeformed_AC : s.undeformed_AC), (viscous ? m_stiffness_damp : m_stiffness));
  for (int i = 0; i < NumDof; i++)
  {
    force[i] = -e.gradient(i);
  }
}

void ShellShellVertexJointCouplingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
  std::vector<VertexHandle> vh = getVertices(s);
  Vec3d A = defoObj().getVertexPosition(vh[0]);
  Vec3d B = defoObj().getVertexPosition(vh[1]);
  Vec3d C = defoObj().getVertexPosition(vh[2]);
  Vec3d D = defoObj().getVertexPosition(vh[3]);
  Vec3d E = defoObj().getVertexPosition(vh[4]);
  
  adreal<NumDof, 1, Scalar> e = adEnergy<1>(*this, A, B, C, D, E, (viscous ? s.damping_undeformed_AB : s.undeformed_AB), (viscous ? s.damping_undeformed_AC : s.undeformed_AC), (viscous ? m_stiffness_damp : m_stiffness));
  for (int i = 0; i < NumDof; i++)
    for (int j = 0; j < NumDof; j++)
    {
      jacobian(i, j) = -e.hessian(i, j);
    }
}

void ShellShellVertexJointCouplingForce::updateStiffness()
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

void ShellShellVertexJointCouplingForce::updateViscousReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.damping_undeformed_AB = s.AB;
    s.damping_undeformed_AC = s.AC;
  }
}

void ShellShellVertexJointCouplingForce::updateProperties()
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
    Vec3d md1 = (D + E - A * 2).normalized();
    Vec3d md2 = (D - E - (D - E).dot(md1) * md1).normalized();
    Vec3d md3 = md1.cross(md2);
    Vec3d AB = B - A;
    Vec3d AC = C - A;
    s.AB = Vec3d(AB.dot(md1), AB.dot(md2), AB.dot(md3));
    s.AC = Vec3d(AC.dot(md1), AC.dot(md2), AC.dot(md3));
  }
}

void ShellShellVertexJointCouplingForce::computeReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.undeformed_AB = s.AB;
    s.undeformed_AC = s.AC;
  }
}

void ShellShellVertexJointCouplingForce::startStep(Scalar time, Scalar timestep)
{
  updateViscousReferenceStrain();
}

void ShellShellVertexJointCouplingForce::endStep(Scalar time, Scalar timestep)
{
  updateStiffness();
}

void ShellShellVertexJointCouplingForce::startIteration(Scalar time, Scalar timestep)
{
  
}

void ShellShellVertexJointCouplingForce::endIteration(Scalar time, Scalar timestep)
{
  updateProperties();
}

