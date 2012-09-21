#include "RodTwistSolidTetCouplingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"

using namespace BASim;

RodTwistSolidTetCouplingForce::RodTwistSolidTetCouplingForce(ElasticRodModel & rod, ElasticSolid & solid, const std::vector<Stencil> & stencils, Scalar stiffness, Scalar stiffness_damp, Scalar timestep) :
  DefoObjForce(rod.getDefoObj(), timestep, "RodTwistSolidTetCouplingForce"),
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
    s.undeformed_delta = 0;
    s.damping_undeformed_delta = 0;
    
    s.dofindices.resize(NumDof);
    
    int dofbase;
    dofbase = defoObj().getPositionDofBase(s.v1);
    s.dofindices[0] = dofbase;
    s.dofindices[1] = dofbase + 1;
    s.dofindices[2] = dofbase + 2;
    dofbase = defoObj().getPositionDofBase(s.v2);
    s.dofindices[3] = dofbase;
    s.dofindices[4] = dofbase + 1;
    s.dofindices[5] = dofbase + 2;
    dofbase = defoObj().getPositionDofBase(defoObj().fromVertex(s.e));
    s.dofindices[6] = dofbase;
    s.dofindices[7] = dofbase + 1;
    s.dofindices[8] = dofbase + 2;
    dofbase = defoObj().getPositionDofBase(defoObj().toVertex(s.e));
    s.dofindices[9] = dofbase;
    s.dofindices[10] = dofbase + 1;
    s.dofindices[11] = dofbase + 2;
    
    dofbase = rod.getEdgeDofBase(s.e);
    s.dofindices[12] = dofbase;
    
    s.delta = 0;  // will be computed by updateProperties() below
    
    m_stencils.push_back(s);
  }

  updateProperties();
  updateStiffness();
  updateViscousReferenceStrain();
  computeReferenceStrain();
}

RodTwistSolidTetCouplingForce::~RodTwistSolidTetCouplingForce()
{
  
}

RodTwistSolidTetCouplingForce::Vector3d RodTwistSolidTetCouplingForce::vec2vector(const Vec3d & input)
{
  Vector3d output;
  output.x() = input.x();
  output.y() = input.y();
  output.z() = input.z();
  return output;
}

template <int DO_HESS>
adreal<RodTwistSolidTetCouplingForce::NumDof, DO_HESS, Scalar> 
RodTwistSolidTetCouplingForce::adEnergy(const RodTwistSolidTetCouplingForce & mn, const Vec3d & A, const Vec3d & B, const Vec3d & C, const Vec3d & D, Scalar theta, Scalar delta, const Vec3d & ref1, const Vec3d & ref2, Scalar undeformed_delta, Scalar stiffness) 
{  
  // typedefs to simplify code below
  typedef adreal<RodTwistSolidTetCouplingForce::NumDof, DO_HESS, Scalar> adrealElast;
  typedef CVec3T<adrealElast> advecElast;
  Mat3T<adrealElast> temp;
  typedef Mat3T<adrealElast> admatElast;

  Vector3d vA = vec2vector(A);
  Vector3d vB = vec2vector(B);
  Vector3d vC = vec2vector(C);
  Vector3d vD = vec2vector(D);
  Vector3d vRef1 = vec2vector(ref1);
  Vector3d vRef2 = vec2vector(ref2);
  
  advecElast p[4];
  adrealElast t;
  set_independent(p[0], vA, 0);
  set_independent(p[1], vB, 3);
  set_independent(p[2], vC, 6);
  set_independent(p[3], vD, 9);
  t.set_independent(theta, 12);

  adrealElast e(0);
  
  adrealElast oldvec_x = cos(delta + t);    // direction of the old delta, in ref frame
  adrealElast oldvec_y = sin(delta + t);
  adrealElast newvec_x = dot((p[0] + p[1]) * 0.5 - p[2], vRef1);   // projection of A-B into the frame plane, in ref frame
  adrealElast newvec_y = dot((p[0] + p[1]) * 0.5 - p[2], vRef2);
  adrealElast newdelta = delta + atan2(oldvec_x * newvec_y - oldvec_y * newvec_x, oldvec_x * newvec_x + oldvec_y * newvec_y);
  
  e = stiffness * ((newdelta - undeformed_delta) * (newdelta - undeformed_delta));

  return e;
}

Scalar RodTwistSolidTetCouplingForce::globalEnergy()
{
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    // energy only include non-viscous forces
    energy += localEnergy(m_stencils[i], false);
  }
  return energy;
}

void RodTwistSolidTetCouplingForce::globalForce(VecXd & force)
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

void RodTwistSolidTetCouplingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
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

Scalar RodTwistSolidTetCouplingForce::localEnergy(Stencil & s, bool viscous)
{
  Vec3d A = defoObj().getVertexPosition(s.v1);
  Vec3d B = defoObj().getVertexPosition(s.v2);
  Vec3d C = defoObj().getVertexPosition(defoObj().fromVertex(s.e));
  Vec3d D = defoObj().getVertexPosition(defoObj().toVertex(s.e));
  Scalar theta = rod().getEdgeTheta(s.e);
  Scalar delta = s.delta;
  
  Vec3d & ref1 = rod().getReferenceDirector1(s.e);
  Vec3d & ref2 = rod().getReferenceDirector2(s.e);
  
  adreal<NumDof, 0, Scalar> e = adEnergy<0>(*this, A, B, C, D, theta, delta, ref1, ref2, (viscous ? s.damping_undeformed_delta : s.undeformed_delta), (viscous ? m_stiffness_damp : m_stiffness));
  Scalar energy = e.value();

  return energy;
}

void RodTwistSolidTetCouplingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
  Vec3d A = defoObj().getVertexPosition(s.v1);
  Vec3d B = defoObj().getVertexPosition(s.v2);
  Vec3d C = defoObj().getVertexPosition(defoObj().fromVertex(s.e));
  Vec3d D = defoObj().getVertexPosition(defoObj().toVertex(s.e));
  Scalar theta = rod().getEdgeTheta(s.e);
  Scalar delta = s.delta;
  
  Vec3d & ref1 = rod().getReferenceDirector1(s.e);
  Vec3d & ref2 = rod().getReferenceDirector2(s.e);
  
  adreal<NumDof, 0, Scalar> e = adEnergy<0>(*this, A, B, C, D, theta, delta, ref1, ref2, (viscous ? s.damping_undeformed_delta : s.undeformed_delta), (viscous ? m_stiffness_damp : m_stiffness));
  for (int i = 0; i < NumDof; i++)
  {
    force[i] = -e.gradient(i);
  }
}

void RodTwistSolidTetCouplingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
  Vec3d A = defoObj().getVertexPosition(s.v1);
  Vec3d B = defoObj().getVertexPosition(s.v2);
  Vec3d C = defoObj().getVertexPosition(defoObj().fromVertex(s.e));
  Vec3d D = defoObj().getVertexPosition(defoObj().toVertex(s.e));
  Scalar theta = rod().getEdgeTheta(s.e);
  Scalar delta = s.delta;
  
  Vec3d & ref1 = rod().getReferenceDirector1(s.e);
  Vec3d & ref2 = rod().getReferenceDirector2(s.e);
  
  adreal<NumDof, 1, Scalar> e = adEnergy<1>(*this, A, B, C, D, theta, delta, ref1, ref2, (viscous ? s.damping_undeformed_delta : s.undeformed_delta), (viscous ? m_stiffness_damp : m_stiffness));
  for (int i = 0; i < NumDof; i++)
    for (int j = 0; j < NumDof; j++)
    {
      jacobian(i, j) = -e.hessian(i, j);
    }
}

void RodTwistSolidTetCouplingForce::updateStiffness()
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

void RodTwistSolidTetCouplingForce::updateViscousReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.damping_undeformed_delta = s.delta;
  }
}

void RodTwistSolidTetCouplingForce::updateProperties()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    
    Vec3d A = defoObj().getVertexPosition(s.v1);
    Vec3d B = defoObj().getVertexPosition(s.v2);
    Vec3d C = defoObj().getVertexPosition(defoObj().fromVertex(s.e));
    Vec3d D = defoObj().getVertexPosition(defoObj().toVertex(s.e));
    Vec3d md1 = rod().getMaterialDirector1(s.e);
    Vec3d md2 = rod().getMaterialDirector2(s.e);
    Vec2d oldvec = Vec2d(cos(s.delta), sin(s.delta)); // direction of the old delta, in mat frame
    Vec2d newvec = Vec2d(((A + B) * 0.5 - C).dot(md1), ((A + B) * 0.5 - C).dot(md2)); // projection of A-B into the frame plane, in mat frame
    s.delta += atan2((oldvec.x() * newvec.y() - oldvec.y() * newvec.x()), oldvec.dot(newvec));
  }
}

void RodTwistSolidTetCouplingForce::computeReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.undeformed_delta = s.delta;
  }
}

void RodTwistSolidTetCouplingForce::startStep(Scalar time, Scalar timestep)
{
  updateViscousReferenceStrain();
}

void RodTwistSolidTetCouplingForce::endStep(Scalar time, Scalar timestep)
{
  updateStiffness();
}

void RodTwistSolidTetCouplingForce::startIteration(Scalar time, Scalar timestep)
{
  
}

void RodTwistSolidTetCouplingForce::endIteration(Scalar time, Scalar timestep)
{
  updateProperties();
}



