#include "RodModelBendingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"

#define EXACT_BENDING_JACOBIAN

using namespace BASim;

RodModelBendingForce::RodModelBendingForce(ElasticRodModel & rod, const std::vector<Stencil> & stencils, Scalar youngs_modulus, Scalar youngs_modulus_damping, Scalar timestep) :
  RodModelForce(rod, timestep, "RodModelStretchingForce"),
  m_stencils(stencils),
  m_youngs_modulus(youngs_modulus),
  m_youngs_modulus_damping(youngs_modulus_damping),
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
      force(m_stencils[i].dofindices[j]) += localforce(j) / timeStep();
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
    Jacobian.add(m_stencils[i].dofindices, m_stencils[i].dofindices, scale / timeStep() * localjacobian);
  }
}

Scalar RodModelBendingForce::localEnergy(Stencil & s, bool viscous)
{
  const Mat2d & B = (viscous ? m_viscous_stiffness[s.v2] : m_stiffness[s.v2]);
  const Scalar len = m_reference_voronoi_length[s.v2];
  const Vec2d & kappa = m_kappa[s.v2];
  const Vec2d & kappaBar = (viscous ? m_damping_undeformed_kappa[s.v2] : m_undeformed_kappa[s.v2]);
  
  return 0.5 / len * (kappa - kappaBar).dot(B * (kappa - kappaBar));
}

void RodModelBendingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
  const Mat2d & B = (viscous ? m_viscous_stiffness[s.v2] : m_stiffness[s.v2]);
  const Scalar len = m_reference_voronoi_length[s.v2];
  const Vec2d & kappa = m_kappa[s.v2];
  const Vec2d & kappaBar = (viscous ? m_damping_undeformed_kappa[s.v2] : m_undeformed_kappa[s.v2]);

  //TODO: cache grad kappa and hess kappa like BASim does
  ElementBiForce gradKappa = computeGradKappa(s);
  force = -1.0 / len * gradKappa * B * (kappa - kappaBar);
}

void RodModelBendingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
  const Mat2d & B = (viscous ? m_viscous_stiffness[s.v2] : m_stiffness[s.v2]);
  const Scalar len = m_reference_voronoi_length[s.v2];
  
  //TODO: cache grad kappa and hess kappa like BASim does
  ElementBiForce gradKappa = computeGradKappa(s);
  symBProduct<11>(jacobian, B, gradKappa);
  jacobian /= -len;
  
  // this part of the computation can be traded off for performance. refer to StrandSim::BendingForce::computeLocalJacobian()
#ifdef EXACT_BENDING_JACOBIAN
  const Vec2d & kappa = m_kappa[s.v2];
  const Vec2d & kappaBar = (viscous ? m_damping_undeformed_kappa[s.v2] : m_undeformed_kappa[s.v2]);

  ElementBiJacobian hessKappa = computeHessKappa(s);
  const Vec2d temp = -(kappa - kappaBar).transpose() * B / len;
  jacobian += temp(0) * hessKappa.first + temp(1) * hessKappa.second;
#endif
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

RodModelBendingForce::ElementBiForce RodModelBendingForce::computeGradKappa(Stencil & s)
{
  // copied from BASim::RodBendingForceSym::getGradKappa()
  ElementBiForce gradKappa = ElementBiForce::Zero();
  
  Scalar norm_e = rod().getEdgeLength(s.e1);
  Scalar norm_f = rod().getEdgeLength(s.e2);
  
  const Vec3d& te = rod().getEdgeTangent(s.e1);
  const Vec3d& tf = rod().getEdgeTangent(s.e2);
  
  const Vec3d& d1e = rod().getMaterialDirector1(s.e1);
  const Vec3d& d2e = rod().getMaterialDirector2(s.e1);
  const Vec3d& d1f = rod().getMaterialDirector1(s.e2);
  const Vec3d& d2f = rod().getMaterialDirector2(s.e2);
  
  Scalar chi = 1.0 + te.dot(tf);
  Vec3d tilde_t = (te + tf) / chi;
  Vec3d tilde_d1 = (d1e + d1f) / chi;
  Vec3d tilde_d2 = (d2e + d2f) / chi;
  
  const Vec2d& kappa = m_kappa[s.v2];
  Scalar kappa1 = kappa(0);
  Scalar kappa2 = kappa(1);
  
  Vec3d Dkappa1De = 1.0 / norm_e * (-kappa1 * tilde_t + tf.cross(tilde_d2));
  Vec3d Dkappa1Df = 1.0 / norm_f * (-kappa1 * tilde_t - te.cross(tilde_d2));
  Vec3d Dkappa2De = 1.0 / norm_e * (-kappa2 * tilde_t - tf.cross(tilde_d1));
  Vec3d Dkappa2Df = 1.0 / norm_f * (-kappa2 * tilde_t + te.cross(tilde_d1));
  
  gradKappa.block<3, 1> (0, 0) = -Dkappa1De;
  gradKappa.block<3, 1> (4, 0) = Dkappa1De - Dkappa1Df;
  gradKappa.block<3, 1> (8, 0) = Dkappa1Df;
  gradKappa.block<3, 1> (0, 1) = -Dkappa2De;
  gradKappa.block<3, 1> (4, 1) = Dkappa2De - Dkappa2Df;
  gradKappa.block<3, 1> (8, 1) = Dkappa2Df;
  
  const Vec3d& kb = rod().getCurvatureBinormal(s.v2);
  
  gradKappa(3, 0) = -0.5 * kb.dot(d1e);
  gradKappa(7, 0) = -0.5 * kb.dot(d1f);
  gradKappa(3, 1) = -0.5 * kb.dot(d2e);
  gradKappa(7, 1) = -0.5 * kb.dot(d2f);
  
  return gradKappa;
}

RodModelBendingForce::ElementBiJacobian RodModelBendingForce::computeHessKappa(Stencil & s)
{
  // copied from BASim::RodBendingForceSym::getHessKappa()
  ElementBiJacobian hessKappa;
  ElementJacobian& DDkappa1 = hessKappa.first;
  ElementJacobian& DDkappa2 = hessKappa.second;
  DDkappa1.setZero();
  DDkappa2.setZero();
  
  const Scalar norm_e = rod().getEdgeLength(s.e1);
  const Scalar norm_f = rod().getEdgeLength(s.e2);
  
  const Scalar norm2_e = square(norm_e);
  const Scalar norm2_f = square(norm_f);
  
  const Vec3d& te = rod().getEdgeTangent(s.e1);
  const Vec3d& tf = rod().getEdgeTangent(s.e2);
  
  const Vec3d& d1e = rod().getMaterialDirector1(s.e1);
  const Vec3d& d2e = rod().getMaterialDirector2(s.e1);
  const Vec3d& d1f = rod().getMaterialDirector1(s.e2);
  const Vec3d& d2f = rod().getMaterialDirector2(s.e2);
  
  const Scalar chi = 1.0 + te.dot(tf);
  const Vec3d tilde_t = (te + tf) / chi;
  const Vec3d tilde_d1 = (d1e + d1f) / chi;
  const Vec3d tilde_d2 = (d2e + d2f) / chi;
  
  const Vec2d& kappa = m_kappa[s.v2];
  const Scalar kappa1 = kappa(0);
  const Scalar kappa2 = kappa(1);
  
  const Vec3d& kb = rod().getCurvatureBinormal(s.v2);
  
  const Mat3d tt_o_tt = outerProd(tilde_t, tilde_t);
  const Mat3d tf_c_d2t_o_tt = outerProd(tf.cross(tilde_d2), tilde_t);
  const Mat3d tt_o_tf_c_d2t = tf_c_d2t_o_tt.transpose();
  const Mat3d kb_o_d2e = outerProd(kb, d2e);
  const Mat3d d2e_o_kb = kb_o_d2e.transpose();
  
  const Mat3d Id = Mat3d::Identity();
  
  const Mat3d D2kappa1De2 = 1.0 / norm2_e * (2 * kappa1 * tt_o_tt - (tf_c_d2t_o_tt + tt_o_tf_c_d2t)) - kappa1 / (chi
                            * norm2_e) * (Id - outerProd(te, te)) + 1.0 / (4.0 * norm2_e) * (kb_o_d2e + d2e_o_kb);
  
  const Mat3d te_c_d2t_o_tt = outerProd(te.cross(tilde_d2), tilde_t);
  const Mat3d tt_o_te_c_d2t = te_c_d2t_o_tt.transpose();
  const Mat3d kb_o_d2f = outerProd(kb, d2f);
  const Mat3d d2f_o_kb = kb_o_d2f.transpose();
  
  const Mat3d D2kappa1Df2 = 1.0 / norm2_f * (2 * kappa1 * tt_o_tt + (te_c_d2t_o_tt + tt_o_te_c_d2t)) - kappa1 / (chi
                            * norm2_f) * (Id - outerProd(tf, tf)) + 1.0 / (4.0 * norm2_f) * (kb_o_d2f + d2f_o_kb);
  
  const Mat3d D2kappa1DeDf = -kappa1 / (chi * norm_e * norm_f) * (Id + outerProd(te, tf)) + 1.0 / (norm_e * norm_f) * (2
                             * kappa1 * tt_o_tt - tf_c_d2t_o_tt + tt_o_te_c_d2t - crossMat(tilde_d2));
  const Mat3d D2kappa1DfDe = D2kappa1DeDf.transpose();
  
  const Mat3d tf_c_d1t_o_tt = outerProd(tf.cross(tilde_d1), tilde_t);
  const Mat3d tt_o_tf_c_d1t = tf_c_d1t_o_tt.transpose();
  const Mat3d kb_o_d1e = outerProd(kb, d1e);
  const Mat3d d1e_o_kb = kb_o_d1e.transpose();
  
  const Mat3d D2kappa2De2 = 1.0 / norm2_e * (2 * kappa2 * tt_o_tt + (tf_c_d1t_o_tt + tt_o_tf_c_d1t)) - kappa2 / (chi
                            * norm2_e) * (Id - outerProd(te, te)) - 1.0 / (4.0 * norm2_e) * (kb_o_d1e + d1e_o_kb);
  
  const Mat3d te_c_d1t_o_tt = outerProd(te.cross(tilde_d1), tilde_t);
  const Mat3d tt_o_te_c_d1t = te_c_d1t_o_tt.transpose();
  const Mat3d kb_o_d1f = outerProd(kb, d1f);
  const Mat3d d1f_o_kb = kb_o_d1f.transpose();
  
  const Mat3d D2kappa2Df2 = 1.0 / norm2_f * (2 * kappa2 * tt_o_tt - (te_c_d1t_o_tt + tt_o_te_c_d1t)) - kappa2 / (chi
                            * norm2_f) * (Id - outerProd(tf, tf)) - 1.0 / (4.0 * norm2_f) * (kb_o_d1f + d1f_o_kb);
  
  const Mat3d D2kappa2DeDf = -kappa2 / (chi * norm_e * norm_f) * (Id + outerProd(te, tf)) + 1.0 / (norm_e * norm_f) * (2
                             * kappa2 * tt_o_tt + tf_c_d1t_o_tt - tt_o_te_c_d1t + crossMat(tilde_d1));
  const Mat3d D2kappa2DfDe = D2kappa2DeDf.transpose();
  
  const Scalar D2kappa1Dthetae2 = -0.5 * kb.dot(d2e);
  const Scalar D2kappa1Dthetaf2 = -0.5 * kb.dot(d2f);
  const Scalar D2kappa2Dthetae2 = 0.5 * kb.dot(d1e);
  const Scalar D2kappa2Dthetaf2 = 0.5 * kb.dot(d1f);
  
  const Vec3d D2kappa1DeDthetae = 1.0 / norm_e * (0.5 * kb.dot(d1e) * tilde_t - 1.0 / chi * tf.cross(d1e));
  const Vec3d D2kappa1DeDthetaf = 1.0 / norm_e * (0.5 * kb.dot(d1f) * tilde_t - 1.0 / chi * tf.cross(d1f));
  const Vec3d D2kappa1DfDthetae = 1.0 / norm_f * (0.5 * kb.dot(d1e) * tilde_t + 1.0 / chi * te.cross(d1e));
  const Vec3d D2kappa1DfDthetaf = 1.0 / norm_f * (0.5 * kb.dot(d1f) * tilde_t + 1.0 / chi * te.cross(d1f));
  const Vec3d D2kappa2DeDthetae = 1.0 / norm_e * (0.5 * kb.dot(d2e) * tilde_t - 1.0 / chi * tf.cross(d2e));
  const Vec3d D2kappa2DeDthetaf = 1.0 / norm_e * (0.5 * kb.dot(d2f) * tilde_t - 1.0 / chi * tf.cross(d2f));
  const Vec3d D2kappa2DfDthetae = 1.0 / norm_f * (0.5 * kb.dot(d2e) * tilde_t + 1.0 / chi * te.cross(d2e));
  const Vec3d D2kappa2DfDthetaf = 1.0 / norm_f * (0.5 * kb.dot(d2f) * tilde_t + 1.0 / chi * te.cross(d2f));
  
  DDkappa1.block<3, 3> (0, 0) = D2kappa1De2;
  DDkappa1.block<3, 3> (0, 4) = -D2kappa1De2 + D2kappa1DeDf;
  DDkappa1.block<3, 3> (4, 0) = -D2kappa1De2 + D2kappa1DfDe;
  DDkappa1.block<3, 3> (4, 4) = D2kappa1De2 - (D2kappa1DeDf + D2kappa1DfDe) + D2kappa1Df2;
  DDkappa1.block<3, 3> (0, 8) = -D2kappa1DeDf;
  DDkappa1.block<3, 3> (8, 0) = -D2kappa1DfDe;
  DDkappa1.block<3, 3> (4, 8) = D2kappa1DeDf - D2kappa1Df2;
  DDkappa1.block<3, 3> (8, 4) = D2kappa1DfDe - D2kappa1Df2;
  DDkappa1.block<3, 3> (8, 8) = D2kappa1Df2;
  DDkappa1(3, 3) = D2kappa1Dthetae2;
  DDkappa1(7, 7) = D2kappa1Dthetaf2;
  DDkappa1.block<3, 1> (0, 3) = -D2kappa1DeDthetae;
  DDkappa1.block<1, 3> (3, 0) = DDkappa1.block<3, 1> (0, 3).transpose();
  DDkappa1.block<3, 1> (4, 3) = D2kappa1DeDthetae - D2kappa1DfDthetae;
  DDkappa1.block<1, 3> (3, 4) = DDkappa1.block<3, 1> (4, 3).transpose();
  DDkappa1.block<3, 1> (8, 3) = D2kappa1DfDthetae;
  DDkappa1.block<1, 3> (3, 8) = DDkappa1.block<3, 1> (8, 3).transpose();
  DDkappa1.block<3, 1> (0, 7) = -D2kappa1DeDthetaf;
  DDkappa1.block<1, 3> (7, 0) = DDkappa1.block<3, 1> (0, 7).transpose();
  DDkappa1.block<3, 1> (4, 7) = D2kappa1DeDthetaf - D2kappa1DfDthetaf;
  DDkappa1.block<1, 3> (7, 4) = DDkappa1.block<3, 1> (4, 7).transpose();
  DDkappa1.block<3, 1> (8, 7) = D2kappa1DfDthetaf;
  DDkappa1.block<1, 3> (7, 8) = DDkappa1.block<3, 1> (8, 7).transpose();
  
  assert(isSymmetric(DDkappa1));
  
  DDkappa2.block<3, 3> (0, 0) = D2kappa2De2;
  DDkappa2.block<3, 3> (0, 4) = -D2kappa2De2 + D2kappa2DeDf;
  DDkappa2.block<3, 3> (4, 0) = -D2kappa2De2 + D2kappa2DfDe;
  DDkappa2.block<3, 3> (4, 4) = D2kappa2De2 - (D2kappa2DeDf + D2kappa2DfDe) + D2kappa2Df2;
  DDkappa2.block<3, 3> (0, 8) = -D2kappa2DeDf;
  DDkappa2.block<3, 3> (8, 0) = -D2kappa2DfDe;
  DDkappa2.block<3, 3> (4, 8) = D2kappa2DeDf - D2kappa2Df2;
  DDkappa2.block<3, 3> (8, 4) = D2kappa2DfDe - D2kappa2Df2;
  DDkappa2.block<3, 3> (8, 8) = D2kappa2Df2;
  DDkappa2(3, 3) = D2kappa2Dthetae2;
  DDkappa2(7, 7) = D2kappa2Dthetaf2;
  DDkappa2.block<3, 1> (0, 3) = -D2kappa2DeDthetae;
  DDkappa2.block<1, 3> (3, 0) = DDkappa2.block<3, 1> (0, 3).transpose();
  DDkappa2.block<3, 1> (4, 3) = D2kappa2DeDthetae - D2kappa2DfDthetae;
  DDkappa2.block<1, 3> (3, 4) = DDkappa2.block<3, 1> (4, 3).transpose();
  DDkappa2.block<3, 1> (8, 3) = D2kappa2DfDthetae;
  DDkappa2.block<1, 3> (3, 8) = DDkappa2.block<3, 1> (8, 3).transpose();
  DDkappa2.block<3, 1> (0, 7) = -D2kappa2DeDthetaf;
  DDkappa2.block<1, 3> (7, 0) = DDkappa2.block<3, 1> (0, 7).transpose();
  DDkappa2.block<3, 1> (4, 7) = D2kappa2DeDthetaf - D2kappa2DfDthetaf;
  DDkappa2.block<1, 3> (7, 4) = DDkappa2.block<3, 1> (4, 7).transpose();
  DDkappa2.block<3, 1> (8, 7) = D2kappa2DfDthetaf;
  DDkappa2.block<1, 3> (7, 8) = DDkappa2.block<3, 1> (8, 7).transpose();
  
  assert(isSymmetric(DDkappa2));
  
  return hessKappa;
}


