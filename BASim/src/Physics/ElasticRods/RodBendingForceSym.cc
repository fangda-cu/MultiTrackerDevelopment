/**
 * \file RodBendingForceSym.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 02/18/2010
 */

#ifdef WETA
#include "RodBendingForceSym.hh"
#include "ElasticRod.hh"
#include "../../Math/Math.hh"
#else
#include "BASim/src/Physics/ElasticRods/RodBendingForceSym.hh"
#include "BASim/src/Math/Math.hh"
#endif

using namespace std;

namespace BASim {

RodBendingForceSym::RodBendingForceSym(ElasticRod& rod, bool vscs, bool runinit)
  : RodForceT<VertexStencil>(rod, "RodBendingForceSym")
{
  if( runinit )
  {
    m_rod.add_property(m_kappa, "material curvature vector");
    m_rod.add_property(m_gradKappaValid, "grad kappa valid", false);
    m_rod.add_property(m_hessKappaValid, "hess kappa valid", false);
    m_rod.add_property(m_gradKappa, "gradient of material curvature vector",MatXd(11,2));
    m_rod.add_property(m_hessKappa, "Hessian of material curvature vector",make_pair<MatXd,MatXd>(MatXd(11,11), MatXd(11,11)));

    if( !vscs ) 
    {
      m_rod.add_property(m_B, "bending metric");
      m_rod.add_property(m_refVertexLength, "Voronoi length of vertex");
      m_rod.add_property(m_kappaBar, "undeformed material curvature vector");
    }
    else 
    {
      m_rod.add_property(m_B, "viscous bending metric");
      m_rod.add_property(m_refVertexLength, "viscous Voronoi length of vertex");
      m_rod.add_property(m_kappaBar, "viscous undeformed material curvature vector");  
    }

    setViscous(vscs);

    updateProperties();
    updateUndeformedStrain();
    updateStiffness();
    updateReferenceDomain();
  }
  else
  {
    m_viscous = vscs;
  }
}

void RodBendingForceSym::updateUndeformedConfiguration(std::vector<Scalar>& vals) 
{
  int i=0;
  
//  std::cout << "kappa\n";
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil, ++i) {
    Vec2d newKappaBar = Vec2d(vals[i*2], vals[i*2+1]);
    //std::cout << newKappaBar <<"\n";
    vertex_handle& vh = m_stencil.handle();
    setKappaBar(vh, newKappaBar);
  }
}

void RodBendingForceSym::reattatchProperties()
{
  m_rod.add_property(m_kappa, "material curvature vector");
  m_rod.add_property(m_gradKappaValid, "grad kappa valid", false);
  m_rod.add_property(m_hessKappaValid, "hess kappa valid", false);
  m_rod.add_property(m_gradKappa, "gradient of material curvature vector",MatXd(11,2));
  m_rod.add_property(m_hessKappa, "Hessian of material curvature vector",make_pair<MatXd,MatXd>(MatXd(11,11), MatXd(11,11)));
  
  if( !m_viscous ) 
  {
    m_rod.add_property(m_B, "bending metric");
    m_rod.add_property(m_refVertexLength, "Voronoi length of vertex");
    m_rod.add_property(m_kappaBar, "undeformed material curvature vector");
  }
  else 
  {
    m_rod.add_property(m_B, "viscous bending metric");
    m_rod.add_property(m_refVertexLength, "viscous Voronoi length of vertex");
    m_rod.add_property(m_kappaBar, "viscous undeformed material curvature vector");  
  }  
}


void RodBendingForceSym::updateProperties()
{
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    const Vec3d& kb = m_rod.getCurvatureBinormal(vh);
    int i = vh.idx();
    const Vec3d& m1e = m_rod.getMaterial1(i - 1);
    const Vec3d& m2e = m_rod.getMaterial2(i - 1);
    const Vec3d& m1f = m_rod.getMaterial1(i);
    const Vec3d& m2f = m_rod.getMaterial2(i);
    setKappa(vh, Vec2d(0.5 * kb.dot(m2e + m2f), -0.5 * kb.dot(m1e + m1f)));
  }

  m_rod.property(m_gradKappaValid) = false;
  m_rod.property(m_hessKappaValid) = false;
}

void RodBendingForceSym::updateStiffness()
{
  Scalar E = m_rod.getYoungsModulus();
  if (viscous()) {
    E = 3 * m_rod.getViscosity() / m_rod.getTimeStep();
  }

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    edge_handle eh0 = m_stencil.inEdge();
    edge_handle eh1 = m_stencil.outEdge();
    Scalar a = (m_rod.radiusA(eh0) + m_rod.radiusA(eh1)) / 2.0;
    Scalar b = (m_rod.radiusB(eh0) + m_rod.radiusB(eh1)) / 2.0;
    Mat2d B(Mat2d::Zero());
    B(0, 0) = E * M_PI * cube(a) * b / 4.0;
    B(1, 1) = E * M_PI * a * cube(b) / 4.0;
    
    // rotate cross section
    Mat2d rot(Mat2d::Zero());
    rot(0,0) = cos(m_rod.baseRotation());
    rot(1,0) = sin(m_rod.baseRotation());
    rot(0,1) = -1*rot(1,0);
    rot(1,1) = rot(0,0);
    
    B = rot*B;
    B = B*rot.transpose();

    setB(vh, B);
  }
}

void RodBendingForceSym::updateUndeformedStrain()
{
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    setKappaBar(vh, getKappa(vh));
  }
}

void RodBendingForceSym::updateReferenceDomain()
{
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    edge_handle eh0 = m_stencil.inEdge();
    edge_handle eh1 = m_stencil.outEdge();
    Scalar len = (m_rod.getEdgeLength(eh0) + m_rod.getEdgeLength(eh1)) / 2.0;
    setRefVertexLength(vh, len);
  }
}

const Mat2d& RodBendingForceSym::getB(const vertex_handle& vh) const
{
  return m_rod.property(m_B)[vh];
}

void RodBendingForceSym::setB(const vertex_handle& vh, const Mat2d& B)
{
  m_rod.property(m_B)[vh] = B;
}

const Vec2d& RodBendingForceSym::getKappa(const vertex_handle& vh) const
{
  return m_rod.property(m_kappa)[vh];
}

void RodBendingForceSym::setKappa(const vertex_handle& vh, const Vec2d& kappa)
{
  m_rod.property(m_kappa)[vh] = kappa;
}

const Vec2d& RodBendingForceSym::getKappaBar(const vertex_handle& vh) const
{
  return m_rod.property(m_kappaBar)[vh];
}

void RodBendingForceSym::setKappaBar(const vertex_handle& vh,
                                     const Vec2d& kappaBar)
{
  m_rod.property(m_kappaBar)[vh] = kappaBar;
}

Scalar RodBendingForceSym::globalEnergy()
{
  if (viscous() && m_rod.getViscosity() == 0.0) return 0.0;

  Scalar energy = 0;
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    energy += localEnergy(vh);
  }
  return energy;
}

Scalar RodBendingForceSym::localEnergy(const vertex_handle& vh)
{
  // Unused? int i = vh.idx();
  Mat2d B = getB(vh);
  Scalar len = getRefVertexLength(vh);

  const Vec2d& kappa = getKappa(vh);
  const Vec2d& kappaBar = getKappaBar(vh);

  return 0.5 / len * (kappa - kappaBar).dot(B * (kappa - kappaBar));
}

void RodBendingForceSym::globalForce(VecXd& force)
{
  if (viscous() && m_rod.getViscosity() == 0.0) return;

  //VecXd force1 = force;

  computeGradKappa();

  VecXd f(11);
  IndexArray indices(11);
  // Unused? unsigned int nv = m_rod.nv();

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    f.setZero();
    localForce(f, vh);
    m_stencil.indices(indices);
    for (int j = 0; j < f.size(); ++j) {
      force(indices[j]) += f(j);
    }
  }

  //std::cout << (viscous() ? "VISCOUS " : "") << "BENDING FORCE\n";
  //std::cout << "norm = " << (force - force1).norm() << " force = " << force-force1 << "\n\n";
}

void RodBendingForceSym::localForce(VecXd& force, const vertex_handle& vh)
{
  // Unused? int i = vh.idx();
  Mat2d B = getB(vh);
  Scalar len = getRefVertexLength(vh);

  const Vec2d& kappa = getKappa(vh);
  const Vec2d& kappaBar = getKappaBar(vh);

  force = -1.0/len * getGradKappa(vh) * B * (kappa - kappaBar);
}

void RodBendingForceSym::globalJacobian(int baseidx, Scalar scale, MatrixBase& J)
{
  if (viscous() && m_rod.getViscosity() == 0.0) return;

  computeGradKappa();
  computeHessKappa();

  MatXd localJ(11, 11);
  IndexArray indices(11);
  iterator end = m_stencil.end();

  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    //localJacobian(localJ, i);

    {
      const Mat2d& B = getB(vh);
      Scalar len = getRefVertexLength(vh);

      const Vec2d& kappa = getKappa(vh);
      const Vec2d& kappaBar = getKappaBar(vh);
      const MatXd& gradKappa = getGradKappa(vh);

      localJ = -1.0 / len * gradKappa * B * gradKappa.transpose();

      const pair<MatXd, MatXd>& hessKappa = getHessKappa(vh);
      Vec2d temp = -1.0 / len * (kappa - kappaBar).transpose() * B;
      localJ += temp(0) * hessKappa.first + temp(1) * hessKappa.second;
    }

    m_stencil.indices(indices);
    for( int i = 0; i < (int) indices.size(); ++i ) indices[i] += baseidx;
    localJ *= scale;
    J.add(indices, indices, localJ);
  }
}

inline void RodBendingForceSym::localJacobian(MatXd& localJ,
                                              const vertex_handle& vh)
{
  Mat2d B = getB(vh);
  Scalar len = getRefVertexLength(vh);

  const Vec2d& kappa = getKappa(vh);
  const Vec2d& kappaBar = getKappaBar(vh);
  const MatXd& gradKappa = getGradKappa(vh);

  localJ = -1.0 / len * gradKappa * B * gradKappa.transpose();

  const pair<MatXd, MatXd>& hessKappa = getHessKappa(vh);
  Vec2d temp = -1.0 / len * (kappa - kappaBar).transpose() * B;
  localJ += temp(0) * hessKappa.first + temp(1) * hessKappa.second;
}

void RodBendingForceSym::globalReverseJacobian(MatrixBase& J) 
{
  if (viscous()) return;

  //return;
  
  computeGradKappa();
  computeHessKappa();

  //MatXd localJ(11, 11);
  IndexArray indices(11);
  iterator end = m_stencil.end();

  uint vid = 1;
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil, ++vid) {
    vertex_handle& vh = m_stencil.handle();
    
    Mat2d B = getB(vh);
    Scalar refLen = getRefVertexLength(vh);
    const Vec2d& kappa = getKappa(vh);
    const Vec2d& kappaBar = getKappaBar(vh);

    VecXd fde(11);
    MatXd fdk(11, 2);
    
    fde = 0.5 / (refLen * refLen) * getGradKappa(vh) * B * (kappa - kappaBar);
    fdk = 1.0 / refLen * getGradKappa(vh) * B; 

    m_stencil.indices(indices);
    
    for(int i=0; i<indices.size(); i++) indices[i] -= 7;
    
//    std::cout << "bend id  " << indices << "\n";
    
    // edge length bar
    if (vid == 1) {
      for(int i=0; i<11; i++) {
        if (indices[i] >= 0)
          J.add(indices[i], (vid - 1) * 4 + 3, fde(i));
      }
    } else {
      for(int i=0; i<11; i++) {
        if (indices[i] >= 0) {
          J.add(indices[i], (vid - 2) * 4 + 3, fde(i));
          J.add(indices[i], (vid - 1) * 4 + 3, fde(i));
        }
      }
    }
    
    // kappa bar
    for(int i=0; i<11; i++) {
      for(int j=0; j<2; j++) {
        if (indices[i] >= 0) 
          J.add(indices[i], (vid - 1) * 4 + j, fdk(i, j));
      }
    }
    
    m_stencil.indices(indices);
    
//    std::cout << "bend id  recover? " << indices << "\n";
    
  }
}

void RodBendingForceSym::updateReverseUndeformedStrain(const VecXd& e)
{
  if (viscous()) return; 
    
  m_stencil = m_stencil.begin();
  iterator end = m_stencil.end();
  
  uint vid = 1;  // eh.id() ?
  
  for (; m_stencil != end; ++m_stencil, ++vid) {
    vertex_handle& vh = m_stencil.handle();
    Vec2d newKappaBar = Vec2d(e( (vid-1) * 4 ), e( (vid-1) * 4 + 1 ));
    
    setKappaBar(vh, newKappaBar);
    
    setRefVertexLength(vh, e( (vid-1) * 4 + 3 ));
    if (vid > 1) {
      setRefVertexLength(vh, e( (vid-2) * 4 + 3 ));
    }
  }  
}
  
const MatXd& RodBendingForceSym::getGradKappa(const vertex_handle& vh) const
{
  return m_rod.property(m_gradKappa)[vh];
}

const pair<MatXd, MatXd>&
RodBendingForceSym::getHessKappa(const vertex_handle& vh) const
{
  return m_rod.property(m_hessKappa)[vh];
}

void RodBendingForceSym::computeGradKappa()
{
  if (m_rod.property(m_gradKappaValid)) return;

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    int i = vh.idx();

    MatXd& gradKappa = m_rod.property(m_gradKappa)[vh];

    // Unused? const Vec3d& e = m_rod.getEdge(i-1);
    // Unused? const Vec3d& f = m_rod.getEdge(i);

    Scalar norm_e = m_rod.getEdgeLength(i-1);
    Scalar norm_f = m_rod.getEdgeLength(i);

    const Vec3d& te = m_rod.getTangent(i-1);
    const Vec3d& tf = m_rod.getTangent(i);

    const Vec3d& d1e = m_rod.getMaterial1(i-1);
    const Vec3d& d2e = m_rod.getMaterial2(i-1);
    const Vec3d& d1f = m_rod.getMaterial1(i);
    const Vec3d& d2f = m_rod.getMaterial2(i);

    Scalar chi = 1.0 + te.dot(tf);
    Vec3d tilde_t = (te + tf) / chi;
    Vec3d tilde_d1 = (d1e + d1f) / chi;
    Vec3d tilde_d2 = (d2e + d2f) / chi;

    const Vec2d& kappa = getKappa(vh);
    Scalar kappa1 = kappa(0);
    Scalar kappa2 = kappa(1);

    Vec3d Dkappa1De = 1.0 / norm_e * (-kappa1 * tilde_t + tf.cross(tilde_d2));
    Vec3d Dkappa1Df = 1.0 / norm_f * (-kappa1 * tilde_t - te.cross(tilde_d2));

    Vec3d Dkappa2De = 1.0 / norm_e * (-kappa2 * tilde_t - tf.cross(tilde_d1));
    Vec3d Dkappa2Df = 1.0 / norm_f * (-kappa2 * tilde_t + te.cross(tilde_d1));

    gradKappa.block<3,1>(0,0) = -Dkappa1De;
    gradKappa.block<3,1>(4,0) =  Dkappa1De - Dkappa1Df;
    gradKappa.block<3,1>(8,0) =              Dkappa1Df;

    gradKappa.block<3,1>(0,1) = -Dkappa2De;
    gradKappa.block<3,1>(4,1) =  Dkappa2De - Dkappa2Df;
    gradKappa.block<3,1>(8,1) =              Dkappa2Df;

    const Vec3d& kb = m_rod.getCurvatureBinormal(i);

    gradKappa(3,0) = -0.5 * kb.dot(d1e);
    gradKappa(7,0) = -0.5 * kb.dot(d1f);
    gradKappa(3,1) = -0.5 * kb.dot(d2e);
    gradKappa(7,1) = -0.5 * kb.dot(d2f);
  }

  m_rod.property(m_gradKappaValid) = true;
}

void RodBendingForceSym::computeHessKappa()
{
  if (m_rod.property(m_hessKappaValid)) return;

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    int i = vh.idx();

    MatXd& DDkappa1 = m_rod.property(m_hessKappa)[vh].first;
    MatXd& DDkappa2 = m_rod.property(m_hessKappa)[vh].second;
    DDkappa1.setZero();
    DDkappa2.setZero();

    // Unused? const Vec3d& e = m_rod.getEdge(i-1);
    // Unused? const Vec3d& f = m_rod.getEdge(i);

    Scalar norm_e = m_rod.getEdgeLength(i-1);
    Scalar norm_f = m_rod.getEdgeLength(i);

    Scalar norm2_e = square(norm_e);
    Scalar norm2_f = square(norm_f);

    const Vec3d& te = m_rod.getTangent(i-1);
    const Vec3d& tf = m_rod.getTangent(i);

    const Vec3d& d1e = m_rod.getMaterial1(i-1);
    const Vec3d& d2e = m_rod.getMaterial2(i-1);
    const Vec3d& d1f = m_rod.getMaterial1(i);
    const Vec3d& d2f = m_rod.getMaterial2(i);

    Scalar chi = 1.0 + te.dot(tf);
    Vec3d tilde_t = (te + tf) / chi;
    Vec3d tilde_d1 = (d1e + d1f) / chi;
    Vec3d tilde_d2 = (d2e + d2f) / chi;

    const Vec2d& kappa = getKappa(vh);
    Scalar kappa1 = kappa(0);
    Scalar kappa2 = kappa(1);

    const Vec3d& kb = m_rod.getCurvatureBinormal(i);

    Mat3d tt_o_tt = outerProd(tilde_t, tilde_t);
    Mat3d tf_c_d2t_o_tt = outerProd(tf.cross(tilde_d2), tilde_t);
    Mat3d tt_o_tf_c_d2t = tf_c_d2t_o_tt.transpose();
    Mat3d kb_o_d2e = outerProd(kb, d2e);
    Mat3d d2e_o_kb = kb_o_d2e.transpose();

    Mat3d Id = Mat3d::Identity();

    Mat3d D2kappa1De2
      = 1.0 / norm2_e * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt - tt_o_tf_c_d2t)
      - kappa1 / (chi * norm2_e) * (Id - outerProd(te, te))
      + 1.0 / (4.0 * norm2_e) * (kb_o_d2e + d2e_o_kb);

    Mat3d te_c_d2t_o_tt = outerProd(te.cross(tilde_d2), tilde_t);
    Mat3d tt_o_te_c_d2t = te_c_d2t_o_tt.transpose();
    Mat3d kb_o_d2f = outerProd(kb, d2f);
    Mat3d d2f_o_kb = kb_o_d2f.transpose();

    Mat3d D2kappa1Df2
      = 1.0 / norm2_f * (2 * kappa1 * tt_o_tt + te_c_d2t_o_tt + tt_o_te_c_d2t)
      - kappa1 / (chi * norm2_f) * (Id - outerProd(tf, tf))
      + 1.0 / (4.0 * norm2_f) * (kb_o_d2f + d2f_o_kb);

    Mat3d D2kappa1DeDf
      = -kappa1/(chi * norm_e * norm_f) * (Id + outerProd(te, tf))
      + 1.0 / (norm_e*norm_f) * (2 * kappa1 * tt_o_tt - tf_c_d2t_o_tt + tt_o_te_c_d2t - crossMat(tilde_d2));
    Mat3d D2kappa1DfDe = D2kappa1DeDf.transpose();

    Mat3d tf_c_d1t_o_tt = outerProd(tf.cross(tilde_d1), tilde_t);
    Mat3d tt_o_tf_c_d1t = tf_c_d1t_o_tt.transpose();
    Mat3d kb_o_d1e = outerProd(kb, d1e);
    Mat3d d1e_o_kb = kb_o_d1e.transpose();

    Mat3d D2kappa2De2
      = 1.0 / norm2_e * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt + tt_o_tf_c_d1t)
      - kappa2 / (chi * norm2_e) * (Id - outerProd(te, te))
      - 1.0 / (4.0 * norm2_e) * (kb_o_d1e + d1e_o_kb);

    Mat3d te_c_d1t_o_tt = outerProd(te.cross(tilde_d1), tilde_t);
    Mat3d tt_o_te_c_d1t = te_c_d1t_o_tt.transpose();
    Mat3d kb_o_d1f = outerProd(kb, d1f);
    Mat3d d1f_o_kb =  kb_o_d1f.transpose();

    Mat3d D2kappa2Df2
      = 1.0 / norm2_f * (2 * kappa2 * tt_o_tt - te_c_d1t_o_tt - tt_o_te_c_d1t)
      - kappa2 / (chi * norm2_f) * (Id - outerProd(tf, tf))
      - 1.0 / (4.0 * norm2_f) * (kb_o_d1f + d1f_o_kb);

    Mat3d D2kappa2DeDf
      = -kappa2/(chi * norm_e * norm_f) * (Id + outerProd(te, tf))
      + 1.0 / (norm_e*norm_f) * (2 * kappa2 * tt_o_tt + tf_c_d1t_o_tt - tt_o_te_c_d1t + crossMat(tilde_d1));
    Mat3d D2kappa2DfDe = D2kappa2DeDf.transpose();

    Scalar D2kappa1Dthetae2 = -0.5 * kb.dot(d2e);
    Scalar D2kappa1Dthetaf2 = -0.5 * kb.dot(d2f);
    Scalar D2kappa2Dthetae2 =  0.5 * kb.dot(d1e);
    Scalar D2kappa2Dthetaf2 =  0.5 * kb.dot(d1f);

    Vec3d D2kappa1DeDthetae
      = 1.0 / norm_e * (0.5 * kb.dot(d1e) * tilde_t - 1.0 / chi * tf.cross(d1e));

    Vec3d D2kappa1DeDthetaf
      = 1.0 / norm_e * (0.5 * kb.dot(d1f) * tilde_t - 1.0 / chi * tf.cross(d1f));

    Vec3d D2kappa1DfDthetae
      = 1.0 / norm_f * (0.5 * kb.dot(d1e) * tilde_t + 1.0 / chi * te.cross(d1e));

    Vec3d D2kappa1DfDthetaf
      = 1.0 / norm_f * (0.5 * kb.dot(d1f) * tilde_t + 1.0 / chi * te.cross(d1f));

    Vec3d D2kappa2DeDthetae
      = 1.0 / norm_e * (0.5 * kb.dot(d2e) * tilde_t - 1.0 / chi * tf.cross(d2e));

    Vec3d D2kappa2DeDthetaf
      = 1.0 / norm_e * (0.5 * kb.dot(d2f) * tilde_t - 1.0 / chi * tf.cross(d2f));

    Vec3d D2kappa2DfDthetae
      = 1.0 / norm_f * (0.5 * kb.dot(d2e) * tilde_t + 1.0 / chi * te.cross(d2e));

    Vec3d D2kappa2DfDthetaf
      = 1.0 / norm_f * (0.5 * kb.dot(d2f) * tilde_t + 1.0 / chi * te.cross(d2f));

    DDkappa1.block<3,3>(0,0) =   D2kappa1De2;
    DDkappa1.block<3,3>(0,4) = - D2kappa1De2 + D2kappa1DeDf;
    DDkappa1.block<3,3>(0,8) =               - D2kappa1DeDf;
    DDkappa1.block<3,3>(4,0) = - D2kappa1De2                + D2kappa1DfDe;
    DDkappa1.block<3,3>(4,4) =   D2kappa1De2 - D2kappa1DeDf - D2kappa1DfDe + D2kappa1Df2;
    DDkappa1.block<3,3>(4,8) =                 D2kappa1DeDf                - D2kappa1Df2;
    DDkappa1.block<3,3>(8,0) =                              - D2kappa1DfDe;
    DDkappa1.block<3,3>(8,4) =                                D2kappa1DfDe - D2kappa1Df2;
    DDkappa1.block<3,3>(8,8) =                                               D2kappa1Df2;

    DDkappa1(3,3) = D2kappa1Dthetae2;
    DDkappa1(7,7) = D2kappa1Dthetaf2;

    DDkappa1.block<3,1>(0,3) = - D2kappa1DeDthetae;
    DDkappa1.block<3,1>(4,3) =   D2kappa1DeDthetae - D2kappa1DfDthetae;
    DDkappa1.block<3,1>(8,3) =                       D2kappa1DfDthetae;
    DDkappa1.block<1,3>(3,0) = DDkappa1.block<3,1>(0,3).transpose();
    DDkappa1.block<1,3>(3,4) = DDkappa1.block<3,1>(4,3).transpose();
    DDkappa1.block<1,3>(3,8) = DDkappa1.block<3,1>(8,3).transpose();

    DDkappa1.block<3,1>(0,7) = - D2kappa1DeDthetaf;
    DDkappa1.block<3,1>(4,7) =   D2kappa1DeDthetaf - D2kappa1DfDthetaf;
    DDkappa1.block<3,1>(8,7) =                       D2kappa1DfDthetaf;
    DDkappa1.block<1,3>(7,0) = DDkappa1.block<3,1>(0,7).transpose();
    DDkappa1.block<1,3>(7,4) = DDkappa1.block<3,1>(4,7).transpose();
    DDkappa1.block<1,3>(7,8) = DDkappa1.block<3,1>(8,7).transpose();

    DDkappa2.block<3,3>(0,0) =   D2kappa2De2;
    DDkappa2.block<3,3>(0,4) = - D2kappa2De2 + D2kappa2DeDf;
    DDkappa2.block<3,3>(0,8) =               - D2kappa2DeDf;
    DDkappa2.block<3,3>(4,0) = - D2kappa2De2                + D2kappa2DfDe;
    DDkappa2.block<3,3>(4,4) =   D2kappa2De2 - D2kappa2DeDf - D2kappa2DfDe + D2kappa2Df2;
    DDkappa2.block<3,3>(4,8) =                 D2kappa2DeDf                - D2kappa2Df2;
    DDkappa2.block<3,3>(8,0) =                              - D2kappa2DfDe;
    DDkappa2.block<3,3>(8,4) =                                D2kappa2DfDe - D2kappa2Df2;
    DDkappa2.block<3,3>(8,8) =                                               D2kappa2Df2;

    DDkappa2(3,3) = D2kappa2Dthetae2;
    DDkappa2(7,7) = D2kappa2Dthetaf2;

    DDkappa2.block<3,1>(0,3) = - D2kappa2DeDthetae;
    DDkappa2.block<3,1>(4,3) =   D2kappa2DeDthetae - D2kappa2DfDthetae;
    DDkappa2.block<3,1>(8,3) =                       D2kappa2DfDthetae;
    DDkappa2.block<1,3>(3,0) = DDkappa2.block<3,1>(0,3).transpose();
    DDkappa2.block<1,3>(3,4) = DDkappa2.block<3,1>(4,3).transpose();
    DDkappa2.block<1,3>(3,8) = DDkappa2.block<3,1>(8,3).transpose();

    DDkappa2.block<3,1>(0,7) = - D2kappa2DeDthetaf;
    DDkappa2.block<3,1>(4,7) =   D2kappa2DeDthetaf - D2kappa2DfDthetaf;
    DDkappa2.block<3,1>(8,7) =                       D2kappa2DfDthetaf;
    DDkappa2.block<1,3>(7,0) = DDkappa2.block<3,1>(0,7).transpose();
    DDkappa2.block<1,3>(7,4) = DDkappa2.block<3,1>(4,7).transpose();
    DDkappa2.block<1,3>(7,8) = DDkappa2.block<3,1>(8,7).transpose();
  }

  m_rod.property(m_hessKappaValid) = true;
}

Scalar RodBendingForceSym::getRefVertexLength(const vertex_handle& vh) const
{
  return m_rod.property(m_refVertexLength)[vh];
}

void RodBendingForceSym::setRefVertexLength(const vertex_handle& vh,
                                            const Scalar& length)
{
  m_rod.property(m_refVertexLength)[vh] = length;
}

} // namespace BASim
