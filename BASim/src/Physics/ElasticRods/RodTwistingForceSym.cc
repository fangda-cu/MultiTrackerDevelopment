/**
 * \file RodTwistingForceSym.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 02/18/2010
 */

#ifdef WETA
#include "RodTwistingForceSym.hh"
#include "ElasticRod.hh"
#include "../../Math/Math.hh"
#include "../../Core/Timer.hh"
#else
#include "BASim/Physics"
#endif

using namespace std;

namespace BASim {

RodTwistingForceSym::RodTwistingForceSym(ElasticRod& rod, bool vscs, bool runinit)
  : RodForceT<VertexStencil>(rod, "RodTwistingForceSym")
{
  if( runinit )
  {
    m_rod.add_property(m_twist, "twist");
    m_rod.add_property(m_gradTwistValid, "grad twist valid", false);
    m_rod.add_property(m_hessTwistValid, "hess twist valid", false);
    m_rod.add_property(m_gradTwist, "gradient of twist", VecXd(11));
    m_rod.add_property(m_hessTwist, "Hessian of twist", MatXd(11, 11));
    
    if( !vscs )
    {
      m_rod.add_property(m_kt, "twist stiffness");
      m_rod.add_property(m_undeformedTwist, "undeformed twist");
      m_rod.add_property(m_refVertexLength, "twisting ref length");
    }
    else
    {
      m_rod.add_property(m_kt, "viscous twist stiffness");
      m_rod.add_property(m_undeformedTwist, "viscous undeformed twist");
      m_rod.add_property(m_refVertexLength, "viscous twisting ref length");
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

void RodTwistingForceSym::updateUndeformedConfiguration(std::vector<Scalar>& vals) {
  /*
  std::cout << "twist origin\n";
  
  {
    iterator end = m_stencil.end();
    for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
      vertex_handle& vh = m_stencil.handle();
      std::cout << getUndeformedTwist(vh) << "\n";
    }  
  }*/

  //std::cout << "twist new\n";
  
  int i=0;
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil, ++i) {
    vertex_handle& vh = m_stencil.handle();
    setUndeformedTwist(vh, vals[i]);
//    std::cout << getUndeformedTwist(vh) << "\n";
  }  
  
}

void RodTwistingForceSym::updateProperties()
{
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    edge_handle eh0 = m_stencil.inEdge();
    edge_handle eh1 = m_stencil.outEdge();
    setTwist(vh, m_rod.getReferenceTwist(vh)
             + m_rod.getTheta(eh1) - m_rod.getTheta(eh0));
  }
  m_rod.property(m_gradTwistValid) = false;
  m_rod.property(m_hessTwistValid) = false;
}

void RodTwistingForceSym::updateStiffness()
{
  Scalar G = m_rod.getShearModulus();
  if (viscous()) {
    G = m_rod.getViscosity() / m_rod.getTimeStep();
  }

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    edge_handle eh0 = m_stencil.inEdge();
    edge_handle eh1 = m_stencil.outEdge();
    Scalar a = (m_rod.radiusA(eh0) + m_rod.radiusA(eh1)) / 2.0;
    Scalar b = (m_rod.radiusB(eh0) + m_rod.radiusB(eh1)) / 2.0;
    setKt(vh, G * M_PI * a * b * (square(a) + square(b)) / 4.0);
    
//    std::cout << G * M_PI * a * b * (square(a) + square(b)) / 4.0 << "\n";
    
  }
}

void RodTwistingForceSym::updateReferenceDomain()
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

// static inline bool stencilFixed(ElasticRod& rod, const vertex_handle& vh)
// {
//   return (rod.vertFixed(i-1) && rod.vertFixed(i) && rod.vertFixed(i+1)
//           && rod.edgeFixed(i-1) && rod.edgeFixed(i))
//     || rod.edgeCut(i-1) || rod.edgeCut(i);
// }

Scalar RodTwistingForceSym::globalEnergy()
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

Scalar RodTwistingForceSym::localEnergy(const vertex_handle& vh)
{
  Scalar kt = getKt(vh);
  Scalar len = getRefVertexLength(vh);
  Scalar undefTwist = getUndeformedTwist(vh);
  Scalar twist = getTwist(vh);

  return kt/(2.0*len) * square(twist - undefTwist);
}

void RodTwistingForceSym::globalForce(VecXd& force)
{
  if (viscous() && m_rod.getViscosity() == 0.0) return;

  //START_TIMER("globalForce");
  computeGradTwist();

  VecXd f(11);
  IndexArray indices(11);

  VecXd force1 = force;
  
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    //if (stencilFixed(m_rod, i)) continue;
    //std::cout << "TWISTING local FORCE\n";

    vertex_handle& vh = m_stencil.handle();
    localForce(f, vh);
    m_stencil.indices(indices);
    for (int j = 0; j < f.size(); ++j) {
      force(indices[j]) += f(j);
    }
  }
  //STOP_TIMER("globalForce");
  
  
  
  //std::cout << "TWISTING FORCE\n";
  //std::cout << force - force1 << "\n\n";
    
}

void RodTwistingForceSym::localForce(VecXd& force, const vertex_handle& vh)
{
  Scalar value = getKt(vh) / getRefVertexLength(vh)
    * (getTwist(vh) - getUndeformedTwist(vh));

  //std::cout << "-value    : " << -value << "\n";
  //std::cout << "GradTwist : " << getGradTwist(vh) << "\n";
  //std::cout << "local for : " << -value * getGradTwist(vh) << "\n";
  
  force = -value * getGradTwist(vh);
}

void RodTwistingForceSym::globalJacobian(int baseidx, Scalar scale, MatrixBase& J)
{
  if (viscous() && m_rod.getViscosity() == 0.0) return;

  //START_TIMER("globalJacobian");
  computeGradTwist();
  computeHessTwist();

  MatXd localJ(11, 11);
  IndexArray indices(11);

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    //if (stencilFixed(m_rod, i)) continue;
    vertex_handle& vh = m_stencil.handle();
    localJ.setZero();
    localJacobian(localJ, vh);
    m_stencil.indices(indices);
    for( int i = 0; i < (int) indices.size(); ++i ) indices[i] += baseidx;
    localJ *= scale;
    J.add(indices, indices, localJ);
  }
  //STOP_TIMER("globalJacobian");
}

inline void RodTwistingForceSym::localJacobian(MatXd& J, const vertex_handle& vh)
{
  Scalar kt = getKt(vh);
  Scalar len = getRefVertexLength(vh);
  Scalar twist = getTwist(vh);
  Scalar undeformedTwist = getUndeformedTwist(vh);

  const VecXd& gradTwist = getGradTwist(vh);
  const MatXd& hessTwist = getHessTwist(vh);

  J = -kt / len * ((twist - undeformedTwist) * hessTwist
                   + gradTwist * gradTwist.transpose());
}

void RodTwistingForceSym::globalReverseJacobian(MatrixBase& J)
{
  if (viscous()) return;

//  return;

  computeGradTwist();
  computeHessTwist(); // do we need hessian?

  //MatXd localJ(11, 11);
  IndexArray indices(11);
  iterator end = m_stencil.end();

  uint vid = 1;
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil, ++vid) {
    vertex_handle& vh = m_stencil.handle();
    
    Scalar m = getTwist(vh);
    Scalar mbar = getUndeformedTwist(vh);
    Scalar kt = getKt(vh);
    Scalar refLen = getRefVertexLength(vh);

    VecXd fde(11);
    VecXd fdm(11);
    
    fde = 0.5 * kt / (refLen * refLen) * (m - mbar) * getGradTwist(vh);
    fdm = 1.0 * kt / refLen * getGradTwist(vh); 

    m_stencil.indices(indices);
    
    for(int i=0; i<indices.size(); i++) indices[i] -= 7;
    
//    std::cout << "twst fde  " << fde << "\n";
//    std::cout << "twst fdm  " << fdm << "\n";
    
    // edge length bar
    if (vid == 1) {
      for(int i=0; i<11; i++) {
        if (indices[i] >= 0) {
          J.add(indices[i], (vid - 1) * 4 + 3, fde(i));
        }
      }
    } else {
      for(int i=0; i<11; i++) {
        if (indices[i] >= 0) {
          J.add(indices[i], (vid - 2) * 4 + 3, fde(i));
          J.add(indices[i], (vid - 1) * 4 + 3, fde(i));
        }
      }
    }
    
    // m bar
    for(int i=0; i<11; i++) {
      if (indices[i] >= 0) {
        J.add(indices[i], (vid - 1) * 4 + 2, fdm(i));
      }
    }
  }  
}

void RodTwistingForceSym::updateReverseUndeformedStrain(const VecXd& e)
{
  if (viscous()) return; 
    
  m_stencil = m_stencil.begin();
  iterator end = m_stencil.end();
  
  uint vid = 1;  // eh.id() ?
  
  for (; m_stencil != end; ++m_stencil, ++vid) {
    vertex_handle& vh = m_stencil.handle();
    setUndeformedTwist(vh, e( (vid-1) * 4 + 2 ));
    
    setRefVertexLength(vh, e( (vid-1) * 4 + 3 ));
    if (vid > 1) {
      setRefVertexLength(vh, e( (vid-2) * 4 + 3 ));
    }
  }  
}

void RodTwistingForceSym::updateUndeformedStrain()
{
  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    setUndeformedTwist(vh, getTwist(vh));
  }
}

Scalar RodTwistingForceSym::getKt(const vertex_handle& vh) const
{
  return m_rod.property(m_kt)[vh];
}

void RodTwistingForceSym::setKt(const vertex_handle& vh, const Scalar& kt)
{
  m_rod.property(m_kt)[vh] = kt;
}

Scalar RodTwistingForceSym::getTwist(const vertex_handle& vh) const
{
  return m_rod.property(m_twist)[vh];
}

void RodTwistingForceSym::setTwist(const vertex_handle& vh, const Scalar& twist)
{
  m_rod.property(m_twist)[vh] = twist;
}

Scalar
RodTwistingForceSym::getUndeformedTwist(const vertex_handle& vh) const
{
  return m_rod.property(m_undeformedTwist)[vh];
}

void RodTwistingForceSym::setUndeformedTwist(const vertex_handle& vh,
                                             const Scalar& undeformedTwist)
{
  m_rod.property(m_undeformedTwist)[vh] = undeformedTwist;
}

Scalar
RodTwistingForceSym::getRefVertexLength(const vertex_handle& vh) const
{
  return m_rod.property(m_refVertexLength)[vh];
}

void RodTwistingForceSym::setRefVertexLength(const vertex_handle& vh,
                                             const Scalar& length)
{
  m_rod.property(m_refVertexLength)[vh] = length;
}

const VecXd& RodTwistingForceSym::getGradTwist(const vertex_handle& vh) const
{
  return m_rod.property(m_gradTwist)[vh];
}

const MatXd& RodTwistingForceSym::getHessTwist(const vertex_handle& vh) const
{
  return m_rod.property(m_hessTwist)[vh];
}

void RodTwistingForceSym::computeGradTwist()
{
  if (m_rod.property(m_gradTwistValid)) return;

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    VecXd& Dtwist = m_rod.property(m_gradTwist)[vh];
    const Vec3d& kb = m_rod.getCurvatureBinormal(vh);
    int i = vh.idx();
    Dtwist.segment<3>(0) = -0.5 / (m_rod.getEdgeLength(i - 1)) * kb;
    Dtwist.segment<3>(8) =  0.5 / (m_rod.getEdgeLength(i)) * kb;
    Dtwist.segment<3>(4) = -(Dtwist.segment<3>(0) + Dtwist.segment<3>(8));
    Dtwist(3) = -1;
    Dtwist(7) = 1;
  }

  m_rod.property(m_gradTwistValid) = true;
}

void RodTwistingForceSym::computeHessTwist()
{
  if (m_rod.property(m_hessTwistValid)) return;

  iterator end = m_stencil.end();
  for (m_stencil = m_stencil.begin(); m_stencil != end; ++m_stencil) {
    vertex_handle& vh = m_stencil.handle();
    MatXd& DDtwist = m_rod.property(m_hessTwist)[vh];
    DDtwist.setZero();
    int i = vh.idx();
    const Vec3d& te = m_rod.getTangent(i-1);
    const Vec3d& tf = m_rod.getTangent(i);
    Scalar norm_e = m_rod.getEdgeLength(i-1);
    Scalar norm_f = m_rod.getEdgeLength(i);
    const Vec3d& kb = m_rod.getCurvatureBinormal(vh);

    Scalar chi = 1 + te.dot(tf);
    Vec3d tilde_t = 1.0 / chi * (te + tf);

    Mat3d D2mDe2 = -0.25 / square(norm_e) * (outerProd(kb, te + tilde_t)
                                             + outerProd(te + tilde_t, kb));
    Mat3d D2mDf2 = -0.25 / square(norm_f) * (outerProd(kb, tf + tilde_t)
                                             + outerProd(tf + tilde_t, kb));
    Mat3d D2mDeDf = 0.5 / (norm_e * norm_f)
      * (2.0 / chi * crossMat(te) - outerProd(kb, tilde_t));
    Mat3d D2mDfDe = D2mDeDf.transpose();

    DDtwist.block<3,3>(0,0) =   D2mDe2;
    DDtwist.block<3,3>(0,4) = - D2mDe2 + D2mDeDf;
    DDtwist.block<3,3>(0,8) =          - D2mDeDf;
    DDtwist.block<3,3>(4,0) = - D2mDe2           + D2mDfDe;
    DDtwist.block<3,3>(4,4) =   D2mDe2 - D2mDeDf - D2mDfDe + D2mDf2;
    DDtwist.block<3,3>(4,8) =            D2mDeDf           - D2mDf2;
    DDtwist.block<3,3>(8,0) =                    - D2mDfDe;
    DDtwist.block<3,3>(8,4) =                      D2mDfDe - D2mDf2;
    DDtwist.block<3,3>(8,8) =                                D2mDf2;
  }

  m_rod.property(m_hessTwistValid) = true;
}

} // namespace BASim