#include "RodModelTwistingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"

using namespace BASim;

RodModelTwistingForce::RodModelTwistingForce(ElasticRodModel & rod, const std::vector<Stencil> & stencils, Scalar shear_modulus, Scalar shear_modulus_damping, Scalar timestep) :
  RodModelForce(rod, timestep, "RodModelStretchingForce"),
  m_stencils(stencils),
  m_shear_modulus(shear_modulus),
  m_shear_modulus_damping(shear_modulus_damping),
  m_stiffness(&rod.getDefoObj()),
  m_viscous_stiffness(&rod.getDefoObj()),
  m_undeformed_twist(&rod.getDefoObj()),
  m_damping_undeformed_twist(&rod.getDefoObj()),
  m_reference_voronoi_length(&rod.getDefoObj()),
  m_twist(&rod.getDefoObj())
{
  updateProperties();
  updateStiffness();
  updateViscousReferenceStrain();
  computeReferenceStrain();
}

RodModelTwistingForce::~RodModelTwistingForce()
{
  
}

Scalar RodModelTwistingForce::globalEnergy()
{
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    // energy only include non-viscous forces
    energy += localEnergy(m_stencils[i], false);
  }
  return energy;
}

void RodModelTwistingForce::globalForce(VecXd & force)
{
  ElementForce localforce;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    
    // non-viscous force
    localForce(localforce, s, false);
    localforce(3) *= (s.e1flip ? -1 : 1);
    localforce(7) *= (s.e2flip ? -1 : 1);
    for (size_t j = 0; j < s.dofindices.size(); j++)
      force(s.dofindices[j]) += localforce(j);
    
    // viscous force
    localForce(localforce, s, true);
    localforce(3) *= (s.e1flip ? -1 : 1);
    localforce(7) *= (s.e2flip ? -1 : 1);
    for (size_t j = 0; j < s.dofindices.size(); j++)
      force(s.dofindices[j]) += localforce(j) / timeStep();
  }
}

void RodModelTwistingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
{
  ElementJacobian localjacobian;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    
    // non-viscous force
    localJacobian(localjacobian, s, false);
    if (s.e1flip) 
    {
      localjacobian.row(3) *= -1;
      localjacobian.col(3) *= -1;
    }
    if (s.e2flip)
    {
      localjacobian.row(7) *= -1;
      localjacobian.col(7) *= -1;
    }
    Jacobian.add(s.dofindices, s.dofindices, scale * localjacobian);
    
    // viscous force
    localJacobian(localjacobian, s, true);
    if (s.e1flip) 
    {
      localjacobian.row(3) *= -1;
      localjacobian.col(3) *= -1;
    }
    if (s.e2flip)
    {
      localjacobian.row(7) *= -1;
      localjacobian.col(7) *= -1;
    }
    Jacobian.add(s.dofindices, s.dofindices, scale / timeStep() * localjacobian);
  }
}

Scalar RodModelTwistingForce::localEnergy(Stencil & s, bool viscous)
{
  Scalar kt = (viscous ? m_viscous_stiffness[s.v2] : m_stiffness[s.v2]);
  Scalar len = m_reference_voronoi_length[s.v2];
  Scalar undefTwist = (viscous ? m_damping_undeformed_twist[s.v2] : m_undeformed_twist[s.v2]);
  Scalar twist = m_twist[s.v2];
  
  return kt / ( 2.0 * len ) * square( twist - undefTwist );
}

void RodModelTwistingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
  Scalar kt = (viscous ? m_viscous_stiffness[s.v2] : m_stiffness[s.v2]);
  Scalar len = m_reference_voronoi_length[s.v2];
  Scalar undefTwist = (viscous ? m_damping_undeformed_twist[s.v2] : m_undeformed_twist[s.v2]);
  Scalar twist = m_twist[s.v2];

  force = -kt / len * (twist - undefTwist) * computeGradTwist(s);
}

void RodModelTwistingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
  Scalar kt = (viscous ? m_viscous_stiffness[s.v2] : m_stiffness[s.v2]);
  Scalar len = m_reference_voronoi_length[s.v2];
  Scalar undefTwist = (viscous ? m_damping_undeformed_twist[s.v2] : m_undeformed_twist[s.v2]);
  Scalar twist = m_twist[s.v2];
  
  ElementForce gradTwist = computeGradTwist(s);
  ElementJacobian hessTwist = computeHessTwist(s);
  
  jacobian = -kt / len * ( ( twist - undefTwist ) * hessTwist + gradTwist * gradTwist.transpose() );
  
  assert( isSymmetric( jacobian ) );
}

void RodModelTwistingForce::updateStiffness()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    
    Vec2d ra = rod().getRadii(s.e1);
    Vec2d rb = rod().getRadii(s.e2);
    
    Scalar a = (ra(0) + rb(0)) / 2.0;
    Scalar b = (ra(1) + rb(1)) / 2.0;

    Scalar k = M_PI * a * b * (square(a) + square(b)) / 4.0;
    
    m_stiffness[s.v2] = k * m_shear_modulus;
    m_viscous_stiffness[s.v2] = k * m_shear_modulus_damping;
  }
}

void RodModelTwistingForce::updateViscousReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    m_damping_undeformed_twist[s.v2] = m_twist[s.v2];
  }
}

void RodModelTwistingForce::updateProperties()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    m_twist[s.v2] = s.referenceTwist + rod().getEdgeTheta(s.e2) * (s.e2flip ? -1 : 1) - rod().getEdgeTheta(s.e1) * (s.e1flip ? -1 : 1);
  }
}

void RodModelTwistingForce::computeReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    m_undeformed_twist[s.v2] = m_twist[s.v2];
    m_reference_voronoi_length[s.v2] = s.voronoiLength;
  }
}

RodModelTwistingForce::ElementForce RodModelTwistingForce::computeGradTwist(Stencil & s)
{
  ElementForce Dtwist = ElementForce::Zero();
  const Vec3d& kb = s.curvatureBinormal;
  Dtwist.segment<3> ( 0 ) = -0.5 / ( rod().getEdgeLength(s.e1) ) * kb;
  Dtwist.segment<3> ( 8 ) = 0.5 / ( rod().getEdgeLength(s.e2) ) * kb;
  Dtwist.segment<3> ( 4 ) = -( Dtwist.segment<3> ( 0 ) + Dtwist.segment<3> ( 8 ) );
  Dtwist( 3 ) = -1;
  Dtwist( 7 ) = 1;
  return Dtwist;
}

RodModelTwistingForce::ElementJacobian RodModelTwistingForce::computeHessTwist(Stencil & s)
{
  ElementJacobian DDtwist;
  DDtwist.setZero();

  const Vec3d  te = rod().getEdgeTangent(s.e1) * (s.e1flip ? -1 : 1);
  const Vec3d  tf = rod().getEdgeTangent(s.e2) * (s.e2flip ? -1 : 1);
  const Scalar norm_e = rod().getEdgeLength(s.e1);
  const Scalar norm_f = rod().getEdgeLength(s.e2);
  const Vec3d& kb = s.curvatureBinormal;
  
  const Scalar chi = 1 + te.dot( tf );
  const Vec3d tilde_t = 1.0 / chi * ( te + tf );
  
  const Mat3d D2mDe2 = -0.25 / square( norm_e ) * ( outerProd( kb, te + tilde_t )
                       + outerProd( te + tilde_t, kb ) );
  assert( isSymmetric( D2mDe2 ) );
  
  const Mat3d D2mDf2 = -0.25 / square( norm_f ) * ( outerProd( kb, tf + tilde_t )
                       + outerProd( tf + tilde_t, kb ) );
  assert( isSymmetric( D2mDf2 ) );
  
  const Mat3d D2mDeDf = 0.5 / ( norm_e * norm_f ) * ( 2.0 / chi * crossMat( te ) - outerProd(
                        kb, tilde_t ) );
  const Mat3d D2mDfDe = D2mDeDf.transpose();
  
  DDtwist.block<3, 3> ( 0, 0 ) = D2mDe2;
  DDtwist.block<3, 3> ( 0, 4 ) = -D2mDe2 + D2mDeDf;
  DDtwist.block<3, 3> ( 4, 0 ) = -D2mDe2 + D2mDfDe;
  DDtwist.block<3, 3> ( 4, 4 ) = D2mDe2 - ( D2mDeDf + D2mDfDe ) + D2mDf2;
  DDtwist.block<3, 3> ( 0, 8 ) = -D2mDeDf;
  DDtwist.block<3, 3> ( 8, 0 ) = -D2mDfDe;
  DDtwist.block<3, 3> ( 8, 4 ) = D2mDfDe - D2mDf2;
  DDtwist.block<3, 3> ( 4, 8 ) = D2mDeDf - D2mDf2;
  DDtwist.block<3, 3> ( 8, 8 ) = D2mDf2;
  
  assert( isSymmetric( DDtwist ) );
  
  return DDtwist;
}
