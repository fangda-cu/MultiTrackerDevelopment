#include "RodModelTwistingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"

using namespace BASim;

RodModelTwistingForce::RodModelTwistingForce(ElasticRodModel & rod, const std::vector<ElasticRodModel::JointStencil> & stencils, Scalar shear_modulus, Scalar shear_modulus_damping, Scalar timestep) :
  RodModelForce(rod, timestep, "RodModelStretchingForce"),
  m_stencils(),
  m_shear_modulus(shear_modulus),
  m_shear_modulus_damping(shear_modulus_damping)
{
  for (size_t i = 0; i < stencils.size(); i++)
  {
    Stencil s(stencils[i]);
    s.stiffness = 0;
    s.viscous_stiffness = 0;
    s.undeformed_twist = 0;
    s.damping_undeformed_twist = 0;
    s.reference_voronoi_length = 0;
    s.twist = 0;
    
    m_stencils.push_back(s);
  }

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
  Scalar kt = (viscous ? s.viscous_stiffness : s.stiffness);
  Scalar len = s.reference_voronoi_length;
  Scalar undefTwist = (viscous ? s.damping_undeformed_twist : s.undeformed_twist);
  Scalar twist = s.twist;
  
  return kt / ( 2.0 * len ) * square( twist - undefTwist );
}

void RodModelTwistingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
  Scalar kt = (viscous ? s.viscous_stiffness : s.stiffness);
  Scalar len = s.reference_voronoi_length;
  Scalar undefTwist = (viscous ? s.damping_undeformed_twist : s.undeformed_twist);
  Scalar twist = s.twist;

  force = -kt / len * (twist - undefTwist) * computeGradTwist(s);
}

void RodModelTwistingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
  Scalar kt = (viscous ? s.viscous_stiffness : s.stiffness);
  Scalar len = s.reference_voronoi_length;
  Scalar undefTwist = (viscous ? s.damping_undeformed_twist : s.undeformed_twist);
  Scalar twist = s.twist;
  
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
    
    s.stiffness = k * m_shear_modulus;
    s.viscous_stiffness = k * m_shear_modulus_damping;
  }
}

void RodModelTwistingForce::updateViscousReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.damping_undeformed_twist = s.twist;
  }
}

void RodModelTwistingForce::updateProperties()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.twist = s.referenceTwist + rod().getEdgeTheta(s.e2) * (s.e2flip ? -1 : 1) - rod().getEdgeTheta(s.e1) * (s.e1flip ? -1 : 1);
  }
}

void RodModelTwistingForce::updateAfterRemesh()
{
  
  const std::vector<ElasticRodModel::JointStencil>& stencils = m_rod.getJointStencils();
  size_t stencilcount = stencils.size();
  if(m_stencils.size() != stencilcount) {
    m_stencils.clear();
    for(size_t i = 0; i < stencilcount; ++i) {
      Stencil stencil(stencils[i]);
      stencil.copyData(stencils[i]);
      m_stencils.push_back(stencil);
    }
  }

  updateProperties();
  updateStiffness();
  updateViscousReferenceStrain();
  computeReferenceStrain();

}

void RodModelTwistingForce::computeReferenceStrain()
{
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    Stencil & s = m_stencils[i];
    s.undeformed_twist = s.twist;
    s.reference_voronoi_length = s.voronoiLength;
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
