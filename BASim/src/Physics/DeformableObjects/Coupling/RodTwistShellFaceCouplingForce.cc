#include "RodTwistShellFaceCouplingForce.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"

using namespace BASim;

RodTwistShellFaceCouplingForce::RodTwistShellFaceCouplingForce(ElasticRodModel & rod, ElasticShell & shell, const std::vector<ElasticRodModel::EdgeStencil> & stencils, Scalar stiffness, Scalar stiffness_damp, Scalar timestep) :
  DefoObjForce(rod.getDefoObj(), timestep, "RodTwistShellFaceCouplingForce"),
  m_stencils(),
  m_stiffness(stiffness),
  m_stiffness_damp(stiffness_damp)
{
  for (size_t i = 0; i < stencils.size(); i++)
  {
    Stencil s(stencils[i]);
    s.stiffness = 0;
    s.viscous_stiffness = 0;
    s.undeformed_length = 0;
    s.damping_undeformed_length = 0;
    
    m_stencils.push_back(s);
  }

  updateProperties();
  updateStiffness();
  updateViscousReferenceStrain();
  computeReferenceStrain();
}

RodTwistShellFaceCouplingForce::~RodTwistShellFaceCouplingForce()
{
  
}

template <int DO_HESS>
adreal<RodTwistShellFaceCouplingForce::NumDof,DO_HESS,Scalar> 
RodTwistShellFaceCouplingForce::adEnergy(const RodTwistShellFaceCouplingForce & mn, const std::vector<Scalar> & deformed, const std::vector<Scalar> & undeformed) 
{  
  
  // typedefs to simplify code below
  typedef adreal<RodTwistShellFaceCouplingForce::NumDof,DO_HESS,Scalar> adrealElast;
  typedef CVec3T<adrealElast> advecElast;
  Mat3T<adrealElast> temp;
  typedef Mat3T<adrealElast> admatElast;
//  
//  Vector3d* s_deformed = (Vector3d*)(&deformed[0]);
//  Vector3d* s_undeformed = (Vector3d*)(&undeformed[0]);
//  
//  // indep variables
//  advecElast   p[4]; // vertex positions
//  set_independent( p[0], s_deformed[0], 0 );
//  set_independent( p[1], s_deformed[1], 3 );
//  set_independent( p[2], s_deformed[2], 6 );    
//  set_independent( p[3], s_deformed[3], 9 );    
//  
//  //dependent variable
  adrealElast e(0);
//  
//  //Compute green strain, following Teran 2003 's formulation  ("Finite Volume Methods for the Simulation of Skeletal Muscle")
//  advecElast ds1 = p[1]-p[0];
//  advecElast ds2 = p[2]-p[0];
//  advecElast ds3 = p[3]-p[0];
//  
//  Vector3d dm1 = s_undeformed[1]-s_undeformed[0];
//  Vector3d dm2 = s_undeformed[2]-s_undeformed[0];
//  Vector3d dm3 = s_undeformed[3]-s_undeformed[0];
//  
//  admatElast ds_mat(ds1[0], ds2[0], ds3[0],
//                    ds1[1], ds2[1], ds3[1],
//                    ds1[2], ds2[2], ds3[2]);
//  
//  admatElast dm_mat(dm1[0], dm2[0], dm3[0],
//                    dm1[1], dm2[1], dm3[1],
//                    dm1[2], dm2[2], dm3[2]);
//  
//  admatElast dm_inv = dm_mat.inverse();
//  
//  // Compute deformation gradient
//  admatElast F = ds_mat * dm_inv;  
//  
//  // Compute green strain
//  admatElast G = F.transpose()*F - admatElast::Identity(); 
//  
//  // Convert from Youngs+Poisson to lame coeffs
//  adrealElast intermediate = Youngs / (1.0+Poisson);
//  adrealElast lambda = intermediate * Poisson / (1.0 - 2.0*Poisson);
//  adrealElast mu = intermediate / 2.0;
//  
//  //Compute stress using linear elasticity
//  admatElast CG = adrealElast(2.0 * mu) * G + adrealElast(lambda) * G.trace() * admatElast::Identity();
//  
//  // Compute elastic potential density: double contraction of stress:strain (expressed as trace of matrix product)
//  e = adrealElast(0.5)*admatElast::doublecontraction(CG, G); 
//  
//  // Scale by (rest) volume 
//  e *= dot(dm1, cross(dm2,dm3))/6.0;
  
  return e;
}

Scalar RodTwistShellFaceCouplingForce::globalEnergy()
{
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    // energy only include non-viscous forces
    energy += localEnergy(m_stencils[i], false);
  }
  return energy;
}

void RodTwistShellFaceCouplingForce::globalForce(VecXd & force)
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

void RodTwistShellFaceCouplingForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
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

Scalar RodTwistShellFaceCouplingForce::localEnergy(Stencil & s, bool viscous)
{
//  Scalar ks = (viscous ? s.viscous_stiffness : s.stiffness);
//  
//  Scalar reflen = (viscous ? s.damping_undeformed_length : s.undeformed_length);
//  Scalar len = rod().getEdgeLength(s.e);
//  
//  return ks / 2.0 * square(len / reflen - 1.0) * reflen;
  return 0;
}

void RodTwistShellFaceCouplingForce::localForce(ElementForce & force, Stencil & s, bool viscous)
{
//  Scalar ks = (viscous ? s.viscous_stiffness : s.stiffness);
//  
//  Scalar reflen = (viscous ? s.damping_undeformed_length : s.undeformed_length);
//  Scalar len = rod().getEdgeLength(s.e);
//  Vec3d tangent = rod().getEdgeTangent(s.e);
//  
//  Vec3d f = ks * (len / reflen - 1.0) * tangent;
//  force.segment<3> (0) = f;
//  force.segment<3> (3) = -f;
}

void RodTwistShellFaceCouplingForce::localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous)
{
//  Scalar ks = (viscous ? s.viscous_stiffness : s.stiffness);
//  
//  Vec3d e = rod().getEdge(s.e);
//  Scalar reflen = (viscous ? s.damping_undeformed_length : s.undeformed_length);
//  Scalar len = rod().getEdgeLength(s.e);
//  Mat3d M = ks * ((1.0 / reflen - 1.0 / len) * Mat3d::Identity() + 1.0 / len * outerProd(e, e) / square(len));
//  
//  //TODO: this is copied from RodStretchingForce::elemengJacobian(). can't this be implemented using blocks?
//  for (int i = 0; i < 3; i++)
//  {
//    for (int j = 0; j < 3; j++)
//    {
//      jacobian(i, j) = jacobian(3 + i, 3 + j) = -M(i, j);
//      jacobian(3 + i, j) = jacobian(i, 3 + j) = M(i, j);
//    }
//  }
//  
//  assert(isSymmetric(jacobian));
}

void RodTwistShellFaceCouplingForce::updateStiffness()
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

void RodTwistShellFaceCouplingForce::updateViscousReferenceStrain()
{
//  for (size_t i = 0; i < m_stencils.size(); i++)
//  {
//    Stencil & s = m_stencils[i];
//    s.damping_undeformed_length = rod().getEdgeLength(s.e);
//  }
}

void RodTwistShellFaceCouplingForce::updateProperties()
{
  
}

void RodTwistShellFaceCouplingForce::computeReferenceStrain()
{
//  for (size_t i = 0; i < m_stencils.size(); i++)
//  {
//    Stencil & s = m_stencils[i];
//    s.undeformed_length = rod().getEdgeLength(s.e);
//  }
}

