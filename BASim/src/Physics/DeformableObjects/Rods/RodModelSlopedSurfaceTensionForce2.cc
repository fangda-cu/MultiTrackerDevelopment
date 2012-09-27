#include "RodModelSlopedSurfaceTensionForce2.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Math/Math.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"

//auto-diff stuff
#include "BASim/src/Math/ADT/adreal.h"
#include "BASim/src/Math/ADT/advec.h"

namespace BASim {
typedef Scalar Real;
typedef CVec3T<Real> Vector3d;
}

using namespace BASim;

RodModelSlopedSurfaceTensionForce2::RodModelSlopedSurfaceTensionForce2(ElasticRodModel & rod, const std::vector<ElasticRodModel::ThreeEdgeStencil> & stencils, Scalar surface_tension_coeff) :
  RodModelForce(rod, 0, "RodModelSlopedSurfaceTensionForce2"),
  m_stencils(),
  m_surface_tension_coeff(surface_tension_coeff)
{
  for (size_t i = 0; i < stencils.size(); i++)
  {
    Stencil s(stencils[i]);
    m_stencils.push_back(s);
  }

}

RodModelSlopedSurfaceTensionForce2::~RodModelSlopedSurfaceTensionForce2  ()
{
  
}

Scalar RodModelSlopedSurfaceTensionForce2::globalEnergy()
{
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    energy += localEnergy(m_stencils[i]);
  }
  return energy;
}

void RodModelSlopedSurfaceTensionForce2::globalForce(VecXd & force)
{
  ElementForce localforce;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    localForce(localforce, m_stencils[i]);
    for (size_t j = 0; j < 12; j++)
      force(m_stencils[i].dofindices[j]) += localforce(j);
  }
}

void RodModelSlopedSurfaceTensionForce2::globalJacobian(Scalar scale, MatrixBase & Jacobian)
{
  ElementJacobian localjacobian;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    localJacobian(localjacobian, m_stencils[i]);
    Jacobian.add(m_stencils[i].dofindices, m_stencils[i].dofindices, scale * localjacobian);
  }
}

template <int DO_HESS>
adreal<12,DO_HESS,Real> STEnergy(const std::vector<Scalar>& deformed, const Vec3d& volumes, bool prevExists, bool nextExists, Scalar surf_coeff) {  

  // typedefs to simplify code below
  typedef adreal<12,DO_HESS,Real> adrealST;
  typedef CVec3T<adrealST> advecST;

  Vector3d* s_deformed = (Vector3d*)(&deformed[0]);

  // indep variables
  advecST   p[4]; // vertex positions
  set_independent( p[0], s_deformed[0], 0 );
  set_independent( p[1], s_deformed[1], 3 );
  set_independent( p[2], s_deformed[2], 6 );    
  set_independent( p[3], s_deformed[3], 9 );    

  //get edge lengths and volumes
  adrealST len1 = len(p[2] - p[1]);
  adrealST rad1 = sqrt(volumes[1] / len1 / M_PI);

  adrealST len0, len2;
  adrealST rad0, rad2;
  
  adrealST radEnd0, radEnd1;
  if(prevExists) {
    len0 = len(p[1] - p[0]);
    rad0 = sqrt(volumes[0] / len0 / M_PI);
    radEnd0 = 0.5*(rad0+rad1);
  }
  else {
    len0 = 0;
    rad0 = 0;
    radEnd0 = rad1;
  }
  
  if(nextExists) {
    len2 = len(p[3] - p[2]);
    rad2 = sqrt(volumes[2] / len2 / M_PI);
    radEnd1 = 0.5*(rad1+rad2);
  }
  else {
    radEnd1 = rad1;
    len2 = 0;
    rad2 = 0;
  }
 
  //use the truncated cone formula
  adrealST e(0);
  e += surf_coeff * M_PI * (radEnd0 + radEnd1) * sqrt(len1*len1 + (radEnd0-radEnd1)*(radEnd0-radEnd1));

  //add end caps
  if(!prevExists)
    e += surf_coeff * M_PI * radEnd0;
  
  if(!nextExists)
    e+= surf_coeff * M_PI * radEnd1;

  return e;
}

void RodModelSlopedSurfaceTensionForce2::collectData(const Stencil& s, std::vector<Scalar>& deformed_data, Vec3d& vols, bool& prevValid, bool& nextValid) {
  
  deformed_data.resize(12);

  prevValid = s.e0.isValid();
  nextValid = s.e2.isValid();

  //left edge
  if(prevValid) {
    vols[0]  = rod().getVolume(s.e0);

    Vec3d x0 = rod().getDefoObj().getVertexPosition(s.v0);
    deformed_data[0] = x0[0]; deformed_data[1] = x0[1]; deformed_data[2] = x0[2];
  }

  //central edge
  vols[1] = rod().getVolume(s.e1);

  Vec3d x1 = rod().getDefoObj().getVertexPosition(s.v1);
  deformed_data[3] = x1[0]; deformed_data[4] = x1[1]; deformed_data[5] = x1[2];

  Vec3d x2 = rod().getDefoObj().getVertexPosition(s.v2);
  deformed_data[6] = x2[0]; deformed_data[7] = x2[1]; deformed_data[8] = x2[2];

  //right edge
  if(nextValid) {
    vols[2] = rod().getVolume(s.e2);

    Vec3d x3 = rod().getDefoObj().getVertexPosition(s.v3);
    deformed_data[9] = x3[0]; deformed_data[10] = x3[1]; deformed_data[11] = x3[2];
  }

}
Scalar RodModelSlopedSurfaceTensionForce2::localEnergy(Stencil & s)
{
 
  std::vector<Scalar> deformed_data;
  bool prevValid, nextValid;
  Vec3d vols;
  collectData(s, deformed_data, vols, prevValid, nextValid);
  
  adreal<12,0,Real> e = STEnergy<0>(deformed_data, vols, prevValid, nextValid, m_surface_tension_coeff);     
  
  return e.value();
}

void RodModelSlopedSurfaceTensionForce2::localForce(ElementForce & force, Stencil & s)
{
  std::vector<Scalar> deformed_data;
  bool prevValid, nextValid;
  Vec3d vols;
  collectData(s, deformed_data, vols, prevValid, nextValid);

  adreal<12,0,Real> e = STEnergy<0>(deformed_data, vols, prevValid, nextValid, m_surface_tension_coeff);     

  for(int i = 0; i < 12; ++i)
    force[i] = -e.gradient(i);
    

}

void RodModelSlopedSurfaceTensionForce2::localJacobian(ElementJacobian & jacobian, Stencil & s)
{
  std::vector<Scalar> deformed_data;
  bool prevValid, nextValid;
  Vec3d vols;
  collectData(s, deformed_data, vols, prevValid, nextValid);

  adreal<12,1,Real> e = STEnergy<1>(deformed_data, vols, prevValid, nextValid, m_surface_tension_coeff);     

  //modify the indices to skip over every 4th one, since those are the 
  for(int i = 0; i < 12; i++) {
    for(int j = 0; j < 12; ++j) {
      jacobian(i,j) = -e.hessian(i,j);
    }
  }

  assert(isSymmetric(jacobian));
}


