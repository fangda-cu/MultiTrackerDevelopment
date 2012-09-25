#include "RodModelSlopedSurfaceTensionForce.hh"
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

RodModelSlopedSurfaceTensionForce::RodModelSlopedSurfaceTensionForce(ElasticRodModel & rod, const std::vector<ElasticRodModel::JointStencil> & stencils, Scalar surface_tension_coeff) :
  RodModelForce(rod, 0, "RodModelSlopedSurfaceTensionForce"),
  m_stencils(),
  m_surface_tension_coeff(surface_tension_coeff)
{
  for (size_t i = 0; i < stencils.size(); i++)
  {
    Stencil s(stencils[i]);
    m_stencils.push_back(s);
  }

}

RodModelSlopedSurfaceTensionForce::~RodModelSlopedSurfaceTensionForce  ()
{
  
}

Scalar RodModelSlopedSurfaceTensionForce::globalEnergy()
{
  Scalar energy = 0;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    energy += localEnergy(m_stencils[i]);
  }
  return energy;
}

void RodModelSlopedSurfaceTensionForce::globalForce(VecXd & force)
{
  ElementForce localforce;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    localForce(localforce, m_stencils[i]);
    for (size_t j = 0; j < m_stencils[i].dofindices.size(); j++)
      force(m_stencils[i].dofindices[j]) += localforce(j);
  }
}

void RodModelSlopedSurfaceTensionForce::globalJacobian(Scalar scale, MatrixBase & Jacobian)
{
  ElementJacobian localjacobian;
  for (size_t i = 0; i < m_stencils.size(); i++)
  {
    localJacobian(localjacobian, m_stencils[i]);
    Jacobian.add(m_stencils[i].dofindices, m_stencils[i].dofindices, scale * localjacobian);
  }
}

template <int DO_HESS>
adreal<9,DO_HESS,Real> STEnergy(const std::vector<Scalar>& deformed, const Vec2d& volumes, Scalar surf_coeff) {  

  // typedefs to simplify code below
  typedef adreal<9,DO_HESS,Real> adrealST;
  typedef CVec3T<adrealST> advecST;

  Vector3d* s_deformed = (Vector3d*)(&deformed[0]);

  // indep variables
  advecST   p[3]; // vertex positions
  set_independent( p[0], s_deformed[0], 0 );
  set_independent( p[1], s_deformed[1], 3 );
  set_independent( p[2], s_deformed[2], 6 );    

  adrealST len0 = len(p[1] - p[0]);
  adrealST len1 = len(p[2] - p[1]);
  adrealST rad0 = sqrt(volumes[0] / len0 / M_PI);
  adrealST rad1 = sqrt(volumes[1] / len1 / M_PI);
  
  adrealST midlen = 0.5*(len0+len1);
  
  //truncated cone formula
  adrealST e(0);
  e += surf_coeff * M_PI * (rad0 + rad1) * sqrt(midlen*midlen + (rad0-rad1)*(rad0-rad1));

  return e;
}


Scalar RodModelSlopedSurfaceTensionForce::localEnergy(Stencil & s)
{
 
  Scalar vol0 = rod().getVolume(s.e1);
  Scalar vol1 = rod().getVolume(s.e2);
  Vec2d vols(vol0, vol1);
  Vec3d x0 = rod().getDefoObj().getVertexPosition(s.v1);
  Vec3d x1 = rod().getDefoObj().getVertexPosition(s.v2);
  Vec3d x2 = rod().getDefoObj().getVertexPosition(s.v3);

  std::vector<Scalar> deformed_data(9);
  deformed_data[0] = x0[0]; deformed_data[1] = x0[1];  deformed_data[2] = x0[2]; 
  deformed_data[3] = x1[0]; deformed_data[4] = x1[1];  deformed_data[5] = x1[2]; 
  deformed_data[6] = x2[0]; deformed_data[7] = x2[1];  deformed_data[8] = x2[2]; 

  adreal<9,0,Real> e = STEnergy<0>(deformed_data, vols, m_surface_tension_coeff);     
  
  return e.value();
}

void RodModelSlopedSurfaceTensionForce::localForce(ElementForce & force, Stencil & s)
{
  Scalar vol0 = rod().getVolume(s.e1);
  Scalar vol1 = rod().getVolume(s.e2);
  Vec2d vols(vol0, vol1);
  Vec3d x0 = rod().getDefoObj().getVertexPosition(s.v1);
  Vec3d x1 = rod().getDefoObj().getVertexPosition(s.v2);
  Vec3d x2 = rod().getDefoObj().getVertexPosition(s.v3);

  std::vector<Scalar> deformed_data(9);
  deformed_data[0] = x0[0]; deformed_data[1] = x0[1];  deformed_data[2] = x0[2]; 
  deformed_data[3] = x1[0]; deformed_data[4] = x1[1];  deformed_data[5] = x1[2]; 
  deformed_data[6] = x2[0]; deformed_data[7] = x2[1];  deformed_data[8] = x2[2]; 

  adreal<9,0,Real> e = STEnergy<0>(deformed_data, vols, m_surface_tension_coeff);     
  for(int i = 0; i < 3; ++i)
    for(int j = 0; j < 3; ++j)
      force[i*4 + j] = -e.gradient(i*3 + j);
    

}

void RodModelSlopedSurfaceTensionForce::localJacobian(ElementJacobian & jacobian, Stencil & s)
{
  Scalar vol0 = rod().getVolume(s.e1);
  Scalar vol1 = rod().getVolume(s.e2);
  Vec2d vols(vol0, vol1);
  Vec3d x0 = rod().getDefoObj().getVertexPosition(s.v1);
  Vec3d x1 = rod().getDefoObj().getVertexPosition(s.v2);
  Vec3d x2 = rod().getDefoObj().getVertexPosition(s.v3);

  std::vector<Scalar> deformed_data(9);
  deformed_data[0] = x0[0]; deformed_data[1] = x0[1];  deformed_data[2] = x0[2]; 
  deformed_data[3] = x1[0]; deformed_data[4] = x1[1];  deformed_data[5] = x1[2]; 
  deformed_data[6] = x2[0]; deformed_data[7] = x2[1];  deformed_data[8] = x2[2]; 

  adreal<9,1,Real> e = STEnergy<1>(deformed_data, vols, m_surface_tension_coeff);     
 
  //modify the indices to skip over every 4th one, since those are the 
  for( uint i = 0; i < 3; i++) {
    for( uint j = 0; j < 3; j++) {
      uint indLeft0 = i*4+j;
      uint indRight0 = i*3+j;
      for(uint k = 0; k < 3; ++k) {
        for(uint m = 0; m < 3; ++m) {
          uint indLeft1 = k*4+m;
          uint indRight1 = k*3+m;
          jacobian(indLeft0,indLeft1) = -e.hessian(indRight0,indRight1);
        }
      }
    }
  }

  assert(isSymmetric(jacobian));
}


