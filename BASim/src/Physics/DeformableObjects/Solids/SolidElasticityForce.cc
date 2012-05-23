/**
 * \file SolidElasticityForce.cc
 *
 * \author batty@cs.columbia.edu
 * \date May 23, 2012
 */

#include "BASim/src/Physics/DeformableObjects/Solids/SolidElasticityForce.hh"
#include "BASim/src/Core/EigenIncludes.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/MatrixBase.hh"

using Eigen::Matrix;

namespace BASim {

SolidElasticityForce::SolidElasticityForce( 
  ElasticSolid& solid, 
  const std::string& name, 
  Scalar Youngs,
  Scalar Poisson,
  Scalar Youngs_damp, 
  Scalar Poisson_damp,
  Scalar timestep)
: ElasticSolidForce(solid, name), m_Poisson(Poisson), m_Youngs(Youngs), m_Youngs_damp(Youngs_damp), m_Poisson_damp(Poisson_damp), m_timestep(timestep)
{
  
}

bool SolidElasticityForce::gatherDOFs(const TetHandle& th, std::vector<Vec3d>& deformed, std::vector<Vec3d>& undeformed_damp, std::vector<Vec3d>& undeformed, std::vector<int>& indices) const {
  assert(deformed.size() == 4);
  
  if(!m_solid.isTetActive(th)) return false;

  //extract the relevant data for the local element
  TetVertexIterator tv_it = m_solid.getDefoObj().tv_iter(th);
  int i = 0;
  for(;tv_it; ++tv_it) {
    const VertexHandle& vh = *tv_it;
    deformed[i] = m_solid.getVertexPosition(vh);
    undeformed[i] = m_solid.getVertexUndeformed(vh);
    undeformed_damp[i] = m_solid.getVertexDampingUndeformed(vh);
    int dofBase = m_solid.getVertexDofBase(vh);
    indices[i*3] = dofBase;
    indices[i*3+1] = dofBase+1;
    indices[i*3+2] = dofBase+2;
    ++i;
  }
  
  return true;
}


template <int DO_HESS>
adreal<NumElasticityDof,DO_HESS,Real> ElasticEnergy(const SolidElasticityForce& mn, const std::vector<Scalar>& deformed, const std::vector<Scalar>& undeformed, Scalar Youngs, Scalar Poisson) {  

  // typedefs to simplify code below
  typedef adreal<NumElasticityDof,DO_HESS,Real> adrealElast;
  typedef CVec3T<adrealElast> advecElast;

  Vector3d* s_deformed = (Vector3d*)(&deformed[0]);
  
  // indep variables
  advecElast   p[4]; // vertex positions
  set_independent( p[0], s_deformed[0], 0 );
  set_independent( p[1], s_deformed[1], 3 );
  set_independent( p[2], s_deformed[2], 6 );    
  set_independent( p[3], s_deformed[3], 9 );    
    
  adrealElast e(0);
  //TODO Fill in elastic force.
  
  return e;
}

Scalar SolidElasticityForce::globalEnergy() const
{
  Scalar energy = 0;
  std::vector<int> indices(12);
  std::vector<Vec3d> deformed(4), undeformed(4), undeformed_damp(4);

  TetIterator tit = m_solid.getDefoObj().tets_begin();
  for (;tit != m_solid.getDefoObj().tets_end(); ++tit) {
    const TetHandle& th = *tit;
    
    gatherDOFs(th, deformed, undeformed_damp, undeformed, indices);
    
    //determine the energy for this element
    energy += elementEnergy(deformed, undeformed, m_Youngs, m_Poisson);
  }
  return energy;
}

void SolidElasticityForce::globalForce( VecXd& force )  const
{

  std::vector<int> indices(12);
  std::vector<Vec3d> deformed(4), undeformed(4), undeformed_damp(4);
  Eigen::Matrix<Scalar, 12, 1> localForce;

  TetIterator tit = m_solid.getDefoObj().tets_begin();
  for (;tit != m_solid.getDefoObj().tets_end(); ++tit) {
    const TetHandle& th = *tit;
   
    bool valid = gatherDOFs(th, deformed, undeformed_damp, undeformed, indices);
    if(!valid) continue;

    if(m_Youngs != 0) {
      elementForce(deformed, undeformed, m_Youngs, m_Poisson, localForce);
      for (unsigned int i = 0; i < indices.size(); ++i)
        force(indices[i]) += localForce(i);
    }
    
    if(m_Youngs_damp != 0) {
      elementForce(deformed, undeformed_damp, m_Youngs_damp, m_Poisson_damp, localForce);
      for (unsigned int i = 0; i < indices.size(); ++i)
        force(indices[i]) += localForce(i) / m_timestep;  //division by timestep to do Viscous Threads-style viscosity/damping
    }
   
  }
  
}

void SolidElasticityForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  std::vector<int> indices(12);
  std::vector<Vec3d> deformed(4), undeformed(4), undeformed_damp(4);
  Eigen::Matrix<Scalar, 12, 12> localMatrix;

  TetIterator tit = m_solid.getDefoObj().tets_begin();
  for (;tit != m_solid.getDefoObj().tets_end(); ++tit) {
    const TetHandle& th = *tit;

    bool valid = gatherDOFs(th, deformed, undeformed_damp, undeformed, indices);
    if(!valid) continue;

    if(m_Youngs != 0) {
      elementJacobian(deformed, undeformed, m_Youngs, m_Poisson, localMatrix);
      Jacobian.add(indices, indices, scale * localMatrix);
    }
    if(m_Youngs_damp != 0) {
      elementJacobian(deformed, undeformed_damp, m_Youngs_damp, m_Poisson_damp, localMatrix);
      Jacobian.add(indices, indices, scale / m_timestep * localMatrix);
    }

  }
  
}


Scalar SolidElasticityForce::elementEnergy(const std::vector<Vec3d>& deformed, const std::vector<Vec3d>& undeformed, 
  Scalar Youngs, Scalar Poisson) const
{
  
  std::vector<Scalar> deformed_data(NumElasticityDof), undeformed_data(NumElasticityDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
    
    undeformed_data[3*i] = undeformed[i][0];
    undeformed_data[3*i+1] = undeformed[i][1];
    undeformed_data[3*i+2] = undeformed[i][2];
  }

  adreal<NumElasticityDof,0,Real> e = ElasticEnergy<0>( *this, deformed_data, undeformed_data, Youngs, Poisson);
  Scalar energy = e.value();

  return energy;
}

void SolidElasticityForce::elementForce(const std::vector<Vec3d>& deformed, const std::vector<Vec3d>& undeformed,
                                    Scalar Youngs, Scalar Poisson, 
                                    Eigen::Matrix<Scalar, 12, 1>& force) const
{
  assert(deformed.size() == 4);
 
  std::vector<Scalar> deformed_data(NumElasticityDof), undeformed_data(NumElasticityDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];

    undeformed_data[3*i] = undeformed[i][0];
    undeformed_data[3*i+1] = undeformed[i][1];
    undeformed_data[3*i+2] = undeformed[i][2];
  }

  //AutoDiff version
  adreal<NumElasticityDof,0,Real> e = ElasticEnergy<0>(*this, deformed_data, undeformed_data, Youngs, Poisson);     
  for( uint i = 0; i < NumElasticityDof; i++ )
  {
    force[i] = -e.gradient(i);
  }
  
}

void SolidElasticityForce::elementJacobian(const std::vector<Vec3d>& deformed, const std::vector<Vec3d>& undeformed,
                                        Scalar Youngs, Scalar Poisson, 
                                        Eigen::Matrix<Scalar,12,12>& jac) const
{
  assert(deformed.size() == 4);

  Eigen::Matrix<Scalar,12,12> jac2;
  std::vector<Scalar> deformed_data(NumElasticityDof);
  std::vector<Scalar> undeformed_data(NumElasticityDof);

  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  
    undeformed_data[3*i] = undeformed[i][0];
    undeformed_data[3*i+1] = undeformed[i][1];
    undeformed_data[3*i+2] = undeformed[i][2];
  }

  jac2.setZero();

  adreal<NumElasticityDof,1,Real> e = ElasticEnergy<1>(*this, deformed_data, undeformed_data, Youngs, Poisson);     
  // insert in the element jacobian matrix
  for( uint i = 0; i < NumElasticityDof; i++ )
  {
    for( uint j = 0; j < NumElasticityDof; j++ )
    {
      jac2(i,j) = -e.hessian(i,j);
    }
  }

  
}



} //namespace BASim