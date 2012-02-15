/**
 * \file ShellVertexPointSpringForce.cc
 *
 * \author batty@cs.columbia.edu
 * \date Dec 30, 2011
 */

#include "BASim/src/Physics/DeformableObjects/Shells/ShellVertexPointSpringForce.hh"
#include "BASim/src/Core/EigenIncludes.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/MatrixBase.hh"

using Eigen::Matrix;

namespace BASim {

ShellVertexPointSpringForce::ShellVertexPointSpringForce( 
  ElasticShell& shell, 
  const std::string& name, 
  Scalar timestep)
: ElasticShellForce(shell, name), m_timestep(timestep)
{
  
}

bool ShellVertexPointSpringForce::gatherDOFs(const VertexHandle& vh, std::vector<Vec3d>& deformed, std::vector<Vec3d>& undef_damp, std::vector<int>& indices) const {
  assert(deformed.size() == NumSpringVerts);
  assert(m_shell.isVertexActive(vh));

  //extract the vertex data
  deformed[0] = m_shell.getVertexPosition(vh);
  undef_damp[0] = m_shell.getVertexDampingUndeformed(vh);
  int dofBase = m_shell.getVertexDofBase(vh);
  indices[0] = dofBase;
  indices[1] = dofBase+1;
  indices[2] = dofBase+2;
  
  return true;
}


template <int DO_HESS>
adreal<NumSpringDofs,DO_HESS,Real> VertPointSpringEnergy(const ShellVertexPointSpringForce& vt, const std::vector<Scalar>& deformed, const Vec3d& position, Scalar strength, Scalar undef_len) {  

  // typedefs to simplify code below
  typedef adreal<NumSpringDofs,DO_HESS,Real> adrealVT;
  typedef CVec3T<adrealVT> advecVT;

  Vector3d* s_deformed = (Vector3d*)(&deformed[0]);
  Vector3d target_position(position[0], position[1], position[2]);
  // indep variables
  advecVT   p[NumSpringVerts]; // vertex positions
  set_independent( p[0], s_deformed[0], 0 );
    
  adrealVT e(0);
  //e += 0.5*strength*sqr(len(p[0] - target_position) - undef_len);  
  e += 0.5*strength*(len(p[0] - target_position) - undef_len)*(len(p[0] - target_position) - undef_len);  

  return e;
}

Scalar ShellVertexPointSpringForce::globalEnergy() const
{
  Scalar energy = 0;
  std::vector<int> indices(NumSpringDofs);
  std::vector<Vec3d> deformed(NumSpringVerts);
  std::vector<Vec3d> undef_damp(NumSpringVerts);

  for (unsigned int i = 0; i < m_vertices.size(); ++i) {
    const VertexHandle& vh = m_vertices[i];

    gatherDOFs(vh, deformed, undef_damp, indices);
    energy += elementEnergy(deformed, m_positions[i], m_stiffnesses[i], m_restlen[i]);
  }
  return energy;
}

void ShellVertexPointSpringForce::globalForce( VecXd& force )  const
{

  std::vector<int> indices(NumSpringDofs);
  std::vector<Vec3d> deformed(NumSpringVerts);
  std::vector<Vec3d> undeformed_damp(NumSpringVerts);
  Eigen::Matrix<Scalar, NumSpringDofs, 1> localForce;

  for (unsigned int f = 0; f < m_vertices.size(); ++f) {
    const VertexHandle& vh = m_vertices[f];
   
    bool valid = gatherDOFs(vh, deformed, undeformed_damp, indices);
    if(!valid) continue;

    //elastic force
    elementForce(deformed, m_positions[f], m_stiffnesses[f], m_restlen[f], localForce);
    for (unsigned int i = 0; i < indices.size(); ++i) {
      force(indices[i]) += localForce(i);
    }
    
    //viscous force
    localForce.setZero();
    Scalar damp_restlen = (undeformed_damp[0]-m_positions[f]).norm();
    elementForce(deformed, m_positions[f], m_damping[f], damp_restlen, localForce);
    for (unsigned int i = 0; i < indices.size(); ++i) {
      force(indices[i]) += (1.0 / m_timestep) * localForce(i);
    }

  }
  
}

void ShellVertexPointSpringForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  std::vector<int> indices(NumSpringDofs);
  std::vector<Vec3d> deformed(NumSpringVerts);
  std::vector<Vec3d> undeformed_damp(NumSpringVerts);

  Eigen::Matrix<Scalar, NumSpringDofs, NumSpringDofs> localMatrix;

  for (unsigned int f = 0; f < m_vertices.size(); ++f) {
    const VertexHandle& vh = m_vertices[f];

    bool valid = gatherDOFs(vh, deformed, undeformed_damp, indices);
    if(!valid) continue;
    
    elementJacobian(deformed, m_positions[f], m_stiffnesses[f], m_restlen[f], localMatrix);
    for (unsigned int i = 0; i < indices.size(); ++i)
      for(unsigned int j = 0; j < indices.size(); ++j)
        Jacobian.add(indices[i], indices[j], scale * localMatrix(i,j));

    //viscosity/damping
    Scalar damp_restlen = (undeformed_damp[0]-m_positions[f]).norm();
    elementJacobian(deformed, m_positions[f], m_damping[f], damp_restlen, localMatrix);
    for (unsigned int i = 0; i < indices.size(); ++i)
      for(unsigned int j = 0; j < indices.size(); ++j)
        Jacobian.add(indices[i], indices[j], scale / m_timestep * localMatrix(i,j));

  }
  
}


Scalar ShellVertexPointSpringForce::elementEnergy(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords, Scalar strength, Scalar restlen) const
{
  
  std::vector<Scalar> deformed_data(NumSpringDofs);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  adreal<NumSpringDofs,0,Real> e = VertPointSpringEnergy<0>( *this, deformed_data, baryCoords, strength, restlen);
  Scalar energy = e.value();

  return energy;
}

void ShellVertexPointSpringForce::elementForce(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords, Scalar strength, Scalar restlen,
                                    Eigen::Matrix<Scalar, NumSpringDofs, 1>& force) const
{
  assert(deformed.size() == NumSpringVerts);
 
  std::vector<Scalar> deformed_data(NumSpringDofs);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  //AutoDiff version
  adreal<NumSpringDofs,0,Real> e = VertPointSpringEnergy<0>(*this, deformed_data, baryCoords, strength, restlen);     
  for( uint i = 0; i < NumSpringDofs; i++ )
  {
    force[i] = -e.gradient(i);
  }
   
}

void ShellVertexPointSpringForce::elementJacobian(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords, Scalar strength, Scalar restlen,
                                       Eigen::Matrix<Scalar,NumSpringDofs,NumSpringDofs>& jac) const
{
  assert(deformed.size() == NumSpringVerts);

  Eigen::Matrix<Scalar,NumSpringDofs,NumSpringDofs> jac2;
  std::vector<Scalar> deformed_data(NumSpringDofs);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  jac2.setZero();

  adreal<NumSpringDofs,1,Real> e = VertPointSpringEnergy<1>(*this, deformed_data, baryCoords, strength, restlen);     
  // insert in the element jacobian matrix
  for( uint i = 0; i < NumSpringDofs; i++ )
  {
    for( uint j = 0; j < NumSpringDofs; j++ )
    {
      jac2(i,j) = -e.hessian(i,j);
    }
  }

}

void ShellVertexPointSpringForce::addSpring(const VertexHandle& vh, const Vec3d& position,
                                          Scalar stiffness, Scalar damping, Scalar restlen) {

  m_vertices.push_back(vh);
  m_positions.push_back(position);
  m_stiffnesses.push_back(stiffness);
  m_damping.push_back(damping);
  m_restlen.push_back(restlen);

}

bool ShellVertexPointSpringForce::hasSpring(const VertexHandle& vh) {
  return std::find(m_vertices.begin(), m_vertices.end(), vh) != m_vertices.end();
}


} //namespace BASim