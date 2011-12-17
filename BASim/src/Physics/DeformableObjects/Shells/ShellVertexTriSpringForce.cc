/**
 * \file ShellVertexTriSpringForce.cc
 *
 * \author batty@cs.columbia.edu
 * \date Nov 22, 2011
 */

#include "BASim/src/Physics/DeformableObjects/Shells/ShellVertexTriSpringForce.hh"
#include "BASim/src/Core/EigenIncludes.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/MatrixBase.hh"

using Eigen::Matrix;

namespace BASim {

ShellVertexTriSpringForce::ShellVertexTriSpringForce( 
  ElasticShell& shell, 
  const std::string& name, 
  Scalar strength )
: ElasticShellForce(shell, name), m_strength(strength)
{
  
}

bool ShellVertexTriSpringForce::gatherDOFs(const FaceHandle& fh, const VertexHandle& vh, std::vector<Vec3d>& deformed, std::vector<int>& indices) const {
  assert(deformed.size() == NumVTVerts);
  
  assert(m_shell.isFaceActive(fh));
  assert(m_shell.isVertexActive(vh));

  //extract the relevant data for the face
  FaceVertexIterator fv_it = m_shell.getDefoObj().fv_iter(fh);
  int i = 0;
  for(;fv_it; ++fv_it) {
    const VertexHandle& facevert = *fv_it;
    deformed[i] = m_shell.getVertexPosition(facevert);
    int dofBase = m_shell.getVertexDofBase(facevert);
    indices[i*3] = dofBase;
    indices[i*3+1] = dofBase+1;
    indices[i*3+2] = dofBase+2;
    ++i;
  }

  //extract the vertex data
  deformed[3] = m_shell.getVertexPosition(vh);
  int dofBase = m_shell.getVertexDofBase(vh);
  indices[9] = dofBase;
  indices[10] = dofBase+1;
  indices[11] = dofBase+2;
  
  return true;
}


template <int DO_HESS>
adreal<NumVTDof,DO_HESS,Real> VTSpringEnergy(const ShellVertexTriSpringForce& vt, const std::vector<Scalar>& deformed, Vec3d baryCoords, Scalar strength) {  

  // typedefs to simplify code below
  typedef adreal<NumVTDof,DO_HESS,Real> adrealVT;
  typedef CVec3T<adrealVT> advecVT;

  Vector3d* s_deformed = (Vector3d*)(&deformed[0]);
  
  // indep variables
  advecVT   p[NumVTVerts]; // vertex positions
  set_independent( p[0], s_deformed[0], 0 );
  set_independent( p[1], s_deformed[1], 3 );
  set_independent( p[2], s_deformed[2], 6 );    
  set_independent( p[3], s_deformed[3], 9 );    
    
  adrealVT e(0);
  e += 0.5f*strength*lenSq( (baryCoords[0]*p[0] + baryCoords[1]*p[1]+baryCoords[2]*p[2]) - p[3]);  

  return e;
}

Scalar ShellVertexTriSpringForce::globalEnergy() const
{
  Scalar energy = 0;
  std::vector<int> indices(NumVTDof);
  std::vector<Vec3d> deformed(NumVTVerts);

  if(m_strength == 0) return 0;

  for (unsigned int i = 0; i < m_faces.size(); ++i) {
    const FaceHandle& fh = m_faces[i];
    const VertexHandle& vh = m_vertices[i];
    
    gatherDOFs(fh, vh, deformed, indices);
    energy += elementEnergy(deformed, m_barycoords[i]);
  }
  return energy;
}

void ShellVertexTriSpringForce::globalForce( VecXd& force )  const
{

  std::vector<int> indices(NumVTDof);
  std::vector<Vec3d> deformed(NumVTVerts);
  Eigen::Matrix<Scalar, NumVTDof, 1> localForce;

  if(m_strength == 0) return;
  
  for (unsigned int f = 0; f < m_faces.size(); ++f) {
    const FaceHandle& fh = m_faces[f];
    const VertexHandle& vh = m_vertices[f];
   
    bool valid = gatherDOFs(fh, vh, deformed, indices);
    if(!valid) continue;

    elementForce(deformed, m_barycoords[f], localForce);
    for (unsigned int i = 0; i < indices.size(); ++i)
      force(indices[i]) += localForce(i);
   
  }
  
}

void ShellVertexTriSpringForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  std::vector<int> indices(NumVTDof);
  std::vector<Vec3d> deformed(NumVTVerts);
  Eigen::Matrix<Scalar, NumVTDof, NumVTDof> localMatrix;

  if(m_strength == 0) return;

  for (unsigned int f = 0; f < m_faces.size(); ++f) {
    const FaceHandle& fh = m_faces[f];
    const VertexHandle& vh = m_vertices[f];

    bool valid = gatherDOFs(fh, vh, deformed, indices);
    if(!valid) continue;

    elementJacobian(deformed, m_barycoords[f], localMatrix);
    for (unsigned int i = 0; i < indices.size(); ++i)
      for(unsigned int j = 0; j < indices.size(); ++j)
        Jacobian.add(indices[i], indices[j], scale * localMatrix(i,j));

  }
  
}


Scalar ShellVertexTriSpringForce::elementEnergy(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords) const
{
  
  std::vector<Scalar> deformed_data(NumVTDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  adreal<NumVTDof,0,Real> e = VTSpringEnergy<0>( *this, deformed_data, baryCoords, m_strength );
  Scalar energy = e.value();

  return energy;
}

void ShellVertexTriSpringForce::elementForce(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords,
                                    Eigen::Matrix<Scalar, NumVTDof, 1>& force) const
{
  assert(deformed.size() == NumVTVerts);
 
  std::vector<Scalar> deformed_data(NumVTDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  //AutoDiff version
  adreal<NumVTDof,0,Real> e = VTSpringEnergy<0>(*this, deformed_data, baryCoords, m_strength);     
  for( uint i = 0; i < NumVTDof; i++ )
  {
    force[i] = -e.gradient(i);
  }
  
 
  
}

void ShellVertexTriSpringForce::elementJacobian(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords,
                                       Eigen::Matrix<Scalar,NumVTDof,NumVTDof>& jac) const
{
  assert(deformed.size() == NumVTVerts);

  Eigen::Matrix<Scalar,NumVTDof,NumVTDof> jac2;
  std::vector<Scalar> deformed_data(NumVTDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  jac2.setZero();

  adreal<NumVTDof,1,Real> e = VTSpringEnergy<1>(*this, deformed_data, baryCoords, m_strength);     
  // insert in the element jacobian matrix
  for( uint i = 0; i < NumVTDof; i++ )
  {
    for( uint j = 0; j < NumVTDof; j++ )
    {
      jac2(i,j) = -e.hessian(i,j);
    }
  }

}

void ShellVertexTriSpringForce::addSpring(const FaceHandle& fh, const VertexHandle& vh, const Vec3d& baryCoords) {
  m_faces.push_back(fh);
  m_vertices.push_back(vh);
  m_barycoords.push_back(baryCoords);
}


} //namespace BASim