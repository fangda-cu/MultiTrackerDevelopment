/**
 * \file ShellVolumeForce.cc
 *
 * \author batty@cs.columbia.edu
 * \date Dec 14, 2011
 */

#include "BASim/src/Physics/DeformableObjects/Shells/ShellVolumeForce.hh"
#include "BASim/src/Core/EigenIncludes.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/MatrixBase.hh"

using Eigen::Matrix;

namespace BASim {

ShellVolumeForce::ShellVolumeForce( 
  ElasticShell& shell, 
  const std::string& name, 
  Scalar strength )
: ElasticShellForce(shell, name), m_strength(strength)
{
  //Get the initial target volume
  Scalar energy = 0;
  std::vector<int> indices(9);
  std::vector<Vec3d> deformed(3);

  m_target_volume = 0;

  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  Scalar volume = 0;
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;

    gatherDOFs(fh, deformed, indices);

    Scalar thickness = m_shell.getThickness(fh);

    //determine the energy for this element
    m_target_volume += elementEnergy(deformed);
  }
}

bool ShellVolumeForce::gatherDOFs(const FaceHandle& fh, std::vector<Vec3d>& deformed, std::vector<int>& indices) const {
  assert(deformed.size() == 3);
  
  if(!m_shell.isFaceActive(fh)) return false;

  //extract the relevant data for the local element
  FaceVertexIterator fv_it = m_shell.getDefoObj().fv_iter(fh);
  int i = 0;
  for(;fv_it; ++fv_it) {
    const VertexHandle& vh = *fv_it;
    deformed[i] = m_shell.getVertexPosition(vh);
    int dofBase = m_shell.getDefoObj().getPositionDofBase(vh);//getVertexDofBase(vh);
    indices[i*3] = dofBase;
    indices[i*3+1] = dofBase+1;
    indices[i*3+2] = dofBase+2;
    ++i;
  }
  
  return true;
}


template <int DO_HESS>
adreal<NumVolDof,DO_HESS,Real> VolumeEnergy(const ShellVolumeForce& mn, const std::vector<Scalar>& deformed, Vec3d ref_point) {  

  // typedefs to simplify code below
  typedef adreal<NumVolDof,DO_HESS,Real> adrealST;
  typedef CVec3T<adrealST> advecST;

  CVec3T<Real> reference_point(ref_point[0], ref_point[1], ref_point[2]);

  Vector3d* s_deformed = (Vector3d*)(&deformed[0]);
  
  // indep variables
  
  advecST   p[3]; // vertex positions
  set_independent( p[0], s_deformed[0], 0 );
  set_independent( p[1], s_deformed[1], 3 );
  set_independent( p[2], s_deformed[2], 6 );    
    
  //energy = signed volume of the tetrahedron
  adrealST e(0);
  e += tripleProd(p[0] - reference_point, p[1]-reference_point, p[2] - reference_point)/6.0;

  return e;
}

void ShellVolumeForce::update() {
    computeRefPoint();
}

Scalar ShellVolumeForce::globalEnergy() const
{
  Scalar energy = 0;
  std::vector<int> indices(9);
  std::vector<Vec3d> deformed(3);

  if(m_strength == 0) return 0;

  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  Scalar volume = 0;
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;
    
    gatherDOFs(fh, deformed, indices);
    
    Scalar thickness = m_shell.getThickness(fh);
    
    //determine the energy for this element
    volume += elementEnergy(deformed);
  }

  //return 0.5 * m_strength * sqr(volume - m_target_volume);
  return 0.5 * m_strength * (volume - m_target_volume)*(volume - m_target_volume);
}

void ShellVolumeForce::globalForce( VecXd& force )  const
{

  std::vector<int> indices(9);
  std::vector<Vec3d> deformed(3);
  Eigen::Matrix<Scalar, 9, 1> localForce;
  
  Scalar volume = 0;
  if(m_strength == 0) return;
  
  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;
   
    bool valid = gatherDOFs(fh, deformed, indices);
    if(!valid) continue;
    
    volume += elementEnergy(deformed);

    elementForce(deformed, localForce);
    for (unsigned int i = 0; i < indices.size(); ++i)
      force(indices[i]) += localForce(i);
  }
  force *= m_strength * (volume - m_target_volume);
  
}

void ShellVolumeForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{

  std::vector<int> indices(9);
  std::vector<Vec3d> deformed(3);
  Eigen::Matrix<Scalar, 9, 9> localMatrix;

  if(m_strength == 0) return;
  
  Scalar volume = 0;
  
  Scalar factor = m_strength * (volume - m_target_volume);

  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;

    bool valid = gatherDOFs(fh, deformed, indices);
    if(!valid) continue;
    
    volume += elementEnergy(deformed);
    elementJacobian(deformed, localMatrix);
    for (unsigned int i = 0; i < indices.size(); ++i)
      for(unsigned int j = 0; j < indices.size(); ++j)
        Jacobian.add(indices[i], indices[j], factor * scale * localMatrix(i,j));

  }
  
  //TODO This Jacobian isn't quite right (lacks the second term), 
  //I guess that makes it a Gauss-Newton approximation?
 
  
}

void ShellVolumeForce::computeRefPoint() {
  VertexIterator vit = m_shell.getDefoObj().vertices_begin();
  Vec3d sum;
  int count = 0;
  for(;vit != m_shell.getDefoObj().vertices_end(); ++vit) {
    VertexHandle v = *vit;
    Vec3d position = m_shell.getVertexPosition(v);
    sum += position;
    ++count;
  }
  m_ref_point = sum / (Scalar)count;

}


Scalar ShellVolumeForce::elementEnergy(const std::vector<Vec3d>& deformed) const
{
  
  std::vector<Scalar> deformed_data(NumVolDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  adreal<NumVolDof,0,Real> e = VolumeEnergy<0>( *this, deformed_data, m_ref_point );
  Scalar energy = e.value();

  return energy;
}

void ShellVolumeForce::elementForce(const std::vector<Vec3d>& deformed,
                                    Eigen::Matrix<Scalar, 9, 1>& force) const
{
  assert(deformed.size() == 3);

  std::vector<Scalar> deformed_data(NumVolDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  //AutoDiff version
  adreal<NumVolDof,0,Real> e = VolumeEnergy<0>(*this, deformed_data, m_ref_point);     
  for( uint i = 0; i < NumVolDof; i++ )
  {
    force[i] = -e.gradient(i);
  }
  
}

void ShellVolumeForce::elementJacobian(const std::vector<Vec3d>& deformed,
                                       Eigen::Matrix<Scalar,9,9>& jac) const
{
  assert(deformed.size() == 3);

  std::vector<Scalar> deformed_data(NumVolDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  jac.setZero();

  adreal<NumVolDof,1,Real> e = VolumeEnergy<1>(*this, deformed_data, m_ref_point);     
  // insert in the element jacobian matrix
  for( uint i = 0; i < NumVolDof; i++ )
  {
    for( uint j = 0; j < NumVolDof; j++ )
    {
      jac(i,j) = -e.hessian(i,j);
    }
  }

}



} //namespace BASim