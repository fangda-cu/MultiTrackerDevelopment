/**
 * \file ShellSurfaceTensionForce.cc
 *
 * \author batty@cs.columbia.edu
 * \date Nov 22, 2011
 */

#include "BASim/src/Physics/DeformableObjects/Shells/ShellSurfaceTensionForce.hh"
#include "BASim/src/Core/EigenIncludes.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/MatrixBase.hh"

using Eigen::Matrix;

namespace BASim {

ShellSurfaceTensionForce::ShellSurfaceTensionForce( 
  ElasticShell& shell, 
  const std::string& name, 
  Scalar surf_coeff )
: ElasticShellForce(shell, name), m_surface_tension_coeff(surf_coeff)
{
  
}

bool ShellSurfaceTensionForce::gatherDOFs(const FaceHandle& fh, std::vector<Vec3d>& deformed, std::vector<int>& indices) const {
  assert(deformed.size() == 3);
  
  if(!m_shell.isFaceActive(fh)) return false;

  //extract the relevant data for the local element
  FaceVertexIterator fv_it = m_shell.getDefoObj().fv_iter(fh);
  int i = 0;
  for(;fv_it; ++fv_it) {
    const VertexHandle& vh = *fv_it;
    deformed[i] = m_shell.getVertexPosition(vh);
    int dofBase = m_shell.getVertexDofBase(vh);
    indices[i*3] = dofBase;
    indices[i*3+1] = dofBase+1;
    indices[i*3+2] = dofBase+2;
    ++i;
  }
  
  return true;
}


template <int DO_HESS>
adreal<NumSTDof,DO_HESS,Real> STEnergy(const ShellSurfaceTensionForce& mn, const std::vector<Scalar>& deformed, Scalar surf_coeff) {  

  // typedefs to simplify code below
  typedef adreal<NumSTDof,DO_HESS,Real> adrealMN;
  typedef CVec3T<adrealMN> advecMN;

  Vector3d* s_deformed = (Vector3d*)(&deformed[0]);
  
  // indep variables
  advecMN   p[3]; // vertex positions
  set_independent( p[0], s_deformed[0], 0 );
  set_independent( p[1], s_deformed[1], 3 );
  set_independent( p[2], s_deformed[2], 6 );    
    
  //energy = 2*coeff*surface_area
  adrealMN e(0);
  e += surf_coeff*len(cross(p[1] - p[0], p[2] - p[0]));  

  return e;
}

Scalar ShellSurfaceTensionForce::globalEnergy() const
{
  Scalar energy = 0;
  std::vector<int> indices(9);
  std::vector<Vec3d> deformed(3);

  if(m_surface_tension_coeff == 0) return 0;

  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;
    
    gatherDOFs(fh, deformed, indices);
    
    Scalar thickness = m_shell.getThickness(fh);
    
    //determine the energy for this element
    energy += elementEnergy(deformed);

  }
  return energy;
}

void ShellSurfaceTensionForce::globalForce( VecXd& force )  const
{

  std::vector<int> indices(9);
  std::vector<Vec3d> deformed(3);
  Eigen::Matrix<Scalar, 9, 1> localForce;

  if(m_surface_tension_coeff == 0) return;
  std::cout << "Computing surf force\n";
  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;
   
    bool valid = gatherDOFs(fh, deformed, indices);
    if(!valid) continue;

    std::cout << "Computing face " << fit->idx() << "\n";

    elementForce(deformed, localForce);
    for (unsigned int i = 0; i < indices.size(); ++i)
      force(indices[i]) += localForce(i);
   
  }
  std::cout << "Done surf force\n";
}

void ShellSurfaceTensionForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  std::vector<int> indices(9);
  std::vector<Vec3d> deformed(3);
  Eigen::Matrix<Scalar, 9, 9> localMatrix;

  if(m_surface_tension_coeff == 0) return;

  std::cout << "Computing surf Jacob\n";

  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;

    bool valid = gatherDOFs(fh, deformed, indices);
    if(!valid) continue;

    elementJacobian(deformed, localMatrix);
    for (unsigned int i = 0; i < indices.size(); ++i)
      for(unsigned int j = 0; j < indices.size(); ++j)
        Jacobian.add(indices[i], indices[j], scale * localMatrix(i,j));

  }
  std::cout << "Done surf Jacob\n";
}


Scalar ShellSurfaceTensionForce::elementEnergy(const std::vector<Vec3d>& deformed) const
{
  
  std::vector<Scalar> deformed_data(NumSTDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  adreal<NumSTDof,0,Real> e = STEnergy<0>( *this, deformed_data, m_surface_tension_coeff );
  Scalar energy = e.value();

  return energy;
}

void ShellSurfaceTensionForce::elementForce(const std::vector<Vec3d>& deformed,
                                    Eigen::Matrix<Scalar, 9, 1>& force) const
{
  assert(deformed.size() == 3);
  std::vector<Scalar> deformed_data(NumSTDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  adreal<NumSTDof,0,Real> e = STEnergy<0>(*this, deformed_data, m_surface_tension_coeff);     
  for( uint i = 0; i < NumSTDof; i++ )
  {
    force[i] = -e.gradient(i);
  }

}

void ShellSurfaceTensionForce::elementJacobian(const std::vector<Vec3d>& deformed,
                                       Eigen::Matrix<Scalar,9,9>& jac) const
{
  assert(deformed.size() == 3);

  std::vector<Scalar> deformed_data(NumSTDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  jac.setZero();

  adreal<NumSTDof,1,Real> e = STEnergy<1>(*this, deformed_data, m_surface_tension_coeff);     
  // insert in the element jacobian matrix
  for( uint i = 0; i < NumSTDof; i++ )
  {
    for( uint j = 0; j < NumSTDof; j++ )
    {
      jac(i,j) = -e.hessian(i,j);
    }
  }
}



} //namespace BASim