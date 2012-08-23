/**
 * \file ShellLinearSurfaceTensionForce.cc
 *
 * \author batty@cs.columbia.edu
 * \date Aug 21, 2012
 */

#include "BASim/src/Physics/DeformableObjects/Shells/ShellLinearSurfaceTensionForce.hh"
#include "BASim/src/Core/EigenIncludes.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"

using Eigen::Matrix;

namespace BASim {

ShellLinearSurfaceTensionForce ::ShellLinearSurfaceTensionForce ( 
  ElasticShell& shell, 
  const std::string& name, 
  Scalar surf_coeff )
: ElasticShellForce(shell, name), m_surface_tension_coeff(surf_coeff)
{
  
}

bool ShellLinearSurfaceTensionForce ::gatherDOFs(const EdgeHandle& eh, std::vector<Vec3d>& deformed, FaceHandle& f1, FaceHandle& f2, std::vector<int>& indices) const {
  assert(deformed.size() == 4);
  
  if(!m_shell.isEdgeActive(eh)) return false;

  //skip neighbourless and non-manifold edges
  int fCount = m_shell.getDefoObj().edgeIncidentFaces(eh);
  if(fCount < 1 || fCount > 2) return false; 

  //extract the relevant data for the local edge
  EdgeVertexIterator ev_it = m_shell.getDefoObj().ev_iter(eh);
  int i = 0;
  for(;ev_it; ++ev_it) {
    const VertexHandle& vh = *ev_it;
    deformed[i] = m_shell.getVertexPosition(vh);
    int dofBase = m_shell.getDefoObj().getPositionDofBase(vh);
    indices[i*3] = dofBase;
    indices[i*3+1] = dofBase+1;
    indices[i*3+2] = dofBase+2;
    ++i;
  }
  
  f2 = FaceHandle(-1); //mark as invalid by default

  EdgeFaceIterator ef_it = m_shell.getDefoObj().ef_iter(eh);
  for(;ef_it; ++ef_it) {
     FaceHandle fh = *ef_it;
     
     if(i==2) f1 = fh;
     else if (i==3) f2 = fh;

     VertexHandle vh;
     bool result = getFaceThirdVertex(m_shell.getDefoObj(), fh, eh, vh);
     deformed[i] = m_shell.getVertexPosition(vh);
     int dofBase = m_shell.getDefoObj().getPositionDofBase(vh);
     indices[i*3] = dofBase;
     indices[i*3+1] = dofBase+1;
     indices[i*3+2] = dofBase+2;
     ++i;
  }
   
  return true;
}


template <int DO_HESS>
adreal<NumLSTDof,DO_HESS,Real> LSTEnergy(const ShellLinearSurfaceTensionForce& mn, 
                                         const std::vector<Scalar>& deformed, 
                                         Scalar surf_coeff, 
                                         Scalar vol1, Scalar vol2) {  

  // typedefs to simplify code below
  typedef adreal<NumLSTDof,DO_HESS,Real> adrealST;
  typedef CVec3T<adrealST> advecST;

  Vector3d* s_deformed = (Vector3d*)(&deformed[0]);
  
  // indep variables
  advecST   p[4]; // vertex positions
  set_independent( p[0], s_deformed[0], 0 );
  set_independent( p[1], s_deformed[1], 3 );
  set_independent( p[2], s_deformed[2], 6 );  
  if(vol2 != 0) {
    set_independent( p[3], s_deformed[3], 9 );
  }

  //edge length
  adrealST le = len(p[1] - p[0]);
  
  //adjacent triangle areas
  adrealST A1 = 0.5*len(cross(p[1] - p[0], p[2] - p[0]));

  adrealST A2 = 0;
  if(vol2 != 0)
    A2 += 0.5*len(cross(p[1] - p[0], p[3] - p[0]));
  
  //thickness values dependent on vertex positions
  adrealST h1 = vol1 / A1;
  adrealST h2 = 0;
  if(vol2 != 0)
    h2 += vol2 / A2;
  
  adrealST e(0);
  if(vol2 != 0) { //both sides
    e += surf_coeff * sqrt( 4.0/9.0*(A1+A2)*(A1+A2) + le*le*(h2-h1)*(h2-h1)/4.0 );
  }
  else {          //only one side
    e += surf_coeff * sqrt( 4.0/9.0*A1*A1 + le*le*h1*h1/4.0 );
  }

  return e;
}

Scalar ShellLinearSurfaceTensionForce ::globalEnergy() const
{
  Scalar energy = 0;
  std::vector<int> indices(12);
  std::vector<Vec3d> deformed(4);

  if(m_surface_tension_coeff == 0) return 0;

  EdgeIterator eit = m_shell.getDefoObj().edges_begin();
  for (;eit != m_shell.getDefoObj().edges_end(); ++eit) {
    const EdgeHandle& eh = *eit;
    
    FaceHandle f1, f2;
    bool valid = gatherDOFs(eh, deformed, f1, f2, indices);
    if(!valid) continue;
    
    if(f2.isValid())
      energy += elementEnergy(deformed, m_shell.getVolume(f1), m_shell.getVolume(f2));
    else
      energy += elementEnergy(deformed, m_shell.getVolume(f1), 0);

  }
  return energy;
}

void ShellLinearSurfaceTensionForce::globalForce( VecXd& force )  const
{

  std::vector<int> indices(12);
  std::vector<Vec3d> deformed(4);
  Eigen::Matrix<Scalar, 12, 1> localForce;

  if(m_surface_tension_coeff == 0) return;
  
  EdgeIterator eit = m_shell.getDefoObj().edges_begin();
  for (;eit != m_shell.getDefoObj().edges_end(); ++eit) {
    const EdgeHandle& eh = *eit;
   
    FaceHandle f1, f2;
    bool valid = gatherDOFs(eh, deformed, f1, f2, indices);
    if(!valid) continue;

    if(f2.isValid())
      elementForce(deformed, m_shell.getVolume(f1), m_shell.getVolume(f2), localForce);
    else 
      elementForce(deformed, m_shell.getVolume(f1), 0, localForce);

    for (int i = 0; i < (f2.isValid()?12:9); ++i)
      force(indices[i]) += localForce(i);
   
  }
  
}

void ShellLinearSurfaceTensionForce ::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  std::vector<int> indices(12);
  std::vector<Vec3d> deformed(4);
  Eigen::Matrix<Scalar, NumLSTDof, NumLSTDof> localMatrix;

  if(m_surface_tension_coeff == 0) return;

  EdgeIterator eit = m_shell.getDefoObj().edges_begin();
  for (;eit != m_shell.getDefoObj().edges_end(); ++eit) {
    const EdgeHandle& eh = *eit;

    FaceHandle f1, f2;
    bool valid = gatherDOFs(eh, deformed, f1, f2, indices);
    if(!valid) continue;

    if(f2.isValid())
      elementJacobian(deformed, m_shell.getVolume(f1), m_shell.getVolume(f2), localMatrix);
    else
      elementJacobian(deformed, m_shell.getVolume(f1), 0, localMatrix);

    for (int i = 0; i < (f2.isValid()?12:9); ++i)
      for(int j = 0; j < (f2.isValid()?12:9); ++j)
        Jacobian.add(indices[i], indices[j], scale * localMatrix(i,j));

  }
  
}


Scalar ShellLinearSurfaceTensionForce ::elementEnergy(const std::vector<Vec3d>& deformed, Scalar vol1, Scalar vol2) const
{
  
  std::vector<Scalar> deformed_data(NumLSTDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  adreal<NumLSTDof,0,Real> e = LSTEnergy<0>( *this, deformed_data, m_surface_tension_coeff, vol1, vol2 );
  Scalar energy = e.value();

  return energy;
}

void ShellLinearSurfaceTensionForce ::elementForce(const std::vector<Vec3d>& deformed, Scalar vol1, Scalar vol2,
                                    Eigen::Matrix<Scalar, NumLSTDof, 1>& force) const
{
  assert(deformed.size() == 4);

  std::vector<Scalar> deformed_data(NumLSTDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  //AutoDiff version
  adreal<NumLSTDof,0,Real> e = LSTEnergy<0>(*this, deformed_data, m_surface_tension_coeff, vol1, vol2);     
  for( int i = 0; i < (vol2 != 0? 12 : 9); i++ )
  {
    force[i] = -e.gradient(i);
  }
  
}

void ShellLinearSurfaceTensionForce ::elementJacobian(const std::vector<Vec3d>& deformed, Scalar vol1, Scalar vol2,
                                       Eigen::Matrix<Scalar,NumLSTDof,NumLSTDof>& jac) const
{
  assert(deformed.size() == 4);

  std::vector<Scalar> deformed_data(NumLSTDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] =deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  jac.setZero();

  adreal<NumLSTDof,1,Real> e = LSTEnergy<1>(*this, deformed_data, m_surface_tension_coeff, vol1, vol2);     
  // insert in the element jacobian matrix
  for( int i = 0; i < (vol2 != 0? 12 : 9); i++ )
  {
    for( int j = 0; j < (vol2 != 0? 12 : 9); j++ )
    {
      jac(i,j) = -e.hessian(i,j);
    }
  }

  


}



} //namespace BASim