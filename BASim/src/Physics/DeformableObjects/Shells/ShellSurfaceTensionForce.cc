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
  Scalar surf_coeff,
  int active_scene)
: ElasticShellForce(shell, name), m_surface_tension_coeff(surf_coeff), m_active_scene(active_scene)
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
    int dofBase = m_shell.getDefoObj().getPositionDofBase(vh);
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
  typedef adreal<NumSTDof,DO_HESS,Real> adrealST;
  typedef CVec3T<adrealST> advecST;

  Vector3d* s_deformed = (Vector3d*)(&deformed[0]);
  
  // indep variables
  advecST   p[3]; // vertex positions
  set_independent( p[0], s_deformed[0], 0 );
  set_independent( p[1], s_deformed[1], 3 );
  set_independent( p[2], s_deformed[2], 6 );    
    
  //energy = 2*coeff*surface_area
  adrealST e(0);
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
    
    //determine the energy for this element
    energy += elementEnergy(deformed);

  }
  return energy;
}

void ShellSurfaceTensionForce::globalForce( VecXd& force )  const
{
//  std::cout << "SF energy = " << globalEnergy() << std::endl;

  std::vector<int> indices(9);
  std::vector<Vec3d> deformed(3);
  Eigen::Matrix<Scalar, 9, 1> localForce;

  if(m_surface_tension_coeff == 0) return;
  
  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;
   
    bool valid = gatherDOFs(fh, deformed, indices);
    if(!valid) continue;

    elementForce(deformed, localForce);
    for (unsigned int i = 0; i < indices.size(); ++i)
      force(indices[i]) += localForce(i);
   
  }
  
  if (m_active_scene == 8)
  {
    // special for the reuleaux tet test: constrain in-plane motion of wall vertices    
    for (FaceIterator fit = m_shell.getDefoObj().faces_begin();fit != m_shell.getDefoObj().faces_end(); ++fit)
    {
      const FaceHandle& fh = *fit;
     
      bool valid = gatherDOFs(fh, deformed, indices);
      if(!valid) continue;

      for (int i = 0; i < 3; i++)
      {
        if (deformed[i].x() <= 0 || deformed[i].x() >= 1 || deformed[i].y() <= 0 || deformed[i].y() >= 1 || deformed[i].z() <= 0 || deformed[i].z() >= 1)
        {
          force[indices[i * 3 + 0]] = 0;
          force[indices[i * 3 + 1]] = 0;
          force[indices[i * 3 + 2]] = 0;
        }
      }
    }
  }
}

void ShellSurfaceTensionForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  std::vector<int> indices(9);
  std::vector<Vec3d> deformed(3);
  Eigen::Matrix<Scalar, 9, 9> localMatrix;

  if(m_surface_tension_coeff == 0) return;

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
 /* 
  std::vector<Scalar> deformed_data(NumSTDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  //AutoDiff version
  adreal<NumSTDof,0,Real> e = STEnergy<0>(*this, deformed_data, m_surface_tension_coeff);     
  for( uint i = 0; i < NumSTDof; i++ )
  {
    force[i] = -e.gradient(i);
  }
  */
  
  //Derived by directly taking the derivative of len(cross(v1,v2))
  Vec3d v01 = deformed[1] - deformed[0];
  Vec3d v20 = deformed[0] - deformed[2];
  Vec3d v12 = deformed[2] - deformed[1];
  Scalar v01norm = v01.norm();
  Scalar v20norm = v20.norm();
  Scalar v12norm = v12.norm();
  
  Vec3d A = v01.cross(-v20);
  Scalar Anorm = A.norm();
  Vec3d mul = A / Anorm;
  Vec3d p2part = m_surface_tension_coeff*v01.cross(mul);
  Vec3d p1part = m_surface_tension_coeff*mul.cross(-v20);
  Vec3d p0part = -(p1part + p2part);

//  // handle degenerate cases
//  Scalar mean_edge_length = (v01norm + v12norm + v20norm) / 3;
//  if (Anorm < 0.01 * mean_edge_length * mean_edge_length)
//  {
//    Scalar r = Anorm / (0.01 * mean_edge_length * mean_edge_length);
//    Scalar t = std::max(0.0, r * 2.0 - 1.0); // t = 0 to 1 for Anorm = 0.5% to 1% of mean edge length squared
//    
//    Vec3d & max_edge = ((v01norm > v12norm && v01norm > v20norm) ? v01 : (v12norm > v20norm ? v12 : v20));
//    Vec3d max_edge_perp;
//    if (fabs(max_edge.x()) < fabs(max_edge.y()) && fabs(max_edge.x()) < fabs(max_edge.z()))
//      max_edge_perp = Vec3d(0.0, -max_edge.z(), max_edge.y()).normalized() * max_edge.norm();
//    else if (fabs(max_edge.y()) < fabs(max_edge.x()) && fabs(max_edge.y()) < fabs(max_edge.z()))
//      max_edge_perp = Vec3d(-max_edge.z(), 0.0, max_edge.x()).normalized() * max_edge.norm();
//    else
//      max_edge_perp = Vec3d(-max_edge.y(), max_edge.x(), 0.0).normalized() * max_edge.norm();
//    
//    Vec3d p2part_limit;
//    Vec3d p1part_limit;
//    Vec3d p0part_limit;
//    
//    if (v01norm > v12norm && v01norm > v20norm)
//    {
//      // v01 is the longest edge
//      p2part_limit = max_edge_perp;
//      p1part_limit = -max_edge_perp * v20norm / (v12norm + v20norm);
//      p0part_limit = -max_edge_perp * v12norm / (v12norm + v20norm);
//    } else if (v12norm > v20norm)
//    {
//      // v12 is the longest edge
//      p2part_limit = -max_edge_perp * v01norm / (v01norm + v20norm);
//      p1part_limit = -max_edge_perp * v20norm / (v01norm + v20norm);
//      p0part_limit = max_edge_perp;
//    } else
//    {
//      // v20 is the longest edge
//      p2part_limit = -max_edge_perp * v01norm / (v01norm + v12norm);
//      p1part_limit = max_edge_perp;
//      p0part_limit = -max_edge_perp * v12norm / (v01norm + v12norm);
//    }
//    
//    p2part = p2part * t + p2part_limit * (1 - t);
//    p1part = p1part * t + p1part_limit * (1 - t);
//    p0part = p0part * t + p0part_limit * (1 - t);
//    
//  }
  
  /*  
  //Alternative implementation
  //using the fact that gradient of triangle area w.r.t p_i is
  //half the vector (p_j-p_k) of the opposing edge, rotated 90 degrees 
  //to point towards p_i, in the triangle plane.
  Vec3d v0 = deformed[1] - deformed[0];
  Vec3d v1 = deformed[2] - deformed[1];
  Vec3d v2 = deformed[0] - deformed[2];
  Vec3d A = v1.cross(v2);
  Vec3d mul = A / A.norm();
  p0part = -m_surface_tension_coeff*mul.cross(v1);
  p1part = -m_surface_tension_coeff*mul.cross(v2);
  p2part = -m_surface_tension_coeff*mul.cross(v0);
  */

  //Assign result back to the force vector
  force[0] = p0part[0];
  force[1] = p0part[1];
  force[2] = p0part[2];
  force[3] = p1part[0];
  force[4] = p1part[1];
  force[5] = p1part[2];
  force[6] = p2part[0];
  force[7] = p2part[1];
  force[8] = p2part[2];
  
}

void ShellSurfaceTensionForce::elementJacobian(const std::vector<Vec3d>& deformed,
                                       Eigen::Matrix<Scalar,9,9>& jac) const
{
  assert(deformed.size() == 3);

  //Eigen::Matrix<Scalar,9,9> jac2;
  //std::vector<Scalar> deformed_data(NumSTDof);
  //for(unsigned int i = 0; i < deformed.size(); ++i) {
  //  deformed_data[3*i] = deformed[i][0];
  //  deformed_data[3*i+1] = deformed[i][1];
  //  deformed_data[3*i+2] = deformed[i][2];
  //}

  //jac2.setZero();

  //adreal<NumSTDof,1,Real> e = STEnergy<1>(*this, deformed_data, m_surface_tension_coeff);     
  //// insert in the element jacobian matrix
  //for( uint i = 0; i < NumSTDof; i++ )
  //{
  //  for( uint j = 0; j < NumSTDof; j++ )
  //  {
  //    jac2(i,j) = -e.hessian(i,j);
  //  }
  //}

  //now construct from the explicitly worked out derivatives
  //(See the PDF by Christopher Batty)
  Vec3d a = deformed[0];
  Vec3d b = deformed[1];
  Vec3d c = deformed[2];

  Vec3d ab = b-a;
  Vec3d ac = c-a;

  //Scalar K = (ab.dot(ab))*(ac.dot(ac)) - sqr(ab.dot(ac));
  Scalar K = (ab.dot(ab))*(ac.dot(ac)) - (ab.dot(ac))*(ab.dot(ac));
  
  Vec3d dkdc = 2*(ab.dot(ab)*ac - ab.dot(ac)*ab);
  Vec3d dkdb = 2*(ac.dot(ac)*ab - ac.dot(ab)*ac);
  Vec3d dkda = -(dkdc+dkdb);

  Eigen::Matrix<Scalar,3,3> d2kdbdc;
  Eigen::Matrix<Scalar,3,3> d2kdcdb;
  Eigen::Matrix<Scalar,3,3> d2kdbdb;
  Eigen::Matrix<Scalar,3,3> d2kdbda;
  Eigen::Matrix<Scalar,3,3> d2kdcdc;
  Eigen::Matrix<Scalar,3,3> d2kdcda;
  Eigen::Matrix<Scalar,3,3> d2kdada;

  for(int i = 0; i < 3; ++i) {
    for(int k = 0; k < 3; ++k) {
      d2kdbdc(i,k) = 2 *( 
        2*(ab(i))*(ac(k))
        -(ac(i))*(ab(k))
        -(i==k?ab.dot(ac):0));
      d2kdbdb(i,k) = 2 *( 
        (i==k?ac.dot(ac):0)
        -(ac(i))*(ac(k))
        );
      d2kdcdc(i,k) = 2 *( 
        (i==k?ab.dot(ab):0)
        -(ab(i))*(ab(k))
        );
    }
  }
  d2kdcdb = d2kdbdc.transpose();
  d2kdbda = -(d2kdbdb + d2kdbdc);
  d2kdcda = -(d2kdcdc + d2kdcdb);
  d2kdada = -(d2kdbda + d2kdcda);
  Scalar factor1 = -1.0 / (8.0*pow(K, 1.5));
  Scalar factor2 = +1.0 / (4.0*sqrt(K));

  Eigen::Matrix<Scalar,3,3> dAdbdc = factor1 * dkdb*dkdc.transpose() + factor2 * d2kdbdc;
  Eigen::Matrix<Scalar,3,3> dAdbdb = factor1 * dkdb*dkdb.transpose() + factor2 * d2kdbdb; 
  Eigen::Matrix<Scalar,3,3> dAdbda = factor1 * dkdb*dkda.transpose() + factor2 * d2kdbda;
  Eigen::Matrix<Scalar,3,3> dAdcdc = factor1 * dkdc*dkdc.transpose() + factor2 * d2kdcdc; 
  Eigen::Matrix<Scalar,3,3> dAdcda = factor1 * dkdc*dkda.transpose() + factor2 * d2kdcda;
  Eigen::Matrix<Scalar,3,3> dAdada = factor1 * dkda*dkda.transpose() + factor2 * d2kdada;

  //assemble into the main matrix
  jac.block(0,0, 3,3) = dAdada;            
  jac.block(3,0, 3,3) = dAdbda;            
  jac.block(6,0, 3,3) = dAdcda;            

  jac.block(0,3, 3,3) = dAdbda.transpose();
  jac.block(3,3, 3,3) = dAdbdb;            
  jac.block(6,3, 3,3) = dAdbdc.transpose();

  jac.block(0,6, 3,3) = dAdcda.transpose();
  jac.block(3,6, 3,3) = dAdbdc;            
  jac.block(6,6, 3,3) = dAdcdc;            
  
  jac *= -m_surface_tension_coeff * 2;

  /*for(int i = 0; i < 9; ++i)
    for(int j = 0; j < 9; ++j)
      if(fabs(jac(i,j) - jac2(i,j)) > 1e-10)
        printf("mismatch\n");*/
}



} //namespace BASim