/**
 * \file ShellVertexTriSpringForce.cc
 *
 * \author batty@cs.columbia.edu
 * \date Nov 22, 2011
 */

#include "BASim/src/Physics/DeformableObjects/Shells/ShellStickyRepulsionForce.hh"
#include "BASim/src/Core/EigenIncludes.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/MatrixBase.hh"

using Eigen::Matrix;

namespace BASim {

ShellStickyRepulsionForce::ShellStickyRepulsionForce( 
  ElasticShell& shell, 
  const std::string& name, 
  Scalar timestep)
: ElasticShellForce(shell, name), m_timestep(timestep)
{
  
}

bool ShellStickyRepulsionForce::gatherDOFs(const FaceHandle& fh, const VertexHandle& vh, std::vector<Vec3d>& deformed, std::vector<Vec3d>& undef_damp, std::vector<int>& indices) const {
  assert(deformed.size() == NumRepulsionVerts );
  
  assert(m_shell.isFaceActive(fh));
  assert(m_shell.isVertexActive(vh));

  //extract the relevant data for the face
  FaceVertexIterator fv_it = m_shell.getDefoObj().fv_iter(fh);
  int i = 0;
  for(;fv_it; ++fv_it) {
    const VertexHandle& facevert = *fv_it;
    deformed[i] = m_shell.getVertexPosition(facevert);
    undef_damp[i] = m_shell.getVertexDampingUndeformed(facevert);
    int dofBase = m_shell.getVertexDofBase(facevert);
    indices[i*3] = dofBase;
    indices[i*3+1] = dofBase+1;
    indices[i*3+2] = dofBase+2;
    ++i;
  }

  //extract the vertex data
  deformed[3] = m_shell.getVertexPosition(vh);
  undef_damp[3] = m_shell.getVertexDampingUndeformed(vh);
  int dofBase = m_shell.getVertexDofBase(vh);
  indices[9] = dofBase;
  indices[10] = dofBase+1;
  indices[11] = dofBase+2;
  
  return true;
}


template <int DO_HESS>
adreal<NumRepulsionDof,DO_HESS,Real> RepulsionEnergy(const ShellStickyRepulsionForce& vt, const std::vector<Scalar>& deformed, Vec3d baryCoords, Scalar strength, Scalar undef_len) {  

  // typedefs to simplify code below
  typedef adreal<NumRepulsionDof,DO_HESS,Real> adrealVT;
  typedef CVec3T<adrealVT> advecVT;

  Vector3d* s_deformed = (Vector3d*)(&deformed[0]);
  
  // indep variables
  advecVT   p[NumRepulsionVerts]; // vertex positions
  set_independent( p[0], s_deformed[0], 0 );
  set_independent( p[1], s_deformed[1], 3 );
  set_independent( p[2], s_deformed[2], 6 );    
  set_independent( p[3], s_deformed[3], 9 );    
    
  adrealVT e(0);
  
  Vector3d normalReal = cross(s_deformed[2] - s_deformed[0], s_deformed[1] - s_deformed[0]); 
  advecVT normalVbl = cross(p[2] - p[0], p[1] - p[0]);
  
  normalize(normalVbl);
  normalize(normalReal);

  advecVT offset = p[3] - (baryCoords[0]*p[0] + baryCoords[1]*p[1] + baryCoords[2]*p[2]);

  if(dot(s_deformed[3] - s_deformed[0], normalReal) > 0) {
    e += 0.5*strength*sqr( dot(offset, normalReal) - undef_len);  
    //e += 0.5*strength*sqr( dot(offset, normalVbl) - undef_len);  
  }
  else {
    e += 0.5*strength*sqr( -dot(offset, normalReal) - undef_len);  
    //e += 0.5*strength*sqr( -dot(offset, normalVbl) - undef_len);  
  }
 
  //tangential energy
  
  //note, this is an approx for speed.
  e += 0.5*strength*lenSq(offset - dot(offset,normalReal)*normalReal); 
  //the true version
  //e += 0.5*strength*lenSq(offset - dot(offset,normalVbl)*normalVbl); 

  return e;
}

Scalar ShellStickyRepulsionForce::globalEnergy() const
{
  Scalar energy = 0;
  std::vector<int> indices(NumRepulsionDof );
  std::vector<Vec3d> deformed(NumRepulsionVerts );
  std::vector<Vec3d> undef_damp(NumRepulsionVerts );

  for (unsigned int i = 0; i < m_faces.size(); ++i) {
    const FaceHandle& fh = m_faces[i];
    const VertexHandle& vh = m_vertices[i];

    gatherDOFs(fh, vh, deformed, undef_damp, indices);
    energy += elementEnergy(deformed, m_barycoords[i], m_stiffnesses[i], m_restlen[i]);
  }
  return energy;
}

void ShellStickyRepulsionForce::globalForce( VecXd& force )  const
{

  std::vector<int> indices(NumRepulsionDof );
  std::vector<Vec3d> deformed(NumRepulsionVerts );
  std::vector<Vec3d> undeformed_damp(NumRepulsionVerts );
  Eigen::Matrix<Scalar, NumRepulsionDof , 1> localForce;

  for (unsigned int f = 0; f < m_faces.size(); ++f) {
    const FaceHandle& fh = m_faces[f];
    const VertexHandle& vh = m_vertices[f];
   
    bool valid = gatherDOFs(fh, vh, deformed, undeformed_damp, indices);
    if(!valid) continue;

    //elastic force
    elementForce(deformed, m_barycoords[f], m_stiffnesses[f], m_restlen[f], localForce);
    for (unsigned int i = 0; i < indices.size(); ++i) {
      force(indices[i]) += localForce(i);
    }

    
    //viscous force
    localForce.setZero();
    Vec3d facePt = undeformed_damp[0]*m_barycoords[f][0] + 
      undeformed_damp[1]*m_barycoords[f][1] +   
      undeformed_damp[2]*m_barycoords[f][2];
    Vec3d normal = (undeformed_damp[2] - undeformed_damp[0]).cross(undeformed_damp[1] - undeformed_damp[0]);
    normal.normalize();
    Scalar damp_restlen = fabs((facePt - undeformed_damp[3]).dot(normal));
    elementForce(deformed, m_barycoords[f], m_damping[f], damp_restlen, localForce);
    for (unsigned int i = 0; i < indices.size(); ++i) {
      force(indices[i]) += (1.0 / m_timestep) * localForce(i);
    }

  }
  
}

void ShellStickyRepulsionForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  std::vector<int> indices(NumRepulsionDof );
  std::vector<Vec3d> deformed(NumRepulsionVerts );
  std::vector<Vec3d> undeformed_damp(NumRepulsionVerts );

  Eigen::Matrix<Scalar, NumRepulsionDof , NumRepulsionDof > localMatrix;

  for (unsigned int f = 0; f < m_faces.size(); ++f) {
    const FaceHandle& fh = m_faces[f];
    const VertexHandle& vh = m_vertices[f];

    bool valid = gatherDOFs(fh, vh, deformed, undeformed_damp, indices);
    if(!valid) continue;
    
    elementJacobian(deformed, m_barycoords[f], m_stiffnesses[f], m_restlen[f], localMatrix);
    for (unsigned int i = 0; i < indices.size(); ++i)
      for(unsigned int j = 0; j < indices.size(); ++j)
        Jacobian.add(indices[i], indices[j], scale * localMatrix(i,j));

    
    Vec3d facePt = undeformed_damp[0]*m_barycoords[f][0] + 
      undeformed_damp[1]*m_barycoords[f][1] +   
      undeformed_damp[2]*m_barycoords[f][2];
    Vec3d normal = (undeformed_damp[2] - undeformed_damp[0]).cross(undeformed_damp[1] - undeformed_damp[0]);
    normal.normalize();
    Scalar damp_restlen = fabs((facePt - undeformed_damp[3]).dot(normal));
    elementJacobian(deformed, m_barycoords[f], m_damping[f], damp_restlen, localMatrix);
    for (unsigned int i = 0; i < indices.size(); ++i)
      for(unsigned int j = 0; j < indices.size(); ++j)
        Jacobian.add(indices[i], indices[j], scale / m_timestep * localMatrix(i,j));

  }
  
}


Scalar ShellStickyRepulsionForce::elementEnergy(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords, Scalar strength, Scalar restlen) const
{
  
  std::vector<Scalar> deformed_data(NumRepulsionDof );
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  adreal<NumRepulsionDof ,0,Real> e = RepulsionEnergy<0>( *this, deformed_data, baryCoords, strength, restlen);
  Scalar energy = e.value();

  return energy;
}

void ShellStickyRepulsionForce::elementForce(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords, Scalar strength, Scalar restlen,
                                    Eigen::Matrix<Scalar, NumRepulsionDof, 1>& force) const
{
  assert(deformed.size() == NumRepulsionVerts );
 
  std::vector<Scalar> deformed_data(NumRepulsionDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  bool useAutodiff = true;
  if(useAutodiff) {
    //AutoDiff version
    adreal<NumRepulsionDof,0,Real> e = RepulsionEnergy<0>(*this, deformed_data, baryCoords, strength, restlen);     
    for( uint i = 0; i < NumRepulsionDof; i++ )
    {
      force[i] = -e.gradient(i);
    }
  }
  else {

    //Explicit version. p0,p1,p2 is the triangle. p3 is the other vertex.
    Eigen::Matrix<Scalar, NumRepulsionDof, 1> force2;

    //Normal to the triangle
    Vec3d normal = (deformed[2] - deformed[0]).cross(deformed[1] - deformed[0]);
    normal.normalize();

    //Offset vector from point on the triangle to the vertex
    Vec3d offset = deformed[3] - (baryCoords[0]*deformed[0] + baryCoords[1]*deformed[1] + baryCoords[2]*deformed[2]);

    //Normal component spring, with non-zero rest length
    Vec3d forceFactor;
    if(offset.dot(normal) > 0)
      forceFactor = strength * (offset.dot(normal) - restlen) * normal;
    else
      forceFactor = -strength * (-offset.dot(normal) - restlen) * normal;
  
    //Tangential zero-rest-length spring
    Vec3d normalComponent = offset.dot(normal)*normal;
    forceFactor += strength * (offset-normalComponent); 

    Vec3d forceV0 = baryCoords[0]*forceFactor;
    Vec3d forceV1 = baryCoords[1]*forceFactor;
    Vec3d forceV2 = baryCoords[2]*forceFactor;
    Vec3d forceV3 = -forceFactor;

    force.block<3,1>(0,0) = forceV0;
    force.block<3,1>(3,0) = forceV1;
    force.block<3,1>(6,0) = forceV2;
    force.block<3,1>(9,0) = forceV3;
  }

 
  /* std::cout << "XXXXX\n";

  std::cout << "Correct force: ";
  for(int i = 0; i < 12; ++i)
  std::cout << force[i] << " " << std::endl;

  std::cout << "Newest force!: ";
  for(int i = 0; i < 12; ++i)
  std::cout << force2[i] << " " << std::endl;

  std::cout << "XXXXX\n";*/
}

void ShellStickyRepulsionForce::elementJacobian(const std::vector<Vec3d>& deformed, const Vec3d& baryCoords, Scalar strength, Scalar restlen,
                                       Eigen::Matrix<Scalar,NumRepulsionDof,NumRepulsionDof>& jac) const
{
  assert(deformed.size() == NumRepulsionVerts );

  Eigen::Matrix<Scalar,NumRepulsionDof,NumRepulsionDof> jac2;
  std::vector<Scalar> deformed_data(NumRepulsionDof);
  for(unsigned int i = 0; i < deformed.size(); ++i) {
    deformed_data[3*i] = deformed[i][0];
    deformed_data[3*i+1] = deformed[i][1];
    deformed_data[3*i+2] = deformed[i][2];
  }

  jac2.setZero();

  adreal<NumRepulsionDof,1,Real> e = RepulsionEnergy<1>(*this, deformed_data, baryCoords, strength, restlen);     
  // insert in the element jacobian matrix
  for( uint i = 0; i < NumRepulsionDof; i++ )
  {
    for( uint j = 0; j < NumRepulsionDof; j++ )
    {
      jac2(i,j) = -e.hessian(i,j);
    }
  }

}

void ShellStickyRepulsionForce::addSpring(const FaceHandle& fh, const VertexHandle& vh, const Vec3d& baryCoords,
                                          Scalar stiffness, Scalar damping, Scalar restlen) {
  m_faces.push_back(fh);
  m_vertices.push_back(vh);
  m_barycoords.push_back(baryCoords);
  m_stiffnesses.push_back(stiffness);
  m_damping.push_back(damping);
  m_restlen.push_back(restlen);

  m_springset.insert(std::make_pair(fh.idx(),vh.idx()));
  
  /*
  Scalar energy = 0;
  std::vector<int> indices(NumRepulsionDof );
  std::vector<Vec3d> deformed(NumRepulsionVerts );
  std::vector<Vec3d> undef_damp(NumRepulsionVerts );

  gatherDOFs(fh, vh, deformed, undef_damp, indices);
  energy += elementEnergy(deformed, baryCoords, stiffness, restlen);
  
  std::cout << "Energy: " << energy << std::endl;
  
  Eigen::Matrix<Scalar, NumRepulsionDof , 1> localForce;
  
  elementForce(deformed, baryCoords, stiffness, restlen, localForce);
  for(int i = 0; i < 12; ++i)
    std::cout << localForce(i) << ", ";
  std::cout << std::endl;
  */
}

void ShellStickyRepulsionForce::clearSprings() {
  m_faces.clear();
  m_vertices.clear();
  m_barycoords.clear();
  m_stiffnesses.clear();
  m_damping.clear();
  m_restlen.clear();
  
  m_springset.clear();
}

void ShellStickyRepulsionForce::clearSprings(VertexHandle& v) {
  bool deletedElt = true;
  while(deletedElt) {
    int target = -1;
    deletedElt = false;
    for(unsigned int i = 0; i < m_vertices.size(); ++i) {
      if(m_vertices[i] == v) {
        target = i;
        deletedElt = true;
        break;
      }
    }

    if(target != -1) {
      m_faces.erase(m_faces.begin() + target);
      m_vertices.erase(m_vertices.begin() + target);
      m_barycoords.erase(m_barycoords.begin() + target);
      m_stiffnesses.erase(m_stiffnesses.begin() + target);
      m_damping.erase(m_damping.begin() + target);
      m_restlen.erase(m_restlen.begin() + target);
    }
  }

  //adjust the springset too
  for (std::set< std::pair<int,int> >::iterator it = m_springset.begin(); it != m_springset.end(); ) {
    std::pair<int,int> springPair = *it;
    if (springPair.second == v.idx()) {
      m_springset.erase(it++);
    }
    else {
      ++it;
    }
  }

}

void ShellStickyRepulsionForce::clearSprings(FaceHandle& f) {
  bool deletedElt = true;
  while(deletedElt) {
    int target = -1;
    deletedElt = false;
    for(unsigned int i = 0; i < m_faces.size(); ++i) {
      if(m_faces[i] == f) {
        target = i;
        deletedElt = true;
        break;
      }
    }

    if(target != -1) {
      m_faces.erase(m_faces.begin() + target);
      m_vertices.erase(m_vertices.begin() + target);
      m_barycoords.erase(m_barycoords.begin() + target);
      m_stiffnesses.erase(m_stiffnesses.begin() + target);
      m_damping.erase(m_damping.begin() + target);
      m_restlen.erase(m_restlen.begin() + target);
    }
  }

  //adjust the springset too
  for (std::set< std::pair<int,int> >::iterator it = m_springset.begin(); it != m_springset.end(); ) {
    std::pair<int,int> springPair = *it;
    if (springPair.first == f.idx()) {
      m_springset.erase(it++);
    }
    else {
      ++it;
    }
  }

}

void ShellStickyRepulsionForce::getSpringLists(std::vector<VertexHandle> &verts, std::vector<FaceHandle>& tris, std::vector<Vec3d>& barycoords) {
  verts = m_vertices;
  tris = m_faces;
  barycoords = m_barycoords;
}

bool ShellStickyRepulsionForce::springExists(const FaceHandle& f, const VertexHandle& v) {
  return m_springset.find(std::make_pair(f.idx(), v.idx())) != m_springset.end();
}

bool ShellStickyRepulsionForce::isVertexInUse(const VertexHandle& vh) {
  return std::find(m_vertices.begin(), m_vertices.end(), vh) != m_vertices.end();
}

bool ShellStickyRepulsionForce::isFaceInUse(const FaceHandle& fh) {
  return std::find(m_faces.begin(), m_faces.end(), fh) != m_faces.end();
}



} //namespace BASim