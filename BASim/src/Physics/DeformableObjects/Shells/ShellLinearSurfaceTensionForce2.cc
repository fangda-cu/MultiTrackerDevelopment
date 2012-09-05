/**
 * \file ShellLinearSurfaceTensionForce2.cc
 *
 * \author batty@cs.columbia.edu
 * \date Aug 21, 2012
 */

#include "BASim/src/Physics/DeformableObjects/Shells/ShellLinearSurfaceTensionForce2.hh"
#include "BASim/src/Core/EigenIncludes.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"

using Eigen::Matrix;

namespace BASim {

ShellLinearSurfaceTensionForce2 ::ShellLinearSurfaceTensionForce2 ( 
  ElasticShell& shell, 
  const std::string& name, 
  Scalar surf_coeff )
: ElasticShellForce(shell, name), m_surface_tension_coeff(surf_coeff)
{
  
}

bool ShellLinearSurfaceTensionForce2 ::gatherDOFs(const FaceHandle& fh, 
        std::vector<Vec3d>& vertices, 
        std::vector<Vec3i>& faces,
        std::vector<Scalar>& volumes,
        Vec3i& incidentFacesPerVertex,
        Vec3i& isBdryVertex,
        std::vector<int>& indices) const {
  
  if(!m_shell.isFaceActive(fh)) return false;

  indices.clear();
  vertices.clear();
  faces.clear();
  volumes.clear();
  incidentFacesPerVertex = Vec3i(0,0,0);

  std::map<int,int> vertex_map; //a unique set of the vertices, and their "local" index

  //extract the relevant data for the central face
  FaceVertexIterator fv_it = m_shell.getDefoObj().fv_iter(fh);
  int i = 0;
  Vec3i centreFace;
  for(;fv_it; ++fv_it) {
    const VertexHandle& vh = *fv_it;
    vertices.push_back(m_shell.getVertexPosition(vh));
    int dofBase = m_shell.getDefoObj().getPositionDofBase(vh);
    indices.push_back(dofBase);
    indices.push_back(dofBase+1);
    indices.push_back(dofBase+2);
    vertex_map[vh.idx()] = i;
    centreFace[i] = i;
    ++i;
  }
  faces.push_back(centreFace);
  volumes.push_back(m_shell.getVolume(fh));
  
  
  //now get the vertices for all the other faces
  FaceVertexIterator fvit = m_shell.getDefoObj().fv_iter(fh);
  int vertNum = 0;
  for(;fvit; ++fvit) { //look at the three vertices of the central face
    VertexHandle vh = *fvit;
    VertexFaceIterator vfit = m_shell.getDefoObj().vf_iter(vh);
    incidentFacesPerVertex[vertNum] = 0;

    //determine if this vertex is on a boundary
    isBdryVertex[vertNum] = 0;
    for(VertexEdgeIterator veit = m_shell.getDefoObj().ve_iter(vh); veit; ++veit) {
      EdgeHandle eh = *veit;
      if(m_shell.getDefoObj().edgeIncidentFaces(eh) == 1) {
        isBdryVertex[vertNum] = 1;
      }
    }

    for(; vfit; ++vfit) { //look at the faces incident on the particular vertex
      ++incidentFacesPerVertex[vertNum];
      FaceHandle fh2 = *vfit;
      volumes.push_back(m_shell.getVolume(fh2));
      FaceVertexIterator fvit2 = m_shell.getDefoObj().fv_iter(fh2);
      Vec3i curFace;
      int inner_ind = 0;
      for(; fvit2; ++fvit2) { //finally, iterate over 3 vertices of that face
        
        const VertexHandle& vh_inner = *fvit2;
        if(vertex_map.find(vh_inner.idx()) != vertex_map.end()) {
          curFace[inner_ind] = vertex_map[vh_inner.idx()];
        }
        else {
          vertices.push_back(m_shell.getVertexPosition(vh_inner));
          int dofBase = m_shell.getDefoObj().getPositionDofBase(vh_inner);
          indices.push_back(dofBase);
          indices.push_back(dofBase+1);
          indices.push_back(dofBase+2);
          curFace[inner_ind] = vertices.size()-1;
          vertex_map[vh_inner.idx()] = vertices.size()-1;
        }
        ++inner_ind;
      }
      faces.push_back(curFace);
     
    }
    ++vertNum;  
  }
  
  assert(vertices.size() == indices.size()/3);

  return true;
}


template <int DO_HESS>
adreal<NumLST2Dof,DO_HESS,Real> LST2Energy(const ShellLinearSurfaceTensionForce2& mn, 
  const std::vector<Scalar>& vertices_flat, 
  const std::vector<Vec3i>& faces,
  const std::vector<Scalar>& volumes,
  const Vec3i& incidentFacesPerMainVertex,
  const Vec3i& isBoundaryVertex,
  Scalar surf_coeff) 
{  

  // typedefs to simplify code below
  typedef adreal<NumLST2Dof,DO_HESS,Real> adrealST;
  typedef CVec3T<adrealST> advecST;

  Vector3d* s_vertices = (Vector3d*)(&vertices_flat[0]);

  assert(vertices_flat.size()/3 < NumLST2VertsMax);

  // independent variables
  advecST   p[NumLST2VertsMax]; // vertex positions 
  for(unsigned int i = 0; i < vertices_flat.size()/3; ++i) {
    set_independent( p[i], s_vertices[i], i*3 );
  }

  //compute face areas from vertex positions
  adrealST faceAreas[3*NumLST2VertsMax];
  for(unsigned int i = 0; i < faces.size(); ++i) {
    Vec3i faceData = faces[i];
    faceAreas[i] = 0.5 * len(cross(p[faceData[1]] - p[faceData[0]], p[faceData[2]] - p[faceData[0]]));
  }

  //work out averaged vertex thicknesses and normals
  adrealST thicknesses[3];
  advecST normals[3];
  int startInd = 1; //first face is the central one, so skip it here
  for(int vertexNumber = 0; vertexNumber < 3; ++vertexNumber) {
    adrealST volAccum(0);
    adrealST areaAccum(0);
    
    if(isBoundaryVertex[vertexNumber] > 0) {
      for(int faceNumber = 0; faceNumber < incidentFacesPerMainVertex[vertexNumber]; ++faceNumber) {
        int curFace = startInd + faceNumber;
        areaAccum += faceAreas[curFace];
        volAccum += volumes[curFace];
      }
      thicknesses[vertexNumber] = volAccum / areaAccum;
    }
    else {
      thicknesses[vertexNumber] = 0;
    }
    startInd += incidentFacesPerMainVertex[vertexNumber];
  }

  //Compute main triangle's face normal
  advecST normal = cross(p[1] - p[0], p[2] - p[0]);
  normalize(normal);

  //Compute effective offset positions for surface vertices based on thicknesses computed above
  advecST vert[3];
  for(int i = 0; i < 3; ++i)
    vert[i] = p[i] + 0.5*thicknesses[i]*normal;//s[i];

  advecST centre = (p[0]+p[1]+p[2])/3.0 + 0.5*(volumes[0] / faceAreas[0]) * normal;

  advecST vertb[3];
  for(int i = 0; i < 3; ++i)
    vertb[i] = p[i] - 0.5*thicknesses[i]*normal;//s[i];

  advecST centreb = (p[0]+p[1]+p[2])/3.0 - 0.5*(volumes[0] / faceAreas[0]) * normal;

  adrealST e(0);

  //normal linear triangle terms
  //need to do both sides in order to ensure cancellation of spurious out of plane forces
  e += 0.5*surf_coeff * len(cross(vert[1] - vert[0], vert[2] - vert[0])); 
  e += 0.5*surf_coeff * len(cross(vertb[1] - vertb[0], vertb[2] - vertb[0]));

  //Assume linear triangles between the vertices and the triangle central peak
  //This doesn't seem to work; instead it partly acts within each triangle to make it
  //near equilateral...
  /*e += 0.5*surf_coeff * len(cross(vert[1] - centre, vert[2] - centre)); 
  e += 0.5*surf_coeff * len(cross(vert[0] - centre, vert[1] - centre)); 
  e += 0.5*surf_coeff * len(cross(vert[2] - centre, vert[0] - centre)); 

  e += 0.5*surf_coeff * len(cross(vertb[1] - centreb, vertb[2] - centreb)); 
  e += 0.5*surf_coeff * len(cross(vertb[0] - centreb, vertb[1] - centreb)); 
  e += 0.5*surf_coeff * len(cross(vertb[2] - centreb, vertb[0] - centreb)); */

  
  //boundary edge terms, if we want to use them
  /*if(isBoundaryVertex[0] && isBoundaryVertex[1]) {
  e += surf_coeff * 0.5*(thicknesses[0]+thicknesses[1]) * len(p[0] - p[1]);
  }
  if(isBoundaryVertex[2] && isBoundaryVertex[1]) {
  e += surf_coeff * 0.5*(thicknesses[2]+thicknesses[1]) * len(p[2] - p[1]);
  }
  if(isBoundaryVertex[2] && isBoundaryVertex[0]) {
  e += surf_coeff * 0.5*(thicknesses[0]+thicknesses[2]) * len(p[2] - p[0]);
  }*/
  

  return e;
}
Scalar ShellLinearSurfaceTensionForce2 ::globalEnergy() const
{
  Scalar energy = 0;
  std::vector<int> indices;
  std::vector<Vec3d> vertices;
  std::vector<Scalar> volumes;
  std::vector<Vec3i> faces;
  Vec3i incidentFacesPerVertex;
  Vec3i isBdryVertex;

  if(m_surface_tension_coeff == 0) return 0;

  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {

    const FaceHandle& fh = *fit;
    
    bool valid = gatherDOFs(fh, vertices, faces, volumes, incidentFacesPerVertex, isBdryVertex, indices);
    if(!valid) continue;
    
    energy += elementEnergy(vertices, faces, volumes, incidentFacesPerVertex, isBdryVertex);

  }
  return energy;
}

void ShellLinearSurfaceTensionForce2::globalForce( VecXd& force )  const
{


  std::cout << "Computing global forces\n";
  std::vector<int> indices;
  std::vector<Vec3d> vertices;
  std::vector<Scalar> volumes;
  std::vector<Vec3i> faces;
  Vec3i incidentFacesPerVertex, isBdryVertex;

  Eigen::Matrix<Scalar, NumLST2Dof, 1> localForce;

  if(m_surface_tension_coeff == 0) return;
  
  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {

    const FaceHandle& fh = *fit;

    bool valid = gatherDOFs(fh, vertices, faces, volumes, incidentFacesPerVertex, isBdryVertex, indices);
    if(!valid) continue;

    elementForce(vertices, faces, volumes, incidentFacesPerVertex, isBdryVertex, localForce);
    for (unsigned int i = 0; i < 3*vertices.size(); ++i)
      force(indices[i]) += localForce(i);
   
  }
  std::cout << "Done global forces\n";
  
}

void ShellLinearSurfaceTensionForce2 ::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  return;

  std::cout << "Computing global Jacobian\n";
  std::vector<int> indices;
  std::vector<Vec3d> vertices;
  std::vector<Scalar> volumes;
  std::vector<Vec3i> faces;
  Vec3i incidentFacesPerVertex, isBdryVertex;

  Eigen::Matrix<Scalar, NumLST2Dof, NumLST2Dof> localMatrix;

  if(m_surface_tension_coeff == 0) return;

  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {

    const FaceHandle& fh = *fit;

    bool valid = gatherDOFs(fh, vertices, faces, volumes, incidentFacesPerVertex, isBdryVertex, indices);
    if(!valid) continue;
    localMatrix.setZero();

    std::cout << "Calling for local Jacobian\n";
    elementJacobian(vertices, faces, volumes, incidentFacesPerVertex, isBdryVertex, localMatrix);
    
    std::cout << "Scattering Jacobian\n";
    for (unsigned int i = 0; i < 3*vertices.size(); ++i)
      for(unsigned int j = 0; j < 3*vertices.size(); ++j)
        Jacobian.add(indices[i], indices[j], scale * localMatrix(i,j));

  }
  std::cout << "Done global Jacobian\n";
  
}


Scalar ShellLinearSurfaceTensionForce2 ::elementEnergy(const std::vector<Vec3d>& vertices, 
  const std::vector<Vec3i>& faces,
  const std::vector<Scalar>& volumes,
  const Vec3i& incidentFacesPerMainVertex,
  const Vec3i& isBdryVertex) const
{
  
  std::vector<Scalar> deformed_data(3*vertices.size());
  for(unsigned int i = 0; i < vertices.size(); ++i) {
    deformed_data[3*i] = vertices[i][0];
    deformed_data[3*i+1] = vertices[i][1];
    deformed_data[3*i+2] = vertices[i][2];
  }

  adreal<NumLST2Dof,0,Real> e = LST2Energy<0>( *this, deformed_data, faces, volumes, incidentFacesPerMainVertex, isBdryVertex, m_surface_tension_coeff );
  Scalar energy = e.value();

  return energy;
}

void ShellLinearSurfaceTensionForce2 ::elementForce(const std::vector<Vec3d>& vertices, 
                                    const std::vector<Vec3i>& faces,
                                    const std::vector<Scalar>& volumes,
                                    const Vec3i& incidentFacesPerMainVertex,
                                    const Vec3i& isBdryVertex, 
                                    Eigen::Matrix<Scalar, NumLST2Dof, 1>& force) const
{
  //std::cout << "Computing element force\n";
  std::vector<Scalar> deformed_data(3*vertices.size());
  for(unsigned int i = 0; i < vertices.size(); ++i) {
    deformed_data[3*i] = vertices[i][0];
    deformed_data[3*i+1] = vertices[i][1];
    deformed_data[3*i+2] = vertices[i][2];
  }
  //std::cout << "Getting energy function.\n";
  //AutoDiff version
  adreal<NumLST2Dof,0,Real> e = LST2Energy<0>(*this, deformed_data, faces, volumes, incidentFacesPerMainVertex, isBdryVertex, m_surface_tension_coeff );     
  
  //std::cout << "Scattering back.\n";
  for(unsigned  int i = 0; i < 3*vertices.size(); i++ )
  {
    force[i] = -e.gradient(i);
  }
  //std::cout << "Done element force\n";
  
}

void ShellLinearSurfaceTensionForce2 ::elementJacobian(const std::vector<Vec3d>& vertices, 
                                                        const std::vector<Vec3i>& faces,
                                                        const std::vector<Scalar>& volumes,
                                                        const Vec3i& incidentFacesPerMainVertex,
                                                        const Vec3i& isBdryVertex, 
                                                        Eigen::Matrix<Scalar,NumLST2Dof,NumLST2Dof>& jac) const
{
  std::cout << "Computing element Jacobian\n";
  std::vector<Scalar> deformed_data(3*vertices.size());
  for(unsigned int i = 0; i < vertices.size(); ++i) {
    deformed_data[3*i] =vertices[i][0];
    deformed_data[3*i+1] = vertices[i][1];
    deformed_data[3*i+2] = vertices[i][2];
  }
  std::cout << "Set it to zero\n";
  jac.setZero();
  
  std::cout << "Call energy function\n";
  adreal<NumLST2Dof,1,Real> e = LST2Energy<1>(*this, deformed_data, faces, volumes, incidentFacesPerMainVertex, isBdryVertex, m_surface_tension_coeff );     
  // insert in the element Jacobian matrix
  std::cout << "Distributing to matrix\n";
  for(unsigned  int i = 0; i < 3*vertices.size(); i++ )
  {
    for(unsigned  int j = 0; j < 3*vertices.size(); j++ )
    {
      jac(i,j) = -e.hessian(i,j);
    }
  }

  


}



} //namespace BASim