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
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"
#include <queue>

using Eigen::Matrix;

namespace BASim {
  
template <class S, class T>
struct less_pair_first
{
  bool operator() (const std::pair<S, T> & x, const std::pair<S, T> & y) const { return x.first < y.first; }
};

template <class S, class T>
struct less_pair_second
{
  bool operator() (const std::pair<S, T> & x, const std::pair<S, T> & y) const { return x.second < y.second; }
};

void outbin(int n)
{
  for (int i = 0; i < 8; i++)
    std::cout << ((n & (1 << i)) ? '1' : '0');
}

ShellVolumeForce::ShellVolumeForce( 
  ElasticShell& shell, 
  const std::string& name, 
  Scalar strength )
: ElasticShellForce(shell, name), m_strength(strength)
{  
  computeRefPoint();
  
  // account for contributions from BB walls
  DeformableObject & obj = m_shell.getDefoObj();
  
  std::vector<VertexHandle> new_vertices;
  std::vector<EdgeHandle> new_edges;
  std::vector<FaceHandle> new_faces;
  
  triangulateBBWalls(new_vertices, new_edges, new_faces);
  
  //Get the initial target volumes
  Scalar energy = 0;
  std::vector<int> indices(9);
  std::vector<Vec3d> deformed(3);

  int maxRegion = 0;
  //count the regions
  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;

    Vec2i labels = m_shell.getFaceLabel(fh);
    maxRegion = max(maxRegion, max(labels[0], labels[1]));
  }

  m_target_volumes.resize(maxRegion+1, 0);

  // compute the volumes due to existing faces
  fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;

    Vec2i labels = m_shell.getFaceLabel(fh);
    gatherDOFs(fh, deformed, indices);

    //determine the energy for this element
    
    //inside face
    if(labels[0] != -1)
      m_target_volumes[labels[0]] += elementEnergy(deformed);
    
    //outside face
    if(labels[1] != -1)
      m_target_volumes[labels[1]] -= elementEnergy(deformed);

  }

  // clean up
  for (size_t i = 0; i < new_faces.size(); i++)
  {
//    FaceVertexIterator fvit = obj.fv_iter(new_faces[i]); assert(fvit);
//    VertexHandle v0 = *fvit; ++fvit; assert(fvit);
//    VertexHandle v1 = *fvit; ++fvit; assert(fvit);
//    VertexHandle v2 = *fvit; ++fvit; assert(!fvit);
//    std::cout << "Deleting face: " << v0.idx() << " " << v1.idx() << " " << v2.idx() << std::endl;
    bool t = obj.deleteFace(new_faces[i], false);
    assert(t);
  }
  for (size_t i = 0; i < new_edges.size(); i++)
  {
//    std::cout << "Deleting edge: " << obj.fromVertex(new_edges[i]).idx() << " " << obj.toVertex(new_edges[i]).idx() << std::endl;
    bool t = obj.deleteEdge(new_edges[i], false);
    assert(t);
  }
  for (size_t i = 0; i < new_vertices.size(); i++)
  {
//    std::cout << "Deleting vertex: " << new_vertices[i].idx() << std::endl;
//    for (VertexEdgeIterator veit = obj.ve_iter(new_vertices[i]); veit; ++veit)
//    {
//      std::cout << "  Adjacent vertex: " << getEdgesOtherVertex(obj, *veit, new_vertices[i]).idx() << std::endl;
//    }
    bool t = obj.deleteVertex(new_vertices[i]);
    assert(t);
  }
  
//  std::cout << "After clean up nv = " << obj.nv() << " ne = " << obj.ne() << " nf = " << obj.nf() << " nt = " << obj.nt() << std::endl;
  
}

bool ShellVolumeForce::gatherDOFs(const FaceHandle& fh, std::vector<Vec3d>& deformed, std::vector<int>& indices) const {
  assert(deformed.size() == 3);
  
//  if(!m_shell.isFaceActive(fh)) return false;

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

int ShellVolumeForce::onBBWall(const Vec3d & pos) const
{
  int walls = 0;
  if (pos.x() < 0 + 1e-6)
    walls |= (1 << 0);
  if (pos.y() < 0 + 1e-6)
    walls |= (1 << 1);
  if (pos.z() < 0 + 1e-6)
    walls |= (1 << 2);
  if (pos.x() > 1 - 1e-6)
    walls |= (1 << 3);
  if (pos.y() > 1 - 1e-6)
    walls |= (1 << 4);
  if (pos.z() > 1 - 1e-6)
    walls |= (1 << 5);
  
  return walls;
}
  
void ShellVolumeForce::triangulateBBWalls(std::vector<VertexHandle> & new_vertices, std::vector<EdgeHandle> & new_edges, std::vector<FaceHandle> & new_faces) const
{
  bool verbose = true;
  
  if (verbose) std::cout << "=========================================================================" << std::endl;
  
  //////////////////////////////////////////////////////////////////////////////////
  // close the phases by triangulating the bounding box walls

  DeformableObject & obj = m_shell.getDefoObj();
  
  if (verbose)
  {
    for (VertexIterator vit = obj.vertices_begin(); vit != obj.vertices_end(); ++vit)
    {
      std::cout << "Vertex " << (*vit).idx() << ": " << obj.getVertexPosition(*vit) << std::endl;
    }
  }
         
  std::vector<std::vector<EdgeHandle> > wall_edges(6);
  EdgeProperty<Vec2i> wall_edge_labels(&obj);
  std::vector<VertexHandle> corners;
  
  // count how many regions there are
  int max_label = 0;
  for (FaceIterator fit = obj.faces_begin(); fit != obj.faces_end(); ++fit)
  {
    Vec2i labels = m_shell.getFaceLabel(*fit);
    
    if (labels.x() >= max_label)
      max_label = labels.x();
    if (labels.y() >= max_label)
      max_label = labels.y();
  }
  
  int nregion = max_label + 1;
  
  // sort the film boundary edges to six walls
  wall_edge_labels.assign(Vec2i(-1, -1));
  Eigen::Matrix<int, Eigen::Dynamic, 1> rcounts;
  
  for (EdgeIterator eit = obj.edges_begin(); eit != obj.edges_end(); ++eit)
  {
    if (obj.isBoundary(*eit))
    {
      int walls0 = onBBWall(m_shell.getVertexPosition(obj.fromVertex(*eit)));
      int walls1 = onBBWall(m_shell.getVertexPosition(obj.toVertex(*eit)));
      assert(walls0 & walls1);
      
      Vec2i edge_label(-1, -1);
      rcounts.setZero(nregion);
      for (EdgeFaceIterator efit = obj.ef_iter(*eit); efit; ++efit)
      {
        Vec2i label = m_shell.getFaceLabel(*efit);
        assert(label.x() >= 0 && label.y() >= 0);
        rcounts[label.x()] += obj.getRelativeOrientation(*efit, *eit);
        rcounts[label.y()] -= obj.getRelativeOrientation(*efit, *eit);
      }
      
      for (int i = 0; i < nregion; i++)
      {
        if (rcounts[i] == 1)
        {
          assert(edge_label.x() == -1);
          edge_label.x() = i; // x component is the region label on the right of the edge
        } else if (rcounts[i] == -1)
        {
          assert(edge_label.y() == -1);
          edge_label.y() = i; // y component is the region label on the left of the edge
        } else
        {
          assert(rcounts[i] == 0);
        }
      }
      
      assert(edge_label.x() >= 0 && edge_label.y() >= 0);
      
      int wall = walls0 & walls1;
      std::vector<int> walls;
      for (int i = 0; i < 6; i++)
        if (wall & (1 << i))
          wall_edges[i].push_back(*eit), walls.push_back(i);
      wall_edge_labels[*eit] = edge_label;
      
      assert(walls.size() == 0 || walls.size() == 1 || walls.size() == 2);
      
    }
  }
  
  if (verbose) 
  {
    std::cout << "Wall edges: " << std::endl;
    for (int i = 0; i < 6; i++)
      for (int j = 0; j < wall_edges[i].size(); j++)
        std::cout << "i = " << i << " j = " << j << " edge: " << obj.fromVertex(wall_edges[i][j]).idx() << " - " << obj.toVertex(wall_edges[i][j]).idx() << " label: " << wall_edge_labels[wall_edges[i][j]] << std::endl;
    
    std::cout << "nv = " << obj.nv() << " ne = " << obj.ne() << " nf = " << obj.nf() << " nt = " << obj.nt() << std::endl;
  }
  
  // temporarily add vertices and edges in order to triangulate BB walls
  corners.push_back(obj.addVertex()); m_shell.setVertexPosition(corners.back(), Vec3d(0, 0, 0));
  corners.push_back(obj.addVertex()); m_shell.setVertexPosition(corners.back(), Vec3d(1, 0, 0));
  corners.push_back(obj.addVertex()); m_shell.setVertexPosition(corners.back(), Vec3d(0, 1, 0));
  corners.push_back(obj.addVertex()); m_shell.setVertexPosition(corners.back(), Vec3d(0, 0, 1));
  corners.push_back(obj.addVertex()); m_shell.setVertexPosition(corners.back(), Vec3d(1, 1, 0));
  corners.push_back(obj.addVertex()); m_shell.setVertexPosition(corners.back(), Vec3d(1, 0, 1));
  corners.push_back(obj.addVertex()); m_shell.setVertexPosition(corners.back(), Vec3d(0, 1, 1));
  corners.push_back(obj.addVertex()); m_shell.setVertexPosition(corners.back(), Vec3d(1, 1, 1));
  std::vector<int> corner_labels(corners.size(), -1);
  
  if (verbose) std::cout << "-----------------------------------------------------------------------" << std::endl;
  
  std::vector<Vec4i> bb_edges(12);  // x = from corner; y = to corner; z = left wall; w = right wall
  bb_edges[ 0].x() = 0;  bb_edges[ 0].y() = 1;  bb_edges[ 0].z() = 1;   bb_edges[ 0].w() = 2; // BB edge: vertices 0, 1
  bb_edges[ 1].x() = 0;  bb_edges[ 1].y() = 2;  bb_edges[ 1].z() = 2;   bb_edges[ 1].w() = 0; // BB edge: vertices 0, 2
  bb_edges[ 2].x() = 0;  bb_edges[ 2].y() = 3;  bb_edges[ 2].z() = 0;   bb_edges[ 2].w() = 1; // BB edge: vertices 0, 3
  bb_edges[ 3].x() = 1;  bb_edges[ 3].y() = 4;  bb_edges[ 3].z() = 3;   bb_edges[ 3].w() = 2; // BB edge: vertices 1, 4
  bb_edges[ 4].x() = 1;  bb_edges[ 4].y() = 5;  bb_edges[ 4].z() = 1;   bb_edges[ 4].w() = 3; // BB edge: vertices 1, 5
  bb_edges[ 5].x() = 2;  bb_edges[ 5].y() = 4;  bb_edges[ 5].z() = 2;   bb_edges[ 5].w() = 4; // BB edge: vertices 2, 4
  bb_edges[ 6].x() = 2;  bb_edges[ 6].y() = 6;  bb_edges[ 6].z() = 4;   bb_edges[ 6].w() = 0; // BB edge: vertices 2, 6
  bb_edges[ 7].x() = 3;  bb_edges[ 7].y() = 5;  bb_edges[ 7].z() = 5;   bb_edges[ 7].w() = 1; // BB edge: vertices 3, 5
  bb_edges[ 8].x() = 3;  bb_edges[ 8].y() = 6;  bb_edges[ 8].z() = 0;   bb_edges[ 8].w() = 5; // BB edge: vertices 3, 6
  bb_edges[ 9].x() = 4;  bb_edges[ 9].y() = 7;  bb_edges[ 9].z() = 3;   bb_edges[ 9].w() = 4; // BB edge: vertices 4, 7
  bb_edges[10].x() = 5;  bb_edges[10].y() = 7;  bb_edges[10].z() = 5;   bb_edges[10].w() = 3; // BB edge: vertices 5, 7
  bb_edges[11].x() = 6;  bb_edges[11].y() = 7;  bb_edges[11].z() = 4;   bb_edges[11].w() = 5; // BB edge: vertices 6, 7
  
  std::vector<Vec3d> wall_normals(6);
  wall_normals[0] = Vec3d(-1, 0, 0);
  wall_normals[1] = Vec3d(0, -1, 0);
  wall_normals[2] = Vec3d(0, 0, -1);
  wall_normals[3] = Vec3d(1, 0, 0);
  wall_normals[4] = Vec3d(0, 1, 0);
  wall_normals[5] = Vec3d(0, 0, 1);
  
  // first find the boundary vertices lying on BB edges
  std::vector<std::vector<VertexHandle> > edge_verts(12);
  for (VertexIterator vit = obj.vertices_begin(); vit != obj.vertices_end(); ++vit)
    if (obj.isBoundary(*vit))
      for (int i = 0; i < 12; i++)
        if ((~onBBWall(m_shell.getVertexPosition(*vit)) & ((1 << bb_edges[i].z()) | (1 << bb_edges[i].w()))) == 0)  // this vertex is on both walls
          edge_verts[i].push_back(*vit);

  // add edges lying on the BB edges to close the regions
  // TODO: the following code can use optimization
  for (int i = 0; i < 12; i++)  // for each BB edge
  {
    VertexHandle corner0 = corners[bb_edges[i].x()];
    VertexHandle corner1 = corners[bb_edges[i].y()];
    int wall0 = bb_edges[i].z();
    int wall1 = bb_edges[i].w();
    
    int edge_mask = ((1 << wall0) | (1 << wall1));
    
    if (verbose) std::cout << "#######################\nedge i = " << i << " from corner " << bb_edges[i].x() << " to corner " << bb_edges[i].y() << " with wall " << wall0 << " on the left and wall " << wall1 << " on the right " << std::endl;
    
    std::vector<std::pair<VertexHandle, Scalar> > evs;  // vertices on this edge, but sortable by distance from vertex corner0
    for (size_t j = 0; j < edge_verts[i].size(); j++)
      evs.push_back(std::pair<VertexHandle, Scalar>(edge_verts[i][j], (m_shell.getVertexPosition(edge_verts[i][j]) - m_shell.getVertexPosition(corner0)).norm()));
    
    if (evs.size() == 0)
    {
      if (verbose) std::cout << "no vertex" << std::endl;
      
      // no vertex on this edge
      if (!findEdge(obj, corner0, corner1).isValid())
      {
        new_edges.push_back(obj.addEdge(corner0, corner1));
        wall_edges[wall0].push_back(new_edges.back());
        wall_edges[wall1].push_back(new_edges.back());
        wall_edge_labels[new_edges.back()] = Vec2i(-1, -1);
        if (verbose) std::cout << "Edge " << obj.fromVertex(new_edges.back()).idx() << " - " << obj.toVertex(new_edges.back()).idx() << " added to walls " << wall0 << " and " << wall1 << " and labeled -1" << std::endl;
      }
    } else
    {
      std::sort(evs.begin(), evs.end(), less_pair_second<VertexHandle, Scalar>());
      std::vector<int> edge_labels(evs.size() + 1, -1);
      for (size_t l = 0; l < evs.size(); l++)
      {
        VertexHandle & v = evs[l].first;

        if (verbose) std::cout << "left wall: " << std::endl;
        int lregion0 = -1;
        int lregion1 = -1;
        int head = -1;
        int tail = -1;
        Vec3d head_vec, tail_vec;
        for (size_t j = 0; j < wall_edges[wall0].size(); j++)
        {
          Vec3d e;
          bool b = false;
          if (obj.fromVertex(wall_edges[wall0][j]) == v)
          {
            if (~onBBWall(m_shell.getVertexPosition(obj.toVertex(wall_edges[wall0][j]))) & edge_mask) // this vertex is not also on this edge
            {
              e = m_shell.getVertexPosition(obj.toVertex(wall_edges[wall0][j])) - m_shell.getVertexPosition(obj.fromVertex(wall_edges[wall0][j]));
              b = true;
            }
          } else if (obj.toVertex(wall_edges[wall0][j]) == v)
          {
            if (~onBBWall(m_shell.getVertexPosition(obj.fromVertex(wall_edges[wall0][j]))) & edge_mask) // this vertex is not also on this edge
            {
              e = m_shell.getVertexPosition(obj.fromVertex(wall_edges[wall0][j])) - m_shell.getVertexPosition(obj.toVertex(wall_edges[wall0][j]));
              b = true;
            }
          }
          
          if (b)
          {
            if (head < 0 || e.cross(head_vec).dot(wall_normals[wall0]) > 0)
            {
              head = j;
              head_vec = e;
            }
            if (tail < 0 || e.cross(head_vec).dot(wall_normals[wall0]) < 0)
            {
              tail = j;
              tail_vec = e;
            }
          }
        }
        
        if (head >= 0 && tail >= 0)
        {
          EdgeHandle headedge = wall_edges[wall0][head];
          EdgeHandle tailedge = wall_edges[wall0][tail];
          lregion0 = (obj.fromVertex(headedge) == v ? wall_edge_labels[headedge].y() : wall_edge_labels[headedge].x());
          lregion1 = (obj.fromVertex(tailedge) == v ? wall_edge_labels[tailedge].x() : wall_edge_labels[tailedge].y());
        }
        
        if (verbose) std::cout << "right wall: " << std::endl;
        int rregion0 = -1;
        int rregion1 = -1;
        head = -1;
        tail = -1;        
        for (size_t j = 0; j < wall_edges[wall1].size(); j++)
        {
          Vec3d e;
          bool b = false;
          if (obj.fromVertex(wall_edges[wall1][j]) == v)
          {
            if (~onBBWall(m_shell.getVertexPosition(obj.toVertex(wall_edges[wall1][j]))) & edge_mask) // this vertex is not also on this edge
            {
              e = m_shell.getVertexPosition(obj.toVertex(wall_edges[wall1][j])) - m_shell.getVertexPosition(obj.fromVertex(wall_edges[wall1][j]));
              b = true;
            }
          } else if (obj.toVertex(wall_edges[wall1][j]) == v)
          {
            if (~onBBWall(m_shell.getVertexPosition(obj.fromVertex(wall_edges[wall1][j]))) & edge_mask) // this vertex is not also on this edge
            {
              e = m_shell.getVertexPosition(obj.fromVertex(wall_edges[wall1][j])) - m_shell.getVertexPosition(obj.toVertex(wall_edges[wall1][j]));
              b = true;
            }
          }
          
          if (b)
          {
            if (head < 0 || e.cross(head_vec).dot(wall_normals[wall1]) < 0)
            {
              head = j;
              head_vec = e;
            }
            if (tail < 0 || e.cross(head_vec).dot(wall_normals[wall1]) > 0)
            {
              tail = j;
              tail_vec = e;
            }
          }
        }
        
        if (head >= 0 && tail >= 0)
        {
          EdgeHandle headedge = wall_edges[wall1][head];
          EdgeHandle tailedge = wall_edges[wall1][tail];
          rregion0 = (obj.fromVertex(headedge) == v ? wall_edge_labels[headedge].x() : wall_edge_labels[headedge].y());
          rregion1 = (obj.fromVertex(tailedge) == v ? wall_edge_labels[tailedge].y() : wall_edge_labels[tailedge].x());
        }

//        assert(!(lregion0 >= 0 && rregion0 >= 0 && lregion0 != rregion0));
//        assert(!(lregion1 >= 0 && rregion1 >= 0 && lregion1 != rregion1));

        if (lregion0 >= 0)
          edge_labels[l] = lregion0;
        if (rregion0 >= 0)
          edge_labels[l] = rregion0;
        if (lregion1 >= 0)
          edge_labels[l + 1] = lregion1;
        if (rregion1 >= 0)
          edge_labels[l + 1] = rregion1;
      }
      
      if (verbose)
      {
        std::cout << "Edge segment vertices: ";
        for (size_t i = 0; i < evs.size(); i++)
        {
          std::cout << evs[i].first.idx() << " ";
        }
        std::cout << std::endl << "Edge segment labels: ";
        for (size_t i = 0; i < evs.size() + 1; i++)
        {
          std::cout << edge_labels[i] << " ";
        }
        std::cout << std::endl;
      }
      
      for (size_t l = 0; l < evs.size() + 1; l++)
      {
        VertexHandle v0 = (l == 0 ?          corner0 : evs[l - 1].first);
        VertexHandle v1 = (l == evs.size() ? corner1 : evs[l].first);
        
        if (!findEdge(obj, v0, v1).isValid())
        {
          int edge_label = edge_labels[l];
          assert(edge_label >= 0);  // if edge_labels[l] is left -1, it must have been because there is an edge between v0 and v1
          
          new_edges.push_back(obj.addEdge(v0, v1));
          
          Vec2i label(edge_label, edge_label);
          wall_edges[wall0].push_back(new_edges.back());
          wall_edges[wall1].push_back(new_edges.back());
          wall_edge_labels[new_edges.back()] = Vec2i(edge_label, edge_label);
          
          if (verbose) std::cout << "Edge " << obj.fromVertex(new_edges.back()).idx() << " - " << obj.toVertex(new_edges.back()).idx() << " added to walls " << wall0 << " and " << wall1 << " and labeled " << edge_label << std::endl;
        }
      }
      
      if (edge_labels.front() >= 0)
      {
        assert(corner_labels[bb_edges[i].x()] <= 0 || corner_labels[bb_edges[i].x()] == edge_labels.front());
        corner_labels[bb_edges[i].x()] = edge_labels.front();
      }
      if (edge_labels.back() >= 0)
      {
        assert(corner_labels[bb_edges[i].y()] <= 0 || corner_labels[bb_edges[i].y()] == edge_labels.back());
        corner_labels[bb_edges[i].y()] = edge_labels.back();
      }
    }
  }
  
  // a few passes to propagate the corner labels through BB edges
  for (int i = 0; i < 8; i++) 
  {
    if (verbose) 
    {
      std::cout << "Propagation iteration i = " << i << std::endl;
      for (int j = 0; j < 8; j++)
      {
        std::cout << "  Vertex " << corners[j].idx() << " label " << corner_labels[j] << " edges: ";
        for (VertexEdgeIterator veit = obj.ve_iter(corners[j]); veit; ++veit)
          std::cout << " vertex " << getEdgesOtherVertex(obj, *veit, corners[j]).idx() << " label " << wall_edge_labels[*veit] << "; ";
        std::cout << std::endl;
      }
    }    
    
    for (int j = 0; j < 8; j++)
    {
      if (corner_labels[j] < 0)
      {
        for (VertexEdgeIterator veit = obj.ve_iter(corners[j]); veit; ++veit)
        {
          if (wall_edge_labels[*veit].x() >= 0)
          {
            assert(wall_edge_labels[*veit].x() == wall_edge_labels[*veit].y());
            corner_labels[j] = wall_edge_labels[*veit].x();
          }
        }
      } else
      {
        for (VertexEdgeIterator veit = obj.ve_iter(corners[j]); veit; ++veit)
        {
          if (wall_edge_labels[*veit].x() >= 0)
          {
            assert(wall_edge_labels[*veit].x() == wall_edge_labels[*veit].y());
            assert(corner_labels[j] == wall_edge_labels[*veit].x());
          } else
          {
            wall_edge_labels[*veit].x() = wall_edge_labels[*veit].y() = corner_labels[j];
          }
        }
      }
    }
  }
  
  for (int i = 0; i < 8; i++)
  {
    assert(corner_labels[i] >= 0);  // cannot allow the case where no phase interface intersects any bounding box wall.
  }
  
  // add one vertex at the center of each wall serving as triangulation pivots. 
  // actually any arbitrary vertex is fine for this task but using a non-existent vertex avoids degenerate triangles (not that they are any problems, just cleaner)
  std::vector<VertexHandle> wall_pivots;
//  wall_pivots.push_back(obj.addVertex()); m_shell.setVertexPosition(wall_pivots.back(), Vec3d(0.0, 0.5, 0.5));
//  wall_pivots.push_back(obj.addVertex()); m_shell.setVertexPosition(wall_pivots.back(), Vec3d(0.5, 0.0, 0.5));
//  wall_pivots.push_back(obj.addVertex()); m_shell.setVertexPosition(wall_pivots.back(), Vec3d(0.5, 0.5, 0.0));
//  wall_pivots.push_back(obj.addVertex()); m_shell.setVertexPosition(wall_pivots.back(), Vec3d(1.0, 0.5, 0.5));
//  wall_pivots.push_back(obj.addVertex()); m_shell.setVertexPosition(wall_pivots.back(), Vec3d(0.5, 1.0, 0.5));
//  wall_pivots.push_back(obj.addVertex()); m_shell.setVertexPosition(wall_pivots.back(), Vec3d(0.5, 0.5, 1.0));
  Scalar pivot_extrusion = 0.0; // non zero values are only for visual debugging
  wall_pivots.push_back(obj.addVertex()); m_shell.setVertexPosition(wall_pivots.back(), Vec3d(-pivot_extrusion, 0.5, 0.5));
  wall_pivots.push_back(obj.addVertex()); m_shell.setVertexPosition(wall_pivots.back(), Vec3d(0.5, -pivot_extrusion, 0.5));
  wall_pivots.push_back(obj.addVertex()); m_shell.setVertexPosition(wall_pivots.back(), Vec3d(0.5, 0.5, -pivot_extrusion));
  wall_pivots.push_back(obj.addVertex()); m_shell.setVertexPosition(wall_pivots.back(), Vec3d(1.0 + pivot_extrusion, 0.5, 0.5));
  wall_pivots.push_back(obj.addVertex()); m_shell.setVertexPosition(wall_pivots.back(), Vec3d(0.5, 1.0 + pivot_extrusion, 0.5));
  wall_pivots.push_back(obj.addVertex()); m_shell.setVertexPosition(wall_pivots.back(), Vec3d(0.5, 0.5, 1.0 + pivot_extrusion));

  // triangulate the walls, using the labels on each boundary edge in the wall
  for (int i = 0; i < 6; i++)
  {
    if (wall_edges[i].size() == 0)
      continue;
    
    VertexHandle wall_pivot = wall_pivots[i];
    
    for (size_t j = 0; j < wall_edges[i].size(); j++)
    {
      EdgeHandle e = wall_edges[i][j];

      if (!findEdge(obj, wall_pivot, obj.fromVertex(e)).isValid())
        new_edges.push_back(obj.addEdge(wall_pivot, obj.fromVertex(e)));
      if (!findEdge(obj, wall_pivot, obj.toVertex(e)).isValid())
        new_edges.push_back(obj.addEdge(wall_pivot, obj.toVertex(e)));
      new_faces.push_back(obj.addFace(wall_pivot, obj.fromVertex(e), obj.toVertex(e)));
      
      std::vector<Vec3d> pos(3);
      pos[0] = obj.getVertexPosition(wall_pivot);
      pos[1] = obj.getVertexPosition(obj.fromVertex(e));
      pos[2] = obj.getVertexPosition(obj.toVertex(e));
            
      if (wall_edge_labels[e].x() == wall_edge_labels[e].y())
      {
        // this is an edge along an edge of BB
        // therefore, the following signed volume cannot be zero
        int sign = ((pos[0] - m_ref_point).cross(pos[1] - m_ref_point).dot(pos[2] - m_ref_point) > 0 ? 1 : -1);        
        m_shell.setFaceLabel(new_faces.back(), sign > 0 ? Vec2i(wall_edge_labels[e].x(), -1) : Vec2i(-1, wall_edge_labels[e].x()));
      } else
      {
        // an edge within the wall
        m_shell.setFaceLabel(new_faces.back(), Vec2i(wall_edge_labels[e].y(), wall_edge_labels[e].x()));
      }
    }
  }
  
  new_vertices.reserve(corners.size() + wall_pivots.size());
  for (size_t i = 0; i < corners.size(); i++)
    new_vertices.push_back(corners[i]);
  for (size_t i = 0; i < wall_pivots.size(); i++)
    new_vertices.push_back(wall_pivots[i]);
    

}

Scalar ShellVolumeForce::globalEnergy() const
{
  // account for the BB walls
  DeformableObject & obj = m_shell.getDefoObj();
  
  std::vector<VertexHandle> new_vertices;
  std::vector<EdgeHandle> new_edges;
  std::vector<FaceHandle> new_faces;
  
  triangulateBBWalls(new_vertices, new_edges, new_faces);
  
  for (size_t i = 0; i < new_vertices.size(); i++)
    obj.setPositionDofBase(new_vertices[i], -10);

  // compute the volumes due to existing faces
  Scalar energy = 0;
  std::vector<int> indices(9);
  std::vector<Vec3d> deformed(3);
  
  if(m_strength == 0) return 0;
  
  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  
  std::vector<Scalar> volumes(m_target_volumes.size(), 0);
  
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;
    
    gatherDOFs(fh, deformed, indices);
    
    Vec2i labels = m_shell.getFaceLabel(fh);
    
    if(labels[0] != -1)
      volumes[labels[0]] += elementEnergy(deformed);
    
    if(labels[1] != -1)
      volumes[labels[1]] -= elementEnergy(deformed);    
  }

  Scalar sum = 0;
  for(unsigned int r = 0; r < volumes.size(); ++r)
    sum += 0.5 * m_strength * (volumes[r] - m_target_volumes[r])*(volumes[r] - m_target_volumes[r]);
  
//  std::cout << "Total energy: " << sum << std::endl;
  
  // clean up
  for (size_t i = 0; i < new_faces.size(); i++)
  {
//    FaceVertexIterator fvit = obj.fv_iter(new_faces[i]); assert(fvit);
//    VertexHandle v0 = *fvit; ++fvit; assert(fvit);
//    VertexHandle v1 = *fvit; ++fvit; assert(fvit);
//    VertexHandle v2 = *fvit; ++fvit; assert(!fvit);
//    std::cout << "Deleting face: " << v0.idx() << " " << v1.idx() << " " << v2.idx() << std::endl;
    bool t = obj.deleteFace(new_faces[i], false);
    assert(t);
  }
  for (size_t i = 0; i < new_edges.size(); i++)
  {
//    std::cout << "Deleting edge: " << obj.fromVertex(new_edges[i]).idx() << " " << obj.toVertex(new_edges[i]).idx() << std::endl;
    bool t = obj.deleteEdge(new_edges[i], false);
    assert(t);
  }
  for (size_t i = 0; i < new_vertices.size(); i++)
  {
//    std::cout << "Deleting vertex: " << new_vertices[i].idx() << std::endl;
    bool t = obj.deleteVertex(new_vertices[i]);
    assert(t);
  }
  
//  std::cout << "After clean up nv = " << obj.nv() << " ne = " << obj.ne() << " nf = " << obj.nf() << " nt = " << obj.nt() << std::endl;
  
  return sum;
}

void ShellVolumeForce::globalForce( VecXd& force )  const
{
  if (m_strength == 0) return;
  
  DeformableObject & obj = m_shell.getDefoObj();
  
//  std::cout << "ndof = " << force.size() << std::endl;
//  
//  std::cout << "Vertex positions:" << std::endl;
//  for (VertexIterator vit = obj.vertices_begin(); vit != obj.vertices_end(); ++vit)
//    std::cout << "Vertex " << (*vit).idx() << " position " << m_shell.getVertexPosition(*vit) << std::endl;
  
  // account for the BB walls  
  std::vector<VertexHandle> new_vertices;
  std::vector<EdgeHandle> new_edges;
  std::vector<FaceHandle> new_faces;
  
  triangulateBBWalls(new_vertices, new_edges, new_faces);
  
  for (size_t i = 0; i < new_vertices.size(); i++)
    obj.setPositionDofBase(new_vertices[i], -10);
  
  // compute the volumes due to existing faces
  std::vector<int> indices(9);
  std::vector<Vec3d> deformed(3);
  Eigen::Matrix<Scalar, 9, 1> localForce;
  
  std::vector<Scalar> volumes(m_target_volumes.size(), 0);
  
  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;
    
    bool valid = gatherDOFs(fh, deformed, indices);
//    if(!valid) continue;
    
    Vec2i labels = m_shell.getFaceLabel(fh);
    
    if(labels[0] != -1)
      volumes[labels[0]] += elementEnergy(deformed);
    if(labels[1] != -1)
      volumes[labels[1]] -= elementEnergy(deformed);
    
  }

  //then compute forces, which relies on the volumes above
  fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;
   
    bool valid = gatherDOFs(fh, deformed, indices);
//    if(!valid) continue;
    
    Vec2i labels = m_shell.getFaceLabel(fh);
    
    elementForce(deformed, localForce);
    
    if(labels[0] != -1) {
      for (unsigned int i = 0; i < indices.size(); ++i)
        if (indices[i] >= 0)
          force(indices[i]) += m_strength * (volumes[labels[0]] - m_target_volumes[labels[0]]) * localForce(i);
    }
    
    if(labels[1] != -1) {
      for (unsigned int i = 0; i < indices.size(); ++i)
        if (indices[i] >= 0)
          force(indices[i]) -= m_strength * (volumes[labels[1]] - m_target_volumes[labels[1]]) * localForce(i);
    }
    
  }
  
  // clean up
  for (size_t i = 0; i < new_faces.size(); i++)
    obj.deleteFace(new_faces[i], false);
  for (size_t i = 0; i < new_edges.size(); i++)
    obj.deleteEdge(new_edges[i], false);
  for (size_t i = 0; i < new_vertices.size(); i++)
    obj.deleteVertex(new_vertices[i]);
  
}

void ShellVolumeForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  if (m_strength == 0) return;
  
  DeformableObject & obj = m_shell.getDefoObj();
  
  // account for the BB walls
  std::vector<VertexHandle> new_vertices;
  std::vector<EdgeHandle> new_edges;
  std::vector<FaceHandle> new_faces;
  
  triangulateBBWalls(new_vertices, new_edges, new_faces);
  
  for (size_t i = 0; i < new_vertices.size(); i++)
    obj.setPositionDofBase(new_vertices[i], -10);
  
  // compute the volumes due to existing faces
  std::vector<int> indices(9);
  std::vector<Vec3d> deformed(3);
  Eigen::Matrix<Scalar, 9, 9> localMatrix;
  Eigen::Matrix<Scalar, 9, 1> localForce;

  std::vector<Scalar> volumes(m_target_volumes.size(), 0);

  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;

    bool valid = gatherDOFs(fh, deformed, indices);
//    if(!valid) continue;

    Vec2i labels = m_shell.getFaceLabel(fh);
    
    if(labels[0] != -1)
      volumes[labels[0]] += elementEnergy(deformed);

    if(labels[1] != -1)
      volumes[labels[1]] -= elementEnergy(deformed);
  }
  
//  //compute the total gradient of volume
//  std::vector<VecXd> tgv(m_target_volumes.size(), VecXd::Zero(obj.ndof()));
//  
//  fit = m_shell.getDefoObj().faces_begin();
//  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
//    const FaceHandle& fh = *fit;
//    
//    bool valid = gatherDOFs(fh, deformed, indices);
////    if(!valid) continue;
//    
//    Vec2i labels = m_shell.getFaceLabel(fh);
//    
//    elementForce(deformed, localForce);
//    
//    if(labels[0] != -1) {
//      for (unsigned int i = 0; i < indices.size(); ++i)
//        if (indices[i] >= 0)
//          tgv[labels[0]](indices[i]) += localForce(i);
//    }
//    
//    if(labels[1] != -1) {
//      for (unsigned int i = 0; i < indices.size(); ++i)
//        if (indices[i] >= 0)
//          tgv[labels[1]](indices[i]) -= localForce(i);
//    }
//    
//  }

  //compute force jacobians, which relies on the volumes above
  fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;

    bool valid = gatherDOFs(fh, deformed, indices);
//    if(!valid) continue;
    
    Vec2i labels = m_shell.getFaceLabel(fh);
    elementJacobian(deformed, localMatrix);
    elementForce(deformed, localForce);
    
//    std::cout << "Face with vertices: ";
//    for (FaceVertexIterator fvit = obj.fv_iter(*fit); fvit; ++fvit)
//      std::cout << (*fvit).idx() << " ";
//    std::cout << " with dofs: ";
//    for (int i = 0; i < 9; i++)
//      std::cout << indices[i] << " ";
//    std::cout << std::endl;
//    std::cout << " local force = " << localForce << std::endl;
    
    if(labels[0] != -1) {
      for (unsigned int i = 0; i < indices.size(); ++i)
        for(unsigned int j = 0; j < indices.size(); ++j)
          if (indices[i] >= 0 && indices[j] >= 0)
            Jacobian.add(indices[i], indices[j], +m_strength * scale * ((volumes[labels[0]] - m_target_volumes[labels[0]]) * localMatrix(i,j)));
//      for (unsigned int i = 0; i < indices.size(); ++i)
//        for (int j = 0; j < obj.ndof(); ++j)
//          if (indices[i] >= 0)
//            Jacobian.add(indices[i], j, -m_strength * scale * tgv[labels[0]][j] * localForce[i]);
    }

    if(labels[1] != -1) {
      for (unsigned int i = 0; i < indices.size(); ++i)
        for(unsigned int j = 0; j < indices.size(); ++j)
          if (indices[i] >= 0 && indices[j] >= 0)
            Jacobian.add(indices[i], indices[j], -m_strength * scale * ((volumes[labels[1]] - m_target_volumes[labels[1]]) * localMatrix(i,j)));
//      for (unsigned int i = 0; i < indices.size(); ++i)
//        for (int j = 0; j < obj.ndof(); ++j)
//          if (indices[i] >= 0)
//            Jacobian.add(indices[i], j, +m_strength * scale * tgv[labels[1]][j] * localForce[i]);
    }
    
  }  
  
  // clean up
  for (size_t i = 0; i < new_faces.size(); i++)
    obj.deleteFace(new_faces[i], false);
  for (size_t i = 0; i < new_edges.size(); i++)
    obj.deleteEdge(new_edges[i], false);
  for (size_t i = 0; i < new_vertices.size(); i++)
    obj.deleteVertex(new_vertices[i]);
  
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