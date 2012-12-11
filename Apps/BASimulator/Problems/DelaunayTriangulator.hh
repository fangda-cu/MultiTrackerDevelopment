/**
 * \file DelaunayTriangulator.hh
 *
 *  3D Delaunay triangulation 
 *
 *  reference:
 *    http://www.cse.iitk.ac.in/users/vision/dipakmj/papers/04276112.pdf
 *
 *  The implementation is highly unoptimized.
 *
 * \author fang@cs.columbia.edu
 * \date Dec 10, 2012
 */

#ifndef DELAUNAY_TRIANGULATOR_HH
#define DELAUNAY_TRIANGULATOR_HH

#include <vector>
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"

using namespace BASim;

namespace DelaunayTriangulator
{
  typedef TopologicalObject Mesh;
  
  class DelaunayTriangulator
  {
  public:
    DelaunayTriangulator(Mesh * mesh, VertexProperty<Vec3d> & positions);
    
  public:
    Mesh * mesh() { return m_mesh; }
    VertexProperty<Vec3d> & positions() { return m_positions; }
    
  public:
    bool insertVertex(Vec3d & v);
    void flip23(TetHandle pabc, TetHandle abcd, FaceHandle abc, TetHandle & pabd, TetHandle & pacd, TetHandle & pbcd);
    void flip32(TetHandle pabc, TetHandle abcd, TetHandle pabd, FaceHandle abc, TetHandle & pacd, TetHandle & pbcd);
    
    bool checkDelaunay();
    
    bool extractVoronoiDiagram(Mesh * tomesh, VertexProperty<Vec3d> & pos);
    
  public:
    static Scalar predicateOriented(const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d & p);
    static Scalar predicateInSphere(const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d & d, const Vec3d & p);
    
    static bool predicateInTetrahedron(const Vec3d & a, const Vec3d & b, const Vec3d & c, const Vec3d & d, const Vec3d & p);
    
    Scalar predicateOriented(VertexHandle a, VertexHandle b, VertexHandle c, VertexHandle p) const;
    Scalar predicateInSphere(VertexHandle a, VertexHandle b, VertexHandle c, VertexHandle d, VertexHandle p) const;
    
    bool predicateInTetrahedron(VertexHandle a, VertexHandle b, VertexHandle c, VertexHandle d, VertexHandle p) const;
    
    Vec3d & pos(VertexHandle v) { return m_positions[v]; }
    const Vec3d & pos(VertexHandle v) const { return m_positions[v]; }
    
  protected:
    Mesh * m_mesh;
    VertexProperty<Vec3d> m_positions;
    
    std::vector<VertexHandle> m_vertices;
    
  };

}

#endif
