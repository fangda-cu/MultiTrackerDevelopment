/**
 * \file TriangleMeshRenderer.cc
 *
 * \author smith@cs.columbia.edu
 * \date 03/15/2010
 */

#include "TriangleMeshRenderer.hh"

#ifdef WETA
#include "Color.hh"
#include "OpenGLDecl.hh"
#else
#include "BASim/src/Render/Color.hh"
#include "BASim/src/Render/OpenGLDecl.hh"
#endif

namespace BASim 
{
  
TriangleMeshRenderer::TriangleMeshRenderer( const TriangleMesh& mesh )
: m_mesh(mesh)
, m_mode(FLAT)
{}
  
void TriangleMeshRenderer::render()
{
  if( m_mode == FLAT )
  {
    glEnable(GL_LIGHTING);

    GLfloat gray[] = {(GLfloat)0.8f,(GLfloat)0.8f,(GLfloat)0.8f,(GLfloat)1.0f}; 
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,gray);
    
    // Render all faces
    glBegin(GL_TRIANGLES);
    //OpenGL::color(Color(255,0,0));
    for( TriangleMesh::face_iter fit = m_mesh.faces_begin(); fit != m_mesh.faces_end(); ++fit )
    {
      std::vector<Vec3d> v;
      for( TriangleMesh::FaceVertexIter fvit = m_mesh.fv_iter(*fit); fvit; ++fvit )
      {
        v.push_back(m_mesh.getVertex(*fvit));
      }
      // Compute a normal for the face
      Vec3d e1 = v[1]-v[0];
      Vec3d e2 = v[2]-v[0];
      Vec3d n = e1.cross(e2);
      if( n.norm() != 0 ) n.normalize();
      glNormal3f((GLfloat)n.x(),(GLfloat)n.y(),(GLfloat)n.z());
      glVertex3f((GLfloat)v[0].x(),(GLfloat)v[0].y(),(GLfloat)v[0].z());
      glVertex3f((GLfloat)v[1].x(),(GLfloat)v[1].y(),(GLfloat)v[1].z());
      glVertex3f((GLfloat)v[2].x(),(GLfloat)v[2].y(),(GLfloat)v[2].z());
    }
    glEnd();
    
    glDisable(GL_LIGHTING);
  }
  else if( m_mode == DBG )
  {
    glDisable(GL_LIGHTING);

    // Render all edges
    glLineWidth(2);
    glBegin(GL_LINES);
    OpenGL::color(Color(0,0,0));
    for( TriangleMesh::edge_iter eit = m_mesh.edges_begin(); eit != m_mesh.edges_end(); ++eit )
    {
      OpenGL::vertex(m_mesh.getVertex(m_mesh.fromVertex(*eit)));
      OpenGL::vertex(m_mesh.getVertex(m_mesh.toVertex(*eit)));
    }
    glEnd();
    
    // Render all faces
    glBegin(GL_TRIANGLES);
    OpenGL::color(Color(255,0,0));
    for( TriangleMesh::face_iter fit = m_mesh.faces_begin(); fit != m_mesh.faces_end(); ++fit )
    {
      for( TriangleMesh::FaceVertexIter fvit = m_mesh.fv_iter(*fit); fvit; ++fvit )
      {
        OpenGL::vertex(m_mesh.getVertex(*fvit));
      }      
    }
    glEnd();
    
    // Render all vertices
    glPointSize(5);
    glBegin(GL_POINTS);
    OpenGL::color(Color(0,0,0));
    for( TriangleMesh::vertex_iter vit = m_mesh.vertices_begin(); vit != m_mesh.vertices_end(); ++vit ) 
      OpenGL::vertex(m_mesh.getVertex(*vit));
    glEnd();
    
    glEnable(GL_LIGHTING);
  }
}

Vec3d TriangleMeshRenderer::calculateObjectCenter()
{
  Vec3d center(0.0,0.0,0.0);

  for( TriangleMesh::vertex_iter vit = m_mesh.vertices_begin(); vit != m_mesh.vertices_end(); ++vit ) 
  {
    center += m_mesh.getVertex(*vit);
  }

  if( m_mesh.nv() != 0 ) center /= ((double)m_mesh.nv());
  
  return center;
}

double TriangleMeshRenderer::calculateObjectBoundingRadius( const Vec3d& center )
{
  Scalar radius = 0.0;
  
  for( TriangleMesh::vertex_iter vit = m_mesh.vertices_begin(); vit != m_mesh.vertices_end(); ++vit )
  {
    radius = std::max(radius, (m_mesh.getVertex(*vit) - center).norm());
  }
  
  return radius;
}

} // namespace BASim
