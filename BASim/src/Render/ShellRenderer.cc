/**
 * \file ShellRenderer.cpp
 *
 * \author batty@cs.columbia.edu
 * \date April 20, 2011
 */

#include "BASim/src/Render/ShellRenderer.hh"
#include "BASim/src/Render/Color.hh"
#include "BASim/src/Render/OpenGLDecl.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShell.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/Math.hh"

namespace BASim 
{

void glVertVec3d(const Vec3d& v) {
  glVertex3f((GLfloat)v.x(), (GLfloat)v.y(), (GLfloat)v.z());
}

void drawThickTri(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2, Scalar thickness) {
  
  Vec3d normal = (v1-v0).cross(v2-v0);
  normal /= normal.norm();

  Vec3d v0_f = v0 + normal*thickness;
  Vec3d v0_b = v0 - normal*thickness;
  Vec3d v1_f = v1 + normal*thickness;
  Vec3d v1_b = v1 - normal*thickness;
  Vec3d v2_f = v2 + normal*thickness;
  Vec3d v2_b = v2 - normal*thickness;

  //glBegin(GL_TRIANGLES);
  glNormal3f((GLfloat)normal.x(),(GLfloat) normal.y(), (GLfloat)normal.z());
  glVertVec3d(v0_f);
  glVertVec3d(v1_f);
  glVertVec3d(v2_f);
  
  glNormal3f((GLfloat)-normal.x(), (GLfloat)-normal.y(), (GLfloat)-normal.z());
  glVertVec3d(v0_b);
  glVertVec3d(v1_b);
  glVertVec3d(v2_b);
  //glEnd();

}

void ShellRenderer::cycleMode() { 
   m_mode = (ShellRenderer::DrawMode) ((m_mode + 1) % 3); 
}

ShellRenderer::ShellRenderer( const ElasticShell& shell )
: m_shell(shell)
, m_mode(FLAT)
{}
  
void ShellRenderer::render()
{

  if( m_mode == FLAT )
  {
    glEnable(GL_LIGHTING);

    GLfloat gray[] = {0.8f,0.8f,0.8f,1.0f}; 
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,gray);

    // Render all faces
    glBegin(GL_TRIANGLES);
    //OpenGL::color(Color(255,0,0));
    const DeformableObject& mesh = m_shell.getDefoObj();
    for( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit )
    {
      std::vector<Vec3d> v;
      for( FaceVertexIterator fvit = mesh.fv_iter(*fit); fvit; ++fvit )
      {
        v.push_back(m_shell.getVertexPosition(*fvit));
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
      //drawThickTri(v[0], v[1], v[2], m_shell.getThickness(*fit));
    }
    glEnd();

    glDisable(GL_LIGHTING);
  }
  else if( m_mode == DBG )
  {
    glDisable(GL_LIGHTING);

    const DeformableObject& mesh = m_shell.getDefoObj();

    // Render all edges
    glLineWidth(2);
    glBegin(GL_LINES);
    OpenGL::color(Color(0,0,0));
    for( EdgeIterator eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit )
    {
      Vec3d p0 = m_shell.getVertexPosition(mesh.fromVertex(*eit));
      Vec3d p1 = m_shell.getVertexPosition(mesh.toVertex(*eit));
      Vec3d dir = (p1-p0);
      //p0 = p0 + 0.05*dir;
      //p1 = p1 - 0.05*dir;
      OpenGL::vertex(p0);
      OpenGL::vertex(p1);
    }
    glEnd();

    // Render all faces
    glBegin(GL_TRIANGLES);
    for( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit )
    {
    
      Vec3d barycentre;
      for( FaceVertexIterator fvit = mesh.fv_iter(*fit); fvit; ++fvit )
      {
        Vec3d pos = m_shell.getVertexPosition(*fvit);
        barycentre += pos;
      }
      barycentre /= 3.0;

      
      Scalar thickness = m_shell.getThickness(*fit);
      int colorVal = (int) (255.0 * thickness / 0.1);
      colorVal = clamp(colorVal, 0, 255);
      OpenGL::color(Color(colorVal,0,0));
      std::vector<Vec3d> points(3);
      int i = 0;
      for( FaceVertexIterator fvit = mesh.fv_iter(*fit); fvit; ++fvit )
      {
        Vec3d pos = m_shell.getVertexPosition(*fvit);
        //pos = pos - 0.05*(pos-barycentre);
        OpenGL::vertex(pos);
        points[i] = pos;
        ++i;
      }      
      
    }
    glEnd();


    // Render all vertices
    glPointSize(5);
    glBegin(GL_POINTS);
    OpenGL::color(Color(0,0,0));
    for( VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit ) {
      Vec3d vertPos = m_shell.getVertexPosition(*vit); 
      OpenGL::vertex(vertPos);
    }
    glEnd();

   /* glBegin(GL_QUADS);
    glVertex3f(-2.0f, -0.2, -2.0f);
    glVertex3f(2.0f, -0.2, -2.0f);
    glVertex3f(2.0f, -0.2, 2.0f);
    glVertex3f(-2.0f, -0.2, 2.0f);
    glEnd();*/

    glEnable(GL_LIGHTING);

    //Draw springs...

  }

}

Vec3d ShellRenderer::calculateObjectCenter()
{
  Vec3d center(0.0,0.0,0.0);
  
  const DeformableObject& mesh = m_shell.getDefoObj();
  for( VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit ) 
  {
    center += m_shell.getVertexPosition(*vit);
  }

  if( mesh.nv() != 0 ) {
    center /= ((double)mesh.nv());
  }
  
  return center;
}

double ShellRenderer::calculateObjectBoundingRadius( const Vec3d& center )
{
  Scalar radius = 0.0;
  
  const DeformableObject& mesh = m_shell.getDefoObj();
  for( VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit )
  {
    radius = std::max(radius, (m_shell.getVertexPosition(*vit) - center).norm());
  }
  
  return radius;
}

} // namespace BASim
