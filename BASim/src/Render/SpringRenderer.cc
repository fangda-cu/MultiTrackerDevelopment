#include "SpringRenderer.hh"

namespace BASim 
{

SpringRenderer::SpringRenderer( const RodRodSpringForce& force )
: m_spring(force)
{
}

SpringRenderer::~SpringRenderer()
{
}

void SpringRenderer::cycleMode() {

}

void SpringRenderer::render()
{
  const Vec3d& vert0 = m_spring.getVertexA();
  const Vec3d& vert1 = m_spring.getVertexB();

  glDisable(GL_LIGHTING);
  
  glLineStipple(1, 0x00FF);
  glEnable(GL_LINE_STIPPLE);

  glLineWidth(3.0);
  glColor3d(1.0,0.0,0.0);
  glBegin(GL_LINES);
  glVertex3d(vert0.x(),vert0.y(),vert0.z());
  glVertex3d(vert1.x(),vert1.y(),vert1.z());
  glEnd();
  
  glDisable(GL_LINE_STIPPLE);
}

Vec3d SpringRenderer::calculateObjectCenter()
{
  const Vec3d& vert0 = m_spring.getVertexA();
  const Vec3d& vert1 = m_spring.getVertexB();
  
  return 0.5*(vert0+vert1);
}

Scalar SpringRenderer::calculateObjectBoundingRadius(const Vec3d& center)
{
  // TODO: implement this
  return 0.0;
}


} // namespace BASim
