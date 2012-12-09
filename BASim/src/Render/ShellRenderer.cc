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
#include "BASim/src/Physics/DeformableObjects/Shells/MNBendingForce.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"
#include "BASim/src/Render/Curvature.hh"

//ElTopo in order to test out the remeshing operations
#include "surftrack.h"


namespace BASim 
{
GLfloat black[] ={0.0f, 0.0f, 0.0f, 1.0f};
GLfloat red[] ={0.0f, 0.0f, 1.0f, 1.0f};

//Color routine to use HSV colors instead of RGB taken
//from http://en.wikipedia.org/wiki/HSL_and_HSV#Converting_to_RGB
// HSV \in (0-360, 0-1, 0-1)
void hsvToRgb(const GLdouble _h, const GLdouble _s, const GLdouble _v,
        GLdouble & r,GLdouble & g,GLdouble & b){
    GLdouble c, m, x, hp;
    GLdouble h = clamp<Scalar>(_h, 0, 360);
    GLdouble s = clamp<Scalar>(_s, 0, 1);
    GLdouble v = clamp<Scalar>(_v, 0, 1);
    c = v * s;
    m = v - c;
    hp = h / 60.;
    x = c * (GLdouble) ( 1 - abs((GLint)hp % 2 - 1));
    assert ( hp < 6. && hp >= 0);
    if ( hp < 1){
        r = c;
        g = x;
        b = 0;
    } else if ( hp < 2){
        r = x;
        g = c;
        b = 0;
    } else if ( hp < 3){
        r = 0;
        g = c;
        b = x;
    } else if ( hp < 4 ){
        r = 0;
        g = x;
        b = c;
    } else if ( hp < 5){
        r = x;
        g = 0;
        b = c;
    } else if ( hp < 6){
        r = c;
        g = 0;
        b = x;
    }
    r += m;
    g += m;
    c += m;
}
void glColorHSV3d(const GLdouble h, const GLdouble s, const GLdouble v){
    GLdouble r, g, b;
    hsvToRgb ( h, s, v, r, g, b);
    glColor3d ( r, g, b);
}

//Draws an arrow whose cone tip has radius equal to base along the z axis,
//from 0 to 1
void glutArrow (GLdouble base){
    //draw from 0 to 1 along the z direction
    glPushMatrix();
        glLineWidth(2.0);
        glBegin(GL_LINES);
            glVertex3d(0.0,0.0,0.0);
            glVertex3d(0.0,0.0,0.8);
        glEnd();
//        gluQuadricNormals(quad, GLU_SMOOTH);
//        gluCylinder(quad, 0.1*base, 0.1*base, 0.8, 6, 1);
        glTranslated(0.0, 0.0, 0.8);
        glutSolidCone(base, 0.2, 6, 3);
    glPopMatrix();
//    gluDeleteQuadric(quad);
}

//Draws an arrow from point a to point b. Base controls the radius of the arrow
//See glutArrow for how base controls the size
void glutDirectedArrow(const Vec3d& a, const Vec3d& b, GLdouble base){
    const Vec3d z =Vec3d(0.0, 0.0, 1.0);
    const Vec3d dir = b-a;
    Scalar eps = 1e-8;

    glPushMatrix();

        glTranslated((GLdouble)a.x(), (GLdouble)a.y(), (GLdouble)a.z());
        Vec3d rotAxis = z.cross(dir);
////        if (rotAxis.isZero(eps)){//this means it is parallel
////            if ( dir.z() < 0){
////                glScaled(0.0, 0.0, -1.0);
////            }
////        } else{
            Scalar theta = angle(z, dir)*180/pi;
            glRotated((GLdouble)theta, (GLdouble)rotAxis.x(), (GLdouble)rotAxis.y(), (GLdouble)rotAxis.z());
////        }
        glScaled(1.0, 1.0, 0.9*(GLdouble)dir.norm());
        glutArrow(base);
    glPopMatrix();
}
void glVertVec3d(const Vec3d& v) {
//    OpenGL::vertex((GLfloat)v.x(), (GLfloat)v.y(), (GLfloat)v.z());
  glVertex3d((GLdouble)v.x(), (GLdouble)v.y(), (GLdouble)v.z());
}
void glNormalVec3d(const Vec3d& v) {
//  OpenGL::normal((GLfloat)v.x(), (GLfloat)v.y(), (GLfloat)v.z());
    glNormal3d((GLdouble)v.x(), (GLdouble)v.y(), (GLdouble)v.z());
}
void drawTri (const Vec3d& a, const Vec3d& b,const Vec3d& c,const Vec3d& n){
  glNormalVec3d(n);
  glVertVec3d(a);
  glVertVec3d(b);
  glVertVec3d(c);
}
void drawStitch (const Vec3d& a, const Vec3d& b,const Vec3d& c,const Vec3d& d){
  Vec3d n = (d-a).cross(b-a);
  n.normalize();

 drawTri(a, b, c, n);
 drawTri(c, d, a, n);
}

void drawThickTri(const Vec3d& v0, const Vec3d& v1, const Vec3d& v2,
        const Scalar & t0, const Scalar & t1, const Scalar & t2,
        const Vec3d& n0, const Vec3d& n1, const Vec3d& n2, const Vec3d& normal) {

  glColor3fv(black);
  Vec3d v0_f = v0 + n0*t0/2.0;
  Vec3d v0_b = v0 - n0*t0/2.0;
  Vec3d v1_f = v1 + n1*t1/2.0;
  Vec3d v1_b = v1 - n1*t1/2.0;
  Vec3d v2_f = v2 + n2*t2/2.0;
  Vec3d v2_b = v2 - n2*t2/2.0;

  //compute effective face normals
  Vec3d f_normal = (v1_f - v0_f).cross(v2_f - v0_f);
  f_normal.normalize();
  Vec3d b_normal = (v1_b - v0_b).cross(v2_b - v0_b);
  b_normal.normalize();

  drawTri (v0_f, v1_f, v2_f, f_normal);
  drawTri (v0_b, v1_b, v2_b, -b_normal);

}

void ShellRenderer::cycleMode() { 
   m_mode = (ShellRenderer::DrawMode) ((m_mode + 1) % 7);
}

ShellRenderer::ShellRenderer( ElasticShell& shell, const Scalar thickness )
: m_shell(shell)
, m_mode(DBG_MULTIPHASE)
, m_refthickness( 2*thickness)
, m_nregion(0)
, m_solid_boundary_visible(true)
{
}


class FaceComp
{
public:
  bool operator() (const std::pair<FaceHandle, Scalar> & f1, const std::pair<FaceHandle, Scalar> & f2) const 
  {
    return f1.second < f2.second;
  }
};
  
bool junction(const DeformableObject & mesh, EdgeHandle e)
{
  int ne = 0;
  for ( EdgeFaceIterator efit = mesh.ef_iter(e); efit; ++efit)
    ne++;
  return ne > 2;
}
  
bool junctionNeighbor(const DeformableObject & mesh, VertexHandle v)
{
  for ( VertexEdgeIterator veit = mesh.ve_iter(v); veit; ++veit)
    if (junction(mesh, *veit))
      return true;
  
  return false;
}

bool junctionNeighbor(const DeformableObject & mesh, EdgeHandle e)
{
  if (junctionNeighbor(mesh, mesh.fromVertex(e)))
    return true;
  
  if (junctionNeighbor(mesh, mesh.toVertex(e)))
    return true;
  
  return false;
}

bool junctionNeighbor(const DeformableObject & mesh, FaceHandle f)
{
  for ( FaceEdgeIterator feit = mesh.fe_iter(f); feit; ++feit)
    if (junctionNeighbor(mesh, *feit))
      return true;

  return false;
}

void ShellRenderer::keyboard(unsigned char key, int x, int y)
{
  if (key == '=' || key == '+')
  {
    m_current_region = (m_current_region + 1) % m_nregion;
    glutPostRedisplay();
  } else if (key == '-' || key == '_')
  {
    m_current_region = (m_current_region + m_nregion - 1) % m_nregion;
    glutPostRedisplay();
  } else if (key == 'l' || key == 'L')
  {
    while ((int)m_region_visible.size() < m_nregion)
      m_region_visible.push_back(true);    
    m_region_visible[m_current_region] = !m_region_visible[m_current_region];
    glutPostRedisplay();
  } else if (key == 'h' || key == 'H')
  {
    while ((int)m_region_visible.size() < m_nregion)
      m_region_visible.push_back(true);
    bool allinvisible = true;
    for (size_t i = 0; i < m_region_visible.size(); i++)
      if (m_region_visible[i])
        allinvisible = false;
    for (size_t i = 0; i < m_region_visible.size(); i++)
      m_region_visible[i] = allinvisible;
    glutPostRedisplay();
  } else if (key == 'b' || key == 'B')
  {
    m_solid_boundary_visible = !m_solid_boundary_visible;
    glutPostRedisplay();
  }
  
}

void ShellRenderer::render()
{
    glPushMatrix();

  if( m_mode == FLAT )
  {
      /*glDisable(GL_LIGHTING);
      MeshCurvature curvature(m_shell.getDefoObj(), m_shell.getVertexPositions());
      curvature.renderCurvatureDirs();*/

    glEnable(GL_LIGHTING);
    //glEnable(GL_COLOR_MATERIAL);

    //Define the hue palette... red for negative thickness

    GLfloat gray[] = {0.8f,0.8f,0.8f,1.0f};
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,gray);

    // Render all faces
    glBegin(GL_TRIANGLES);
    //OpenGL::color(Color(255,0,0));
    DeformableObject& mesh = m_shell.getDefoObj();

    for( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit )
    {
      std::vector<Vec3d> v;
      //Scalar t = m_shell.getThickness(*fit);
      //Scalar hue;
      //Scalar sat;
      //if ( t < 0){//use only red
      //    hue = 0;
      //    sat = fabs(t/m_refthickness - 0.3);

      //    //use the red if t < 0
      //}else{
      //    hue = 240;
      //    sat = t/m_refthickness + 0.3;
      //}
      //glColorHSV3d(hue,sat , 1);
      glColorMaterial (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
      for( FaceVertexIterator fvit = mesh.fv_iter(*fit); fvit; ++fvit )
      {
        v.push_back(m_shell.getVertexPosition(*fvit));
      }

      // Compute a normal for the face
      Vec3d e1 = v[1]-v[0];
      Vec3d e2 = v[2]-v[0];
      Vec3d n = e1.cross(e2);
      if( n.norm() != 0 ) n.normalize();
      drawTri(v[0], v[1], v[2], n);
      //drawThickTri(v[0], v[1], v[2], m_shell.getThickness(*fit));
    }
    glEnd();

    //glDisable(GL_COLOR_MATERIAL);
    glDisable(GL_LIGHTING);
  }
  else if( m_mode == DBG )
  {

//      renderVelocity();

    

    glDisable(GL_LIGHTING);
    DeformableObject& mesh = m_shell.getDefoObj();

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    FaceProperty<Vec3d> faceNormals(&m_shell.getDefoObj());
    m_shell.getFaceNormals(faceNormals);

  

    // Render all edges
    glLineWidth(2);
    glBegin(GL_LINES);
    OpenGL::color(Color(0,0,0));
    
    for( EdgeIterator eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit )
    {
      EdgeHandle eh = *eit;
      Vec3d p0 = m_shell.getVertexPosition(mesh.fromVertex(*eit));
      Vec3d p1 = m_shell.getVertexPosition(mesh.toVertex(*eit));
      Vec3d dir = (p1-p0);
      //p0 = p0 + 0.05*dir;
      //p1 = p1 - 0.05*dir;
      if ( m_shell.shouldFracture(*eit) ){
          OpenGL::color(Color(1.0, 1.0, 0.0));
      } else if (mesh.isBoundary(*eit)){
          OpenGL::color(Color(0.0, 1.0, 0.0));
      }
      else {
        OpenGL::color(Color(0.0,0.0,0.0));
      }
      OpenGL::vertex(p0);
      OpenGL::vertex(p1);
      
      /*
      //Draw edge avg normal vector
      //get adjacent faces
      EdgeFaceIterator efit = m_shell.getDefoObj().ef_iter(eh);
      Vec3d normal1 = faceNormals[*efit]; ++efit;
      Vec3d normal2;
      if(efit)
          normal2 = faceNormals[*efit];
      else
          normal2 = normal1;

      Vec3d avgNormal = 0.5*(normal1+normal2);
      Scalar magnitude = avgNormal.norm();
      avgNormal /= magnitude;
      
      Vec3d midpoint = 0.5*(p0+p1);
      Vec3d endpoint = midpoint + avgNormal*(p1-p0).norm();
      OpenGL::color(Color(0,255,0));
      OpenGL::vertex(midpoint);
      OpenGL::vertex(endpoint);
      
      std::vector<ElasticShellForce*> forces = m_shell.getForces();
      for(unsigned int f = 0; f < forces.size(); ++f) {
          MNBendingForce* ptr_specific = dynamic_cast<MNBendingForce*>(forces[f]);
          if( ptr_specific ){
             Vec3d normalDir = ptr_specific->getEdgeNormal(eh);
             midpoint = 0.5*(p0+p1);
             endpoint = midpoint + normalDir * (p1-p0).norm();
             OpenGL::color(Color(0,0,255));
             OpenGL::vertex(midpoint);
             OpenGL::vertex(endpoint);
          }
      }
      */
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
      int colorVal = (int) (255.0 * thickness/ 0.25); //rescale
      //int colorVal = (int) (255.0 * (thickness - 0.0025) / 0.0025); //test
      //colorVal = clamp(colorVal, 0, 255);
      colorVal = 255;
      OpenGL::color(Color(colorVal,0,0, 128));
      std::vector<Vec3d> points(3);
      int i = 0;
      for( FaceVertexIterator fvit = mesh.fv_iter(*fit); fvit; ++fvit )
      {
        Vec3d pos = m_shell.getVertexPosition(*fvit);
        //pos = pos - 0.02*(pos-barycentre);
        OpenGL::vertex(pos);
        points[i] = pos;
        ++i;
      }      
      
    }
    
    glEnd();
    
    glColor3f(0.0, 0.0, 0.0);
    /*std::cout << "Calling curvature\n";
    MeshCurvature curvature(m_shell.getDefoObj(), m_shell.getVertexPositions());
    curvature.renderCurvatureDirs();
    std::cout << "Done curvature";*/


     //Render all vertices
    glPointSize(4);
    glBegin(GL_POINTS);
    
    for( VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit ) {
      Vec3d vertPos = m_shell.getVertexPosition(*vit); 
      VertexHandle vh = *vit;
      if(!m_shell.getDefoObj().isConstrained(vh)) {
        OpenGL::color(Color(0,0,0));
      }
      else {
        OpenGL::color(Color(0,255,0));
      }

      OpenGL::vertex(vertPos);
    }
    glEnd();
    

    //Draw collision springs
    std::vector<Vec3d> starts, ends;
    m_shell.getSpringList(starts, ends);
    glLineWidth(5);
    
    glBegin(GL_LINES);
    glColor3f(0.0, 1.0, 0.0);
     for(unsigned int i = 0; i < starts.size(); ++i) {
      OpenGL::vertex(starts[i]);
      OpenGL::vertex(ends[i]);
    }
    glEnd();

    //Vec3d spherePos;
    //Scalar sphereRad;
    //m_shell.getCollisionSphere(spherePos, sphereRad);
    //glPointSize(20);
    //glColor3f(1,0,0);
    //glBegin(GL_POINTS);
    //glVertex3f(spherePos[0], spherePos[1], spherePos[2]);
    //glEnd();

    glPointSize(10);
    glBegin(GL_POINTS);
    glColor3f(0,0,1);
    for(unsigned int i = 0; i < starts.size(); ++i) {
      OpenGL::vertex(starts[i]);
    }
    glEnd();
    glPointSize(10);
    
    glBegin(GL_POINTS);
    glColor3f(0,1,1);
    for(unsigned int i = 0; i < ends.size(); ++i) {
      OpenGL::vertex(ends[i]);
    }
    glEnd();



   /* glBegin(GL_QUADS);
    glVertex3f(-2.0f, -0.2, -2.0f);
    glVertex3f(2.0f, -0.2, -2.0f);
    glVertex3f(2.0f, -0.2, 2.0f);
    glVertex3f(-2.0f, -0.2, 2.0f);
    glEnd();*/




    
    glEnable(GL_LIGHTING);

  }
  else if( m_mode == DBG_BUBBLE || m_mode == DBG_JUNCTION )
  {
//      renderVelocity();
    
    glDisable(GL_LIGHTING);
    DeformableObject& mesh = m_shell.getDefoObj();
    
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    
    FaceProperty<Vec3d> faceNormals(&m_shell.getDefoObj());
    m_shell.getFaceNormals(faceNormals);
    
    // Render all edges
    glLineWidth(2);
    glBegin(GL_LINES);
    OpenGL::color(Color(0,0,0));
    
    for( EdgeIterator eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit )
    {
      EdgeHandle eh = *eit;
      
      if (m_mode == DBG_JUNCTION)
        if (!junctionNeighbor(mesh, eh))
          continue;
      
      Vec3d p0 = m_shell.getVertexPosition(mesh.fromVertex(*eit));
      Vec3d p1 = m_shell.getVertexPosition(mesh.toVertex(*eit));
      Vec3d dir = (p1-p0);
      //p0 = p0 + 0.05*dir;
      //p1 = p1 - 0.05*dir;
      if ( m_shell.shouldFracture(*eit) ){
        OpenGL::color(Color(1.0, 1.0, 0.0));
      } else if (mesh.isBoundary(*eit)){
        OpenGL::color(Color(0.0, 1.0, 0.0));
      }
      else {
        OpenGL::color(Color(0.0,0.0,0.0));
      }
      
      if (m_mode == DBG_JUNCTION)
        if (junction(mesh, eh))
          OpenGL::color(Color(0.5, 1.0, 0.0));
      
      OpenGL::vertex(p0);
      OpenGL::vertex(p1);      
    }
    glEnd();
    
    // stats on total number of labels
    int maxlabel = -1;
    for (FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit)
    {
      FaceHandle f = *fit;
      Vec2i regions = m_shell.getFaceLabel(f);
      if (regions.x() >= maxlabel)
        maxlabel = regions.x();
      if (regions.y() >= maxlabel)
        maxlabel = regions.y();
    }
    m_nregion = maxlabel + 1;
    
    // render face labels
    if (true)
    {
      // generate a list of colors for all the labels present
      std::vector<Vec3d> labelcolors;
      labelcolors.push_back(Vec3d(0, 0, 0));
      for (int i = 0; i < maxlabel + 1; i++)
      {
        float r, g, b;
        float t = (float)i / (maxlabel == 0 ? 1 : maxlabel);
        if (t < 0.5)
        {
          r = 1 - t * 1.5;
          g = t * 1.5;
          b = 0;
        } else
        {
          r = 0;
          g = 1 - (t - 0.5) * 1.5;
          b = (t - 0.5) * 1.5;
        }
        labelcolors.push_back(Vec3d(r, g, b));
      }
      
      glBegin(GL_LINES);
      for( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit )
      {
        FaceHandle f = *fit;
        
        if (m_mode == DBG_JUNCTION)
          if (!junctionNeighbor(mesh, f))
            continue;
        
        Vec2i regions = m_shell.getFaceLabel(f);
        FaceVertexIterator fvit = mesh.fv_iter(f); assert(fvit);
        Vec3d p0 = m_shell.getVertexPosition(*fvit);  ++fvit;   assert(fvit);
        Vec3d p1 = m_shell.getVertexPosition(*fvit);  ++fvit;   assert(fvit);
        Vec3d p2 = m_shell.getVertexPosition(*fvit);  ++fvit;   assert(!fvit);
        
        Vec3d c = (p0 + p1 + p2) / 3;
        Vec3d n = (p1 - p0).cross(p2 - p0).normalized();

        Vec3d color0 = labelcolors[std::max(0, regions.x() + 1)];
        glColor3f(color0.x(), color0.y(), color0.z());
        
        if (regions.x() >= 0)
        {
          OpenGL::vertex(c);
          OpenGL::vertex(Vec3d(c - n * 0.05));
        }
        
        Vec3d color1 = labelcolors[std::max(0, regions.y() + 1)];
        glColor3f(color1.x(), color1.y(), color1.z());
                
        if (regions.y() >= 0)
        {
          OpenGL::vertex(c);
          OpenGL::vertex(Vec3d(c + n * 0.05));          
        }
      }
      glEnd();
    }
    
    // render back to front semi transparent
    float mv[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, mv);
    Vec3d view_vec(mv[2], mv[6], mv[10]);  // assuming ModelView matrix contains only translation, rotation and uniform scaling
    
    std::vector<std::pair<FaceHandle, Scalar> > sorted_faces;
    for( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit )
    {
      if (m_mode == DBG_JUNCTION)
        if (!junctionNeighbor(mesh, *fit))
          continue;
      
      Vec3d barycentre;
      for( FaceVertexIterator fvit = mesh.fv_iter(*fit); fvit; ++fvit )
      {
        Vec3d pos = m_shell.getVertexPosition(*fvit);
        barycentre += pos;
      }
      barycentre /= 3.0;
      
      Scalar depth = barycentre.dot(view_vec);
      sorted_faces.push_back(std::pair<FaceHandle, Scalar>(*fit, depth));
    }
    
    FaceComp fc;
    std::sort(sorted_faces.begin(), sorted_faces.end(), fc);
    
    // Render all faces
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);
    glBegin(GL_TRIANGLES);
    
    for (size_t i = 0; i < sorted_faces.size(); i++)
    {
      Vec3d barycentre;
      for( FaceVertexIterator fvit = mesh.fv_iter(sorted_faces[i].first); fvit; ++fvit )
      {
        Vec3d pos = m_shell.getVertexPosition(*fvit);
        barycentre += pos;
      }
      barycentre /= 3.0;

      Scalar alpha = 0.25;
      OpenGL::color(Color(1.0,0.0,0.0,alpha));
      for( FaceVertexIterator fvit = mesh.fv_iter(sorted_faces[i].first); fvit; ++fvit )
      {
        Vec3d pos = m_shell.getVertexPosition(*fvit);
        pos += (barycentre - pos) * 0.1;
        OpenGL::vertex(pos);
      }      
    }
    
    glEnd();
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
    
    glColor3f(0.0, 0.0, 0.0);
    /*std::cout << "Calling curvature\n";
     MeshCurvature curvature(m_shell.getDefoObj(), m_shell.getVertexPositions());
     curvature.renderCurvatureDirs();
     std::cout << "Done curvature";*/
    
    
    //Render all vertices
    glPointSize(4);
    glBegin(GL_POINTS);
    
    for( VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit ) 
    {
      Vec3d vertPos = m_shell.getVertexPosition(*vit); 
      VertexHandle vh = *vit;
      
      if (m_mode == DBG_JUNCTION)
        if (!junctionNeighbor(mesh, vh))
          continue;
      
      if(!m_shell.getDefoObj().isConstrained(vh)) {
        OpenGL::color(Color(0,0,0));
      }
      else {
        OpenGL::color(Color(0,255,0));
      }
      
      OpenGL::vertex(vertPos);
    }
    glEnd();
    
    //Draw collision springs
    std::vector<Vec3d> starts, ends;
    m_shell.getSpringList(starts, ends);
    glLineWidth(5);
    
    glBegin(GL_LINES);
    glColor3f(0.0, 1.0, 0.0);
    for(unsigned int i = 0; i < starts.size(); ++i) {
      OpenGL::vertex(starts[i]);
      OpenGL::vertex(ends[i]);
    }
    glEnd();
    
    //Vec3d spherePos;
    //Scalar sphereRad;
    //m_shell.getCollisionSphere(spherePos, sphereRad);
    //glPointSize(20);
    //glColor3f(1,0,0);
    //glBegin(GL_POINTS);
    //glVertex3f(spherePos[0], spherePos[1], spherePos[2]);
    //glEnd();
    
    glPointSize(10);
    glBegin(GL_POINTS);
    glColor3f(0,0,1);
    for(unsigned int i = 0; i < starts.size(); ++i) {
      OpenGL::vertex(starts[i]);
    }
    glEnd();
    glPointSize(10);
    
    glBegin(GL_POINTS);
    glColor3f(0,1,1);
    for(unsigned int i = 0; i < ends.size(); ++i) {
      OpenGL::vertex(ends[i]);
    }
    glEnd();
    
    glEnable(GL_LIGHTING);
    
  }
  else if( m_mode == DBG_MULTIPHASE )
  {
    glDisable(GL_LIGHTING);
    DeformableObject& mesh = m_shell.getDefoObj();
    
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    //glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    
    FaceProperty<Vec3d> faceNormals(&m_shell.getDefoObj());
    m_shell.getFaceNormals(faceNormals);
    
    // Render all edges
    glLineWidth(1);
    glBegin(GL_LINES);
    OpenGL::color(Color(0,0,0));
    
    for( EdgeIterator eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit )
    {
      EdgeHandle eh = *eit;
      Vec3d p0 = m_shell.getVertexPosition(mesh.fromVertex(*eit));
      Vec3d p1 = m_shell.getVertexPosition(mesh.toVertex(*eit));

      if (!m_solid_boundary_visible)
      {
        bool visible = false;
        for (EdgeFaceIterator efit = mesh.ef_iter(eh); efit; ++efit)
        {
          Vec2i labels = m_shell.getFaceLabel(*efit);
          if (labels.x() >= 0 && labels.y() >= 0)
            visible = true;
        }

        bool xb = ((p0.x() < 1e-4 && p1.x() < 1e-4) || (p0.x() > 1 - 1e-4 && p1.x() > 1 - 1e-4));
        bool yb = ((p0.y() < 1e-4 && p1.y() < 1e-4) || (p0.y() > 1 - 1e-4 && p1.y() > 1 - 1e-4));
        bool zb = ((p0.z() < 1e-4 && p1.z() < 1e-4) || (p0.z() > 1 - 1e-4 && p1.z() > 1 - 1e-4));
        
        bool cubeedge = ((xb && yb) || (xb && zb) || (yb && zb));
        
        if (!visible && !cubeedge)
          continue;
      }
      
      Vec3d dir = (p1-p0);
      //p0 = p0 + 0.05*dir;
      //p1 = p1 - 0.05*dir;
      if ( m_shell.shouldFracture(*eit) ){
        OpenGL::color(Color(1.0, 1.0, 0.0));
      } else if (mesh.isBoundary(*eit)){
        OpenGL::color(Color(0.0, 1.0, 0.0));
      }
      else {
        OpenGL::color(Color(0.0,0.0,0.0));
      }
      
      if (m_mode == DBG_JUNCTION)
        if (junction(mesh, eh))
          OpenGL::color(Color(0.5, 1.0, 0.0));
      
      OpenGL::vertex(p0);
      OpenGL::vertex(p1);      
    }
    glEnd();
    
    // stats on total number of labels
    int maxlabel = -1;
    for (FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit)
    {
      FaceHandle f = *fit;
      Vec2i regions = m_shell.getFaceLabel(f);
      if (regions.x() >= maxlabel)
        maxlabel = regions.x();
      if (regions.y() >= maxlabel)
        maxlabel = regions.y();
    }
    m_nregion = maxlabel + 1;
    
    while (m_region_visible.size() < m_nregion)
      m_region_visible.push_back(true);
    
    // generate a list of colors for all the labels present
    std::vector<Vec3d> labelcolors;
    labelcolors.push_back(Vec3d(0, 0, 0));
    for (int i = 0; i < maxlabel + 1; i++)
    {
      float r, g, b;
      float t = (float)i / (maxlabel == 0 ? 1 : maxlabel);
      if (t < 0.5)
      {
        r = 1 - t * 1.5;
        g = t * 1.5;
        b = 0;
      } else
      {
        r = 0;
        g = 1 - (t - 0.5) * 1.5;
        b = (t - 0.5) * 1.5;
      }
      labelcolors.push_back(Vec3d(r, g, b));
    }

    // render face labels
    if (false)
    {      
      glBegin(GL_LINES);
      for( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit )
      {
        FaceHandle f = *fit;
        
        if (m_mode == DBG_JUNCTION)
          if (!junctionNeighbor(mesh, f))
            continue;
        
        Vec2i regions = m_shell.getFaceLabel(f);
        FaceVertexIterator fvit = mesh.fv_iter(f); assert(fvit);
        Vec3d p0 = m_shell.getVertexPosition(*fvit);  ++fvit;   assert(fvit);
        Vec3d p1 = m_shell.getVertexPosition(*fvit);  ++fvit;   assert(fvit);
        Vec3d p2 = m_shell.getVertexPosition(*fvit);  ++fvit;   assert(!fvit);
        
        Vec3d c = (p0 + p1 + p2) / 3;
        Vec3d n = (p1 - p0).cross(p2 - p0).normalized();
        
        Vec3d color0 = labelcolors[std::max(0, regions.x() + 1)];
        glColor3f(color0.x(), color0.y(), color0.z());
        
        if (regions.x() >= 0)
        {
          OpenGL::vertex(c);
          OpenGL::vertex(Vec3d(c - n * 0.05));
        }
        
        Vec3d color1 = labelcolors[std::max(0, regions.y() + 1)];
        glColor3f(color1.x(), color1.y(), color1.z());
        
        if (regions.y() >= 0)
        {
          OpenGL::vertex(c);
          OpenGL::vertex(Vec3d(c + n * 0.05));          
        }
      }
      glEnd();
    }
    
    // render back to front semi transparent
    float mv[16];
    glGetFloatv(GL_MODELVIEW_MATRIX, mv);
    Vec3d view_vec(mv[2], mv[6], mv[10]);  // assuming ModelView matrix contains only translation, rotation and uniform scaling
    
    for (size_t i = 0; i < m_region_visible.size(); i++)
      std::cout << "region " << i << " " << (m_region_visible[i] ? "visible" : "invisible") << std::endl;
    std::cout << "current region = " << m_current_region << "/" << m_nregion << std::endl;
    
    std::vector<std::pair<FaceHandle, Scalar> > sorted_faces;
    for( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit )
    {
      if (m_mode == DBG_JUNCTION)
        if (!junctionNeighbor(mesh, *fit))
          continue;
      
      Vec3d barycentre;
      for( FaceVertexIterator fvit = mesh.fv_iter(*fit); fvit; ++fvit )
      {
        Vec3d pos = m_shell.getVertexPosition(*fvit);
        barycentre += pos;
      }
      barycentre /= 3.0;
      
      Scalar depth = barycentre.dot(view_vec);
      sorted_faces.push_back(std::pair<FaceHandle, Scalar>(*fit, depth));
    }
    
    FaceComp fc;
    std::sort(sorted_faces.begin(), sorted_faces.end(), fc);
    
    // Render all faces
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);
    glBegin(GL_TRIANGLES);
    
    for (size_t i = 0; i < sorted_faces.size(); i++)
    {
      Vec3d barycentre;
      for( FaceVertexIterator fvit = mesh.fv_iter(sorted_faces[i].first); fvit; ++fvit )
      {
        Vec3d pos = m_shell.getVertexPosition(*fvit);
        barycentre += pos;
      }
      barycentre /= 3.0;
      
      Vec2i regions = m_shell.getFaceLabel(sorted_faces[i].first);
      if (!m_solid_boundary_visible)
        if (regions.x() < 0 || regions.y() < 0)
          continue;
      
      Vec3d color0 = labelcolors[regions.x() + 1];
      Vec3d color1 = labelcolors[regions.y() + 1];
      Vec3d color_combined = Vec3d::Zero();
      Scalar alpha = 0.05;
      int visible_count = 0;
      if (regions.x() >= 0 && m_region_visible[regions.x()])
      {
        color_combined += color0;
        visible_count++;
        alpha += 0.25;
      }
      if (regions.y() >= 0 && m_region_visible[regions.y()])
      {
        color_combined += color1;
        visible_count++;
        alpha += 0.25;
      }
      if (visible_count == 0)
        color_combined = Vec3d(0.0,0.0,0.0);
      else
        color_combined /= visible_count;
      
      OpenGL::color(Color(color_combined.x(), color_combined.y(), color_combined.z() ,alpha));
      for( FaceVertexIterator fvit = mesh.fv_iter(sorted_faces[i].first); fvit; ++fvit )
      {
        Vec3d pos = m_shell.getVertexPosition(*fvit);
//        pos += (barycentre - pos) * 0.1;
        OpenGL::vertex(pos);
      }      
    }
    
    glEnd();
    glDisable(GL_BLEND);
    glDepthMask(GL_TRUE);
    
    glColor3f(0.0, 0.0, 0.0);
    /*std::cout << "Calling curvature\n";
     MeshCurvature curvature(m_shell.getDefoObj(), m_shell.getVertexPositions());
     curvature.renderCurvatureDirs();
     std::cout << "Done curvature";*/
    
    
    //Render all vertices
    glPointSize(4);
    glBegin(GL_POINTS);
    
    for( VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit ) 
    {
      Vec3d vertPos = m_shell.getVertexPosition(*vit); 
      VertexHandle vh = *vit;

      if (!m_solid_boundary_visible)
      {
        bool visible = false;
        for (VertexFaceIterator vfit = mesh.vf_iter(vh); vfit; ++vfit)
        {
          Vec2i labels = m_shell.getFaceLabel(*vfit);
          if (labels.x() >= 0 && labels.y() >= 0)
            visible = true;
        }
        
        if (!visible)
          continue;
      }
            
      if(!m_shell.getDefoObj().isConstrained(vh)) {
        OpenGL::color(Color(0,0,0));
      }
      else {
        OpenGL::color(Color(0,255,0));
      }
      
      OpenGL::vertex(vertPos);
    }
    glEnd();
    
    //Draw collision springs
    std::vector<Vec3d> starts, ends;
    m_shell.getSpringList(starts, ends);
    glLineWidth(5);
    
    glBegin(GL_LINES);
    glColor3f(0.0, 1.0, 0.0);
    for(unsigned int i = 0; i < starts.size(); ++i) {
      OpenGL::vertex(starts[i]);
      OpenGL::vertex(ends[i]);
    }
    glEnd();
    
    //Vec3d spherePos;
    //Scalar sphereRad;
    //m_shell.getCollisionSphere(spherePos, sphereRad);
    //glPointSize(20);
    //glColor3f(1,0,0);
    //glBegin(GL_POINTS);
    //glVertex3f(spherePos[0], spherePos[1], spherePos[2]);
    //glEnd();
    
    glPointSize(10);
    glBegin(GL_POINTS);
    glColor3f(0,0,1);
    for(unsigned int i = 0; i < starts.size(); ++i) {
      OpenGL::vertex(starts[i]);
    }
    glEnd();
    glPointSize(10);
    
    glBegin(GL_POINTS);
    glColor3f(0,1,1);
    for(unsigned int i = 0; i < ends.size(); ++i) {
      OpenGL::vertex(ends[i]);
    }
    glEnd();
    
    glEnable(GL_LIGHTING);
    
  }
/*else if (m_mode == VOLUMETRIC){
//      glDisable(GL_LIGHTING);
      glDisable(GL_LIGHTING);

      DeformableObject& mesh = m_shell.getDefoObj();

      glEnable(GL_LIGHTING);

      glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
      GLfloat blue[] = {0.8f,0.8f,0.8f,1.0f};
      glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE,blue);

      FaceProperty<Vec3d> faceNormals (&m_shell.getDefoObj());
      VertexProperty<Vec3d> vertexNormals(&m_shell.getDefoObj());

      m_shell.getFaceNormals(faceNormals);
      m_shell.getVertexNormals(vertexNormals);

         // Render all faces
         glBegin(GL_TRIANGLES);
         for( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit )
         {
           Vec3dArray vs;
           Vec3dArray ns;
           ScalarArray ts;

           for( FaceVertexIterator fvit = mesh.fv_iter(*fit); fvit; ++fvit)
           {
             vs.push_back(m_shell.getVertexPosition(*fvit));
             ns.push_back(vertexNormals[*fvit]);
             ts.push_back(m_shell.getThickness(*fvit));
           }

           for (FaceEdgeIterator feit = mesh.fe_iter(*fit); feit; ++feit){
               if (mesh.edgeIncidentFaces(*feit)==1){ //it is a boundary
                   Vec3d n = faceNormals[*fit];

                   int orient = mesh.getRelativeOrientation(*fit, *feit);

                   Vec3d from, to;
                   Scalar thickFrom, thickTo;
                   if (orient == 1){
                       from = m_shell.getVertexPosition(mesh.fromVertex(*feit));
                       to = m_shell.getVertexPosition(mesh.toVertex(*feit));
                       thickFrom = m_shell.getThickness(mesh.fromVertex(*feit));
                       thickTo = m_shell.getThickness(mesh.toVertex(*feit));
                   } else {
                       to = m_shell.getVertexPosition(mesh.fromVertex(*feit));
                       from  = m_shell.getVertexPosition(mesh.toVertex(*feit));
                       thickTo = m_shell.getThickness(mesh.fromVertex(*feit));
                       thickFrom = m_shell.getThickness(mesh.toVertex(*feit));
                   }

                   drawStitch(from - n * thickFrom / 2.0, from + n * thickFrom / 2.0,
                           to + n * thickTo / 2.0, to - n * thickTo/ 2.0);
               }
           }

           drawThickTri(vs[0], vs[1], vs[2],
                   ts[0], ts[1], ts[2],
                   ns[0], ns[1], ns[2],
                   faceNormals[*fit]);

         }
         //         for( EdgeIterator eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit )
         //         {
         //           Vec3d p0 = m_shell.getVertexPosition(mesh.fromVertex(*eit));
         //           Vec3d p1 = m_shell.getVertexPosition(mesh.toVertex(*eit));
         //           glutDirectedArrow(p0, p1, (p1-p0).norm()*0.1);
         //
         //
         //         }
         glEnd();

         glDisable(GL_LIGHTING);

  }*/
  else if (m_mode == VOLUMETRIC){

    DeformableObject& mesh = m_shell.getDefoObj();

    glEnable(GL_LIGHTING);

    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    GLfloat blue[] = {0.8f,0.8f,0.8f,1.0f};
    glMaterialfv(GL_FRONT_AND_BACK,GL_AMBIENT_AND_DIFFUSE, blue);

    FaceProperty<Vec3d> faceNormals (&m_shell.getDefoObj());
    VertexProperty<Vec3d> vertexNormals(&m_shell.getDefoObj());

    m_shell.getFaceNormals(faceNormals);
    m_shell.getVertexNormals(vertexNormals);

    // Render all faces
    glBegin(GL_TRIANGLES);
    for( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit )
    {
      Vec3dArray vs;
      Vec3dArray ns;
      ScalarArray ts;

      for( FaceVertexIterator fvit = mesh.fv_iter(*fit); fvit; ++fvit)
      {
        vs.push_back(m_shell.getVertexPosition(*fvit));
        ns.push_back(vertexNormals[*fvit]);
        
        //determine if this is an edge vertex
        bool insideVert = true;
        for(VertexEdgeIterator veit = mesh.ve_iter(*fvit); veit; ++veit) {
          EdgeHandle eh = *veit;
          if(mesh.edgeIncidentFaces(eh) != 2)
            insideVert = false;
        }
        if(insideVert)
          ts.push_back(m_shell.getThickness(*fvit));
        else
          ts.push_back(0);

      }

      /*
      for (FaceEdgeIterator feit = mesh.fe_iter(*fit); feit; ++feit){
        
        if (mesh.edgeIncidentFaces(*feit)==1){ //it is a boundary
          Vec3d n = faceNormals[*fit];

          int orient = mesh.getRelativeOrientation(*fit, *feit);

          Vec3d from, to;
          Scalar thickFrom, thickTo;
          if (orient == 1){
            from = m_shell.getVertexPosition(mesh.fromVertex(*feit));
            to = m_shell.getVertexPosition(mesh.toVertex(*feit));
            thickFrom = m_shell.getThickness(mesh.fromVertex(*feit));
            thickTo = m_shell.getThickness(mesh.toVertex(*feit));
          } else {
            to = m_shell.getVertexPosition(mesh.fromVertex(*feit));
            from  = m_shell.getVertexPosition(mesh.toVertex(*feit));
            thickTo = m_shell.getThickness(mesh.fromVertex(*feit));
            thickFrom = m_shell.getThickness(mesh.toVertex(*feit));
          }

          drawStitch(from - n * thickFrom / 2.0, from + n * thickFrom / 2.0,
            to + n * thickTo / 2.0, to - n * thickTo/ 2.0);
        }
      }
      */
      
     /* std::cout << "Face thickness: " << m_shell.getThickness(*fit) << std::endl;
      std::cout << "Face area: " << m_shell.getArea(*fit) << std::endl;
      std::cout << "Side length: " << (vs[1] - vs[0]).norm() << std::endl;*/

      drawThickTri(vs[0], vs[1], vs[2],
        ts[0], ts[1], ts[2],
        ns[0], ns[1], ns[2],
        faceNormals[*fit]);

    }
    //         for( EdgeIterator eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit )
    //         {
    //           Vec3d p0 = m_shell.getVertexPosition(mesh.fromVertex(*eit));
    //           Vec3d p1 = m_shell.getVertexPosition(mesh.toVertex(*eit));
    //           glutDirectedArrow(p0, p1, (p1-p0).norm()*0.1);
    //
    //
    //         }
    glEnd();

    glDisable(GL_LIGHTING);

  }
  glPopMatrix();

}
void ShellRenderer::renderVelocity(){
    //Draw velocity field as arrows
    const DeformableObject& mesh = m_shell.getDefoObj();
    glPushMatrix();
    for ( VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit){
        Vec3d p0 = m_shell.getVertexPosition(*vit);
        Vec3d v0 = m_shell.getVertexVelocity(*vit);
        glutDirectedArrow(p0, p0 + v0, v0.norm()*0.1);
    }
    glPopMatrix();
}
void ShellRenderer::renderEdges(){
//    draw edges as arrows!
    const DeformableObject& mesh = m_shell.getDefoObj();
     glPushMatrix();
     for( EdgeIterator eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit )
     {
       Vec3d p0 = m_shell.getVertexPosition(mesh.fromVertex(*eit));
       Vec3d p1 = m_shell.getVertexPosition(mesh.toVertex(*eit));
       glutDirectedArrow(p0, p1, (p1-p0).norm()*0.1);

     }
     glPopMatrix();
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

