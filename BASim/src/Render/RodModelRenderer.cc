//
//  RodModelRenderer.cc
//  BASim
//
//  Created by Fang Da (fang@cs.columbia.edu) on 5/13/12.
//  Copyright (c) 2012 Columbia University. All rights reserved.
//


#include "RodModelRenderer.hh"
#include "RenderUtils.hh"

namespace BASim
{
  RodModelRenderer::RodModelRenderer(ElasticRodModel & rod ) :
    m_rod(rod), 
    m_mode(DBG), 
    m_drawRod(true), 
    m_drawMaterial(false),
    m_drawRootMaterial(false), 
    m_drawReference(false), 
    m_drawArrows(false), 
    m_drawVelocity(false)
  {
    // material colors
    m_palette.push_back( Color( 255, 125, 75 ) );
    m_palette.push_back( Color( 55, 125, 255 ) );
    
    // reference frame colors
    m_palette.push_back( Color( 200, 200, 200 ) );
    m_palette.push_back( Color( 50, 50, 50 ) );
    
    //rod color for smooth mode
    m_simpleRod.push_back( Color( 50, 50, 255 ) );
  }
  
  void RodModelRenderer::cycleMode() 
  { 
    m_mode = (RodModelRenderer::DrawMode) ((m_mode + 1) % 4); 
  }
  
  void RodModelRenderer::render()
  {
    if ( m_drawRod )
    {
      if ( m_mode == VOLUMETRIC )
        drawSmoothRod();
      else if ( m_mode == FLAT || m_mode == DBG )
        drawDebugRod();
    }
    
    if ( m_drawMaterial )
      drawMaterialFrame();
    if ( m_drawRootMaterial )
      drawRootMaterialFrame();
    if ( m_drawReference )
      drawReferenceFrame();    
    if ( m_drawVelocity )
      drawVelocityVector();
  }
  
  void RodModelRenderer::drawDebugRod()
  {    
    glDisable(GL_LIGHTING);
    
    DeformableObject & obj = m_rod.getDefoObj();
    glLineWidth(5);
    glBegin(GL_LINES);
    OpenGL::color(m_simpleRod[0]);
    for (EdgeIterator eit = obj.edges_begin(); eit != obj.edges_end(); ++eit)
    {
      if (m_rod.isEdgeActive(*eit))
      {
        EdgeVertexIterator evit = m_rod.getDefoObj().ev_iter( *eit );
        Vec3d v1 = m_rod.getDefoObj().getVertexPosition(*evit); ++evit;
        Vec3d v2 = m_rod.getDefoObj().getVertexPosition(*evit); ++evit;
        
        OpenGL::vertex(v1);
        OpenGL::vertex(v2);      
      }
    }        
    glEnd();

    /*
    glPointSize(10);
    glBegin(GL_POINTS);
    OpenGL::color(Color(255, 0, 0));
    for (EdgeIterator eit = obj.edges_begin(); eit != obj.edges_end(); ++eit)
    {
      if (m_rod.isEdgeActive(*eit))
      {
        EdgeVertexIterator evit = m_rod.getDefoObj().ev_iter( *eit );
        Vec3d v1 = m_rod.getDefoObj().getVertexPosition(*evit); ++evit;
        Vec3d v2 = m_rod.getDefoObj().getVertexPosition(*evit); ++evit;

        OpenGL::vertex(v1);
        OpenGL::vertex(v2);      
      }
    }        
    glEnd();
    */

    // render material frames
    glLineWidth(1);
    glBegin(GL_LINES);
    for (EdgeIterator eit = obj.edges_begin(); eit != obj.edges_end(); ++eit)
    {
      if (m_rod.isEdgeActive(*eit))
      {
        EdgeVertexIterator evit = m_rod.getDefoObj().ev_iter( *eit );
        Vec3d v1 = m_rod.getDefoObj().getVertexPosition(*evit); ++evit;
        Vec3d v2 = m_rod.getDefoObj().getVertexPosition(*evit); ++evit;
        Vec3d vcenter = (v1 + v2) / 2;
        Vec3d md1 = m_rod.getMaterialDirector1(*eit);
        Vec3d md2 = m_rod.getMaterialDirector2(*eit);
        Vec2d radii = m_rod.getRadii(*eit);
        Scalar r = (radii.x() + radii.y()) / 2;
        OpenGL::color(m_palette[0]);
        OpenGL::vertex(vcenter);
        OpenGL::vertex(Vec3d(vcenter + md1 * r));
        OpenGL::color(m_palette[1]);
        OpenGL::vertex(vcenter);
        OpenGL::vertex(Vec3d(vcenter + md2 * r));
      }
    }    
    glEnd();
  }
  
  void RodModelRenderer::drawSmoothRod()
  {
    const int N = 6;
    DeformableObject & obj = m_rod.getDefoObj();
    
    // precomputation: average edge material frames to vertices. averaged material frame directors
    //  are also scaled by radii already
    VertexProperty<Vec3d> vmd1(&obj);
    VertexProperty<Vec3d> vmd2(&obj);
    VertexProperty<int> edge_count(&obj);
    vmd1.assign(Vec3d(0, 0, 0));
    vmd2.assign(Vec3d(0, 0, 0));
    edge_count.assign(0);
    for (EdgeIterator eit = obj.edges_begin(); eit != obj.edges_end(); ++eit)
    {
      if (m_rod.isEdgeActive(*eit))
      {
        EdgeVertexIterator evit = m_rod.getDefoObj().ev_iter(*eit);
        VertexHandle v1 = *evit; ++evit;
        VertexHandle v2 = *evit; ++evit;
        Vec2d radii = m_rod.getRadii(*eit);
        Vec3d emd1 = m_rod.getMaterialDirector1(*eit) * radii(0);
        Vec3d emd2 = m_rod.getMaterialDirector2(*eit) * radii(1);
        vmd1[v1] += emd1;
        vmd2[v1] += emd2;
        vmd1[v2] += emd1;
        vmd2[v2] += emd2;
        edge_count[v1]++;
        edge_count[v2]++;
      }
    }
    
    // render all rod edges, extruding radially from each vertex by averaged material frame directors
    glEnable(GL_COLOR_MATERIAL);
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
    glShadeModel(GL_SMOOTH);

    glEnable(GL_LIGHTING);
    
    glBegin(GL_TRIANGLES);
    OpenGL::color(m_simpleRod[0]);
    
    for (EdgeIterator eit = obj.edges_begin(); eit != obj.edges_end(); ++eit)
    {
      if (m_rod.isEdgeActive(*eit))
      {
        EdgeVertexIterator evit = m_rod.getDefoObj().ev_iter(*eit);
        VertexHandle v1 = *evit; ++evit;
        VertexHandle v2 = *evit; ++evit;
        Vec3d pos1 = m_rod.getDefoObj().getVertexPosition(v1);
        Vec3d pos2 = m_rod.getDefoObj().getVertexPosition(v2);

        Vec3d n11 = vmd1[v1] / edge_count[v1];  // v1 is incident to an active count so edge_count[v1] can't be zero
        Vec3d n12 = vmd2[v1] / edge_count[v1];  // v1 is incident to an active count so edge_count[v1] can't be zero
        Vec3d n21 = vmd1[v2] / edge_count[v2];
        Vec3d n22 = vmd2[v2] / edge_count[v2];

        //Assuming edge radius drops to zero at end points, e.g. for liquid surface tension.
        if(edge_count[v1] == 1) { n11 *= 0; n12 *= 0; }
        if(edge_count[v2] == 1) { n21 *= 0; n22 *= 0; }

        // render a (possibly twisted) prism
        for (int i = 0; i < N; i++)
        {
          Scalar thistheta = (i) * 2 * M_PI / N;
          Scalar nexttheta = (i + 1) * 2 * M_PI / N;
          Scalar c1 = cos(thistheta);
          Scalar s1 = sin(thistheta);
          Scalar c2 = cos(nexttheta);
          Scalar s2 = sin(nexttheta);

          Vec3d pn11 = n11 * c1 + n12 * s1;
          Vec3d pn12 = n21 * c1 + n22 * s1;
          Vec3d pn21 = n11 * c2 + n12 * s2;
          Vec3d pn22 = n21 * c2 + n22 * s2;
          
          Vec3d p11 = pos1 + pn11;
          Vec3d p12 = pos2 + pn12;
          Vec3d p21 = pos1 + pn21;
          Vec3d p22 = pos2 + pn22;
          
          pn11.normalize();
          pn12.normalize();
          pn21.normalize();
          pn22.normalize();
          
          OpenGL::normal(pn11);
          OpenGL::vertex(p11);
          OpenGL::normal(pn12);
          OpenGL::vertex(p12);
          OpenGL::normal(pn22);
          OpenGL::vertex(p22);
          
          OpenGL::normal(pn22);
          OpenGL::vertex(p22);
          OpenGL::normal(pn21);
          OpenGL::vertex(p21);
          OpenGL::normal(pn11);
          OpenGL::vertex(p11);
        }
      }
    }    
    
    glEnd();

    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
  }
  
  Vec3d RodModelRenderer::calculateObjectCenter()
  {
    Vec3d center = Vec3d::Zero();
    
    for (VertexIterator vit = m_rod.getDefoObj().vertices_begin(); vit != m_rod.getDefoObj().vertices_end(); ++vit)
    {
      center += m_rod.getDefoObj().getVertexPosition(*vit);
    }
    
    center /= m_rod.getDefoObj().nv();
    
    return center;
  }
  
  Scalar RodModelRenderer::calculateObjectBoundingRadius( const Vec3d& center )
  {
    Scalar radius = 0.0;
    
    for (VertexIterator vit = m_rod.getDefoObj().vertices_begin(); vit != m_rod.getDefoObj().vertices_end(); ++vit)
    {
      radius = std::max( radius, ( m_rod.getDefoObj().getVertexPosition(*vit) - center ).norm() );
    }
    
    return radius;
  }
}
