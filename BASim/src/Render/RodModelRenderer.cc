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
    m_mode(SIMPLE), 
    m_drawRod(true), 
    m_drawMaterial(false),
    m_drawRootMaterial(true), 
    m_drawReference(false), 
    m_scaleToRadius(true),
    m_drawArrows(false), 
    m_drawVelocity(false)
  {
    // material colors
    m_palette.push_back( Color( 255, 125, 75 ) );
    m_palette.push_back( Color( 55, 125, 255 ) );
    
    // reference frame colors
    m_palette.push_back( Color( 200, 200, 200 ) );
    m_palette.push_back( Color( 50, 50, 50 ) );
    
    // curvature binormal color
    
    // velocity color
    m_palette.push_back( Color( 255, 0, 0 ) );
    
    //root and tip color for simple mode
    m_simpleRod.push_back( Color( 0, 0, 0 ) );
    m_simpleRod.push_back( Color( 155, 200, 100 ) );
  }
  
  void RodModelRenderer::cycleMode() 
  { 
    m_mode = (RodModelRenderer::DrawMode) ((m_mode + 1) % 3); 
  }
  
  void RodModelRenderer::render()
  {
    if ( m_drawRod )
    {
      if ( m_mode == SMOOTH )
        drawSmoothRod();
      else if ( m_mode == SIMPLE )
        drawSimpleRod();
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
  
  void RodModelRenderer::drawSimpleRod()
  {
    glEnable( GL_COLOR_MATERIAL);
    glColorMaterial( GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE );
    
    glBegin( GL_LINES);
    
    DeformableObject & obj = m_rod.getDefoObj();
    for (EdgeIterator eit = obj.edges_begin(); eit != obj.edges_end(); ++eit)
    {
      EdgeVertexIterator evit = m_rod.getDefoObj().ev_iter( *eit );
      
      do
      {
        Vec3d x = m_rod.getDefoObj().getVertexPosition( *evit );
        OpenGL::color( m_simpleRod[0] );
        OpenGL::vertex( x );
        ++evit;
        
      } while (evit);
    }    
    
    glEnd();
    glDisable( GL_COLOR_MATERIAL );
  }
  
  void RodModelRenderer::drawSmoothRod()
  {
    assert(!"Not implemented");
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
