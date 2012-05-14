//
//  RodModelRenderer.hh
//  BASim
//
//  Created by Fang Da (fang@cs.columbia.edu) on 5/13/12.
//  Copyright (c) 2012 Columbia University. All rights reserved.
//

#ifndef RODMODELRENDERER_HH
#define RODMODELRENDERER_HH

#include "BASim/src/Physics/DeformableObjects/Rods/ElasticRodModel.hh"
#include "RenderBase.hh"
#include "Color.hh"
#include "OpenGLDecl.hh"

namespace BASim 
{  
  /** Class that implements OpenGL rendering for ElasticRodModel (code mostly adapted 
   *  from RodRenderer, but adopting ShellRenderer's draw modes in order to keep in
   *  sync with it when used together). */
  class RodModelRenderer : public RenderBase
  {
  public:    
    enum DrawMode { NONE, DBG, FLAT, VOLUMETRIC };
        
  public:
    explicit RodModelRenderer(ElasticRodModel & rod);
    
    void render();
    
    DrawMode getMode() const { return m_mode; }
    void setMode(DrawMode mode) { m_mode = mode; }
    
    void cycleMode();
    
    bool & drawRod() { return m_drawRod; }
    bool & drawMaterial() { return m_drawMaterial; }
    bool & drawRootMaterial() { return m_drawRootMaterial; }
    bool & drawReference() { return m_drawReference; }
    bool & drawArrows() { return m_drawArrows; }
    bool & drawVelocity() { return m_drawVelocity; }
    
    virtual Vec3d calculateObjectCenter();
    virtual Scalar calculateObjectBoundingRadius(const Vec3d& center);
    
  protected:
    void drawDebugRod();
    void drawSmoothRod();
    
    void drawMaterialFrame() { assert(!"Not implemented"); }
    void drawRootMaterialFrame() { assert(!"Not implemented"); }
    void drawReferenceFrame() { assert(!"Not implemented"); }
    void drawVelocityVector() { assert(!"Not implemented"); }

  protected:
    ElasticRodModel & m_rod;
    DrawMode m_mode;
    
    bool m_drawRod;
    bool m_drawMaterial;
    bool m_drawRootMaterial;
    bool m_drawReference;
    bool m_drawArrows;
    bool m_drawVelocity;
    
    std::vector<Color> m_palette;
    std::vector<Color> m_simpleRod;
  };
  
} // namespace BASim


#endif
