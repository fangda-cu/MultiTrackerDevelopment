/**
 * \file DefoObjRenderer.hh
 *
 * \author batty@cs.columbia.edu
 * \date April 20, 2011
 */

#ifndef SHELLRENDERER_H
#define SHELLRENDERER_H

#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShell.hh"
#include "BASim/src/Render/RenderBase.hh"

namespace BASim {

  /** Class that implements OpenGL rendering for elastic shells. Mimics Breannan's triangle mesh code. */
  class ShellRenderer : public RenderBase
  {
  public:
  
    enum DrawMode { NONE, DBG, DBG_BUBBLE, DBG_JUNCTION, DBG_MULTIPHASE, FLAT, VOLUMETRIC };

    ShellRenderer( ElasticShell& shell, const Scalar thickness = 1.0 );
    
    void render();
    void renderEdges();
    void renderVelocity();
    void cycleMode();
    DrawMode getMode() const { return m_mode; }
    void setMode(DrawMode mode) { m_mode = mode; }
    
    virtual Vec3d calculateObjectCenter();
    virtual Scalar calculateObjectBoundingRadius(const Vec3d& center);
    
    void keyboard(unsigned char key, int x, int y);
    
  protected:
    ElasticShell& m_shell;
    DrawMode m_mode;
    const Scalar m_refthickness;
    int m_current_region;
    int m_nregion;
    std::vector<bool> m_region_visible;
    
  };
  
} // namespace BASim

#endif //DEFOOBJRENDERER_H
