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
  
    enum DrawMode { DBG, FLAT, NONE, VOLUMETRIC };

    ShellRenderer( const ElasticShell& shell, const Scalar thickness = 1.0 );
    
    void render();
    void renderEdges();
    void cycleMode();
    DrawMode getMode() const { return m_mode; }
    void setMode(DrawMode mode) { m_mode = mode; }
    
    virtual Vec3d calculateObjectCenter();
    virtual Scalar calculateObjectBoundingRadius(const Vec3d& center);
    
  protected:
    const ElasticShell& m_shell;
    DrawMode m_mode;
    const Scalar m_refthickness;
  };
  
} // namespace BASim

#endif //DEFOOBJRENDERER_H
