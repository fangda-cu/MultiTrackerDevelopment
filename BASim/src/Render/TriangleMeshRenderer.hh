/**
 * \file TriangleMeshRenderer.hh
 *
 * \author smith@cs.columbia.edu
 * \date 03/15/2010
 */

#ifndef TRIANGLEMESHRENDERER_HH
#define TRIANGLEMESHRENDERER_HH

#ifdef WETA
#include "../Core/TriangleMesh.hh"
#include "../Render/RenderBase.hh"
#else
#include "BASim/src/Core/TriangleMesh.hh"
#include "BASim/src/Render/RenderBase.hh"
#endif

namespace BASim {

  /** Class that implements OpenGL rendering for triangle meshes. */
  class TriangleMeshRenderer : public RenderBase
  {
  public:
  
    enum DrawMode { DBG, FLAT, NONE };

    explicit TriangleMeshRenderer( const TriangleMesh& mesh );
    
    void render();
    
    DrawMode getMode() const { return m_mode; }
    void setMode(DrawMode mode) { m_mode = mode; }
    
    virtual Vec3d calculateObjectCenter();
    virtual Scalar calculateObjectBoundingRadius(const Vec3d& center);
    
  protected:
    const TriangleMesh& m_mesh;
    DrawMode m_mode;
  };
  
} // namespace BASim

#endif // TRIANGLEMESHRENDERER_HH
