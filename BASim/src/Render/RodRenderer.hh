/**
 * \file RodRenderer.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/30/2009
 */

#ifndef RODRENDERER_HH
#define RODRENDERER_HH

#include "../Core/EigenIncludes.hh"
#include "RenderBase.hh"
#include "../Physics/ElasticRods/ElasticRod.hh"
#include "../Physics/ElasticRods/RodTube.hh"
#include "Color.hh"
#include "OpenGLDecl.hh"

namespace BASim {

/** Class that implements OpenGL rendering for rods. */
class RodRenderer : public RenderBase
{
public:

  enum DrawMode { SIMPLE, SMOOTH, NONE };

  explicit RodRenderer(ElasticRod& rod);

  void render();

  DrawMode getMode() const { return m_mode; }
  void setMode(DrawMode mode) { m_mode = mode; }

  bool& drawRod() { return m_drawRod; }
  bool& drawMaterial() { return m_drawMaterial; }
  bool& drawRootMaterial() { return m_drawRootMaterial; }
  bool& drawReference() { return m_drawReference; }
  bool& scaleToRadius() { return m_scaleToRadius; }
  bool& drawArrows() { return m_drawArrows; }
  bool& drawVelocity() { return m_drawVelocity; }
  bool& drawResponse() { return m_drawResponse; }

  virtual Vec3d calculateObjectCenter();
  virtual Scalar calculateObjectBoundingRadius(const Vec3d& center);

  void drawSmoothPartialRod( const int i_startVertex, const int i_endVertex, const Vec3d i_color );

  void setColorInSimpleMode( const Color &i_root, const Color &i_tip)
  {
      m_simpleRod.clear();
      m_simpleRod.push_back( i_root );
      m_simpleRod.push_back( i_tip );
  }

protected:

  void drawSimpleRod();
  void drawSmoothRod();

  void drawMaterialFrame();
  void drawRootMaterialFrame();
  void drawReferenceFrame();

  void drawVelocityVector();
  void drawResponseVector();

  ElasticRod& m_rod;
  RodTube m_tube;

  DrawMode m_mode;

  bool m_drawRod;
  bool m_drawMaterial;
  bool m_drawRootMaterial;
  bool m_drawReference;
  bool m_scaleToRadius;
  bool m_drawArrows;
  bool m_drawVelocity;
  bool m_drawResponse;

  std::vector<Color> m_palette;
  std::vector<Color> m_simpleRod;
};

} // namespace BASim

#endif // RODRENDERER_HH
