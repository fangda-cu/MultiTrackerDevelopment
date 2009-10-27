/**
 * \file RodRenderer.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/30/2009
 */

#ifndef RODRENDERER_HH
#define RODRENDERER_HH

namespace BASim {

/** Class that implements OpenGL rendering for rods. */
class RodRenderer : public RenderBase
{
public:

  enum DrawMode { SIMPLE, SMOOTH, NONE };

  RodRenderer(ElasticRod& rod);

  void render();

  DrawMode getMode() const { return m_mode; }
  void setMode(DrawMode mode) { m_mode = mode; }

  bool& drawMaterial() { return m_drawMaterial; }
  bool& drawReference() { return m_drawReference; }
  bool& scaleToRadius() { return m_scaleToRadius; }
  bool& drawArrows() { return m_drawArrows; }

protected:

  void drawSimpleRod();
  void drawSmoothRod();

  void drawMaterialFrame();
  void drawReferenceFrame();

  ElasticRod& m_rod;
  RodTube m_tube;

  DrawMode m_mode;

  bool m_drawMaterial;
  bool m_drawReference;
  bool m_scaleToRadius;
  bool m_drawArrows;

  std::vector<Color> m_palette;
};

} // namespace BASim

#endif // RODRENDERER_HH
