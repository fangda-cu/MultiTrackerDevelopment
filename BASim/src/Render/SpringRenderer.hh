/**
 * \file SpringRenderer.hh
 *
 * \author smith@cs.columbia.edu
 * \date 07/06/2010
 */

#ifndef SPRINGRENDERER_HH
#define SPRINGRENDERER_HH

#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Render/OpenGLHeaders.hh"

#include "BASim/src/Render/RenderBase.hh"

#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"
#include "BASim/src/Physics/ElasticRods/RodRodSpringForce.hh"

namespace BASim {

class SpringRenderer : public RenderBase
{
public:

  SpringRenderer( const RodRodSpringForce& force );

  virtual ~SpringRenderer();

  virtual void render();

  virtual Vec3d calculateObjectCenter();
  virtual Scalar calculateObjectBoundingRadius(const Vec3d& center);
  void cycleMode();

protected:

  const RodRodSpringForce& m_spring;

};

} // namespace BASim

#endif // SPRINGRENDERER_HH
