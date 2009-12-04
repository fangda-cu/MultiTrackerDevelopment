/**
 * \file RodUtils.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/03/2009
 */

#ifndef RODUTILS_HH
#define RODUTILS_HH

namespace BASim {

struct RodOptions {
  int numVertices;
  Scalar density;
  Scalar radiusA;
  Scalar radiusB;
  Scalar YoungsModulus;
  Scalar ShearModulus;
  bool anisotropic;
  bool elastic;
  bool viscous;
  bool quasistatic;
  bool inextensible;
  ElasticRod::RefFrameType refFrame;

  RodOptions()
    : numVertices(3)
    , density(1)
    , radiusA(0.1)
    , radiusB(0.1)
    , YoungsModulus(1000)
    , ShearModulus(1000)
    , anisotropic(true)
    , elastic(true)
    , viscous(false)
    , quasistatic(false)
    , inextensible(false)
    , refFrame(ElasticRod::TimeParallel)
  {}
};

inline ElasticRod* setupRod(const RodOptions& opts,
                            const std::vector<Vec3d>& initialPosition,
                            const std::vector<Vec3d>& undeformedPosition)
{
  assert(opts.numVertices == (int) initialPosition.size());
  assert(opts.numVertices == (int) undeformedPosition.size());

  ElasticRod* rod = NULL;
  if (opts.anisotropic) rod = new AnisotropicRod(opts.numVertices);
  else rod = new AnisotropicRod(opts.numVertices);

  rod->setRadius(opts.radiusA, opts.radiusB);
  rod->setDensity(opts.density);
  rod->setYoungsModulus(opts.YoungsModulus);
  rod->setShearModulus(opts.ShearModulus);
  rod->setQuasistatic(opts.quasistatic);
  rod->setRefFrameType(opts.refFrame);

  // set up using undeformed positions
  for (int i = 0; i < rod->nv(); ++i)
    rod->setVertex(i, undeformedPosition[i]);
  rod->setup();

  // update to initial positions
  for (int i = 0; i < rod->nv(); ++i)
  {
    rod->setVertex(i, initialPosition[i]);

    // Setup initial collision structures storing previous 
    // positions of the rod
    Vec3d v = rod->getVertex(i);
    rod->getStartPositions()[i] = v;
    rod->getEndPositions()[i] = v;
    rod->getVelocities()[i] = Vec3d(0.0,0.0,0.0);
  }
  rod->updateProperties();

  return rod;
}

inline Vec3d calculateObjectCenter(const ElasticRod& rod)
{
  Vec3d center = Vec3d::Zero();

  ElasticRod::vertex_iter vit, end = rod.vertices_end();
  for (vit = rod.vertices_begin(); vit != end; ++vit) {
    center += rod.getVertex(*vit);
  }

  center /= rod.nv();

  return center;
}

inline Scalar calculateObjectBoundingRadius(const ElasticRod& rod,
                                            const Vec3d& center)
{
  Scalar radius = 0.0;

  ElasticRod::vertex_iter vit, end = rod.vertices_end();
  for (vit = rod.vertices_begin(); vit != end; ++vit) {
    radius = std::max(radius, (rod.getVertex(*vit) - center).norm());
  }

  return radius;
}

} // namespace BASim

#endif // RODUTILS_HH
