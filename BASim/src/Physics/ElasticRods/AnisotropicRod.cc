/**
 * \file AnisotropicRod.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#include "AnisotropicRod.hh"
#include "RodStretchingForce.hh"
#include "RodTwistingForce.hh"
#include "RodBendingForce.hh"
#include "RodAnisoForce.hh"

namespace BASim {

AnisotropicRod::AnisotropicRod(int numVertices)
  : ElasticRod(numVertices)
{
  setRefFrameType(TimeParallel);
}

void AnisotropicRod::setup()
{
  ElasticRod::setup();

  addForce(new RodStretchingForce(*this));
  addForce(new RodTwistingForce(*this));
  if (refFrameType() == TimeParallel) addForce(new RodBendingForce(*this));
  else addForce(new RodAnisoForce(*this));
}

} // namespace BASim
