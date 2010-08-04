/**
 * \file AnisotropicRod.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#include "AnisotropicRod.hh"
#include "RodStretchingForce.hh"
#include "RodTwistingForce.hh"
#include "RodTwistingForceSym.hh"
#include "RodBendingForce.hh"
#include "RodBendingForceSym.hh"
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

  // Add elastic forces
  addForce(new RodStretchingForce(*this));
  addForce(new RodTwistingForceSym(*this));
  if (refFrameType() == TimeParallel) addForce(new RodBendingForceSym(*this));
  else addForce(new RodAnisoForce(*this));

  // Add viscous forces
  addForce(new RodStretchingForce(*this,true));
  addForce(new RodBendingForceSym(*this,true));
  addForce(new RodTwistingForceSym(*this,true));
}

} // namespace BASim
