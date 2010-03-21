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

  addForce(new RodStretchingForce(*this));
  //addForce(new RodTwistingForce(*this));
  addForce(new RodTwistingForceSym(*this));
  if (refFrameType() == TimeParallel) {
    //addForce(new RodBendingForce(*this));
    addForce(new RodBendingForceSym(*this));
  } else addForce(new RodAnisoForce(*this));

  if (getViscosity() != 0.0) {
    RodStretchingForce* stretching = new RodStretchingForce(*this);
    RodBendingForceSym* bending = new RodBendingForceSym(*this);
    RodTwistingForceSym* twisting = new RodTwistingForceSym(*this);

    stretching->setViscous(true);
    bending->setViscous(true);
    twisting->setViscous(true);

    addForce(stretching);
    addForce(bending);
    addForce(twisting);
  }
}

} // namespace BASim
