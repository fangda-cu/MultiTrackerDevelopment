/**
 * \file BentTwisting.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/09/2009
 */

#ifndef BENTTWISTING_HH
#define BENTTWISTING_HH

#include "ProblemBase.hh"

/**
 * Test problem for the rods.
 */

class BentTwisting : public Problem
{
public:

  BentTwisting();

protected:

  void Setup();
  void AtEachTimestep();

  ElasticRod* rod;
  RodTimeStepper* stepper;

  Scalar m_maxTwist;
  Scalar m_twistRate;
  Scalar m_currentTwist;
};

#endif // BENTTWISTING_HH
