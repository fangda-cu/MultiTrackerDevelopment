/**
 * \file SHO.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/09/2009
 */

#ifndef SHO_HH
#define SHO_HH

#include "ProblemBase.hh"

/**
 * Test problem for the rods.
 */

class SHO : public Problem
{
public:

  SHO();
  ~SHO();

protected:

  void Setup();
  void AtEachTimestep();

  ElasticRod* rod;
  RodTimeStepper* stepper;
};

#endif // SHO_HH
