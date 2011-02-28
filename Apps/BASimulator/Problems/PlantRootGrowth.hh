/**
 * \file PlantRootGrowth.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 05/13/2010
 */

#ifndef PLANTROOTGROWTH_HH
#define PLANTROOTGROWTH_HH

#include "ProblemBase.hh"

/**
 * Test problem for the rods.
 */

class PlantRootGrowth : public Problem
{
public:

  PlantRootGrowth();
  virtual ~PlantRootGrowth();

protected:

  void Setup();
  void AtEachTimestep();

  ElasticRod* rod;
  RodTimeStepper* stepper;
  RodGelForce* gelForce;

  Scalar rootLength;
  Scalar gelStiffness;
};

#endif // PLANTROOTGROWTH_HH
