/**
 * \file BentTwisting.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/09/2009
 */

#ifndef BENTTWISTING_HH
#define BENTTWISTING_HH

#include "ProblemBase.hh"

#ifdef WETA
#include <weta/Wfigaro/Physics/ElasticRods/TopologicalObjectSerializer.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodStateBackup.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodTimeStepper.hh>
#else
#include "BASim/src/Physics/ElasticRods/TopologicalObjectSerializer.hh"
#include "BASim/src/Physics/ElasticRods/RodStateBackup.hh"
#endif


/**
 * Test problem for the rods.
 */

class BentTwisting : public Problem
{
public:

  BentTwisting();
  virtual ~BentTwisting();

  virtual void serialize( std::ofstream& of );
  virtual void resumeFromfile( std::ifstream& ifs );

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
