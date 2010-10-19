/**
 * \file MicrotuboliDNA.hh
 *
 * \author smith@cs.columbia.edu
 * \date 08/06/2010
 */

// This example requires a general sparse solver; for now assume Pardiso.
// TODO: Get CG working again.
#ifdef HAVE_PARDISO

#ifndef MICROTUBOLIDNA_HH
#define MICROTUBOLIDNA_HH

#include "ProblemBase.hh"
#include "BASim/src/Physics/ElasticRods/RodRodSpringForce.hh"

#include "BASim/src/Render/SpringRenderer.hh"

class MicrotuboliDNA : public Problem
{
public:

  MicrotuboliDNA();
  virtual ~MicrotuboliDNA();

protected:

  void Setup();
  void AtEachTimestep();

private:

  ElasticRod* createRootedLowerRod( const Vec3d& rootlocation, const RodOptions& opts );
  ElasticRod* createRootedUpperRod( const Vec3d& rootlocation, const RodOptions& opts );
  
  std::vector<ElasticRod*> m_rods;
  std::vector<RodTimeStepper*> m_steppers;
  MultipleRodTimeStepper* m_multiple_stepper;
  
  std::vector<RodRodSpringForce*> m_rod_rod_springs;

  double m_current_rad;
};

#endif

#endif // MICROTUBOLIDNA_HH
