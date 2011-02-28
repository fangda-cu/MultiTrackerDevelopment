/**
 * \file ViscosityTests.hh
 *
 * \author smith@cs.columbia.edu
 * \date 07/04/2010
 */

#ifndef VISCOSITYTESTS_HH
#define VISCOSITYTESTS_HH

#include "ProblemBase.hh"

class ViscosityTests : public Problem
{
public:
  
  ViscosityTests();
  virtual ~ViscosityTests();
  
protected:

  void setupStretchingViscosityTest();
  void setupBendingViscosityTest();
  
  void Setup();

private:
  std::vector<ElasticRod*> m_rods;
  std::vector<RodTimeStepper*> m_steppers;
};

#endif 

