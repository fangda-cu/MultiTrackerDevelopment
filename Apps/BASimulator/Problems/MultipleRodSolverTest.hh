/**
 * \file MultipleRodSolverTest.hh
 *
 * \author smith@cs.columbia.edu
 * \date 05/29/2010
 */

#ifndef MULTIPLERODSOLVERTEST_HH
#define MULTIPLERODSOLVERTEST_HH

#include "ProblemBase.hh"
#include "BASim/src/Physics/ElasticRods/RodRodSpringForce.hh"
//#include "BASim/src/Physics/ElasticRods/BARodStepper.hh"
//#include "BASim/src/Physics/ElasticRods/RodForce.hh"
//#include <vector>
//#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"
//#include "BASim/src/Core/TriangleMesh.hh"
//#include "BASim/src/IO/ObjParser.hh"

#include "BASim/src/Render/SpringRenderer.hh"

class MultipleRodSolverTest : public Problem
{
public:

  MultipleRodSolverTest();
  virtual ~MultipleRodSolverTest();

  void freeRodsSetup();
  void freeRodsAtEachTimestep();
  
  void rodsWithFixedEndSetup();
  void rodsWithFixedEndAtEachTimestep();

  void rodSpinSetup();
  void rodSpinAtEachTimestep();

  void rodSpringSetup();
  void rodSpringAtEachTimestep();

protected:

  void Setup();
  void AtEachTimestep();

private:

  std::vector<ElasticRod*> m_rods;
  std::vector<RodTimeStepper*> m_steppers;
  MultipleRodTimeStepper* m_multiple_stepper;

  double m_current_rad;
};

#endif // MULTIPLERODSOLVERTEST_HH
