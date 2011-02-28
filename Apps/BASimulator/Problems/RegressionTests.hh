/**
 * \file RegressionTests.hh
 *
 * \author smith@cs.columbia.edu
 * \date 04/30/2010
 */

#ifndef REGRESSIONTESTS_HH
#define REGRESSIONTESTS_HH

#include "ProblemBase.hh"

//#include "BASim/src/Core/TriangleMesh.hh"
//#include "BASim/src/IO/ObjParser.hh"
#include "BASim/src/Physics/ElasticRods/BridsonStepper.hh"
#include "BASim/src/Math/SolverUtils.hh"

class RegressionTests : public Problem
{
public:
  
  RegressionTests();
  virtual ~RegressionTests();
  
protected:

  /////////////////////////////////////////////////////////////
  // TESTS WITH UNRESOLVED BUGS

  // A rod that begins in a planar configuration buckles out of
  // plane with no twisting applied.
  void generateOutOfPlaneBuckling();
  void atEachTimestepOutOfPlaneBuckling();

  /////////////////////////////////////////////////////////////
  // TESTS WITH KNOWN ISSUES
  
  // FAILS FOR LAPACK. KNOWN BUG IN QUASISTATIC FRAME TREATMENT
  //   FOR ISOTROPIC RODS.
  // A rod with hair parameters deformed into a sawtooth.
  void generateStiffBendingNonFixed();
  void atEachTimestepStiffBendingNonFixed();

  // FAILS FOR CG. KNOWN ISSUE WITH NONSYMMETRIC MATRIX
  //   WHEN DEGREES OF FREEDOM FIXED.
  // A rod with hair parameters deformed into a sawtooth.
  // The first two vertices and first edge of this rod
  // are fixed.
  void generateStiffBendingFixed();
  void atEachTimestepStiffBendingFixed();



  /////////////////////////////////////////////////////////////
  // TESTS THAT HAVE BEEN FIXED
  
  // Used to fail with LAPACK (still have cg fixed dof problem). 
  // A rod with hair parameters with the first two vertices and first edge
  // fixed. This fixed edge is pulled with constant velocity.
  void generatePullTest();
  void atEachTimestepPullTest();
  
  // Used to fail with LAPACK.
  // Numeric instabilities in parallel transport were exacerbated by
  // rotations causing parallel edges to not be exactly parallel. 
  void generateTwoVertexBendingIsometry();
  void atEachTimestepGenerateTwoVertexBendingIsometry();
  
  
  
  /////////////////////////////////////////////////////////////
  // EXPECTED BEHAVIOR?
  
  // Fails with LAPACK, not with MKL.
  // Rod with an initial configuration that is very strained.
  // Lots of stored energy... possibly 'whip' or aliasing?
  void veryStrainedRodTest();
  void atEachTimestepVeryStrainedRodTest();

  
  void Setup();
  void AtEachTimestep();

private:
  std::vector<ElasticRod*> m_rods;
  std::vector<RodTimeStepper*> m_steppers;
};

#endif 

