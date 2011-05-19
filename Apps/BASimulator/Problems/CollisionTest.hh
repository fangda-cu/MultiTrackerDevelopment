///**
// * \file CollisionTest.hh
// *
// * \author smith@cs.columbia.edu
// * \date 02/17/2010
// */
//
//#ifndef COLLISIONTEST_HH
//#define COLLISIONTEST_HH
//
//#include "ProblemBase.hh"
//#include "BASim/src/Physics/ElasticRods/BARodStepper.hh"
//#include "BASim/src/Physics/ElasticRods/RodForce.hh"
//#include <vector>
//#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"
//#include "BASim/src/Core/TriangleMesh.hh"
//#include "BASim/src/IO/ObjParser.hh"
//
//// TODO: check that all of the other tests place triangulated objects into the vector that gets deleted
//
///**
// * A number of goofy tests to see if the collision system works.
// */
//class CollisionTest : public Problem
//{
//public:
//  
//  CollisionTest();
//  virtual ~CollisionTest();
//  
//protected:
//  
//  void Setup();
//  void AtEachTimestep();
//  
//private:
//  
//  
//  /////////////////////////////////////////////////////////////////////////////
//  // Tests of full collision system
//  
//  // Just tests if BARodStepper can load triangle meshes.
//  void genTorusTest();
//  // Test in which nothing happens, and nothing should happen.
//  void genNullTest();
//  // A clump of hair is pulled through a torus.
//  void genTorusPull(); 
//  void genTorusPullAtEachTimestep();
//
//  /////////////////////////////////////////////////////////////////////////////
//  // Rod-rod penalty tests
//
//  void genRodRodPenaltyBasic00();
//  void genRodRodPenaltyBasic01();
//  void genRodRodPenaltyCorner00();
//  void genRodRodPenaltyFixed00();
//  void genRodRodPenaltyFixed01();
//  void genRodRodPenaltyDifferentRadii00();
//  void genRodRodMultipleContact00(); // Scale this example up more, 
//
//  /////////////////////////////////////////////////////////////////////////////
//  // Rod-object penalty tests
//  
//  void genRodObjectPenaltyBasic00();
//  void genRodObjectPenaltyBasic01(); // This example has some wierdness with collisions speeding things up ...
//  // Drop a bunch of rods on a sphere
//  void genRodObjectPenaltyBasic02();
//  // Drop a bunch of rods into a box
//  void genRodObjectPenaltyBasic03();
//  
//  /////////////////////////////////////////////////////////////////////////////
//  // Edge-edge iterative inelastic impulse tests
//  
//  void genEdgeEdgeImpulseBasic00(); // Causes lots of funniness with root finding code for far apart geometry. 
//  void genEdgeEdgeImpulseDifferentMass00();
//  
//  
//  /////////////////////////////////////////////////////////////////////////////
//  // Vertex-face iterative inelastic impulse tests
//  
//  void genVertexFaceImpulseBasic00();
//  // Rod hits a sphere
//  void genVertexFaceImpulseBasic01();
//
//  
//  std::vector<ElasticRod*> m_rods;
//  std::vector<TriangleMesh*> m_tri_objs;
//  std::vector<RodTimeStepper*> m_controllers;
//  BARodStepper* m_br_stepper;
//  
//  std::vector<ScriptingController*> m_scripting_controllers;
//};
//
//#endif // COLLISIONTEST_HH




