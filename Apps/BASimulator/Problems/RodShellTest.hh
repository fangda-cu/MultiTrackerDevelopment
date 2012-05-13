//
//  RodShellTest.hh
//  BASim
//
//  Created by Fang Da (fang@cs.columbia.edu) on 5/12/12.
//  Copyright (c) 2012 Columbia University. All rights reserved.
//

#ifndef BASim_RodShellTest_hh
#define BASim_RodShellTest_hh

#include "ProblemBase.hh"

#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShell.hh"
#include "BASim/src/Physics/DeformableObjects/Rods/ElasticRodModel.hh"
#include "BASim/src/Physics/DeformableObjects/DefoObjTimeStepper.hh"

/**
 * Test problem for the rod + shell integration
 */

class RodShellTest : public Problem
{
public:
  RodShellTest();
  virtual ~RodShellTest();
  
//  virtual void serialize(std::ofstream & of);
//  virtual void resumeFromfile(std::ifstream & ifs);
  
public:
  void setupScene1(); // shell test: vertical flat sheet
  void setupScene2(); // rod test: bent twisting
  
protected:
  void Setup();
  void AtEachTimestep();
  
  DeformableObject * obj;
  ElasticRodModel * rod;
  ElasticShell * shell;

  DefoObjTimeStepper * stepper;
  
  Scalar m_timestep;
  
  int m_active_scene;

};




#endif
