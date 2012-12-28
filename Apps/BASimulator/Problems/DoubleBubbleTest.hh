/**
 * \file DoubleBubbleTest.hh
 *
 * \author fang@cs.columbia.edu
 * \date Nov 20, 2012
 */

#ifndef DOUBLEBUBBLETEST_HH
#define DOUBLEBUBBLETEST_HH

#include "ProblemBase.hh"

#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShell.hh"
#include "BASim/src/Physics/DeformableObjects/DefoObjTimeStepper.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellVolumeForce.hh"

class DoubleBubbleTest : public Problem
{
public:
  DoubleBubbleTest();
  virtual ~DoubleBubbleTest();

  virtual void serialize( std::ofstream& of ) { assert(!"Not implemented"); }
  virtual void resumeFromfile( std::ifstream& ifs ) { assert(!"Not implemented"); }

protected:
  void Setup();
  void AtEachTimestep();
  void AfterStep();
  
  DeformableObject * shellObj;
  ElasticShell * shell;
  DefoObjTimeStepper * stepper;

  Scalar m_timestep;
//  Scalar m_initial_thickness;
  int m_active_scene;

//  int m_s4_nbubble;
  
  std::vector<VertexHandle> triangulation_added_vertices;
  std::vector<EdgeHandle>   triangulation_added_edges;
  std::vector<FaceHandle>   triangulation_added_faces;
  ShellVolumeForce * svf;
  
public:
  void setupScene1(); // VIIM test: single film in cube
  void setupScene2(); // T1 transition
//  void setupScene3(); // double bubble collision
//  void setupScene4(); // n bubble collision
  void setupScene5(); // VIIM figure 5
  void setupScene6(); // VIIM multiphase cube test
  

};

#endif // DOUBLEBUBBLETEST_HH

