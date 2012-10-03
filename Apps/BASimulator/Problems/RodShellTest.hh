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
#include "BASim/src/Physics/DeformableObjects/Solids/ElasticSolid.hh"
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
  void setupScene3(); // simple rod shell test: scene 1 with a rod rib
  void setupScene4(); // umbrella opening
  void setupScene5(); // car sunshade folding
  void setupScene6(); // shell contraction
  void setupScene7(); // collapsible tunnel folding
  void setupScene8(); // collapsible tunnel opening
  void setupScene9(); // balls in bag
  void setupScene10();// Houdini exporter/loader test
  void setupScene11();// Solid test scene
  void setupScene12();// Rod surface tension test
  void setupScene13();// rod twist shell face coupling test
  void setupScene14();// rod twist solid tet coupling test
  void setupScene15();// non manifold coupling between rod shell and solid

protected:
  void Setup();
  void AtEachTimestep();
  
  DeformableObject * obj;
  ElasticRodModel * rod;
  ElasticShell * shell;
  ElasticSolid * solid;

  DefoObjTimeStepper * stepper;
  
  Scalar m_timestep;
  
  int m_active_scene;
  Scalar m_time;

  //This lets us choose between:
  //  1) Setting the radii uniformly as usual (default); or 
  //  2) Setting irregular radii specifically in the scene setup functions, using m_radii_list
  bool m_rod_radii_assigned;
  EdgeProperty<Vec2d>* m_radii_list;

  ////////////////////////
  // scene specific data
  // scene 5:
  VertexHandle m_s5_l1;
  VertexHandle m_s5_l2;
  VertexHandle m_s5_r1;
  VertexHandle m_s5_r2;

  EdgeHandle m_s5_le;
  EdgeHandle m_s5_re;

  // scene 7:
  std::vector<VertexHandle> m_s7_rightring;
  std::vector<VertexHandle> m_s7_leftring;
  Scalar m_s7_rot_rate;
  Scalar m_s7_vel;
  Scalar m_s7_tmax;
  
  // scene 9:
  std::vector<FaceHandle> m_s9_ball_faces;
  
  // scene 13:
  std::vector<EdgeHandle> m_s13_rod_edges;
  
  // scene 14:
  std::vector<EdgeHandle> m_s14_rod_edges;

  // scene 15:
  std::vector<VertexHandle> m_s15_vertices;
  std::vector<EdgeHandle> m_s15_edges;
  std::vector<FaceHandle> m_s15_faces;
  std::vector<TetHandle> m_s15_tets;
};



#endif
