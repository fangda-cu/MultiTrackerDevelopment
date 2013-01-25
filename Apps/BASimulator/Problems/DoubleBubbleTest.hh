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
#include "ElTopo/eltopo3d/surftrack.h"

class Recording
{
public:
  Recording() : m_recording_name("rec"), m_current_frame(0), m_current_step(0), m_recording(false) { }
  
  void setRecordingName(const std::string & name) { m_recording_name = name; }
  void setCurrentFrame(int frame) { m_current_frame = frame; }
  
  void recordSurfTrack(ElTopo::SurfTrack & st);
  
  void turnOnRecording() { m_recording = true; }
  void turnOffRecording() { m_recording = false; }
  
public:
  static void writeSurfTrack(std::ostream & os, ElTopo::SurfTrack & st);
  static void readSurfTrack(std::istream & is, ElTopo::SurfTrack & st);
  
protected:
  std::string m_recording_name;
  int m_current_frame;  // frame number
  int m_current_step;   // step number within the frame
  
  bool m_recording;
  
  std::ofstream m_of;
  std::ifstream m_if;
};

extern Recording g_recording;

class DoubleBubbleTest : public Problem, public ElasticShell::SteppingCallback
{
public:
  DoubleBubbleTest();
  virtual ~DoubleBubbleTest();

  virtual void serialize( std::ofstream& of );
  virtual void resumeFromfile( std::ifstream& ifs );

  void beforeEndStep();
  
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
  
  int onBBWall(const Vec3d & pos) const;
  void updateBBWallConstraints();
    
  ElTopo::SurfTrack * mesh2surftrack();
  void surftrack2mesh(ElTopo::SurfTrack & st);
  
public:
  void setupScene1(); // VIIM test: single film in cube
  void setupScene2(); // T1 transition
//  void setupScene3(); // double bubble collision
//  void setupScene4(); // n bubble collision
  void setupScene5(); // VIIM figure 17
  void setupScene6(); // VIIM multiphase cube test (figure 24)
    
  void setupScene7(); // Enright test with a sphere
  void setupScene8(); // Reauleux tetrahedron test VIIM figure 18
  void setupScene9(); // Normal motion VIIM figure 20
  void setupScene10(); // Normal motion with cyclic speed relation VIIM figure 21
  

  void s7_enright_velocity(double t, const Vec3d & pos, Vec3d & out);
  
  void createIcoSphere(DeformableObject & mesh, Vec3d & center, Scalar r, int subdivision, std::vector<VertexHandle> & vertList, std::vector<FaceHandle> & faceList, VertexProperty<Vec3d> & positions);

  int m_nregion;

};

#endif // DOUBLEBUBBLETEST_HH

