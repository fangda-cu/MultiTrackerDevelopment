/**
 * \file HairyBall.hh
 *
 * \author smith@cs.columbia.edu
 * \date 05/11/2010
 */

#ifndef HAIRYBALL_WETA_HH
#define HAIRYBALL_WETA_HH

#include "ProblemBase.hh"

#ifdef WETA
#include <weta/Wfigaro/Physics/ElasticRods/BARodStepper.hh>
#include <weta/Wfigaro/Math/SolverUtils.hh>
#include <weta/Wfigaro/Physics/ElasticRods/ParallelStepper.hh>
#include <weta/Wfigaro/Physics/ElasticRods/AdaptiveBinaryStepper.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodForce.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodStretchingForce.hh>
#include <weta/Wfigaro/Core/TriangleMesh.hh>
#include <weta/Wfigaro/Core/ScriptingController.hh>
#include <weta/Wfigaro/IO/ObjParser.hh>
#include <weta/Wfigaro/Physics/ElasticRods/TopologicalObjectSerializer.hh>
#else
#include "BASim/src/Physics/ElasticRods/BARodStepper.hh"
#include "BASim/src/Physics/ElasticRods/ParallelStepper.hh"
#include "BASim/src/Physics/ElasticRods/AdaptiveBinaryStepper.hh"
#include "BASim/src/Physics/ElasticRods/RodForce.hh"
#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"
#include "BASim/src/Core/TriangleMesh.hh"
#include "BASim/src/Core/ScriptingController.hh"
#include "BASim/src/IO/ObjParser.hh"
#include "BASim/src/Physics/ElasticRods/TopologicalObjectSerializer.hh"
#endif

#include <vector>
#include <fstream>
#include <iomanip>

class SphereRotatorW : public ScriptingController
{
public:

  SphereRotatorW( TriangleMesh& mesh, std::vector<ElasticRod*>& rods, double time, double dt );
 ~SphereRotatorW();

  bool execute();
  
  //void setTime( double time );
  //void setDt( double dt );
  
private:

  //double getTime();
  //double getDt();
  
  void transformTriangleObject( TriangleMesh& tri_mesh, const Mat3d& transformation );
  void transformRodRoot( ElasticRod* rod, const Mat3d& transformation );
  float roundto(float val, float prec);

  TriangleMesh& m_mesh;
  std::vector<ElasticRod*>& m_rods;
  //double m_time;
  //double m_dt;
  FILE * m_pFile;
};

class HairyBallWeta : public Problem
{
public:
  
  HairyBallWeta();
  virtual ~HairyBallWeta();
  
protected:
  
  void Setup();
  void AtEachTimestep();
  //void AtEachTimestepOLD();
  
  virtual void serialize( std::ofstream& of );
  virtual void resumeFromfile( std::ifstream& ifs );
  
private:
  
  void loadHandCodedVectors( std::vector<Vec3d>& normals );
  void loadNormalsFile( const std::string& filename, std::vector<Vec3d>& normals );
  
  // These assume the objects are cenetered at the origin. Easy to relax this, I just haven't.
  void transformTriangleObject( TriangleMesh* tri_mesh, const Mat3d& transformation );
  void transformRodRoot( ElasticRod* rod, const Mat3d& transformation );

  void generateStraightHair( const Vec3d& initnorm, const Vec3d& startpoint, const double& dL, const int& nv, std::vector<Vec3d>& vertices );
  void generateCurlyHair( const Vec3d& initnorm, const Vec3d& startpoint, const double& dL, const int& nv, std::vector<Vec3d>& vertices );

  std::vector<TriangleMesh*> m_tri_meshes;
  std::vector<ElasticRod*> m_rods;
  std::vector<ScriptingController*> m_scripting_controllers;
  std::vector<RodTimeStepper*> m_steppers;

  ParallelStepper* m_pr_stepper;
  BARodStepper* m_br_stepper;
  
  SphereRotatorW* m_rotator;
};

#endif // HAIRYBALL_HH
