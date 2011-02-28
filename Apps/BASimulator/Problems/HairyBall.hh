/**
 * \file HairyBall.hh
 *
 * \author smith@cs.columbia.edu
 * \date 05/11/2010
 */

#ifndef HAIRYBALL_HH
#define HAIRYBALL_HH

#include "ProblemBase.hh"
#include "BASim/src/Physics/ElasticRods/BridsonStepper.hh"
#include "BASim/src/Physics/ElasticRods/ParallelStepper.hh"
#include "BASim/src/Physics/ElasticRods/AdaptiveBinaryStepper.hh"
#include "BASim/src/Physics/ElasticRods/RodForce.hh"
#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"
#include "BASim/src/Core/TriangleMesh.hh"
#include "BASim/src/Core/ScriptingController.hh"
#include "BASim/src/IO/ObjParser.hh"
#include "BASim/src/Physics/ElasticRods/TopologicalObjectSerializer.hh"

#include <vector>
#include <fstream>
#include <iomanip>

class SphereRotator : public ScriptingController
{
public:

  SphereRotator( TriangleMesh& mesh, std::vector<ElasticRod*>& rods, double time, double dt );

  bool execute();
  
  //void setTime( double time );
  //void setDt( double dt );
  
private:

  //double getTime();
  //double getDt();
  
  void transformTriangleObject( TriangleMesh& tri_mesh, const Mat3d& transformation );
  void transformRodRoot( ElasticRod* rod, const Mat3d& transformation );

  TriangleMesh& m_mesh;
  std::vector<ElasticRod*>& m_rods;
  //double m_time;
  //double m_dt;
};

class HairyBall : public Problem
{
public:
  
  HairyBall();
  virtual ~HairyBall();
  
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
  BridsonStepper* m_br_stepper;
  
  SphereRotator* m_rotator;
};

#endif // HAIRYBALL_HH
