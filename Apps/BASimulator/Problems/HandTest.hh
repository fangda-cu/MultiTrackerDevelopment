/**
 * \file HandTest.hh
 *
 * \author 
 * \date 
 */

#ifndef HandTest_HH
#define HandTest_HH

#include "ProblemBase.hh"
#include "BASim/src/Core/TriangleMesh.hh"
#include "BASim/src/IO/ObjParser.hh"

#include "BASim/src/Physics/ElasticRods/BridsonStepper.hh"
#include "BASim/src/Core/ScriptingController.hh"

class HandController : public ScriptingController
{
public:

  HandController( TriangleMesh& mesh, double time, double dt, std::vector<Vec3d> &scripted_translation, std::vector<double> &scripted_time );

  bool execute();
  
  void setTime( double time );
  void setDt( double dt );
  
private:

  double getTime();
  double getDt();
  
  TriangleMesh& m_mesh;
//  std::vector<ElasticRod*>& m_rods;
  double m_time;
  double m_dt;

	// scripted motion
  std::vector<Vec3d> hand_translation;
  std::vector<double> hand_time;

};


class HandTest : public Problem
{
public:

  HandTest();
  virtual ~HandTest();

protected:

  void Setup();
  void AtEachTimestep();

private:

	void setupHand();
	
  std::vector<ElasticRod*> m_rods;
//  std::vector<RodTimeStepper*> m_steppers;
  std::vector<TriangleMesh*> m_tri_objs;
  std::vector<RodTimeStepper*> m_controllers;
  BridsonStepper* m_br_stepper;
  std::vector<ScriptingController*> m_scripting_controllers;

  TriangleMesh* hand_mesh;
  
  HandController* hand_controller;
 
};

#endif // HandTest_HH
