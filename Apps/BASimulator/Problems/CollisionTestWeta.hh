/**
 * \file CollisionTestTwo.hh
 *
 * \author smith@cs.columbia.edu
 * \date 04/27/2010
 */

#ifndef COLLISIONTEST_WETA_HH
#define COLLISIONTEST_WETA_HH

#include "ProblemBase.hh"
#include "BASim/src/Physics/ElasticRods/BridsonStepper.hh"
#include "BASim/src/Physics/ElasticRods/RodForce.hh"
#include <vector>
#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"
#include "BASim/src/Core/TriangleMesh.hh"
#include "BASim/src/IO/ObjParser.hh"
#include "BASim/src/Core/ScriptingController.hh"
#include <stdio.h>

class CollisionTestWeta : public Problem
{
public:
  
  CollisionTestWeta();
  virtual ~CollisionTestWeta();
  
protected:
  
  void Setup();
  void AtEachTimestep();
  
private:

  /////////////////////////////////////////////////////////////////////////////
  // ROD-OBJECT SNAGGING

  void RodObjectSnaggingSetup();
  void RodObjectSnaggingAtEachTimestep();
  void MovingSphereSetup();
  void RodObjectSnaggingMovingObjectSetup();
  void RodObjectSnaggingMovingObjectAtEachTimestep();

  void RodObjectSnaggingVertFaceSetup();

  void RodObjectSnaggingVertFaceTwoSetup();

  /////////////////////////////////////////////////////////////////////////////
  // ROD-ROD SNAGGING

  void RodRodSnaggingSetup();
  void RodRodSnaggingAtEachTimestep();

  void RodRodSnaggingTwoSetup();
  void RodRodSnaggingTwoAtEachTimestep();

  void RodRodSnaggingThreeSetup();
  void RodRodSnaggingThreeAtEachTimestep();

  void RodFixedRodSnaggingSetup();
  void RodFixedRodSnaggingAtEachTimestep();

  void RodRodSnaggingSmallSetup();
  void RodRodSnaggingSmallAtEachTimestep();

  void RodRodSnaggingDifferentSizeSetup();

  std::vector<ElasticRod*> m_rods;
  std::vector<TriangleMesh*> m_tri_objs;
  std::vector<RodTimeStepper*> m_controllers;
  BridsonStepper* m_br_stepper;
  
  std::vector<ScriptingController*> m_scripting_controllers;
};

class RodObjSnaggingController : public ScriptingController
{
public:
  RodObjSnaggingController( ElasticRod& rod, double time, double dt, FILE * pfile );
  bool execute();

private:
  ElasticRod& m_rod;
  FILE * m_pFile;
};

class RodRodSnagController : public ScriptingController
{
public:
  RodRodSnagController( ElasticRod& rod, double time, double dt );
  bool execute();

private:
  ElasticRod& m_rod;
};

class RodFixedRodSnagScriptingController : public ScriptingController
{
public:
  RodFixedRodSnagScriptingController( ElasticRod& rod, double time, double dt );
  bool execute();

private:
  ElasticRod& m_rod;
};

class ObjTranslator : public ScriptingController
{
public:
  ObjTranslator( TriangleMesh& mesh, double time, double dt );
  bool execute();

private:
  TriangleMesh& m_mesh;
};

class ObjTranslatorWeta : public ScriptingController
{
public:
  ObjTranslatorWeta( ElasticRod& rod, TriangleMesh& mesh, double time, double dt,  FILE * pfile );
  bool execute();

private:
  TriangleMesh& m_mesh;
  ElasticRod& m_rod;
  FILE * m_pFile;
};

#endif // COLLISIONTEST_HH
