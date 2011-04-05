/*
 * TestCase6.hh
 *
 *  Created on: 6/04/2011
 *      Author: jaubry
 */

#ifndef TESTCASE6_HH_
#define TESTCASE6_HH_

#include "ProblemBase.hh"

#ifdef WETA
#include <weta/Wfigaro/Physics/ElasticRods/BridsonStepper.hh>
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
#include "BASim/src/Physics/ElasticRods/BridsonStepper.hh"
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

namespace BASim
{

class TestCase6: public BASim::Problem
{
public:
    TestCase6();
    virtual ~TestCase6();

    virtual void Setup();
    virtual void AtEachTimestep();
    virtual void AfterLoad();
    virtual void AfterStep();

private:
    std::vector<TriangleMesh*> m_tri_meshes;
    std::vector<ElasticRod*> m_rods;
    std::vector<ScriptingController*> m_scripting_controllers;
    std::vector<RodTimeStepper*> m_steppers;

    ParallelStepper* m_pr_stepper;
    BridsonStepper* m_br_stepper;

};

}

#endif /* TESTCASE6_HH_ */
