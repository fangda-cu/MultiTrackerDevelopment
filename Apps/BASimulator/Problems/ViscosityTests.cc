#include "ViscosityTests.hh"

#include <iostream>
#include <fstream>

using namespace std;

ViscosityTests::ViscosityTests()
: Problem("Viscosity Tests", "A collection of tests to ensure viscosity functions correctly")
, m_rods()
, m_steppers()
{
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();
 
  // Values for the default options
  GetStringOpt("integrator") = "implicit";
  GetVecOpt("gravity") = Vec3d(0,0,0);
  GetScalarOpt("dt") = 0.0001;
  GetScalarOpt("mass-damping") = 0;
  GetScalarOpt("viscosity") = 0;
  GetScalarOpt("radius-scale") = 10.0;

  AddOption("ViscosityTest", "Which viscosity test to execute.", 0);
}

ViscosityTests::~ViscosityTests()
{
  for( int i = 0; i < (int) m_rods.size(); ++i )
  {
    assert( m_rods[i] != NULL );
    delete m_rods[i];
    m_rods[i] = NULL;
  }

  for( int i = 0; i < (int) m_steppers.size(); ++i )
  {
    assert( m_steppers[i] != NULL );
    delete m_steppers[i];
    m_steppers[i] = NULL;
  }
}

void ViscosityTests::setupStretchingViscosityTest()
{
  GetScalarOpt("dt") = 0.000001;
  
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.viscosity = 0.0;
  opts.numVertices = 3;

  std::vector<Vec3d> undeformed;
  undeformed.push_back(Vec3d(-1.0,0.0,0.0));
  undeformed.push_back(Vec3d( 0.0,0.0,0.0));
  undeformed.push_back(Vec3d( 1.0,0.0,0.0));
  
  std::vector<Vec3d> deformed;
  deformed.push_back(Vec3d(-1.0,0.0,0.0));
  deformed.push_back(Vec3d( 0.0,0.0,0.0));
  deformed.push_back(Vec3d( 1.5,0.0,0.0));
  
  ElasticRod* newrod;
  RodTimeStepper* newstepper;

  for( int i = 0; i < 20; ++i )
  {
    newrod = setupRod(opts,deformed,undeformed);
    newrod->getBoundaryCondition()->setDesiredVertexPosition(0,newrod->getVertex(0));
    newrod->getBoundaryCondition()->setDesiredVertexPosition(1,newrod->getVertex(1));
    newrod->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);  
    m_rods.push_back(newrod);
    
    newstepper = getRodTimeStepper(*newrod);
    m_steppers.push_back(newstepper);
    
    opts.viscosity += 500.0;
    for( int i = 0; i < (int) deformed.size(); ++i ) deformed[i] -= Vec3d(0.0,0.5,0.0);
  }

  for( int i = 0; i < (int) m_rods.size(); ++i )     m_world->addObject(m_rods[i]);
  for( int i = 0; i < (int) m_steppers.size(); ++i ) m_world->addController(m_steppers[i]);
  
  for( int i = 0; i < (int) m_rods.size(); ++i )
  {
    ObjPropHandle<Scalar> phndl;
    m_rods[i]->property_handle<Scalar>(phndl,"dynamic viscosity");
    std::cout << "Rod " << i << " viscosity " << m_rods[i]->property(phndl) << std::endl;
  }  
}

void ViscosityTests::setupBendingViscosityTest()
{
  GetScalarOpt("dt") = 0.0001;

  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.viscosity = 0.0;
  opts.numVertices = 5;
  
  std::vector<Vec3d> undeformed;
  undeformed.push_back(Vec3d(0.0,0.0,0.0));
  for( int i = 1; i < opts.numVertices; ++i ) undeformed.push_back(undeformed[undeformed.size()-1]+Vec3d(1.0,0.0,0.0));
  
  std::vector<Vec3d> deformed;
  deformed.push_back(Vec3d(0.0,0.0,0.0));
  deformed.push_back(Vec3d(1.0,0.0,0.0));
  for( int i = 2; i < opts.numVertices; ++i ) deformed.push_back(deformed[deformed.size()-1]+((i%2==0)?Vec3d(1/sqrt(2.0),1/sqrt(2.0),0.0):Vec3d(1/sqrt(2.0),-1/sqrt(2.0),0.0)));
  
  ElasticRod* newrod;
  RodTimeStepper* newstepper;
  
  for( int i = 0; i < 10; ++i )
  {
    newrod = setupRod(opts,deformed,undeformed);
    newrod->getBoundaryCondition()->setDesiredVertexPosition(0,newrod->getVertex(0));
    newrod->getBoundaryCondition()->setDesiredVertexPosition(1,newrod->getVertex(1));
    newrod->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);  
    m_rods.push_back(newrod);
    
    newstepper = getRodTimeStepper(*newrod);
    m_steppers.push_back(newstepper);
    
    opts.viscosity += 0.5e6;
    for( int i = 0; i < (int) deformed.size(); ++i ) deformed[i] -= Vec3d(0.0,1.0,0.0);
  }
  
  for( int i = 0; i < (int) m_rods.size(); ++i )     m_world->addObject(m_rods[i]);
  for( int i = 0; i < (int) m_steppers.size(); ++i ) m_world->addController(m_steppers[i]);
  
  for( int i = 0; i < (int) m_rods.size(); ++i )
  {
    ObjPropHandle<Scalar> phndl;
    m_rods[i]->property_handle<Scalar>(phndl,"dynamic viscosity");
    std::cout << "Rod " << i << " viscosity " << m_rods[i]->property(phndl) << std::endl;
  }
}

void ViscosityTests::Setup()
{
  switch (GetIntOpt("ViscosityTest")) 
  {
    case 0:
    {
      std::cout << "Executing test: Stretching Viscosity Test" << std::endl;
      setupStretchingViscosityTest();
      break;
    }
    case 1:
    {
      std::cout << "Executing test: Bending Viscosity Test" << std::endl;
      setupBendingViscosityTest();
      break;
    }
    default:
    {
      std::cout << "UNKNOWN TEST CASE: " << GetIntOpt("ViscosityTest") << std::endl;
      exit(1);
      break;
    }
  }
}

