/**
 * \file MultipleRodSolverTest.cc
 *
 * \author smith@cs.columbia.edu
 * \date 05/29/2010
 */

#include "MultipleRodSolverTest.hh"

MultipleRodSolverTest::MultipleRodSolverTest()
: Problem("MultipleRodSolverTest", "MultipleRodSolverTest.")
, m_rods()
, m_steppers()
, m_multiple_stepper()
, m_current_rad(0.0)
{
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();
  
  // Set Defaults for built-in options
  GetStringOpt("integrator") = "implicit";
  GetVecOpt("gravity") = Vec3d(0,0,0);
  GetScalarOpt("dt") = 0.01;
  GetIntOpt("iterations") = 1000; // Maximum number of implicit solver itations
  GetIntOpt("nv") = 3;
  //GetScalarOpt("mass-damping") = 0;
  //GetScalarOpt("viscosity") = 0;

  //AddOption("RegressionTest", "Which regression test to execute.", 0);
}

MultipleRodSolverTest::~MultipleRodSolverTest()
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
  
  if( m_multiple_stepper != NULL )
  {
    delete m_multiple_stepper;
    m_multiple_stepper = NULL;
  }
}

void MultipleRodSolverTest::freeRodsSetup()
{
  loadDynamicsProps();
  std::cout << "TOP ROW:    One big system." << std::endl;
  std::cout << "BOTTOM ROW: Multiple isolated systems." << std::endl;
  
  RodOptions opts;
  getRodOptions(opts);
  opts.quasistatic = false;
  opts.radiusScale = 10.0;
  opts.YoungsModulus *= 1.0e-3;
  opts.ShearModulus *= 1.0e-3;
  
  std::vector<Vec3d> undeformed;
  undeformed.push_back(Vec3d(-1.0,0.0,0.0));
  undeformed.push_back(Vec3d( 0.0,0.0,0.0));
  undeformed.push_back(Vec3d( 1.0,0.0,0.0));

  std::vector<Vec3d> deformed;
  deformed.push_back(Vec3d(-1.0/sqrt(2.0),1.0/sqrt(2.0),0.0));
  deformed.push_back(Vec3d( 0.0,0.0,0.0));
  deformed.push_back(Vec3d(1.0/sqrt(2.0),1.0/sqrt(2.0),0.0));




  // Create the first rod
  ElasticRod* newrod = setupRod(opts, deformed, undeformed);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);

  // Create the second rod
  for( int i = 0; i < (int) deformed.size(); ++i ) deformed[i] += Vec3d(2.1,0.0,0.0);
  newrod = setupRod(opts, deformed, undeformed);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);

  // Create the third rod
  for( int i = 0; i < (int) deformed.size(); ++i ) deformed[i] += Vec3d(2.1,0.0,0.0);
  newrod = setupRod(opts, deformed, undeformed);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);


  // Create the multiple rod time stepper
  m_multiple_stepper = getMultipleRodTimeStepper();
  m_world->addController(m_multiple_stepper);

  // Add all rods to the rod time stepper
  for( int i = 0; i < (int) m_rods.size(); ++i ) m_multiple_stepper->addRod(m_rods[i]);




  // Do the same rods but in individual solvers
  for( int i = 0; i < (int) deformed.size(); ++i ) deformed[i] -= Vec3d(4.2,0.0,0.0);
  for( int i = 0; i < (int) deformed.size(); ++i ) deformed[i] -= Vec3d(0.0,2.1,0.0);
  newrod = setupRod(opts, deformed, undeformed);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);
  RodTimeStepper* stepper = getRodTimeStepper(*newrod);
  m_steppers.push_back(stepper);
  m_world->addController(stepper);

  for( int i = 0; i < (int) deformed.size(); ++i ) deformed[i] += Vec3d(2.1,0.0,0.0);
  newrod = setupRod(opts, deformed, undeformed);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);
  stepper = getRodTimeStepper(*newrod);
  m_steppers.push_back(stepper);
  m_world->addController(stepper);
  
  for( int i = 0; i < (int) deformed.size(); ++i ) deformed[i] += Vec3d(2.1,0.0,0.0);
  newrod = setupRod(opts, deformed, undeformed);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);
  stepper = getRodTimeStepper(*newrod);
  m_steppers.push_back(stepper);
  m_world->addController(stepper);
}

void MultipleRodSolverTest::freeRodsAtEachTimestep()
{
}


void MultipleRodSolverTest::rodsWithFixedEndSetup()
{
  GetVecOpt("gravity") = Vec3d(0.0,-981.0,0.0);
  //setGravity(GetVecOpt("gravity"));
  
  loadDynamicsProps();
  std::cout << "TOP ROW:    One big system." << std::endl;
  std::cout << "BOTTOM ROW: Multiple isolated systems." << std::endl;
  
  RodOptions opts;
  getRodOptions(opts);
  opts.quasistatic = true;
  opts.radiusScale = 10.0;
  opts.YoungsModulus *= 1.0e-1;
  opts.ShearModulus *= 1.0e-1;
  opts.numVertices = 30;

  double L = 3.0;
  double dL = L/((double)(opts.numVertices-1));
  Vec3d dx = dL*Vec3d(1.0,0.0,0.0);
  std::vector<Vec3d> undeformed;
  undeformed.push_back(Vec3d::Zero());
  for( int i = 0; i < opts.numVertices-1; ++i ) undeformed.push_back(undeformed.back()+dx);
  
  std::vector<Vec3d> deformed;
  deformed = undeformed;
  
  // Add the first rod
  ElasticRod* newrod = setupRod(opts, deformed, undeformed);
  RodBoundaryCondition* boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);
  
  // Add the second rod
  for( int i = 0; i < opts.numVertices; ++i ) deformed[i] += Vec3d(0.0,0.0,-1.0);
  newrod = setupRod(opts, deformed, undeformed);
  boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);
  
  // Add the third rod
  for( int i = 0; i < opts.numVertices; ++i ) deformed[i] += Vec3d(0.0,0.0,-1.0);
  newrod = setupRod(opts, deformed, undeformed);
  boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);
  
  // Add the fourth rod
  for( int i = 0; i < opts.numVertices; ++i ) deformed[i] += Vec3d(0.0,0.0,-1.0);
  newrod = setupRod(opts, deformed, undeformed);
  boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);
  
  // Create the multiple rod time stepper
  m_multiple_stepper = getMultipleRodTimeStepper();
  m_world->addController(m_multiple_stepper);
  
  // Add all rods to the rod time stepper
  for( int i = 0; i < (int) m_rods.size(); ++i ) m_multiple_stepper->addRod(m_rods[i]);  
  
  
  
  
  // Create a number of rods with their own solvers
  for( int i = 0; i < opts.numVertices; ++i ) deformed[i] += Vec3d(0.0,-1.0,0.0);
  for( int i = 0; i < opts.numVertices; ++i ) deformed[i] += Vec3d(0.0,0.0,3.0);
  newrod = setupRod(opts, deformed, undeformed);
  boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);
  RodTimeStepper* stepper = getRodTimeStepper(*newrod);
  m_steppers.push_back(stepper);
  m_world->addController(stepper);

  for( int i = 0; i < opts.numVertices; ++i ) deformed[i] += Vec3d(0.0,0.0,-1.0);
  newrod = setupRod(opts, deformed, undeformed);
  boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);
  stepper = getRodTimeStepper(*newrod);
  m_steppers.push_back(stepper);
  m_world->addController(stepper);
  
  for( int i = 0; i < opts.numVertices; ++i ) deformed[i] += Vec3d(0.0,0.0,-1.0);
  newrod = setupRod(opts, deformed, undeformed);
  boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);
  stepper = getRodTimeStepper(*newrod);
  m_steppers.push_back(stepper);
  m_world->addController(stepper);
  
  for( int i = 0; i < opts.numVertices; ++i ) deformed[i] += Vec3d(0.0,0.0,-1.0);
  newrod = setupRod(opts, deformed, undeformed);
  boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  m_rods.push_back(newrod);
  m_world->addObject(newrod);
  stepper = getRodTimeStepper(*newrod);
  m_steppers.push_back(stepper);
  m_world->addController(stepper);
}

void MultipleRodSolverTest::rodsWithFixedEndAtEachTimestep()
{
}

void MultipleRodSolverTest::rodSpinSetup()
{
  GetVecOpt("gravity") = Vec3d(0.0,-981.0,0.0);
  //setGravity(GetVecOpt("gravity"));
  
  loadDynamicsProps();
  std::cout << "TOP ROW:    One big system." << std::endl;
  std::cout << "BOTTOM ROW: Multiple isolated systems." << std::endl;

  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;
  opts.numVertices = 20;
  
  // 12 inches == 60.96 cm
  double L = 60.96/3.0;
  
  double dx = 0.3;
  double dy = 0.3;
  double dh = L/((double)opts.numVertices-1);

  int numrods = 8;
  double radius = 3.0;
  double dtheta = 2.0*pi/((double)numrods);
  double theta = 0.0;
  for( int i = 0; i < numrods; ++i )
  {
    Vec3d baseposition(radius*cos(theta),0.0,radius*sin(theta));
    std::vector<Vec3d> vertices;
    for( int k = 0; k < opts.numVertices; ++k )
    {
      Vec3d vert = (baseposition+Vec3d(0,-k*dh,0));
      //Vec3d vert = (baseposition+Vec3d(sin(-k*dh),-k*dh,cos(-k*dh)));
      vertices.push_back(vert);
    }
    assert( opts.numVertices == (int) vertices.size() );
    
    ElasticRod* newrod = setupRod(opts, vertices, vertices);
    RodBoundaryCondition* boundary = newrod->getBoundaryCondition();
    boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
    boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
    boundary->setDesiredEdgeAngle(0,0.0);
    
    m_rods.push_back(newrod);
    m_world->addObject(newrod);
    
    //RodTimeStepper* tstep = getRodTimeStepper(*newrod);
    //m_world->addController(tstep);
    //m_steppers.push_back(tstep);
    theta += dtheta;
  }  
  
  
  
  // Create the multiple rod time stepper
  m_multiple_stepper = getMultipleRodTimeStepper();
  m_world->addController(m_multiple_stepper);
  
  // Add all rods to the rod time stepper
  for( int i = 0; i < (int) m_rods.size(); ++i ) m_multiple_stepper->addRod(m_rods[i]);
  
  
  theta = 0.0;
  for( int i = 0; i < numrods; ++i )
  {
    Vec3d baseposition(radius*cos(theta),0.0,radius*sin(theta));
    baseposition += 6.0*radius*Vec3d(1.0,0.0,0.0);
    std::vector<Vec3d> vertices;
    for( int k = 0; k < opts.numVertices; ++k )
    {
      Vec3d vert = (baseposition+Vec3d(0,-k*dh,0));
      //Vec3d vert = (baseposition+Vec3d(sin(-k*dh),-k*dh,cos(-k*dh)));
      vertices.push_back(vert);
    }
    assert( opts.numVertices == (int) vertices.size() );
    
    ElasticRod* newrod = setupRod(opts, vertices, vertices);
    RodBoundaryCondition* boundary = newrod->getBoundaryCondition();
    boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
    boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
    boundary->setDesiredEdgeAngle(0,0.0);
    
    m_rods.push_back(newrod);
    m_world->addObject(newrod);
    
    RodTimeStepper* tstep = getRodTimeStepper(*newrod);
    m_world->addController(tstep);
    m_steppers.push_back(tstep);
    theta += dtheta;
  }  
}


void MultipleRodSolverTest::rodSpinAtEachTimestep()
{
  double radspersec = 4.0*pi;
  double drad = GetScalarOpt("dt")*radspersec;
  
  //m_current_rad += drad; 
  
  for( int i = 0; i < (int) m_rods.size()/2; ++ i )
  {
    Vec3d v0 = m_rods[i]->getVertex(0);
    Vec3d v1 = m_rods[i]->getVertex(1);

    //Mat3d rotz;
    //rotz << cos(drad),-sin(drad), 0.0,
    //        sin(drad), cos(drad), 0.0,
    //        0.0,       0.0,       1.0;
    
    Mat3d roty;
    roty << cos(drad),0,sin(drad),
            0,1,0,
            -sin(drad),0,cos(drad);

    v0 = roty*v0;
    v1 = roty*v1;
    
    m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(0,v0);
    m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(1,v1);
  }
  
  double radius = 3.0;
  Vec3d translation = 6.0*radius*Vec3d(1.0,0.0,0.0);

  for( int i = (int) m_rods.size()/2; i < (int) m_rods.size(); ++i )
  {
    Vec3d v0 = m_rods[i]->getVertex(0)-translation;
    Vec3d v1 = m_rods[i]->getVertex(1)-translation;
    
    //Mat3d rotz;
    //rotz << cos(drad),-sin(drad), 0.0,
    //        sin(drad), cos(drad), 0.0,
    //        0.0,       0.0,       1.0;
    
    Mat3d roty;
    roty << cos(drad),0,sin(drad),
    0,1,0,
    -sin(drad),0,cos(drad);
    
    v0 = roty*v0+translation;
    v1 = roty*v1+translation;
    
    m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(0,v0);
    m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(1,v1);
  }
  
}

void MultipleRodSolverTest::rodSpringSetup()
{
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.quasistatic = false;
  opts.radiusScale = 10.0;
  opts.numVertices = 50;
  
  
  double L = 3.0;
  double dL = L/((double)(opts.numVertices-1));
  Vec3d dy = dL*Vec3d(0.0,1.0,0.0);
  std::vector<Vec3d> undeformed;
  undeformed.push_back(Vec3d::Zero());
  for( int i = 0; i < opts.numVertices-1; ++i ) undeformed.push_back(undeformed.back()-dy);
  
  std::vector<Vec3d> deformed;
  deformed = undeformed;
  
  // Add the first rod
  ElasticRod* newrod = setupRod(opts, deformed, undeformed);
  RodBoundaryCondition* boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  boundary->setDesiredVertexPosition(newrod->nv()-1, newrod->getVertex(newrod->nv()-1));
  boundary->setDesiredVertexPosition(newrod->nv()-2, newrod->getVertex(newrod->nv()-2));
  boundary->setDesiredEdgeAngle(newrod->ne()-1,0.0);  
  m_rods.push_back(newrod);
  m_world->addObject(newrod);
  
  // Add the second rod
  double dist = 2.0;
  for( int i = 0; i < (int) deformed.size(); ++i ) deformed[i] += Vec3d(dist,0.0,0.0);
  newrod = setupRod(opts, deformed, undeformed);
  boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  boundary->setDesiredVertexPosition(newrod->nv()-1, newrod->getVertex(newrod->nv()-1));
  boundary->setDesiredVertexPosition(newrod->nv()-2, newrod->getVertex(newrod->nv()-2));
  boundary->setDesiredEdgeAngle(newrod->ne()-1,0.0);  
  m_rods.push_back(newrod);
  m_world->addObject(newrod);
  
  
  // Create the multiple rod time stepper
  m_multiple_stepper = getMultipleRodTimeStepper();
  m_world->addController(m_multiple_stepper);
  
  // Add all rods to the rod time stepper
  for( int i = 0; i < (int) m_rods.size(); ++i ) m_multiple_stepper->addRod(m_rods[i]);
  
  
  RodRodSpringForce* rodrodspring = new RodRodSpringForce( m_rods, 0, 1, opts.numVertices/2, opts.numVertices/2, 10000.0, 1.0 );
  
  SpringRenderer* springRenderer = new SpringRenderer(*rodrodspring);
  m_world->addRenderer(springRenderer);
  
  m_multiple_stepper->addRodRodExternalForce(rodrodspring);
}

void MultipleRodSolverTest::rodSpringAtEachTimestep()
{
}


void MultipleRodSolverTest::Setup()
{
  //freeRodsSetup();
  //rodsWithFixedEndSetup();
  //rodSpinSetup();
  rodSpringSetup();
}

void MultipleRodSolverTest::AtEachTimestep()
{
  //freeRodsAtEachTimestep();
  //rodsWithFixedEndAtEachTimestep();
  //rodSpinAtEachTimestep();
  rodSpringAtEachTimestep();
}









