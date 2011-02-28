#include "RegressionTests.hh"

#include <iostream>
#include <fstream>

using namespace std;

RegressionTests::RegressionTests()
: Problem("Regression Tests", "A collection of tests for problems that prove/proved difficult")
, m_rods()
, m_steppers()
{
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();
 
  AddOption("RegressionTest", "Which regression test to execute.", 0);
}

RegressionTests::~RegressionTests()
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



/////////////////////////////////////////////////////////////
// TESTS WITH UNRESOLVED BUGS

void RegressionTests::generateOutOfPlaneBuckling()
{
  // Values for the default options
  GetStringOpt("integrator") = "implicit";
  GetVecOpt("gravity") = Vec3d(0,0,0);
  GetScalarOpt("dt") = 0.001;
  // Maximum number of implicit solver itations
  GetIntOpt("iterations") = 1000; 
  GetIntOpt("nv") = 23;
  GetScalarOpt("mass-damping") = 0;
  GetScalarOpt("viscosity") = 0;
  GetScalarOpt("radius-scale") = 10.0;
  
  GetBoolOpt("quasistatic") = false;
  
  
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  //opts.radiusScale = GetScalarOpt("radius-scale");
  
  // 12 inches == 60.96 cm
  double L = 60.96/3.0;  
  double dh = L/((double)opts.numVertices-1);
  
  Vec3d baseposition(0.0,0.0,0.0);
  
  std::vector<Vec3d> deformed;
  deformed.push_back(baseposition);
  for( int k = 0; k < opts.numVertices-1; ++k )
  {
    Vec3d vert = deformed.back();
    if( k%2 == 0 ) vert += Vec3d(1.0/sqrt(2.0),1.0/sqrt(2.0),0.0);
    else vert += Vec3d(1.0/sqrt(2.0),-1.0/sqrt(2.0),0.0);
    deformed.push_back(vert);
  }
  
  std::vector<Vec3d> undeformed;
  undeformed.push_back(baseposition);
  for( int k = 0; k < opts.numVertices-1; ++k )
  {
    Vec3d vert = undeformed.back();
    vert += Vec3d(1.0,0.0,0.0);
    undeformed.push_back(vert);
  }
  
  
  ElasticRod* newrod = setupRod(opts,deformed,undeformed);
  
  newrod->getBoundaryCondition()->setDesiredVertexPosition(0,newrod->getVertex(0));
  newrod->getBoundaryCondition()->setDesiredVertexPosition(1,newrod->getVertex(1));
  newrod->getBoundaryCondition()->setDesiredVertexPosition(9,newrod->getVertex(9));
  newrod->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
  
  m_world->addObject(newrod);
  m_rods.push_back(newrod);
  
  RodTimeStepper* newstepper = getRodTimeStepper(*newrod);
  m_world->addController(newstepper);
  
  std::cout << "Solver Options: " << std::endl;
  std::cout << "   integrator: " << GetStringOpt("integrator") << std::endl;
  std::cout << "   dt:         " << GetScalarOpt("dt") << std::endl;
  std::cout << "   iterations: " << GetIntOpt("iterations") << std::endl;
  std::cout << "   solver: " << SolverUtils::instance()->getSolverName() << std::endl;
  
  std::cout << "\"World\" Options: " << std::endl;
  std::cout << "   mass-damping: " << GetScalarOpt("mass-damping") << std::endl;
  
  std::cout << "Rod Options: " << std::endl;
  std::cout << "   numVertices:   " << opts.numVertices << std::endl;
  std::cout << "   density:       " << opts.density << std::endl;
  std::cout << "   radiusA:       " << opts.radiusA << std::endl;
  std::cout << "   radiusB:       " << opts.radiusB << std::endl;
  std::cout << "   radiusScale:   " << opts.radiusScale << std::endl;
  std::cout << "   YoungsModulus: " << opts.YoungsModulus << std::endl;
  std::cout << "   ShearModulus:  " << opts.ShearModulus << std::endl;
  std::cout << "   viscosity:     " << opts.viscosity << std::endl;
  std::cout << "   anisotropic:   " << opts.anisotropic << std::endl;
  std::cout << "   elastic:       " << opts.elastic << std::endl;
  std::cout << "   quasistatic:   " << opts.quasistatic << std::endl;
  std::cout << "   inextensible:  " << opts.inextensible << std::endl;
  std::cout << "   reframe:       ";
  switch(opts.refFrame) 
  {
    case ElasticRod::TimeParallel:
      std::cout << "TimeParallel" << std::endl;
      break;
    case ElasticRod::SpaceParallel:
      std::cout << "SpaceParallel" << std::endl;
      break;
    default:
      std::cout << "WARNING NONHANDLED CASE" << std::endl;
      break;
  }  
}

void RegressionTests::atEachTimestepOutOfPlaneBuckling()
{
}

/////////////////////////////////////////////////////////////
// TESTS WITH KNOWN ISSUES

void RegressionTests::generateStiffBendingNonFixed()
{
  GetStringOpt("integrator") = "implicit";
  GetIntOpt("nv") = 13;
  GetIntOpt("iterations") = 100; //1000;
  GetScalarOpt("mass-damping") = 0;
  GetScalarOpt("viscosity") = 0;
  GetVecOpt("gravity") = Vec3d(0.0,0.0,0.0);
  GetScalarOpt("dt") = 0.1;
  
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;
  

  std::vector<Vec3d> undeformed;
  for( int i = 0; i < opts.numVertices; ++i ) undeformed.push_back(Vec3d((double) i,0.0,0.0));
  
  Vec3d up(1.0/sqrt(2.0),1.0/sqrt(2.0),0.0);
  Vec3d down(1.0/sqrt(2.0),-1.0/sqrt(2.0),0.0);
  std::vector<Vec3d> vertices;
  vertices.push_back(Vec3d(0.0,0.0,0.0));
  for( int i = 1; i < opts.numVertices; ++i ) 
  {
    if( i%2 == 1 ) vertices.push_back( vertices[i-1]+up );
    else vertices.push_back( vertices[i-1]+down );
  }

  ElasticRod* newrod = setupRod(opts, vertices, undeformed);
  m_world->addObject(newrod);
  m_rods.push_back(newrod);
  
  RodTimeStepper* newstepper = getRodTimeStepper(*newrod);
  m_world->addController(newstepper);
  m_steppers.push_back(newstepper);

  std::cout << "Solver Options: " << std::endl;
  std::cout << "   integrator: " << GetStringOpt("integrator") << std::endl;
  std::cout << "   dt:         " << GetScalarOpt("dt") << std::endl;
  std::cout << "   iterations: " << GetIntOpt("iterations") << std::endl;
  std::cout << "   solver: " << SolverUtils::instance()->getSolverName() << std::endl;
  
  std::cout << "\"World\" Options: " << std::endl;
  std::cout << "   mass-damping: " << GetScalarOpt("mass-damping") << std::endl;

  std::cout << "Rod Options: " << std::endl;
  std::cout << "   numVertices:   " << opts.numVertices << std::endl;
  std::cout << "   density:       " << opts.density << std::endl;
  std::cout << "   radiusA:       " << opts.radiusA << std::endl;
  std::cout << "   radiusB:       " << opts.radiusB << std::endl;
  std::cout << "   radiusScale:   " << opts.radiusScale << std::endl;
  std::cout << "   YoungsModulus: " << opts.YoungsModulus << std::endl;
  std::cout << "   ShearModulus:  " << opts.ShearModulus << std::endl;
  std::cout << "   viscosity:     " << opts.viscosity << std::endl;
  std::cout << "   anisotropic:   " << opts.anisotropic << std::endl;
  std::cout << "   elastic:       " << opts.elastic << std::endl;
  std::cout << "   quasistatic:   " << opts.quasistatic << std::endl;
  std::cout << "   inextensible:  " << opts.inextensible << std::endl;
  std::cout << "   reframe:       ";
  switch(opts.refFrame) 
  {
    case ElasticRod::TimeParallel:
      std::cout << "TimeParallel" << std::endl;
      break;
    case ElasticRod::SpaceParallel:
      std::cout << "SpaceParallel" << std::endl;
      break;
    default:
      std::cout << "WARNING NONHANDLED CASE" << std::endl;
      break;
  }
}

void RegressionTests::atEachTimestepStiffBendingNonFixed()
{
  // No scripting
}


void RegressionTests::generateStiffBendingFixed()
{  
  GetStringOpt("integrator") = "implicit";
  GetIntOpt("nv") = 13;
  GetIntOpt("iterations") = 100; //1000;
  GetScalarOpt("mass-damping") = 0;
  GetScalarOpt("viscosity") = 0;
  GetVecOpt("gravity") = Vec3d(0.0,0.0,0.0);
  GetScalarOpt("dt") = 0.1;
  
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;
  
  
  std::vector<Vec3d> undeformed;
  for( int i = 0; i < opts.numVertices; ++i ) undeformed.push_back(Vec3d((double) i,0.0,0.0));
  
  Vec3d up(1.0/sqrt(2.0),1.0/sqrt(2.0),0.0);
  Vec3d down(1.0/sqrt(2.0),-1.0/sqrt(2.0),0.0);
  std::vector<Vec3d> vertices;
  vertices.push_back(Vec3d(0.0,0.0,0.0));
  for( int i = 1; i < opts.numVertices; ++i ) 
  {
    if( i%2 == 1 ) vertices.push_back( vertices[i-1]+up );
    else vertices.push_back( vertices[i-1]+down );
  }
  
  ElasticRod* newrod = setupRod(opts, vertices, undeformed);
  RodBoundaryCondition* boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0);
  m_world->addObject(newrod);
  m_rods.push_back(newrod);
  
  RodTimeStepper* newstepper = getRodTimeStepper(*newrod);
  m_world->addController(newstepper);
  m_steppers.push_back(newstepper);
  
  std::cout << "Solver Options: " << std::endl;
  std::cout << "   integrator: " << GetStringOpt("integrator") << std::endl;
  std::cout << "   dt:         " << GetScalarOpt("dt") << std::endl;
  std::cout << "   iterations: " << GetIntOpt("iterations") << std::endl;
  std::cout << "   solver: " << SolverUtils::instance()->getSolverName() << std::endl;
  
  std::cout << "\"World\" Options: " << std::endl;
  std::cout << "   mass-damping: " << GetScalarOpt("mass-damping") << std::endl;
  
  std::cout << "Rod Options: " << std::endl;
  std::cout << "   numVertices:   " << opts.numVertices << std::endl;
  std::cout << "   density:       " << opts.density << std::endl;
  std::cout << "   radiusA:       " << opts.radiusA << std::endl;
  std::cout << "   radiusB:       " << opts.radiusB << std::endl;
  std::cout << "   radiusScale:   " << opts.radiusScale << std::endl;
  std::cout << "   YoungsModulus: " << opts.YoungsModulus << std::endl;
  std::cout << "   ShearModulus:  " << opts.ShearModulus << std::endl;
  std::cout << "   viscosity:     " << opts.viscosity << std::endl;
  std::cout << "   anisotropic:   " << opts.anisotropic << std::endl;
  std::cout << "   elastic:       " << opts.elastic << std::endl;
  std::cout << "   quasistatic:   " << opts.quasistatic << std::endl;
  std::cout << "   inextensible:  " << opts.inextensible << std::endl;
  std::cout << "   reframe:       ";
  switch(opts.refFrame) 
  {
    case ElasticRod::TimeParallel:
      std::cout << "TimeParallel" << std::endl;
      break;
    case ElasticRod::SpaceParallel:
      std::cout << "SpaceParallel" << std::endl;
      break;
    default:
      std::cout << "WARNING NONHANDLED CASE" << std::endl;
      break;
  }
}

void RegressionTests::atEachTimestepStiffBendingFixed()
{
  // No scripting
}




/////////////////////////////////////////////////////////////
// TESTS THAT HAVE BEEN FIXED

void RegressionTests::generateTwoVertexBendingIsometry()
{
  loadDynamicsProps();
  
  // Ensure that gravity is disabled
  this->setGravity(Vec3d(0.0,0.0,0.0));
  
  // Use a small timestep, large number of iterations to help the Newton Solver
  this->setDt(0.0001);
  GetIntOpt("iterations") = 1000;
  
  RodOptions opts;
  getRodOptions(opts);
  // Scale the radius of the rod for rendering/contact
  opts.radiusScale = 10.0;
  
  // Simplest rod with bending has 3 verts
  opts.numVertices = 3;
  
  // A normal in the -z direction
  Vec3d nhat(0,0,-1);
  
  // Rotate the normal (Rotation matrix generated in Mathematica)
  Mat3d rotmat;
  rotmat << 0.921443, 0.0518094, -0.385043, 0.0518094, 0.965831, 0.253942, 0.385043, -0.253942, 0.887274;
  assert( approxEq(rotmat.determinant(),1.0,1.0e-6) );
  nhat = rotmat*nhat; // Comment out this line to create a working example
  assert( approxEq(nhat.norm(),1.0,1.0e-6) );
  
  // Generate a naturally straight rod with three vertices
  std::vector<Vec3d> vertices;
  vertices.push_back(nhat);
  vertices.push_back(nhat+nhat);
  vertices.push_back(nhat+nhat+nhat);
  assert( vertices.size() == (size_t) opts.numVertices );
  
  ElasticRod* newrod = setupRod(opts, vertices, vertices);
  RodBoundaryCondition* boundary = newrod->getBoundaryCondition();
  m_rods.push_back(newrod);
  
  // Fix the first two vertices
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  // Fix the angle of the first material frame
  boundary->setDesiredEdgeAngle(0,0.0);
  
  // Add the rod to the world for rendering, etc
  m_world->addObject(newrod);
  
  // Create a timestepper to evolve the rod in time
  RodTimeStepper* tstep = getRodTimeStepper(*newrod);
  m_world->addController(tstep);
  m_steppers.push_back(tstep);
}

void RegressionTests::atEachTimestepGenerateTwoVertexBendingIsometry()
{
}

void RegressionTests::generatePullTest()
{
  // Set Defaults for built-in options
  GetStringOpt("integrator") = "implicit";
  GetVecOpt("gravity") = Vec3d(0,-981,0);
  GetScalarOpt("dt") = 0.1;
  // Maximum number of implicit solver itations
  GetIntOpt("iterations") = 1000; 
  GetIntOpt("nv") = 40;
  GetScalarOpt("mass-damping") = 0;
  GetScalarOpt("viscosity") = 0;
  
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;

  // 12 inches == 60.96 cm
  double L = 60.96/3.0;  
  double dh = L/((double)opts.numVertices-1);

  Vec3d baseposition(0.0,2.0,0.0);
  std::vector<Vec3d> vertices;

  for( int k = 0; k < opts.numVertices; ++k )
  {
    // Swap these lines to switch from straight to helical undeformed
    Vec3d vert = (baseposition+Vec3d(0,-k*dh,0));
    //Vec3d vert = (baseposition+Vec3d(sin(-k*dh),-k*dh,cos(-k*dh)));
    vertices.push_back(vert);
  }
  
  ElasticRod* newrod = setupRod(opts, vertices, vertices);
  RodBoundaryCondition* boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  
  m_world->addObject(newrod);
  m_rods.push_back(newrod);
  
  RodTimeStepper* newstepper = getRodTimeStepper(*newrod);
  m_world->addController(newstepper);
  m_steppers.push_back(newstepper);
}

void RegressionTests::atEachTimestepPullTest()
{
  // cm / second
  double pull_rate = 5.0;
  
  if( getTime() <= 7.5 )
  {  
    Vec3d dx = Vec3d(pull_rate*GetScalarOpt("dt"),0.0,0.0);
    for( int i = 0; i < (int) m_rods.size(); ++i ) 
    {
      m_rods[i]->setVertex(0,m_rods[i]->getVertex(0)+dx);
      m_rods[i]->setVertex(1,m_rods[i]->getVertex(1)+dx);
      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(0,m_rods[i]->getVertex(0));
      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(1,m_rods[i]->getVertex(1));
      m_rods[i]->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
    }
  }
}



/////////////////////////////////////////////////////////////
// EXPECTED BEHAVIOR?

void RegressionTests::veryStrainedRodTest()
{
  GetVecOpt("gravity") = Vec3d(0,-981,0); 
  
	GetIntOpt("nv") = 60;
	GetIntOpt("iterations") = 1000;
  
	GetStringOpt("integrator") = "implicit"; 
  GetScalarOpt("dt") = 0.01;
  
  GetBoolOpt("quasistatic") = true;
	GetScalarOpt("density") = 1.0;
  
	GetScalarOpt("youngs-modulus") = 1e6;
	GetScalarOpt("shear-modulus") = 1e6;
	GetScalarOpt("major-radius") = 0.05;
	GetScalarOpt("minor-radius") = 0.05;

  GetScalarOpt("viscosity") = 0.0;
  
  
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  
  
  //Scalar radius = GetScalarOpt("major-radius");
  //Scalar rho = GetScalarOpt("density");
  Scalar m_rodLength = 1.5;
  int nv = opts.numVertices;
  
  Scalar edgeLen = m_rodLength; //min(m_fallHeight.getValue() / 20.0, m_rodLength);
  
	std::vector<Vec3d> initial_vertices;
	std::vector<Vec3d> undeformed_vertices;
  
	for (int i = 0; i < opts.numVertices; ++i) {
		undeformed_vertices.push_back(Vec3d(0, 20 + (Scalar)(opts.numVertices-1 - i) * edgeLen, 0));
	}
	
	initial_vertices = undeformed_vertices;
  
  for (int i = 4; i < nv/2; ++i) {
  	Vec3d v0 = initial_vertices[i-1];
  	Vec3d v = initial_vertices[i];
  	double r = 1.5;
  	v = Vec3d((Scalar)(i-4) *1.5, 20 + (Scalar)(nv-1 - i) * edgeLen + (Scalar)(i-4) * (Scalar)(i-4) * 1, 0);
  	
  	v = (v - v0) / ((v - v0).norm())*r + v0;
  	
  	initial_vertices[i] = v;
  }
  for (int i = nv/2; i < nv; ++i) {
  	Vec3d v0 = initial_vertices[i-1];
  	Vec3d v = initial_vertices[i];
  	double r = 1.5;
  	v = Vec3d((Scalar)(i-4) *1.5, 20 + (Scalar)(nv-1 - i) * edgeLen + (Scalar)(i-4) * (Scalar)(i-4) * 1 - (Scalar)(i-4) * (Scalar)(i-4) * 1.3, 0);
  	
  	v = (v - v0) / ((v - v0).norm())*r + v0;
  	
  	initial_vertices[i] = v;
  }
  
  ElasticRod* m_rod = setupRod(opts, initial_vertices, undeformed_vertices);
  m_rods.push_back(m_rod);
  
  
	RodTimeStepper* m_stepper = getRodTimeStepper(*m_rod);
  m_steppers.push_back(m_stepper);
  
	RodBoundaryCondition* boundary = m_rod->getBoundaryCondition();
	
	boundary->setDesiredVertexPosition(0, m_rod->getVertex(0));
	boundary->setDesiredVertexPosition(1, m_rod->getVertex(1));
	boundary->setDesiredEdgeAngle(0, m_rod->getTheta(0));
  
	m_world->addObject(m_rod);
	m_world->addController(m_stepper);	
  
}

void RegressionTests::atEachTimestepVeryStrainedRodTest()
{
}


void RegressionTests::Setup()
{
  switch (GetIntOpt("RegressionTest")) 
  {
    case 0:
    {
      std::cout << "Executing test: generateStiffBendingNonFixed() [KNOWN BUG, UNRESOLVED]" << std::endl;
      generateStiffBendingNonFixed();
      break;
    }
    case 1:
    {
      std::cout << "Executing test: generateStiffBendingFixed() [KNOWN BUG, UNRESOLVED FOR CG]" << std::endl;
      generateStiffBendingFixed();
      break;
    }
    case 2:
    {
      std::cout << "Executing test: generatePullTest() [RESOLVED]" << std::endl;
      generatePullTest();
      break;
    }
    case 3:
    {
      std::cout << "Executing test: generateTwoVertexBendingIsometry() [RESOLVED]" << std::endl;
      generateTwoVertexBendingIsometry();
      break;
    }
    case 4:
    {
      std::cout << "Executing test: veryStrainedRodTest() [EXPECTED BEHAVIOR?]" << std::endl;
      veryStrainedRodTest();
      break;
    }
    case 5:
    {
      std::cout << "Executing test: generateOutOfPlaneBuckling() [UNRESOLVED]" << std::endl;
      generateOutOfPlaneBuckling();
      break;
    }
    default:
    {
      std::cout << "UNKOWN TEST CASE: " << GetIntOpt("RegressionTest") << std::endl;
      assert( false );
      break;
    }
  }
}


void RegressionTests::AtEachTimestep()
{
  switch (GetIntOpt("RegressionTest")) 
  {
    case 0:
    {
      atEachTimestepStiffBendingNonFixed();
      break;
    }
    case 1:
    {
      atEachTimestepStiffBendingNonFixed();
      break;
    }
    case 2:
    {
      atEachTimestepPullTest();
      break;
    }
    case 3:
    {
      atEachTimestepGenerateTwoVertexBendingIsometry();
      break;
    }
    case 4:
    {
      atEachTimestepVeryStrainedRodTest();
      break;
    }
    case 5:
    {
      atEachTimestepOutOfPlaneBuckling();
      break;
    }
    default:
    {
      std::cout << "UNKOWN TEST CASE: " << GetIntOpt("RegressionTest") << std::endl;
      assert( false );
      break;
    }
  }
}




















//void ImplicitEulerTests::curlSupportTest()
//{
//  GetStringOpt("integrator") = "implicit";
//  GetIntOpt("nv") = 55;
//  GetIntOpt("iterations") = 1000;
//  GetScalarOpt("mass-damping") = 0.01;
//  //GetScalarOpt("viscosity") = 0;
//  GetVecOpt("gravity") = Vec3d(0.0,-981,0.0);
//  GetScalarOpt("dt") = 0.01;
//
//  loadDynamicsProps();
//
//  RodOptions opts;
//  getRodOptions(opts);
//  opts.radiusScale = 10.0;
//
//  // 12 inches == 60.96 cm
//  double L = 60.96/3.0;
//  double dL = L/((double)(opts.numVertices-1));
//
//  // "Radii" of the helix
//  double hrA = 60.96/40.0;
//  double hrB = 60.96/40.0;
//  
//  ///////////////// ISOTROPIC ROD
//  std::vector<Vec3d> undeformed;
//  for( int i = 0; i < opts.numVertices; ++i ) undeformed.push_back(Vec3d(hrA*cos((double)-i*dL),(double)-i*dL,hrB*sin((double)-i*dL)));
//
//  // vertices, undeformed
//  ElasticRod* newrod = setupRod(opts, undeformed, undeformed);
//  RodBoundaryCondition* boundary = newrod->getBoundaryCondition();
//  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
//  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
//  boundary->setDesiredEdgeAngle(0,0);
//  m_world->addObject(newrod);
//  
//  RodTimeStepper* newstepper = getRodTimeStepper(*newrod);
//  m_world->addController(newstepper);
//  
//  
//  ///////////////// SLIGHTLY ANISOTROPIC ROD
//  opts.radiusB = 2.0*opts.radiusA;
//  
//  undeformed.clear();
//  for( int i = 0; i < opts.numVertices; ++i ) undeformed.push_back(Vec3d(5.0+hrA*cos((double)-i*dL),(double)-i*dL,hrB*sin((double)-i*dL)));
//  
//  // vertices, undeformed
//  newrod = setupRod(opts, undeformed, undeformed);
//  boundary = newrod->getBoundaryCondition();
//  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
//  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
//  boundary->setDesiredEdgeAngle(0,0);
//  m_world->addObject(newrod);
//  
//  newstepper = getRodTimeStepper(*newrod);
//  m_world->addController(newstepper);
//
//
//  ///////////////// MORE ANISOTROPIC ROD
//  opts.radiusB = 4.0*opts.radiusA;
//  
//  undeformed.clear();
//  for( int i = 0; i < opts.numVertices; ++i ) undeformed.push_back(Vec3d(10.0+hrA*cos((double)-i*dL),(double)-i*dL,hrB*sin((double)-i*dL)));
//  
//  // vertices, undeformed
//  newrod = setupRod(opts, undeformed, undeformed);
//  boundary = newrod->getBoundaryCondition();
//  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
//  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
//  boundary->setDesiredEdgeAngle(0,0);
//  m_world->addObject(newrod);
//  
//  newstepper = getRodTimeStepper(*newrod);
//  m_world->addController(newstepper);
//  
//  
//  ///////////////// EVEN MORE ANISOTROPIC ROD
//  opts.radiusB = 8.0*opts.radiusA;
//  
//  undeformed.clear();
//  for( int i = 0; i < opts.numVertices; ++i ) undeformed.push_back(Vec3d(15.0+hrA*cos((double)-i*dL),(double)-i*dL,hrB*sin((double)-i*dL)));
//  
//  // vertices, undeformed
//  newrod = setupRod(opts, undeformed, undeformed);
//  boundary = newrod->getBoundaryCondition();
//  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
//  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
//  boundary->setDesiredEdgeAngle(0,0);
//  m_world->addObject(newrod);
//  
//  newstepper = getRodTimeStepper(*newrod);
//  m_world->addController(newstepper);
//}



