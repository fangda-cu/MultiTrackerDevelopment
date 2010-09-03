/**
 * \file MicrotuboliDNA.cc
 *
 * \author smith@cs.columbia.edu
 * \date 08/06/2010
 */

#ifdef HAVE_PARDISO

#include "MicrotuboliDNA.hh"

MicrotuboliDNA::MicrotuboliDNA()
: Problem("Microtuboli DNA", "Simulation of Microtuboli DNA")
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
  GetIntOpt("nv") = 50; // Number of vertices per rod
  GetScalarOpt("mass-damping") = 0;
  GetScalarOpt("viscosity") = 0;
  GetBoolOpt("quasistatic") = true; // Quasistatic or dynamic material frame. Dynamic is buggy as of 08/06/2010.
  GetScalarOpt("major-radius") = 1.0/200.0;
  GetScalarOpt("minor-radius") = 1.0/200.0;
  GetScalarOpt("radius-scale") = 1.0; // Scales radius for rendering
  
  // Add new options
  AddOption("rod-length", "length of each rod", 3.0/4.0);
  AddOption("square-side-length", "length and width of each square", 1.0);
  AddOption("plate-distance", "distance between the parallel plates", 1.0);
  
  AddOption("hook-spring-stiffness", "stiffness of Hook springs", 1.0);
  AddOption("hook-spring-rest-length", "rest length of Hook springs", 0.5);
}

MicrotuboliDNA::~MicrotuboliDNA()
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



ElasticRod* MicrotuboliDNA::createRootedLowerRod( const Vec3d& rootlocation, const RodOptions& opts )
{
  double L = GetScalarOpt("rod-length");
  double dL = L/((double)(opts.numVertices-1));

  // Create the undeformed configuration
  Vec3d dy = dL*Vec3d(0.0,1.0,0.0);
  std::vector<Vec3d> undeformed;
  undeformed.push_back(Vec3d::Zero());
  for( int i = 0; i < opts.numVertices-1; ++i ) undeformed.push_back(undeformed.back()+dy);

  // Create the deformed configuration
  std::vector<Vec3d> deformed;
  deformed.push_back(rootlocation);
  for( int i = 0; i < opts.numVertices-1; ++i ) deformed.push_back(deformed.back()+dy);
  
  // Create the rod
  ElasticRod* newrod = setupRod(opts, deformed, undeformed);
  RodBoundaryCondition* boundary = newrod->getBoundaryCondition();

  // Set the boundary conditions
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
//  boundary->setDesiredVertexPosition(newrod->nv()-1, newrod->getVertex(newrod->nv()-1));
//  boundary->setDesiredVertexPosition(newrod->nv()-2, newrod->getVertex(newrod->nv()-2));
//  boundary->setDesiredEdgeAngle(newrod->ne()-1,0.0);
  
  return newrod;
}

ElasticRod* MicrotuboliDNA::createRootedUpperRod( const Vec3d& rootlocation, const RodOptions& opts )
{
  double L = GetScalarOpt("rod-length");
  double dL = L/((double)(opts.numVertices-1));
  
  // Create the undeformed configuration
  Vec3d dy = dL*Vec3d(0.0,1.0,0.0);
  std::vector<Vec3d> undeformed;
  undeformed.push_back(Vec3d::Zero());
  for( int i = 0; i < opts.numVertices-1; ++i ) undeformed.push_back(undeformed.back()+dy);
  
  // Create the deformed configuration
  std::vector<Vec3d> deformed;
  deformed.push_back(rootlocation);
  for( int i = 0; i < opts.numVertices-1; ++i ) deformed.push_back(deformed.back()-dy);
  
  // Create the rod
  ElasticRod* newrod = setupRod(opts, deformed, undeformed);
  RodBoundaryCondition* boundary = newrod->getBoundaryCondition();
  
  // Set the boundary conditions
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  //  boundary->setDesiredVertexPosition(newrod->nv()-1, newrod->getVertex(newrod->nv()-1));
  //  boundary->setDesiredVertexPosition(newrod->nv()-2, newrod->getVertex(newrod->nv()-2));
  //  boundary->setDesiredEdgeAngle(newrod->ne()-1,0.0);
  
  return newrod;
}


void MicrotuboliDNA::Setup()
{
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);

  // Create the upper-plate rods
  m_rods.push_back( createRootedLowerRod( Vec3d(0.0,0.0,0.0), opts ) );
  m_rods.push_back( createRootedLowerRod( Vec3d(GetScalarOpt("square-side-length"),0.0,0.0), opts ) );
  m_rods.push_back( createRootedLowerRod( Vec3d(GetScalarOpt("square-side-length"),0.0,-GetScalarOpt("square-side-length")), opts ) );
  m_rods.push_back( createRootedLowerRod( Vec3d(0.0,0.0,-GetScalarOpt("square-side-length")), opts ) );

  // Create the lower-plate rods
  m_rods.push_back( createRootedUpperRod( Vec3d(0.5*GetScalarOpt("square-side-length"),GetScalarOpt("plate-distance"),-0.5*GetScalarOpt("square-side-length")), opts ) );
  m_rods.push_back( createRootedUpperRod( Vec3d(GetScalarOpt("square-side-length")+0.5*GetScalarOpt("square-side-length"),GetScalarOpt("plate-distance"),-0.5*GetScalarOpt("square-side-length")), opts ) );
  m_rods.push_back( createRootedUpperRod( Vec3d(GetScalarOpt("square-side-length")+0.5*GetScalarOpt("square-side-length"),GetScalarOpt("plate-distance"),-GetScalarOpt("square-side-length")-0.5*GetScalarOpt("square-side-length")), opts ) );
  m_rods.push_back( createRootedUpperRod( Vec3d(0.5*GetScalarOpt("square-side-length"),GetScalarOpt("plate-distance"),-GetScalarOpt("square-side-length")-0.5*GetScalarOpt("square-side-length")), opts ) );


  // Create the multiple rod time stepper
  m_multiple_stepper = getMultipleRodTimeStepper();
  m_world->addController(m_multiple_stepper);

  // Add all rods to the world, to the rod time stepper
  for( size_t i = 0; i < m_rods.size(); ++i ) m_world->addObject(m_rods[i]);
  for( int i = 0; i < (int) m_rods.size(); ++i ) m_multiple_stepper->addRod(m_rods[i]);

  // The springs won't be EXACTLY 2/3 unless we choose nv correctly
  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 0, 1, 2*opts.numVertices/3, 2*opts.numVertices/3, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );
  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 1, 2, 2*opts.numVertices/3, 2*opts.numVertices/3, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );
  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 2, 3, 2*opts.numVertices/3, 2*opts.numVertices/3, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );
  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 3, 0, 2*opts.numVertices/3, 2*opts.numVertices/3, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );

  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 4, 5, 2*opts.numVertices/3, 2*opts.numVertices/3, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );
  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 5, 6, 2*opts.numVertices/3, 2*opts.numVertices/3, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );
  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 6, 7, 2*opts.numVertices/3, 2*opts.numVertices/3, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );
  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 7, 4, 2*opts.numVertices/3, 2*opts.numVertices/3, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );

  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 0, 4, opts.numVertices/2, opts.numVertices/2, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );
  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 1, 4, opts.numVertices/2, opts.numVertices/2, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );
  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 2, 4, opts.numVertices/2, opts.numVertices/2, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );
  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 3, 4, opts.numVertices/2, opts.numVertices/2, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );

  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 2, 5, opts.numVertices/2, opts.numVertices/2, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );
  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 2, 6, opts.numVertices/2, opts.numVertices/2, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );
  m_rod_rod_springs.push_back( new RodRodSpringForce( m_rods, 2, 7, opts.numVertices/2, opts.numVertices/2, GetScalarOpt("hook-spring-stiffness"), GetScalarOpt("hook-spring-rest-length") ) );

  for( size_t i = 0; i < m_rod_rod_springs.size(); ++i ) m_world->addRenderer(new SpringRenderer(*m_rod_rod_springs[i]));
  for( size_t i = 0; i < m_rod_rod_springs.size(); ++i ) m_multiple_stepper->addRodRodExternalForce(m_rod_rod_springs[i]);
}

void MicrotuboliDNA::AtEachTimestep()
{
}

#endif




