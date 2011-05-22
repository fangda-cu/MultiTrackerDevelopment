/**
 * \file CollisionTestTwo.cc
 *
 * \author smith@cs.columbia.edu
 * \date 04/27/2010
 */

#include "CollisionTestTwo.hh"

CollisionTestTwo::CollisionTestTwo()
: Problem("Tricky Collisions", "A number of tests highlighting problems with naive collision response.")
, m_rods()
, m_tri_objs()
, m_controllers()
, m_br_stepper(NULL)
, m_scripting_controllers()
{
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();
  
  // Set Defaults for built-in options
  GetStringOpt("integrator") = "implicit";
  GetVecOpt("gravity") = Vec3d(0,-981,0);
  GetScalarOpt("dt") = 0.01;
  GetIntOpt("iterations") = 1000; // Maximum number of implicit solver itations
  GetIntOpt("nv") = 40;
  GetScalarOpt("mass-damping") = 0;
  GetScalarOpt("viscosity") = 0;

  AddOption("RegressionTest", "Which regression test to execute.", 0);
}

CollisionTestTwo::~CollisionTestTwo()
{
  for( int i = 0; i < (int) m_rods.size(); ++i )
  {
    assert( m_rods[i] != NULL );
    delete m_rods[i];
    m_rods[i] = NULL;
  }
  
  for( int i = 0; i < (int) m_tri_objs.size(); ++i )
  {
    assert( m_tri_objs[i] != NULL );
    delete m_tri_objs[i];
    m_tri_objs[i] = NULL;
  }
  
  for( int i = 0; i < (int) m_controllers.size(); ++i )
  {
    assert( m_controllers[i] != NULL );
    delete m_controllers[i];
    m_controllers[i] = NULL;
  }
  
  for( int i = 0; i < (int) m_scripting_controllers.size(); ++i )
  {
    assert( m_scripting_controllers[i] != NULL );
    delete m_scripting_controllers[i];
    m_scripting_controllers[i] = NULL;
  }
  
  if( m_br_stepper != NULL )
  {
    delete m_br_stepper;
    m_br_stepper = NULL;
  }
}

void CollisionTestTwo::Setup()
{
  switch (GetIntOpt("RegressionTest")) 
  {
    case 0:
    {
      std::cout << "Executing test: RodObjectSnaggingSetup()" << std::endl;
      RodObjectSnaggingSetup();
      break;
    }
    case 1:
    {
      std::cout << "Executing test: RodObjectSnaggingMovingObjectSetup()" << std::endl;
      RodObjectSnaggingMovingObjectSetup();
      break;
    }
    case 2:
    {
      std::cout << "Executing test: RodObjectSnaggingVertFaceSetup()" << std::endl;
      RodObjectSnaggingVertFaceSetup();
      break;
    }
    case 3:
    {
      std::cout << "Executing test: RodObjectSnaggingVertFaceTwoSetup()" << std::endl;
      RodObjectSnaggingVertFaceTwoSetup();
      break;
    }
    case 4:
    {
      std::cout << "Executing test: RodRodSnaggingSetup()" << std::endl;
      RodRodSnaggingSetup();
      break;
    }
    case 5:
    {
      std::cout << "Executing test: RodRodSnaggingTwoSetup()" << std::endl;
      RodRodSnaggingTwoSetup();
      break;
    }
    case 6:
    {
      std::cout << "Executing test: RodRodSnaggingThreeSetup()" << std::endl;
      RodRodSnaggingThreeSetup();
      break;
    }
    case 7:
    {
      std::cout << "Executing test: RodFixedRodSnaggingSetup()" << std::endl;
      RodFixedRodSnaggingSetup();
      break;
    }
    case 8:
    {
      std::cout << "Executing test: RodRodSnaggingSmallSetup()" << std::endl;
      RodRodSnaggingSmallSetup();
      break;
    }
    case 9:
    {
      std::cout << "Executing test: RodRodSnaggingDifferentSizeSetup()" << std::endl;
      RodRodSnaggingDifferentSizeSetup();
      break;
    }
    default:
    {
      std::cerr << "Invalid RegressionTest setting, exiting." << std::endl;
      exit(0);
      break;
    }
  }
}

void CollisionTestTwo::AtEachTimestep()
{
}


///////////////////////////////////////////////////////////////////////////////
// ROD OBJECT SNAGGING

RodObjectSnaggingController::RodObjectSnaggingController( ElasticRod& rod, double time, double dt )
: ScriptingController(time,dt)
, m_rod(rod)
{}
  
bool RodObjectSnaggingController::execute()
{
  // cm / second
  double pull_rate = 5.0;
  
  if( ScriptingController::getTime() <= 25.0 )
  {  
    Vec3d dx = Vec3d(pull_rate*ScriptingController::getDt(),0.0,0.0);

    m_rod.setVertex(0,m_rod.getVertex(0)+dx);
    m_rod.setVertex(1,m_rod.getVertex(1)+dx);
    m_rod.getBoundaryCondition()->setDesiredVertexPosition(0,m_rod.getVertex(0));
    m_rod.getBoundaryCondition()->setDesiredVertexPosition(1,m_rod.getVertex(1));
    m_rod.getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
  }
  
  if( ScriptingController::getTime() > 7.5 ) exit(0);  

  return true;
}

void CollisionTestTwo::RodObjectSnaggingSetup()
{
  loadDynamicsProps();

  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;
  
  // 12 inches == 60.96 cm
  double L = 60.96/3.0;
  
  double dx = 0.3;
  double dy = 0.3;
  double dh = L/((double)opts.numVertices-1);
  
  int numx = 1;
  int numy = 1;
  for ( int i = 0; i < numx; ++i ) for( int j = 0; j < numy; ++j )
  {
    Vec3d baseposition(i*dx-dx*(numx/2),2.0,j*dy-dy*(numy/2));
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
    
    RodTimeStepper* stepper = getRodTimeStepper(*newrod);
    m_controllers.push_back(stepper);
  }
  
  ObjParser objparser;
  TriangleMesh* tri_mesh = new TriangleMesh();
  objparser.loadTriangularMesh( "assets/Square.obj", *tri_mesh );
  for( TriangleMesh::vertex_iter itr = tri_mesh->vertices_begin(); itr != tri_mesh->vertices_end(); ++itr )
  {
    tri_mesh->getVertex(*itr) -= Vec3d(0.0,5.0,0.0);
  }
  m_world->addObject(tri_mesh);
  m_tri_objs.push_back(tri_mesh);

  m_scripting_controllers.push_back(new RodObjectSnaggingController(*m_rods[0],getTime(),GetScalarOpt("dt")));

  m_br_stepper = new BARodStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

void CollisionTestTwo::RodObjectSnaggingAtEachTimestep()
{
}

void CollisionTestTwo::RodObjectSnaggingMovingObjectSetup()
{
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;
  
  // 12 inches == 60.96 cm
  double L = 60.96/3.0;
  
  double dx = 0.3;
  double dy = 0.3;
  double dh = L/((double)opts.numVertices-1);
  
  int numx = 1;
  int numy = 1;
  for ( int i = 0; i < numx; ++i ) for( int j = 0; j < numy; ++j )
  {
    Vec3d baseposition(i*dx-dx*(numx/2),2.0,j*dy-dy*(numy/2));
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
    
    RodTimeStepper* stepper = getRodTimeStepper(*newrod);
    m_controllers.push_back(stepper);
  }
  
  ObjParser objparser;
  TriangleMesh* tri_mesh = new TriangleMesh();
  objparser.loadTriangularMesh( "assets/Square.obj", *tri_mesh );
  for( TriangleMesh::vertex_iter itr = tri_mesh->vertices_begin(); itr != tri_mesh->vertices_end(); ++itr )
  {
    tri_mesh->getVertex(*itr) -= Vec3d(0.0,5.0,0.0);
  }
  m_world->addObject(tri_mesh);
  m_tri_objs.push_back(tri_mesh);

  m_scripting_controllers.push_back(new ObjectTranslator(*tri_mesh,getTime(),GetScalarOpt("dt")));
  
  m_br_stepper = new BARodStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

ObjectTranslator::ObjectTranslator( TriangleMesh& mesh, double time, double dt )
: ScriptingController(time,dt)
, m_mesh(mesh)
{}

bool ObjectTranslator::execute()
{
  double distance = 20.0;
  double transittime = 10.0;
  
  double speed = distance/transittime;
  double translatedist = speed*getDt();
  
  for( TriangleMesh::vertex_iter itr = m_mesh.vertices_begin(); itr != m_mesh.vertices_end(); ++itr )
  {
    m_mesh.getVertex(*itr) += Vec3d(-translatedist,0.0,0.0);
  }
  
  if( getTime() > 7.5 ) exit(0);
  
  return true;
}

void CollisionTestTwo::RodObjectSnaggingMovingObjectAtEachTimestep()
{
}

void CollisionTestTwo::RodObjectSnaggingVertFaceSetup()
{
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;
  
  // 12 inches == 60.96 cm
  double L = 60.96/3.0;
  
  double dx = 0.3;
  double dy = 0.3;
  double dh = L/((double)opts.numVertices-1);
  
  int numx = 1;
  int numy = 1;
  for ( int i = 0; i < numx; ++i ) for( int j = 0; j < numy; ++j )
  {
    Vec3d baseposition(i*dx-dx*(numx/2)+6.0,13.0,-1.5);
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
    
    RodTimeStepper* stepper = getRodTimeStepper(*newrod);
    m_controllers.push_back(stepper);
  }
  
  ObjParser objparser;
  TriangleMesh* tri_mesh = new TriangleMesh();
  objparser.loadTriangularMesh( "assets/Square.obj", *tri_mesh );
  for( TriangleMesh::vertex_iter itr = tri_mesh->vertices_begin(); itr != tri_mesh->vertices_end(); ++itr )
  {
    tri_mesh->getVertex(*itr) -= Vec3d(0.0,5.0,0.0);
  }
  m_world->addObject(tri_mesh);
  m_tri_objs.push_back(tri_mesh);
  
  m_scripting_controllers.push_back(new ObjectTranslator(*tri_mesh,getTime(),GetScalarOpt("dt")));
  
  m_br_stepper = new BARodStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

void CollisionTestTwo::RodObjectSnaggingVertFaceTwoSetup()
{
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;
  
  // 12 inches == 60.96 cm
  double L = 60.96/3.0;
  
  double dx = 0.3;
  double dy = 0.3;
  double dh = L/((double)opts.numVertices-1);
  
  int numx = 1;
  int numy = 1;
  for ( int i = 0; i < numx; ++i ) for( int j = 0; j < numy; ++j )
  {
    Vec3d baseposition(i*dx-dx*(numx/2)+6.0,13.0,-1.5);
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
    
    RodTimeStepper* stepper = getRodTimeStepper(*newrod);
    m_controllers.push_back(stepper);
  }

  ObjParser objparser;
  TriangleMesh* tri_mesh = new TriangleMesh();
  objparser.loadTriangularMesh( "assets/Square.obj", *tri_mesh );
  for( TriangleMesh::vertex_iter itr = tri_mesh->vertices_begin(); itr != tri_mesh->vertices_end(); ++itr )
  {
    tri_mesh->getVertex(*itr) -= Vec3d(0.0,5.0,0.0);
  }
  m_world->addObject(tri_mesh);
  m_tri_objs.push_back(tri_mesh);
  
  m_scripting_controllers.push_back(new RodObjectSnaggingController(*m_rods[0],getTime(),GetScalarOpt("dt")));

  m_br_stepper = new BARodStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}


/////////////////////////////////////////////////////////////////////////////
// ROD-ROD SNAGGING

RodRodSnaggingController::RodRodSnaggingController( ElasticRod& rod, double time, double dt )
: ScriptingController(time,dt)
, m_rod(rod)
{}

bool RodRodSnaggingController::execute()
{
  // cm / second
  double pull_rate = 5.0;
  
  if( ScriptingController::getTime() <= 25.0 )
  {  
    Vec3d dx = Vec3d(0.0,0.0,-pull_rate*ScriptingController::getDt());
    
    m_rod.setVertex(0,m_rod.getVertex(0)+dx);
    m_rod.setVertex(1,m_rod.getVertex(1)+dx);
    m_rod.getBoundaryCondition()->setDesiredVertexPosition(0,m_rod.getVertex(0));
    m_rod.getBoundaryCondition()->setDesiredVertexPosition(1,m_rod.getVertex(1));
    m_rod.getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
  }
  
  if( ScriptingController::getTime() > 7.5 ) exit(0);  
  
  return true;
}

void CollisionTestTwo::RodRodSnaggingSetup()
{
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;
  
  // 12 inches == 60.96 cm
  double L = 60.96/3.0;
  
  double dx = 0.3;
  double dy = 0.3;
  double dh = L/((double)opts.numVertices-1);
  
  // Create a rod that we will pull
  Vec3d baseposition(0.0,2.0,0.0);
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
  m_controllers.push_back(getRodTimeStepper(*newrod));

  
  // Create a rod for the pulled rod to collide against
  baseposition = Vec3d(-L/2.0,-5.0,-4.0);
  vertices.clear();
  for( int k = 0; k < opts.numVertices; ++k )
  {
    Vec3d vert = (baseposition+Vec3d(k*dh,0.0,0.0));
    //Vec3d vert = (baseposition+Vec3d(sin(-k*dh),-k*dh,cos(-k*dh)));
    vertices.push_back(vert);
  }
  assert( opts.numVertices == (int) vertices.size() );

  newrod = setupRod(opts, vertices, vertices);
  boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  boundary->setDesiredVertexPosition(newrod->nv()-1, newrod->getVertex(newrod->nv()-1));
  boundary->setDesiredVertexPosition(newrod->nv()-2, newrod->getVertex(newrod->nv()-2));
  boundary->setDesiredEdgeAngle(newrod->ne()-1,0.0);
  
  m_rods.push_back(newrod);
  m_controllers.push_back(getRodTimeStepper(*newrod));
  
  for( int i = 0; i < (int) m_rods.size(); ++i ) m_world->addObject(m_rods[i]);
  
  m_scripting_controllers.push_back(new RodRodSnaggingController(*m_rods[0],getTime(),GetScalarOpt("dt")));
  
  m_br_stepper = new BARodStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

void CollisionTestTwo::RodRodSnaggingAtEachTimestep()
{
//  // cm / second
//  double pull_rate = 5.0;
//  
//  if( getTime() <= 25.0 )
//  {  
//    Vec3d dz = Vec3d(0.0,0.0,-pull_rate*GetScalarOpt("dt"));
//    for( int i = 0; i < (int) 1; ++i ) 
//    {
//      m_rods[i]->setVertex(0,m_rods[i]->getVertex(0)+dz);
//      m_rods[i]->setVertex(1,m_rods[i]->getVertex(1)+dz);
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(0,m_rods[i]->getVertex(0));
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(1,m_rods[i]->getVertex(1));
//      m_rods[i]->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
//    }
//  }
//  
//  if( getTime() > 7.5 ) exit(0);  
}

void CollisionTestTwo::RodRodSnaggingTwoSetup()
{
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;
  
  // 12 inches == 60.96 cm
  double L = 60.96/3.0;
  
  double dx = 0.3;
  double dy = 0.3;
  double dh = L/((double)opts.numVertices-1);
  
  // Create a rod that we will pull
  Vec3d baseposition(-L/2.0+dh/2.0,2.0,0.0);
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
  m_controllers.push_back(getRodTimeStepper(*newrod));
  
  
  // Create a rod for the pulled rod to collide against
  baseposition = Vec3d(-L/2.0,-5.0,-4.0);
  vertices.clear();
  for( int k = 0; k < opts.numVertices; ++k )
  {
    Vec3d vert = (baseposition+Vec3d(k*dh,0.0,0.0));
    //Vec3d vert = (baseposition+Vec3d(sin(-k*dh),-k*dh,cos(-k*dh)));
    vertices.push_back(vert);
  }
  assert( opts.numVertices == (int) vertices.size() );
  
  newrod = setupRod(opts, vertices, vertices);
  boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  boundary->setDesiredVertexPosition(newrod->nv()-1, newrod->getVertex(newrod->nv()-1));
  boundary->setDesiredVertexPosition(newrod->nv()-2, newrod->getVertex(newrod->nv()-2));
  boundary->setDesiredEdgeAngle(newrod->ne()-1,0.0);
  
  m_rods.push_back(newrod);
  m_controllers.push_back(getRodTimeStepper(*newrod));
  
  for( int i = 0; i < (int) m_rods.size(); ++i ) m_world->addObject(m_rods[i]);

  m_scripting_controllers.push_back(new RodRodSnaggingController(*m_rods[0],getTime(),GetScalarOpt("dt")));

  m_br_stepper = new BARodStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

void CollisionTestTwo::RodRodSnaggingTwoAtEachTimestep()
{
//  // cm / second
//  double pull_rate = 5.0;
//  
//  if( getTime() <= 25.0 )
//  {  
//    Vec3d dz = Vec3d(0.0,0.0,-pull_rate*GetScalarOpt("dt"));
//    for( int i = 0; i < (int) 1; ++i ) 
//    {
//      m_rods[i]->setVertex(0,m_rods[i]->getVertex(0)+dz);
//      m_rods[i]->setVertex(1,m_rods[i]->getVertex(1)+dz);
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(0,m_rods[i]->getVertex(0));
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(1,m_rods[i]->getVertex(1));
//      m_rods[i]->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
//    }
//  }
//  
//  if( getTime() > 7.5 ) exit(0);  
}

void CollisionTestTwo::RodRodSnaggingThreeSetup()
{
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;
  
  // 12 inches == 60.96 cm
  double L = 60.96/3.0;
  
  double dx = 0.3;
  double dy = 0.3;
  double dh = L/((double)opts.numVertices-1);
  
  // Create a rod that we will pull
  Vec3d baseposition(-L/2.0+1.5*dh,2.0,0.0);
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
  m_controllers.push_back(getRodTimeStepper(*newrod));
  
  
  // Create a rod for the pulled rod to collide against
  baseposition = Vec3d(-L/2.0,-5.0,-4.0);
  vertices.clear();
  for( int k = 0; k < opts.numVertices; ++k )
  {
    Vec3d vert = (baseposition+Vec3d(k*dh,0.0,0.0));
    //Vec3d vert = (baseposition+Vec3d(sin(-k*dh),-k*dh,cos(-k*dh)));
    vertices.push_back(vert);
  }
  assert( opts.numVertices == (int) vertices.size() );
  
  newrod = setupRod(opts, vertices, vertices);
  boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  boundary->setDesiredVertexPosition(newrod->nv()-1, newrod->getVertex(newrod->nv()-1));
  boundary->setDesiredVertexPosition(newrod->nv()-2, newrod->getVertex(newrod->nv()-2));
  boundary->setDesiredEdgeAngle(newrod->ne()-1,0.0);
  
  m_rods.push_back(newrod);
  m_controllers.push_back(getRodTimeStepper(*newrod));
  
  for( int i = 0; i < (int) m_rods.size(); ++i ) m_world->addObject(m_rods[i]);
  
  m_scripting_controllers.push_back(new RodRodSnaggingController(*m_rods[0],getTime(),GetScalarOpt("dt")));

  m_br_stepper = new BARodStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

void CollisionTestTwo::RodRodSnaggingThreeAtEachTimestep()
{
//  // cm / second
//  double pull_rate = 5.0;
//  
//  if( getTime() <= 25.0 )
//  {  
//    Vec3d dz = Vec3d(0.0,0.0,-pull_rate*GetScalarOpt("dt"));
//    for( int i = 0; i < (int) 1; ++i ) 
//    {
//      m_rods[i]->setVertex(0,m_rods[i]->getVertex(0)+dz);
//      m_rods[i]->setVertex(1,m_rods[i]->getVertex(1)+dz);
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(0,m_rods[i]->getVertex(0));
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(1,m_rods[i]->getVertex(1));
//      m_rods[i]->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
//    }
//  }
//  
//  if( getTime() > 7.5 ) exit(0);  
}

RodFixedRodSnaggingScriptingController::RodFixedRodSnaggingScriptingController( ElasticRod& rod, double time, double dt )
: ScriptingController(time,dt)
, m_rod(rod)
{}

bool RodFixedRodSnaggingScriptingController::execute()
{
  // cm / second
  double pull_rate = 5.0;
  
  if( ScriptingController::getTime() <= 25.0 )
  {  
      Vec3d dz = Vec3d(0.0,0.0,-pull_rate*ScriptingController::getDt());
      m_rod.setVertex(0,m_rod.getVertex(0)+dz);
      m_rod.setVertex(1,m_rod.getVertex(1)+dz);
      m_rod.getBoundaryCondition()->setDesiredVertexPosition(0,m_rod.getVertex(0));
      m_rod.getBoundaryCondition()->setDesiredVertexPosition(1,m_rod.getVertex(1));
      m_rod.getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
  }  
  if( ScriptingController::getTime() > 7.5 ) exit(0);  
  
  return true;
}

void CollisionTestTwo::RodFixedRodSnaggingSetup()
{
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;
  
  // 12 inches == 60.96 cm
  double L = 60.96/3.0;
  
  double dx = 0.3;
  double dy = 0.3;
  double dh = L/((double)opts.numVertices-1);
  
  // Create a rod that we will pull
  Vec3d baseposition(0.0,2.0,0.0);
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
  m_controllers.push_back(getRodTimeStepper(*newrod));
  
  
  // Create a rod for the pulled rod to collide against
  baseposition = Vec3d(-L/2.0,-5.0,-4.0);
  vertices.clear();
  for( int k = 0; k < opts.numVertices; ++k )
  {
    Vec3d vert = (baseposition+Vec3d(k*dh,0.0,0.0));
    //Vec3d vert = (baseposition+Vec3d(sin(-k*dh),-k*dh,cos(-k*dh)));
    vertices.push_back(vert);
  }
  assert( opts.numVertices == (int) vertices.size() );
  
  newrod = setupRod(opts, vertices, vertices);
  boundary = newrod->getBoundaryCondition();
  for( int i = 0; i < newrod->nv(); ++i ) boundary->setDesiredVertexPosition(i, newrod->getVertex(i));
  for( int i = 0; i < newrod->ne(); ++i ) boundary->setDesiredEdgeAngle(i,0.0);
  
  m_rods.push_back(newrod);
  m_controllers.push_back(getRodTimeStepper(*newrod));
  
  for( int i = 0; i < (int) m_rods.size(); ++i ) m_world->addObject(m_rods[i]);
  
  m_scripting_controllers.push_back(new RodFixedRodSnaggingScriptingController(*m_rods[0],getTime(),GetScalarOpt("dt")));
  
  m_br_stepper = new BARodStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

void CollisionTestTwo::RodFixedRodSnaggingAtEachTimestep()
{
//  std::cout << "Currently broken" << std::endl;
//  exit(1);
//
//  // cm / second
//  double pull_rate = 5.0;
//  
//  if( getTime() <= 25.0 )
//  {  
//    Vec3d dz = Vec3d(0.0,0.0,-pull_rate*GetScalarOpt("dt"));
//    for( int i = 0; i < (int) 1; ++i ) 
//    {
//      m_rods[i]->setVertex(0,m_rods[i]->getVertex(0)+dz);
//      m_rods[i]->setVertex(1,m_rods[i]->getVertex(1)+dz);
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(0,m_rods[i]->getVertex(0));
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(1,m_rods[i]->getVertex(1));
//      m_rods[i]->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
//    }
//  }
//  
//  if( getTime() > 7.5 ) exit(0);  
}














void CollisionTestTwo::RodRodSnaggingSmallSetup()
{
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;
  
  opts.numVertices = 8;

  // 12 inches == 60.96 cm
  double L = 60.96/3.0;
  
  double dx = 0.3;
  double dy = 0.3;
  double dh = L/((double)opts.numVertices-1);
  
  
  // Create a rod that we will pull
  Vec3d baseposition(0.0,2.0+0.5*dh,0.0);
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
  m_controllers.push_back(getRodTimeStepper(*newrod));
  
  
  // Create a rod for the pulled rod to collide against
  baseposition = Vec3d(-L/2.0,-5.0,-4.0);
  vertices.clear();
  for( int k = 0; k < opts.numVertices; ++k )
  {
    Vec3d vert = (baseposition+Vec3d(k*dh,0.0,0.0));
    //Vec3d vert = (baseposition+Vec3d(sin(-k*dh),-k*dh,cos(-k*dh)));
    vertices.push_back(vert);
  }
  assert( opts.numVertices == (int) vertices.size() );
  
  newrod = setupRod(opts, vertices, vertices);
  boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  boundary->setDesiredVertexPosition(newrod->nv()-1, newrod->getVertex(newrod->nv()-1));
  boundary->setDesiredVertexPosition(newrod->nv()-2, newrod->getVertex(newrod->nv()-2));
  boundary->setDesiredEdgeAngle(newrod->ne()-1,0.0);
  
  m_rods.push_back(newrod);
  m_controllers.push_back(getRodTimeStepper(*newrod));
  
  for( int i = 0; i < (int) m_rods.size(); ++i ) m_world->addObject(m_rods[i]);
  
  m_scripting_controllers.push_back(new RodFixedRodSnaggingScriptingController(*m_rods[0],getTime(),GetScalarOpt("dt")));
  
  m_br_stepper = new BARodStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

void CollisionTestTwo::RodRodSnaggingSmallAtEachTimestep()
{
//  std::cout << "Broken" << std::endl;
//  exit(1);
//
//  // cm / second
//  double pull_rate = 5.0;
//  
//  if( getTime() <= 25.0 )
//  {  
//    Vec3d dz = Vec3d(0.0,0.0,-pull_rate*GetScalarOpt("dt"));
//    for( int i = 0; i < (int) 1; ++i ) 
//    {
//      m_rods[i]->setVertex(0,m_rods[i]->getVertex(0)+dz);
//      m_rods[i]->setVertex(1,m_rods[i]->getVertex(1)+dz);
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(0,m_rods[i]->getVertex(0));
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(1,m_rods[i]->getVertex(1));
//      m_rods[i]->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
//    }
//  }
//  
//  if( getTime() > 7.5 ) exit(0);  
}

void CollisionTestTwo::RodRodSnaggingDifferentSizeSetup()
{
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.radiusScale = 10.0;
  
  opts.numVertices = 8;
  
  // 12 inches == 60.96 cm
  double L = 60.96/3.0;
  
  double dx = 0.3;
  double dy = 0.3;
  double dh = L/((double)opts.numVertices-1);
  
  
  // Create a rod that we will pull
  Vec3d baseposition(0.0,2.0+0.5*dh,0.0);
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
  m_controllers.push_back(getRodTimeStepper(*newrod));

  // Double the number of verts in the other rod
  opts.numVertices = 16;
  dh = L/((double)opts.numVertices-1);
  
  // Create a rod for the pulled rod to collide against
  baseposition = Vec3d(-L/2.0,-5.0,-4.0);
  vertices.clear();
  for( int k = 0; k < opts.numVertices; ++k )
  {
    Vec3d vert = (baseposition+Vec3d(k*dh,0.0,0.0));
    //Vec3d vert = (baseposition+Vec3d(sin(-k*dh),-k*dh,cos(-k*dh)));
    vertices.push_back(vert);
  }
  assert( opts.numVertices == (int) vertices.size() );
  
  newrod = setupRod(opts, vertices, vertices);
  boundary = newrod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
  boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
  boundary->setDesiredVertexPosition(newrod->nv()-1, newrod->getVertex(newrod->nv()-1));
  boundary->setDesiredVertexPosition(newrod->nv()-2, newrod->getVertex(newrod->nv()-2));
  boundary->setDesiredEdgeAngle(newrod->ne()-1,0.0);
  
  m_rods.push_back(newrod);
  m_controllers.push_back(getRodTimeStepper(*newrod));
  
  for( int i = 0; i < (int) m_rods.size(); ++i ) m_world->addObject(m_rods[i]);
  
  m_scripting_controllers.push_back(new RodFixedRodSnaggingScriptingController(*m_rods[0],getTime(),GetScalarOpt("dt")));
  
  m_br_stepper = new BARodStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}








