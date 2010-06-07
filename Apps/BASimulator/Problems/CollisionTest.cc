/**
 * \file CollisionTest.cc
 *
 * \author smith@cs.columbia.edu
 * \date 02/17/2010
 */

#include "CollisionTest.hh"

CollisionTest::CollisionTest()
:Problem("Collision Test", "Test cases to exercise the velocity-filter collision response.")
,m_rods()
,m_tri_objs()
,m_br_stepper(NULL)
{
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();
  
  // Set Defaults for built-in options
  GetStringOpt("integrator") = "implicit";
  GetVecOpt("gravity") = Vec3d(0,0,0);
  GetScalarOpt("dt") = 0.0001;
  // Maximum number of implicit solver itations
  GetIntOpt("iterations") = 500; 
  //GetScalarOpt("major-radius") = 0.1;
  //GetScalarOpt("minor-radius") = 0.1;
  GetIntOpt("nv") = 2;
  //GetScalarOpt("shear-modulus") = 100;
  //GetScalarOpt("youngs-modulus") = 500;
  GetScalarOpt("mass-damping") = 0.02;
  //GetScalarOpt("viscosity") = 0;
}

CollisionTest::~CollisionTest()
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
  
  if( m_br_stepper != NULL )
  {
    delete m_br_stepper;
    m_br_stepper = NULL;
  }
}



/////////////////////////////////////////////////////////////////////////////
// Tests of full collision system

void CollisionTest::genTorusTest()
{
  loadDynamicsProps();  
  RodOptions opts;
  getRodOptions(opts);

  ObjParser objparser;
  TriangleMesh* tri_mesh = new TriangleMesh();
  objparser.loadTriangularMesh( "assets/TriangulatedTorus.obj", *tri_mesh );
  m_world->addObject(tri_mesh);
  m_br_stepper->addTriangleMesh(tri_mesh);
  
  m_br_stepper->prepareForExecution();
}

void CollisionTest::genNullTest()
{
  loadDynamicsProps();
  RodOptions opts;
  getRodOptions(opts);
  //opts.ShearModulus = GetScalarOpt("shear-modulus");
  //opts.YoungsModulus = GetScalarOpt("youngs-modulus");
  //opts.viscosity = GetScalarOpt("viscosity");
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 40;

  m_br_stepper->enablePenaltyImpulses();
  m_br_stepper->enableIterativeInelasticImpulses();
  m_br_stepper->setEdgeEdgePenalty(200.0);
  m_br_stepper->setVertexFacePenalty(200.0);
  
  double L = 20.0;
  
  double dx = 0.3;
  double dy = 0.3;
  double dh = L/((double)opts.numVertices-1);
  
  for( int i = -1; i < 2; ++i ) for( int j = -1; j < 2; ++j )
  {
    Vec3d baseposition(i*dx,0.0,j*dy);
    std::vector<Vec3d> vertices;
    for( int k = 0; k < opts.numVertices; ++k )
    {
      Vec3d vert = (baseposition+Vec3d(0,-k*dh,0));
      vertices.push_back(vert);
    }
    
    assert( opts.numVertices == (int) vertices.size() );
    ElasticRod* newrod = setupRod(opts, vertices, vertices);
    m_rods.push_back(newrod);
    RodBoundaryCondition* boundary = newrod->getBoundaryCondition();
    boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
    
    m_br_stepper->addRod(newrod);
    
    m_world->addObject(newrod);
  }

  ObjParser objparser;
  TriangleMesh* tri_mesh = new TriangleMesh();
  objparser.loadTriangularMesh( "assets/TriangulatedTorus.obj", *tri_mesh );
  for( TriangleMesh::vertex_iter itr = tri_mesh->vertices_begin(); itr != tri_mesh->vertices_end(); ++itr )
  {
    tri_mesh->getVertex(*itr) -= Vec3d(0.0,5.0,0.0);
  }
  m_world->addObject(tri_mesh);
  m_br_stepper->addTriangleMesh(tri_mesh);

  m_br_stepper->prepareForExecution();
}

void CollisionTest::genTorusPull()
{
  loadDynamicsProps();
  RodOptions opts;
  getRodOptions(opts);
  opts.numVertices = 40;
  // Scale the radius so we can see it when rendering, but not change dynamics
  opts.radiusScale = 20.0;
  
  // 9.81 m/s^2 = 981 cm/s^2
  m_br_stepper->setGravity(Vec3d(0,-981,0));

  m_br_stepper->enablePenaltyImpulses();
  m_br_stepper->enableIterativeInelasticImpulses();
  m_br_stepper->setEdgeEdgePenalty(200.0);
  m_br_stepper->setVertexFacePenalty(200.0);
  

  // 12 inches == 60.96 cm
  double L = 60.96;

  double dx = 0.3;
  double dy = 0.3;
  double dh = L/((double)opts.numVertices-1);

  int numx = 6;
  int numy = 6;
  for ( int i = 0; i < numx; ++i ) for( int j = 0; j < numy; ++j )
  {
    Vec3d baseposition(i*dx-dx*(numx/2),0.0,j*dy-dy*(numy/2));
    std::vector<Vec3d> vertices;
    for( int k = 0; k < opts.numVertices; ++k )
    {
      Vec3d vert = (baseposition+Vec3d(0,-k*dh,0));
      vertices.push_back(vert);
    }
    assert( opts.numVertices == (int) vertices.size() );

    ElasticRod* newrod = setupRod(opts, vertices, vertices);
    RodBoundaryCondition* boundary = newrod->getBoundaryCondition();
    boundary->setDesiredVertexPosition(0, newrod->getVertex(0));

    m_rods.push_back(newrod);
    m_br_stepper->addRod(newrod);
    m_world->addObject(newrod);
  }

  ObjParser objparser;
  TriangleMesh* tri_mesh = new TriangleMesh();
  objparser.loadTriangularMesh( "assets/TriangulatedTorus.obj", *tri_mesh );
  for( TriangleMesh::vertex_iter itr = tri_mesh->vertices_begin(); itr != tri_mesh->vertices_end(); ++itr )
  {
    tri_mesh->getVertex(*itr) -= Vec3d(0.0,5.0,0.0);
  }
  m_world->addObject(tri_mesh);
  m_br_stepper->addTriangleMesh(tri_mesh);
  m_tri_objs.push_back(tri_mesh);

  m_br_stepper->prepareForExecution();
}

void CollisionTest::genTorusPullAtEachTimestep()
{
  // cm / second
  double pull_rate = 5.0;

  if( getTime() <= 25.0 )
  {  
    Vec3d dx = Vec3d(pull_rate*GetScalarOpt("dt"),0.0,0.0);
    for( int i = 0; i < (int) m_rods.size(); ++i ) 
    {
      m_rods[i]->setVertex(0,m_rods[i]->getVertex(0)+dx);
      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(0,m_rods[i]->getVertex(0));
    }
  }
  
  if( getTime() > 50.0 ) exit(0);
}



/////////////////////////////////////////////////////////////////////////////
// Rod-rod penalty tests

// Two perependicular segments colliding
void CollisionTest::genRodRodPenaltyBasic00()
{
  m_br_stepper->enablePenaltyImpulses();
  m_br_stepper->disableIterativeInelasticImpulses();
  m_br_stepper->setEdgeEdgePenalty(200.0);
  m_br_stepper->setVertexFacePenalty(200.0);
  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 2;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);
  
  // Generate first rod
  std::vector<Vec3d> vertices00;
  vertices00.push_back(Vec3d(-1.0,0.0,0.0));
  vertices00.push_back(Vec3d( 1.0,0.0,0.0));
  
  ElasticRod* newrod = setupRod(opts, vertices00, vertices00);  
  m_rods.push_back(newrod);

  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);

  // Generate second rod
  std::vector<Vec3d> vertices01;
  vertices01.push_back(Vec3d(0.0,-1.0,-2.0));
  vertices01.push_back(Vec3d(0.0,1.0,-2.0));

  newrod = setupRod(opts, vertices01, vertices01);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,0.2));
  newrod->setVelocity(1, Vec3d(0.0,0.0,0.2));

  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  m_br_stepper->prepareForExecution();
}

// Multiple perependicular segments colliding
void CollisionTest::genRodRodPenaltyBasic01()
{
  m_br_stepper->enablePenaltyImpulses();
  m_br_stepper->disableIterativeInelasticImpulses();
  m_br_stepper->setEdgeEdgePenalty(200.0);
  m_br_stepper->setVertexFacePenalty(200.0);
  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 2;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);
  
  // Generate first rod
  std::vector<Vec3d> vertices00;
  vertices00.push_back(Vec3d(-1.0,0.0,0.0));
  vertices00.push_back(Vec3d( 1.0,0.0,0.0));
  
  ElasticRod* newrod = setupRod(opts, vertices00, vertices00);
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  // Generate second rod
  std::vector<Vec3d> vertices01;
  vertices01.push_back(Vec3d(0.0,-1.0,-2.0));
  vertices01.push_back(Vec3d(0.0,1.0,-2.0));
  
  newrod = setupRod(opts, vertices01, vertices01);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,0.2));
  newrod->setVelocity(1, Vec3d(0.0,0.0,0.2));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  // Generate third rod
  std::vector<Vec3d> vertices02;
  vertices02.push_back(Vec3d(0.0,-1.0,2.0));
  vertices02.push_back(Vec3d(0.0,1.0,2.0));
  
  newrod = setupRod(opts, vertices02, vertices02);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,-0.2));
  newrod->setVelocity(1, Vec3d(0.0,0.0,-0.2));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  // Generate fourth rod
  std::vector<Vec3d> vertices03;
  vertices03.push_back(Vec3d(1,-1.0,2.0));
  vertices03.push_back(Vec3d(1,1.0,2.0));
  
  newrod = setupRod(opts, vertices03, vertices03);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,-0.2));
  newrod->setVelocity(1, Vec3d(0.0,0.0,-0.2));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);

  // Generate fifth rod
  std::vector<Vec3d> vertices05;
  vertices05.push_back(Vec3d(-1,-1,0));
  vertices05.push_back(Vec3d(1,-1,0));
  
  newrod = setupRod(opts, vertices05, vertices05);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,-0.2));
  newrod->setVelocity(1, Vec3d(0.0,0.0,-0.2));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);

  m_br_stepper->prepareForExecution();
}

// "Corner" cases (mostly in detection, but this in turn influences response)
void CollisionTest::genRodRodPenaltyCorner00()
{
  m_br_stepper->enablePenaltyImpulses();
  m_br_stepper->disableIterativeInelasticImpulses();
  m_br_stepper->setEdgeEdgePenalty(200.0);
  m_br_stepper->setVertexFacePenalty(200.0);
  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 2;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);

  // Generate first rod
  std::vector<Vec3d> vertices00;
  vertices00.push_back(Vec3d(-1.0,0.0,0.0));
  vertices00.push_back(Vec3d( 1.0,0.0,0.0));
  
  ElasticRod* newrod = setupRod(opts, vertices00, vertices00);
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);

  // Generate second rod
  std::vector<Vec3d> vertices01;
  vertices01.push_back(Vec3d(-1.0,0.0,-3.0));
  vertices01.push_back(Vec3d(1.0,0.0,-3.0));
  
  newrod = setupRod(opts, vertices01, vertices01);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,0.2));
  newrod->setVelocity(1, Vec3d(0.0,0.0,0.2));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);

  // Generate third rod
  std::vector<Vec3d> vertices02;
  vertices02.push_back(Vec3d(-1.0,1.0,-3.0));
  vertices02.push_back(Vec3d(1.0,1.0,-3.0));
  
  newrod = setupRod(opts, vertices02, vertices02);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,0.0));
  newrod->setVelocity(1, Vec3d(0.0,0.0,0.0));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);

  // Generate fourth rod
  std::vector<Vec3d> vertices03;
  vertices03.push_back(Vec3d(0.0,1.0,1.0));
  vertices03.push_back(Vec3d(0.0,1.0,-1.0));
  
  newrod = setupRod(opts, vertices03, vertices03);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,-0.2));
  newrod->setVelocity(1, Vec3d(0.0,0.0,-0.2));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);

  // Generate fifth rod
  std::vector<Vec3d> vertices04;
  vertices04.push_back(Vec3d(-1.0,-1.0,1.0));
  vertices04.push_back(Vec3d(-1.0,-1.0,-1.0));
  
  newrod = setupRod(opts, vertices04, vertices04);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,0.0));
  newrod->setVelocity(1, Vec3d(0.0,0.0,0.0));
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);

  // Generate sixth rod
  std::vector<Vec3d> vertices05;
  vertices05.push_back(Vec3d(0.0,-1.0,1.6));
  vertices05.push_back(Vec3d(0.0,-1.0,-0.4));
  
  newrod = setupRod(opts, vertices05, vertices05);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(-0.2,0.0,0.0));
  newrod->setVelocity(1, Vec3d(-0.2,0.0,0.0));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  m_br_stepper->prepareForExecution();
}

// Fixed rod collision tests
void CollisionTest::genRodRodPenaltyFixed00()
{
  m_br_stepper->enablePenaltyImpulses();
  m_br_stepper->disableIterativeInelasticImpulses();
  m_br_stepper->setEdgeEdgePenalty(200.0);
  m_br_stepper->setVertexFacePenalty(200.0);
  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 2;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);

  // Generate first rod
  std::vector<Vec3d> vertices00;
  vertices00.push_back(Vec3d(-1.0,0.0,0.0));
  vertices00.push_back(Vec3d( 1.0,0.0,0.0));
  
  ElasticRod* newrod = setupRod(opts, vertices00, vertices00);
  newrod->getBoundaryCondition()->setDesiredVertexPosition(0, newrod->getVertex(0));
  newrod->getBoundaryCondition()->setDesiredVertexPosition(1, newrod->getVertex(1));
  m_rods.push_back(newrod);

  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);

  // Generate second rod
  std::vector<Vec3d> vertices01;
  vertices01.push_back(Vec3d(0.0,0.0,-2.0));
  vertices01.push_back(Vec3d(0.0,2.0,-2.0));
  
  newrod = setupRod(opts, vertices01, vertices01);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,0.2));
  newrod->setVelocity(1, Vec3d(0.0,0.0,0.2));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  // Generate third rod
  std::vector<Vec3d> vertices02;
  vertices02.push_back(Vec3d(0.0,0.0,2.0));
  vertices02.push_back(Vec3d(0.0,-2.0,2.0));
  
  newrod = setupRod(opts, vertices02, vertices02);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,-0.2));
  newrod->setVelocity(1, Vec3d(0.0,0.0,-0.2));
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Generate third rod
  std::vector<Vec3d> vertices03;
  vertices03.push_back(Vec3d(-3.0,0.0,2.0));
  vertices03.push_back(Vec3d(-3.0,-2.0,2.0));
  
  newrod = setupRod(opts, vertices03, vertices03);
  m_rods.push_back(newrod);
  newrod->getBoundaryCondition()->setDesiredVertexPosition(0, newrod->getVertex(0));
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);

  // Generate fourth rod
  std::vector<Vec3d> vertices04;
  vertices04.push_back(Vec3d(-4.0,-1.0,-2.0));
  vertices04.push_back(Vec3d(-2.0,-1.0,-2.0));
  
  newrod = setupRod(opts, vertices04, vertices04);
  newrod->setVelocity(0, Vec3d(0.0,0.0,0.2));
  newrod->setVelocity(1, Vec3d(0.0,0.0,0.2));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  m_br_stepper->prepareForExecution();
}

// Swinging rod collides with a fixed rod
void CollisionTest::genRodRodPenaltyFixed01()
{
  m_br_stepper->setGravity(Vec3d(0,-19.81,0));
  
  m_br_stepper->enablePenaltyImpulses();
  m_br_stepper->disableIterativeInelasticImpulses();
  m_br_stepper->setEdgeEdgePenalty(100.0);
  m_br_stepper->setVertexFacePenalty(200.0);

  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 2;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);

  // Generate rod with many vertices
  opts.numVertices = 20;
  std::vector<Vec3d> vertices00;
  
  double length = 4.0;
  for( int i = 0; i < opts.numVertices; ++i )
  {
    vertices00.push_back(Vec3d(0.0,0.0,i*(length/(opts.numVertices-1))));
  }

  ElasticRod* newrod = setupRod(opts, vertices00, vertices00);
  newrod->getBoundaryCondition()->setDesiredVertexPosition(0, newrod->getVertex(0));

  std::vector<RodForce*> forces = newrod->getForces();
  
  for( int i = 0; i < (int) forces.size(); ++i )
  {
    if( forces[i]->getName() == "RodStretchingForce" )
    {
      for( ElasticRod::edge_iter itr = newrod->edges_begin(); itr != newrod->edges_end(); ++itr )
      {
        ((RodStretchingForce*)forces[i])->setKs(*itr,250.0);
      }
      //std::cout << ((RodStretchingForce*) forces[i])->globalEnergy() << std::endl;
    }
  }
  
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Generate second rod
  opts.numVertices = 2;
  std::vector<Vec3d> vertices01;
  vertices01.push_back(Vec3d(-1.0,-1.0,1.0));
  vertices01.push_back(Vec3d(1.0,-1.0,1.0));
  
  newrod = setupRod(opts, vertices01, vertices01);
  newrod->getBoundaryCondition()->setDesiredVertexPosition(0, newrod->getVertex(0));
  newrod->getBoundaryCondition()->setDesiredVertexPosition(1, newrod->getVertex(1));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  m_br_stepper->prepareForExecution();
}

// Rod of different radii colliding
void CollisionTest::genRodRodPenaltyDifferentRadii00()
{
  m_br_stepper->enablePenaltyImpulses();
  m_br_stepper->disableIterativeInelasticImpulses();
  m_br_stepper->setEdgeEdgePenalty(200.0);
  m_br_stepper->setVertexFacePenalty(200.0);
  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 2;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);
  
  // Generate first rod
  opts.radiusA = opts.radiusB = 0.3;
  std::vector<Vec3d> vertices00;
  vertices00.push_back(Vec3d(-1.0,0.0,0.0));
  vertices00.push_back(Vec3d( 1.0,0.0,0.0));
  
  ElasticRod* newrod = setupRod(opts, vertices00, vertices00);
  newrod->getBoundaryCondition()->setDesiredVertexPosition(0, newrod->getVertex(0));
  newrod->getBoundaryCondition()->setDesiredVertexPosition(1, newrod->getVertex(1));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Generate second rod
  opts.radiusA = opts.radiusB = 0.1;
  std::vector<Vec3d> vertices01;
  vertices01.push_back(Vec3d(0.0,-1.0,-2.0));
  vertices01.push_back(Vec3d(0.0,3.0,-2.0));
  
  newrod = setupRod(opts, vertices01, vertices01);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,0.2));
  newrod->setVelocity(1, Vec3d(0.0,0.0,0.2));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  m_br_stepper->prepareForExecution();
}

// Stress the code with lots of multiple contact by dropping many rods
void CollisionTest::genRodRodMultipleContact00()
{
  m_br_stepper->setGravity(Vec3d(0,-9.81,0));
  
  m_br_stepper->enablePenaltyImpulses();
  m_br_stepper->disableIterativeInelasticImpulses();
  m_br_stepper->setEdgeEdgePenalty(200.0);
  m_br_stepper->setVertexFacePenalty(200.0);

  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 2;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);


  int num_rods_cylinder = 150;
  double dtheta = 2.0*pi/((double)num_rods_cylinder);
  for( int i = 0; i < num_rods_cylinder; ++i )
  {
    Vec3d n(cos(i*dtheta),0.0,sin(i*dtheta));
    n.normalize();
    
    double start = 2.0*opts.radiusA/sqrt(2.0*(1.0-cos(dtheta)));
    
    std::vector<Vec3d> vertices;
    vertices.push_back(0.0*n);
    vertices.push_back((1.0+start)*n+Vec3d(0.0,2.0,0.0));
    ElasticRod* newrod = setupRod(opts,vertices,vertices);
    newrod = setupRod(opts,vertices,vertices);
    newrod->getBoundaryCondition()->setDesiredVertexPosition(0, newrod->getVertex(0));
    newrod->getBoundaryCondition()->setDesiredVertexPosition(1, newrod->getVertex(1));
    m_br_stepper->addRod(newrod);  
    m_world->addObject(newrod);
  }
  
  
  int num_falling_rods = 64;
  dtheta = 2.0*pi/16.0;
  for( int i = 0; i < num_falling_rods; ++i )
  {
    Vec3d n(cos(i*dtheta),0.0,sin(i*dtheta));
    std::vector<Vec3d> vertices;
    vertices.push_back(0.0*n+Vec3d(0.0,1.5+i*0.2,0.0));
    vertices.push_back(2.0*n+Vec3d(0.0,1.5+i*0.2,0.0));
    ElasticRod* newrod = setupRod(opts,vertices,vertices);
    newrod = setupRod(opts,vertices,vertices);
    m_br_stepper->addRod(newrod);  
    m_world->addObject(newrod);
  }
  
  m_br_stepper->prepareForExecution();
}

/////////////////////////////////////////////////////////////////////////////
// Rod-object penalty tests

void CollisionTest::genRodObjectPenaltyBasic00()
{
  m_br_stepper->setGravity(Vec3d(0,0,0));
  
  m_br_stepper->enablePenaltyImpulses();
  m_br_stepper->disableIterativeInelasticImpulses();
  m_br_stepper->setEdgeEdgePenalty(200.0);
  m_br_stepper->setVertexFacePenalty(200.0);
  
  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 2;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);


  // Create a triangle
  TriangleMesh* testmesh = new TriangleMesh();
  TriangleMesh::vertex_handle vh0 = testmesh->addVertex();
  TriangleMesh::vertex_handle vh1 = testmesh->addVertex();
  TriangleMesh::vertex_handle vh2 = testmesh->addVertex();

  TriangleMesh::edge_handle eh0 = testmesh->addEdge(vh0,vh1);
  TriangleMesh::edge_handle eh1 = testmesh->addEdge(vh1,vh2);
  TriangleMesh::edge_handle eh2 = testmesh->addEdge(vh2,vh0);
  
  TriangleMesh::face_handle fh0 = testmesh->addFace(vh0,vh1,vh2);
  
  testmesh->getVertex(vh0) = Vec3d(-1,0,-1);
  testmesh->getVertex(vh1) = Vec3d(1,0,-1);
  testmesh->getVertex(vh2) = Vec3d(0,0,1);
  
  m_br_stepper->addTriangleMesh(testmesh);
  m_world->addObject(testmesh);

  
  // Create a rod to hit the face
  std::vector<Vec3d> vertices;
  vertices.push_back(Vec3d(0.0,3.0,-0.5));
  vertices.push_back(Vec3d(0.0,1.0,-0.5));
  ElasticRod* newrod = setupRod(opts,vertices,vertices);
  newrod = setupRod(opts,vertices,vertices);
  newrod->setVelocity(0,Vec3d(0.0,-0.5,0.0));
  newrod->setVelocity(1,Vec3d(0.0,-0.5,0.0));
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Create a second rod to hit the face
  vertices.clear();
  vertices.push_back(Vec3d(0.5,3.0,-0.5));
  vertices.push_back(Vec3d(0.5,1.0,-0.5));
  newrod = setupRod(opts,vertices,vertices);
  newrod = setupRod(opts,vertices,vertices);
  newrod->setVelocity(0,Vec3d(0.0,-0.5,0.0));
  newrod->setVelocity(1,Vec3d(0.0,-0.5,0.0));
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Create a third rod to hit the face
  vertices.clear();
  vertices.push_back(Vec3d(-0.5,3.0,-0.5));
  vertices.push_back(Vec3d(-0.5,1.0,-0.5));
  newrod = setupRod(opts,vertices,vertices);
  newrod = setupRod(opts,vertices,vertices);
  newrod->setVelocity(0,Vec3d(0.0,-0.5,0.0));
  newrod->setVelocity(1,Vec3d(0.0,-0.5,0.0));
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Create a fourth rod to hit the face
  vertices.clear();
  vertices.push_back(Vec3d(0.0,3.0,0.5));
  vertices.push_back(Vec3d(0.0,1.0,0.5));
  newrod = setupRod(opts,vertices,vertices);
  newrod = setupRod(opts,vertices,vertices);
  newrod->setVelocity(0,Vec3d(0.0,-0.5,0.0));
  newrod->setVelocity(1,Vec3d(0.0,-0.5,0.0));
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);  
  
  // Create three rods to hit the edges of the face
  vertices.clear();
  vertices.push_back(Vec3d(0.0,-1.0,-2.0));
  vertices.push_back(Vec3d(0.0,1.0,-2.0));
  newrod = setupRod(opts,vertices,vertices);
  newrod->setVelocity(0,Vec3d(0.0,0.0,0.5));
  newrod->setVelocity(1,Vec3d(0.0,0.0,0.5));
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  vertices.clear();
  vertices.push_back(Vec3d(2.0,-1.0,-0.25));
  vertices.push_back(Vec3d(2.0,1.0,-0.25));
  newrod = setupRod(opts,vertices,vertices);
  newrod->setVelocity(0,Vec3d(-0.4,0.0,0.0));
  newrod->setVelocity(1,Vec3d(-0.4,0.0,0.0));
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  vertices.clear();
  vertices.push_back(Vec3d(-2.0,-1.0,-0.25));
  vertices.push_back(Vec3d(-2.0,1.0,-0.25));
  newrod = setupRod(opts,vertices,vertices);
  newrod->setVelocity(0,Vec3d(0.4,0.0,0.0));
  newrod->setVelocity(1,Vec3d(0.4,0.0,0.0));
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  m_br_stepper->prepareForExecution();
}

void CollisionTest::genRodObjectPenaltyBasic01()
{
  m_br_stepper->setGravity(Vec3d(0,0,0));
  
  m_br_stepper->enablePenaltyImpulses();
  m_br_stepper->disableIterativeInelasticImpulses();
  m_br_stepper->setEdgeEdgePenalty(200.0);
  m_br_stepper->setVertexFacePenalty(200.0);
  
  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 2;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);
  
  
  // Create a triangle
  TriangleMesh* testmesh = new TriangleMesh();
  TriangleMesh::vertex_handle vh0 = testmesh->addVertex();
  TriangleMesh::vertex_handle vh1 = testmesh->addVertex();
  TriangleMesh::vertex_handle vh2 = testmesh->addVertex();
  
  TriangleMesh::edge_handle eh0 = testmesh->addEdge(vh0,vh1);
  TriangleMesh::edge_handle eh1 = testmesh->addEdge(vh1,vh2);
  TriangleMesh::edge_handle eh2 = testmesh->addEdge(vh2,vh0);
  
  TriangleMesh::face_handle fh0 = testmesh->addFace(vh0,vh1,vh2);
  
  testmesh->getVertex(vh0) = Vec3d(-1,0,-1);
  testmesh->getVertex(vh1) = Vec3d(1,0,-1);
  testmesh->getVertex(vh2) = Vec3d(0,0,1);
  
  m_br_stepper->addTriangleMesh(testmesh);
  m_world->addObject(testmesh);
  
  // Create a rod to hit the triangles
  std::vector<Vec3d> vertices;
  vertices.push_back(Vec3d(0.0,1.5,-0.5));
  vertices.push_back(Vec3d(0.0,0.5,-0.5));
  ElasticRod* newrod = setupRod(opts,vertices,vertices);
  newrod = setupRod(opts,vertices,vertices);
  newrod->setVelocity(0,Vec3d(0.0,-0.5,0.0));
  newrod->setVelocity(1,Vec3d(0.0,-0.5,0.0));
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Create a second triangle
  testmesh = new TriangleMesh();
  vh0 = testmesh->addVertex();
  vh1 = testmesh->addVertex();
  vh2 = testmesh->addVertex();
  
  eh0 = testmesh->addEdge(vh0,vh1);
  eh1 = testmesh->addEdge(vh1,vh2);
  eh2 = testmesh->addEdge(vh2,vh0);
  
  fh0 = testmesh->addFace(vh0,vh1,vh2);
  
  testmesh->getVertex(vh0) = Vec3d(-1,2,-1);
  testmesh->getVertex(vh1) = Vec3d(1,2,-1);
  testmesh->getVertex(vh2) = Vec3d(0,2,1);
  
  m_br_stepper->addTriangleMesh(testmesh);
  m_world->addObject(testmesh);
  
  // Create a third triangle
  testmesh = new TriangleMesh();
  vh0 = testmesh->addVertex();
  vh1 = testmesh->addVertex();
  vh2 = testmesh->addVertex();
  
  eh0 = testmesh->addEdge(vh0,vh1);
  eh1 = testmesh->addEdge(vh1,vh2);
  eh2 = testmesh->addEdge(vh2,vh0);
  
  fh0 = testmesh->addFace(vh0,vh1,vh2);
  
  testmesh->getVertex(vh0) = Vec3d(2,-1,-1);
  testmesh->getVertex(vh1) = Vec3d(2,-1,1);
  testmesh->getVertex(vh2) = Vec3d(2,1,0);
  
  m_br_stepper->addTriangleMesh(testmesh);
  m_world->addObject(testmesh);
  
  m_br_stepper->prepareForExecution();
}  

void CollisionTest::genRodObjectPenaltyBasic02()
{
  m_br_stepper->setGravity(Vec3d(0,-0.5,0));
  
  m_br_stepper->enablePenaltyImpulses();
  m_br_stepper->disableIterativeInelasticImpulses();
  m_br_stepper->setEdgeEdgePenalty(200.0);
  m_br_stepper->setVertexFacePenalty(200.0);
  
  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 20;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);

  ObjParser objparser;
  TriangleMesh* tri_mesh = new TriangleMesh();
  objparser.loadTriangularMesh( "TriangulatedSphere.obj", *tri_mesh );
  m_world->addObject(tri_mesh);
  m_br_stepper->addTriangleMesh(tri_mesh);
  
  int numrods = 14;
  double dx = 6/((double)numrods-1.0);
  for( int j = 0; j < numrods; ++j )
  {
    // Create a rod to hit the sphere
    std::vector<Vec3d> vertices;
    for( int i = 0; i < opts.numVertices; ++i ) vertices.push_back(Vec3d(-3+j*dx,3.0,-2.0+0.4*i));

    ElasticRod* newrod = setupRod(opts,vertices,vertices);
    newrod = setupRod(opts,vertices,vertices);
    m_br_stepper->addRod(newrod);
    m_world->addObject(newrod);
  }
  
  
  m_br_stepper->prepareForExecution();
}

void CollisionTest::genRodObjectPenaltyBasic03()
{
  m_br_stepper->setGravity(Vec3d(0,-0.05,0));
  
  m_br_stepper->enablePenaltyImpulses();
  m_br_stepper->disableIterativeInelasticImpulses();
  m_br_stepper->setEdgeEdgePenalty(200.0);
  m_br_stepper->setVertexFacePenalty(200.0);
  
  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 2;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);

  ObjParser objparser;
  TriangleMesh* tri_mesh = new TriangleMesh();
  objparser.loadTriangularMesh( "TriangulatedBoxContainer.obj", *tri_mesh );
  m_world->addObject(tri_mesh);
  m_br_stepper->addTriangleMesh(tri_mesh);
  
  int num_falling_rods = 64;
  double dtheta = 2.0*pi/16.0;
  for( int i = 0; i < num_falling_rods; ++i )
  {
    Vec3d n(cos(i*dtheta),0.0,sin(i*dtheta));
    std::vector<Vec3d> vertices;
    vertices.push_back(0.0*n+Vec3d(0.0,1.5+i*0.2,0.0));
    vertices.push_back(2.0*n+Vec3d(0.0,1.5+i*0.2,0.0));
    ElasticRod* newrod = setupRod(opts,vertices,vertices);
    newrod = setupRod(opts,vertices,vertices);
    m_br_stepper->addRod(newrod);  
    m_world->addObject(newrod);
  }
  
  m_br_stepper->prepareForExecution();
}


/////////////////////////////////////////////////////////////////////////////
// Edge-edge iterative inelastic impulse tests

void CollisionTest::genEdgeEdgeImpulseBasic00()
{
  m_br_stepper->setGravity(Vec3d(0,0,0));
  
  m_br_stepper->disablePenaltyImpulses();
  m_br_stepper->enableIterativeInelasticImpulses();
  
  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 2;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);

  // Generate first rod
  std::vector<Vec3d> vertices00;
  vertices00.push_back(Vec3d(-1.0,0.0,0.0));
  vertices00.push_back(Vec3d( 1.0,0.0,0.0));
  
  ElasticRod* newrod = setupRod(opts, vertices00, vertices00);
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  // Generate second rod
  std::vector<Vec3d> vertices01;
  vertices01.push_back(Vec3d(0.0,-1.0,-2.0));
  vertices01.push_back(Vec3d(0.0,1.0,-2.0));
  
  newrod = setupRod(opts, vertices01, vertices01);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,0.2));
  newrod->setVelocity(1, Vec3d(0.0,0.0,0.2));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  
  // Generate third rod
  vertices01.clear();
  vertices01.push_back(Vec3d(3.0,0.0,-1.0));
  vertices01.push_back(Vec3d(3.0,0.0,1.0));
  
  newrod = setupRod(opts, vertices01, vertices01);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,0.0));
  newrod->setVelocity(1, Vec3d(0.0,0.0,0.0));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  // Generate fourth rod
  vertices01.clear();
  vertices01.push_back(Vec3d(2.0,2.0,2.5));
  vertices01.push_back(Vec3d(4.0,2.0,2.5));
  
  newrod = setupRod(opts, vertices01, vertices01);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,-2.0,-2.0));
  newrod->setVelocity(1, Vec3d(0.0,-2.0,-2.0));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);



  // Generate fifth rod
  vertices01.clear();
  vertices01.push_back(Vec3d(-3.0,0.0,-1.0));
  vertices01.push_back(Vec3d(-3.0,0.0,1.0));
  
  newrod = setupRod(opts, vertices01, vertices01);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,0.0,0.0));
  newrod->setVelocity(1, Vec3d(0.0,0.0,0.0));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  // Generate sixth rod
  vertices01.clear();
  vertices01.push_back(Vec3d(-2.0,2.0,2.5));
  vertices01.push_back(Vec3d(-4.0,2.0,2.5));
  
  newrod = setupRod(opts, vertices01, vertices01);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,-2.0,-2.0));
  newrod->setVelocity(1, Vec3d(0.0,-2.0,-2.0));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);



  // Generate seventh rod
  vertices01.clear();
  vertices01.push_back(Vec3d(-2.0,4.0,-1.0));
  vertices01.push_back(Vec3d(-2.0,4.0,1.0));
  
  newrod = setupRod(opts, vertices01, vertices01);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,10.0,0.0));
  newrod->setVelocity(1, Vec3d(0.0,10.0,0.0));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  // Generate eigth rod
  vertices01.clear();
  vertices01.push_back(Vec3d(-3.0,6.0,0.0));
  vertices01.push_back(Vec3d(-1.0,6.0,0.0));
  
  newrod = setupRod(opts, vertices01, vertices01);
  m_rods.push_back(newrod);
  newrod->setVelocity(0, Vec3d(0.0,-10.0,0.0));
  newrod->setVelocity(1, Vec3d(0.0,-10.0,0.0));
  
  m_br_stepper->addRod(newrod);  
  m_world->addObject(newrod);
  
  
  m_br_stepper->testCoplanarityTime();
  m_br_stepper->prepareForExecution();
}

void CollisionTest::genEdgeEdgeImpulseDifferentMass00()
{
  m_br_stepper->setGravity(Vec3d(0,0,0));
  
  m_br_stepper->disablePenaltyImpulses();
  m_br_stepper->enableIterativeInelasticImpulses();
  
  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 2;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);

  // Generate first fixed rod
  std::vector<Vec3d> vertices00;
  vertices00.push_back(Vec3d(-1.0,0.0,0.0));
  vertices00.push_back(Vec3d( 1.0,0.0,0.0));
  
  ElasticRod* newrod = setupRod(opts, vertices00, vertices00);
  newrod->getBoundaryCondition()->setDesiredVertexPosition(0, newrod->getVertex(0));
  newrod->getBoundaryCondition()->setDesiredVertexPosition(1, newrod->getVertex(1));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Generate free rod to collide with fixed
  vertices00.clear();
  vertices00.push_back(Vec3d(0.0,1.0,-1.0));
  vertices00.push_back(Vec3d(0.0,1.0,1.0));
  
  newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(0.0,-10.0,0.0));
  newrod->setVelocity(1, Vec3d(0.0,-10.0,0.0));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
 
  
  
  // Generate a massive rod to hit a light rod
  vertices00.clear();
  vertices00.push_back(Vec3d(5.0,1.0,-1.0));
  vertices00.push_back(Vec3d(5.0,1.0,1.0));
  
  newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(0.0,-2.0,0.0));
  newrod->setVelocity(1, Vec3d(0.0,-2.0,0.0));
  m_rods.push_back(newrod);
  
  newrod->setVertexMass(0,10.0); newrod->setVertexMass(1,10.0);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Generate a light rod to be hit by a massive rod
  vertices00.clear();
  vertices00.push_back(Vec3d(4.0,0.0,0.0));
  vertices00.push_back(Vec3d(6.0,0.0,0.0));
  
  newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(0.0,0.0,0.0));
  newrod->setVelocity(1, Vec3d(0.0,0.0,0.0));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);

  
  
  // Generate a light rod to hit a massive rod
  vertices00.clear();
  vertices00.push_back(Vec3d(-5.0,1.0,-1.0));
  vertices00.push_back(Vec3d(-5.0,1.0,1.0));
  
  newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(0.0,-2.0,0.0));
  newrod->setVelocity(1, Vec3d(0.0,-2.0,0.0));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Generate a massive rod to be hit by a light rod
  vertices00.clear();
  vertices00.push_back(Vec3d(-6.0,0.0,0.0));
  vertices00.push_back(Vec3d(-4.0,0.0,0.0));
  
  newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(0.0,0.0,0.0));
  newrod->setVelocity(1, Vec3d(0.0,0.0,0.0));
  m_rods.push_back(newrod);
  
  newrod->setVertexMass(0,10.0); newrod->setVertexMass(1,10.0);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  m_br_stepper->prepareForExecution();
}


/////////////////////////////////////////////////////////////////////////////
// Vertex-face iterative inelastic impulse tests

void CollisionTest::genVertexFaceImpulseBasic00()
{
  m_br_stepper->setGravity(Vec3d(0,0,0));
  m_br_stepper->disablePenaltyImpulses();
  m_br_stepper->enableIterativeInelasticImpulses();
  
  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 2;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);


  // Create a fixed triangle
  TriangleMesh* testmesh = new TriangleMesh();
  TriangleMesh::vertex_handle vh0 = testmesh->addVertex();
  TriangleMesh::vertex_handle vh1 = testmesh->addVertex();
  TriangleMesh::vertex_handle vh2 = testmesh->addVertex();
  
  TriangleMesh::edge_handle eh0 = testmesh->addEdge(vh0,vh1);
  TriangleMesh::edge_handle eh1 = testmesh->addEdge(vh1,vh2);
  TriangleMesh::edge_handle eh2 = testmesh->addEdge(vh2,vh0);
  
  TriangleMesh::face_handle fh0 = testmesh->addFace(vh0,vh1,vh2);
  
  testmesh->getVertex(vh0) = Vec3d(-1,0,-1);
  testmesh->getVertex(vh1) = Vec3d(1,0,-1);
  testmesh->getVertex(vh2) = Vec3d(0,0,1);
  
  m_br_stepper->addTriangleMesh(testmesh);
  m_world->addObject(testmesh);


  // Generate a rod to collide with the triangle
  std::vector<Vec3d> vertices00;
  vertices00.push_back(Vec3d(0.0,2.05,0.0));
  vertices00.push_back(Vec3d(0.0,3.05,0.0));
  
  ElasticRod* newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(0.0,-15.0,0.0));
  newrod->setVelocity(1, Vec3d(0.0,-15.0,0.0));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);

  // Generate another rod to collide with the triangle
  vertices00.clear();
  vertices00.push_back(Vec3d(0.0,2.0,-0.3));
  vertices00.push_back(Vec3d(0.0,3.0,-0.3));
  
  newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(0.0,-1.0,0.0));
  newrod->setVelocity(1, Vec3d(0.0,-1.0,0.0));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Generate yet another rod to collide with the triangle
  vertices00.clear();
  vertices00.push_back(Vec3d(0.0,2.0,-0.6));
  vertices00.push_back(Vec3d(0.0,3.0,-0.6));
  
  newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(0.0,-0.5,0.0));
  newrod->setVelocity(1, Vec3d(0.0,-0.5,0.0));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Generate yet yet another rod to collide with the triangle
  vertices00.clear();
  vertices00.push_back(Vec3d(-0.3,2.0,-0.6));
  vertices00.push_back(Vec3d(-0.3,3.0,-0.6));
  
  newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(0.0,-0.25,0.0));
  newrod->setVelocity(1, Vec3d(0.0,-0.25,0.0));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Generate yet yet yet another rod to collide with the triangle
  vertices00.clear();
  vertices00.push_back(Vec3d(0.3,2.0,-0.6));
  vertices00.push_back(Vec3d(0.3,3.0,-0.6));
  
  newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(0.0,-0.25,0.0));
  newrod->setVelocity(1, Vec3d(0.0,-0.25,0.0));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);

  // Generate yet yet yet yet another rod to collide with the triangle
  vertices00.clear();
  vertices00.push_back(Vec3d(0.3,-2.0,-0.9));
  vertices00.push_back(Vec3d(0.3,-3.0,-0.9));
  
  newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(0.0,0.6,0.0));
  newrod->setVelocity(1, Vec3d(0.0,0.6,0.0));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Generate yet yet yet yet yet another rod to collide with the triangle
  vertices00.clear();
  vertices00.push_back(Vec3d(-0.3,-2.0,-0.9));
  vertices00.push_back(Vec3d(-0.3,-3.0,-0.9));
  
  newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(0.0,0.6,0.0));
  newrod->setVelocity(1, Vec3d(0.0,0.6,0.0));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Generate a rod to collide with the tirangle edge-edge 
  vertices00.clear();
  vertices00.push_back(Vec3d(-2.0,-1.0,0.5));
  vertices00.push_back(Vec3d(-2.0, 1.0,0.5));
  
  newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(1.5,0.0,0.0));
  newrod->setVelocity(1, Vec3d(1.5,0.0,0.0));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  // Generate another rod to collide with the tirangle edge-edge 
  vertices00.clear();
  vertices00.push_back(Vec3d(2.0,-1.0,0.5));
  vertices00.push_back(Vec3d(2.0, 1.0,0.5));
  
  newrod = setupRod(opts, vertices00, vertices00);
  newrod->setVelocity(0, Vec3d(-1.5,0.0,0.0));
  newrod->setVelocity(1, Vec3d(-1.5,0.0,0.0));
  m_rods.push_back(newrod);
  
  m_br_stepper->addRod(newrod);
  m_world->addObject(newrod);
  
  
  m_br_stepper->prepareForExecution();
}

void CollisionTest::genVertexFaceImpulseBasic01()
{
  m_br_stepper->setGravity(Vec3d(0,-0.5,0));
  m_br_stepper->disablePenaltyImpulses();
  m_br_stepper->enableIterativeInelasticImpulses();
  
  RodOptions opts;
  opts.radiusA = opts.radiusB = 0.1;
  opts.numVertices = 20;
  opts.ShearModulus = 100.0;
  opts.YoungsModulus = 500.0;
  opts.viscosity = 0.0;
  GetScalarOpt("mass-damping") = 0.0;
  //GetScalarOpt("dt") = 0.01;
  setDt(0.01);
  m_br_stepper->setDt(0.01);
  
  
  ObjParser objparser;
  TriangleMesh* tri_mesh = new TriangleMesh();
  objparser.loadTriangularMesh( "TriangulatedSphere.obj", *tri_mesh );
  m_world->addObject(tri_mesh);
  m_br_stepper->addTriangleMesh(tri_mesh);
  
  int numrods = 14;
  double dx = 6/((double)numrods-1.0);
  for( int j = 0; j < numrods; ++j )
  {
    // Create a rod to hit the sphere
    std::vector<Vec3d> vertices;
    for( int i = 0; i < opts.numVertices; ++i ) vertices.push_back(Vec3d(-3+j*dx,3.0,-2.0+0.4*i));
    
    ElasticRod* newrod = setupRod(opts,vertices,vertices);
    newrod = setupRod(opts,vertices,vertices);
    m_br_stepper->addRod(newrod);
    m_world->addObject(newrod);
  }
  
  m_br_stepper->prepareForExecution();
}

void CollisionTest::Setup()
{
  loadDynamicsProps();
  
  RodTimeStepper::Method integrator;
  if( GetStringOpt("integrator") == "implicit" ) integrator = RodTimeStepper::IMPL_EULER;
  else assert( false );
  
  m_br_stepper = new BridsonStepper( RodTimeStepper::IMPL_EULER, GetIntOpt("iterations"), GetScalarOpt("dt"), GetScalarOpt("mass-damping"), GetVecOpt("gravity") );

  //genRodRodPenaltyBasic00();
  //genRodRodPenaltyBasic01();
  //genRodRodPenaltyCorner00();
  //genRodRodPenaltyFixed00();
  //genRodRodPenaltyFixed01();
  //genRodRodPenaltyDifferentRadii00();
  //genRodRodMultipleContact00();
  
  //genRodObjectPenaltyBasic00();
  //genRodObjectPenaltyBasic01();
  //genRodObjectPenaltyBasic02();
  //genRodObjectPenaltyBasic03();
  
  //genEdgeEdgeImpulseBasic00();
  //genEdgeEdgeImpulseDifferentMass00();

  //genVertexFaceImpulseBasic00();
  //genVertexFaceImpulseBasic01();

  //genTorusTest();
  //genNullTest();
  genTorusPull();

  m_world->addController(m_br_stepper);
}

void CollisionTest::AtEachTimestep()
{
  genTorusPullAtEachTimestep();
}





///////////////////////////////////////////////////////////////////////////////
// Troublesome examples that I'm not sure work correctly! 

//m_br_stepper->setGravity(Vec3d(0,0,0));
//m_br_stepper->disablePenaltyImpulses();
//m_br_stepper->enableIterativeInelasticImpulses();
//
//RodOptions opts;
//opts.ShearModulus = GetScalarOpt("shear-modulus");
//opts.YoungsModulus = GetScalarOpt("youngs-modulus");
//opts.numVertices = 2;
//opts.radiusA = opts.radiusB = 0.1;
//
//
//// Create a fixed triangle
//ScriptedTriangleMesh* testmesh = new ScriptedTriangleMesh();
//ScriptedTriangleMesh::vertex_handle vh0 = testmesh->addVertex();
//ScriptedTriangleMesh::vertex_handle vh1 = testmesh->addVertex();
//ScriptedTriangleMesh::vertex_handle vh2 = testmesh->addVertex();
//
//ScriptedTriangleMesh::edge_handle eh0 = testmesh->addEdge(vh0,vh1);
//ScriptedTriangleMesh::edge_handle eh1 = testmesh->addEdge(vh1,vh2);
//ScriptedTriangleMesh::edge_handle eh2 = testmesh->addEdge(vh2,vh0);
//
//ScriptedTriangleMesh::face_handle fh0 = testmesh->addFace(vh0,vh1,vh2);
//
//testmesh->getVertex(vh0) = Vec3d(-1,0,-1);
//testmesh->getVertex(vh1) = Vec3d(1,0,-1);
//testmesh->getVertex(vh2) = Vec3d(0,0,1);
//
//m_br_stepper->addSriptedTriangleMesh(testmesh);
//m_world->addObject(testmesh);
//
//
//// Generate a rod to collide with the triangle
//std::vector<Vec3d> vertices00;
//vertices00.push_back(Vec3d(0.0,2.0,0.0));
//vertices00.push_back(Vec3d(0.0,3.0,0.0));
//
//ElasticRod* newrod = setupRod(opts, vertices00, vertices00);
//newrod->setVelocity(0, Vec3d(0.0,-1.0,0.0));
//newrod->setVelocity(1, Vec3d(0.0,-1.0,0.0));
//m_rods.push_back(newrod);
//
//m_br_stepper->addRod(newrod);
//m_world->addObject(newrod);

///////////////////////////////////////////////////////////////////////////////






