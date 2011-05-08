/**
 * \file CollisionTestWeta.cc
 *
 * \author smith@cs.columbia.edu
 * \author adapted to weta needs by showard@wetafx.co.nz
 * \date 04/27/2010  and  feb 17 20111
 */

#include "CollisionTestWeta.hh"
#include <stdio.h>

CollisionTestWeta::CollisionTestWeta()
: Problem("Collision Tests Weta", "A number of tests highlighting problems with naive collision response.")
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

CollisionTestWeta::~CollisionTestWeta()
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

void CollisionTestWeta::Setup()
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
    case 10:
    {
      std::cout << "Executing test: MovingSphereSetup()" << std::endl;
      MovingSphereSetup();
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

void CollisionTestWeta::AtEachTimestep()
{
}


///////////////////////////////////////////////////////////////////////////////
// ROD OBJECT SNAGGING

RodObjSnaggingController::RodObjSnaggingController( ElasticRod& rod, double time, double dt, FILE * pFile )
: ScriptingController(time,dt)
, m_rod(rod)
, m_pFile(pFile)
{}
  
bool RodObjSnaggingController::execute()
{
  // cm / second
  double pull_rate = 15.0;
  
  if( ScriptingController::getTime() <= 25.0 )
  {  
    Vec3d dx = Vec3d(pull_rate*ScriptingController::getDt(),0.0,0.0);

    m_rod.setVertex(0,m_rod.getVertex(0)+dx);
    m_rod.setVertex(1,m_rod.getVertex(1)+dx);
    m_rod.getBoundaryCondition()->setDesiredVertexPosition(0,m_rod.getVertex(0));
    m_rod.getBoundaryCondition()->setDesiredVertexPosition(1,m_rod.getVertex(1));
    m_rod.getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
    //prints sepcified vertices to file
    //fprintf (m_pFile, "%f\n",m_rod.getVertex(0));
    Vec3d test0=m_rod.getVertex(0);
    Vec3d test1=m_rod.getVertex(1);
    Vec3d test10=m_rod.getVertex(10);
    Vec3d test19=m_rod.getVertex(19);
    Vec3d test39=m_rod.getVertex(39);
    //std::cout << "rod vertex 0  position "<<test0(0)<< std::endl;
    //std::cout << "rod vertex 0  position "<< m_rod.getVertex(0)<< std::endl;
    //std::cout << "rod vertex 1  position "<< m_rod.getVertex(1)<< std::endl;
   // std::cout << "rod vertex 10  position "<< m_rod.getVertex(10)<< std::endl;
   // std::cout << "rod vertex 19  position "<< m_rod.getVertex(19)<< std::endl;
    //fprintf(m_pFile, "this is a test");
    float v0_0= roundto( test0(0), 10000.0 );
    float v0_1= roundto( test0(1), 10000.0 );
    float v0_2= roundto( test0(2), 10000.0 );
    float v1_0= roundto( test1(0), 10000.0 );
    float v1_1= roundto( test1(1), 10000.0 );
    float v1_2= roundto( test1(2), 10000.0 );
    float v10_0= roundto( test10(0), 10000.0 );
    float v10_1= roundto( test10(1), 10000.0 );
    float v10_2= roundto( test10(2), 10000.0 );
    float v19_0= roundto( test19(0), 10000.0 );
    float v19_1= roundto( test19(1), 10000.0 );
    float v19_2= roundto( test19(2), 10000.0 );
    float v39_0= roundto( test39(0), 10000.0 );
    float v39_1= roundto( test39(1), 10000.0 );
    float v39_2= roundto( test39(2), 10000.0 );
    fprintf(m_pFile, "%9.4f %9.4f %9.4f\n",  v0_0, v0_1, v0_2);
    fprintf(m_pFile, "%9.4f %9.4f %9.4f\n",  v1_0, v1_1, v1_2);
    fprintf(m_pFile, "%9.4f %9.4f %9.4f\n", v10_0, v10_1, v10_2);
    fprintf(m_pFile, "%9.4f %9.4f %9.4f\n", v19_0, v19_1, v19_2);
    fprintf(m_pFile, "%9.4f %9.4f %9.4f\n", v39_0, v39_1, v39_2);
    //sleep(1);
  }
  
  if( ScriptingController::getTime() >2.0 )
  {
      //fprintf( m_pFile, "closing file" );
      fclose (m_pFile);
      exit(0); 
  }
  return true;
}

float  RodObjSnaggingController::roundto(float val, float prec) { return floor( prec * val + 0.5 ) / prec; }


void CollisionTestWeta::RodObjectSnaggingSetup()
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
  FILE * pFile;

  pFile=fopen("./testing/testBASimulator/run7_t0_test.txt","w");
   if ( pFile == NULL )
  {
       std::cerr << "Couldn't open the file!" << std::endl;
   }

  m_scripting_controllers.push_back(new RodObjSnaggingController(*m_rods[0],getTime(),GetScalarOpt("dt"), pFile));

  m_br_stepper = new BridsonStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);


}

void CollisionTestWeta::RodObjectSnaggingAtEachTimestep()
{
}

void CollisionTestWeta::RodObjectSnaggingMovingObjectSetup()
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

  m_scripting_controllers.push_back(new ObjTranslator(*tri_mesh,getTime(),GetScalarOpt("dt")));
  
  m_br_stepper = new BridsonStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
  sleep(1);
}

ObjTranslator::ObjTranslator( TriangleMesh& mesh, double time, double dt )
: ScriptingController(time,dt)
, m_mesh(mesh)
{}

bool ObjTranslator::execute()
{
  double distance = 8000.0;
  double transittime = 10.0;  //was 10.0
  
  double speed = distance/transittime;
  double translatedist = speed*getDt();
  
  for( TriangleMesh::vertex_iter itr = m_mesh.vertices_begin(); itr != m_mesh.vertices_end(); ++itr )
  {
    m_mesh.getVertex(*itr) += Vec3d(-translatedist,0.0,0.0);
    //std::cout << "vertex   position "<< m_mesh.getVertex(*itr)<< std::endl;
  }
   sleep(2);
  if( getTime() > 10.5 ) exit(0);
  
  return true;
}

ObjTranslatorWeta::ObjTranslatorWeta( ElasticRod& rod,TriangleMesh& mesh, double time, double dt, FILE * pFile )
: ScriptingController(time,dt)
, m_mesh(mesh)
, m_rod(rod)
, m_pFile(pFile)
{}

bool ObjTranslatorWeta::execute()
{
  double distance = 8000.0;
  double transittime = 10.0;  //was 10.0
  
  double speed = distance/transittime;
  double translatedist = speed*getDt();
  // cm / second
  double pull_rate = 0.0;

  for( TriangleMesh::vertex_iter itr = m_mesh.vertices_begin(); itr != m_mesh.vertices_end(); ++itr )
  {
    m_mesh.getVertex(*itr) += Vec3d(-translatedist,0.0,0.0);
    //std::cout << "vertex   position "<< m_mesh.getVertex(*itr)<< std::endl;
  }

  if( ScriptingController::getTime() <= 25.0 )
  {  
    Vec3d dx = Vec3d(pull_rate*ScriptingController::getDt(),0.0,0.0);

    m_rod.setVertex(0,m_rod.getVertex(0)+dx);
    m_rod.setVertex(1,m_rod.getVertex(1)+dx);
    m_rod.getBoundaryCondition()->setDesiredVertexPosition(0,m_rod.getVertex(0));
    m_rod.getBoundaryCondition()->setDesiredVertexPosition(1,m_rod.getVertex(1));
    m_rod.getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
    //prints sepcified vertices to file
    //fprintf (m_pFile, "%f\n",m_rod.getVertex(0));
    Vec3d test0=m_rod.getVertex(0);
    Vec3d test1=m_rod.getVertex(1);
    Vec3d test10=m_rod.getVertex(10);
    Vec3d test19=m_rod.getVertex(19);
    Vec3d test39=m_rod.getVertex(39);
    //std::cout << "rod vertex 0  position "<<test0(0)<< std::endl;
    //std::cout << "rod vertex 0  position "<< m_rod.getVertex(0)<< std::endl;
    //std::cout << "rod vertex 1  position "<< m_rod.getVertex(1)<< std::endl;
   // std::cout << "rod vertex 10  position "<< m_rod.getVertex(10)<< std::endl;
   // std::cout << "rod vertex 19  position "<< m_rod.getVertex(19)<< std::endl;
    //fprintf(m_pFile, "this is a test");
    float v0_0= roundto( test0(0), 10000.0 );
    float v0_1= roundto( test0(1), 10000.0 );
    float v0_2= roundto( test0(2), 10000.0 );
    float v1_0= roundto( test1(0), 10000.0 );
    float v1_1= roundto( test1(1), 10000.0 );
    float v1_2= roundto( test1(2), 10000.0 );
    float v10_0= roundto( test10(0), 10000.0 );
    float v10_1= roundto( test10(1), 10000.0 );
    float v10_2= roundto( test10(2), 10000.0 );
    float v19_0= roundto( test19(0), 10000.0 );
    float v19_1= roundto( test19(1), 10000.0 );
    float v19_2= roundto( test19(2), 10000.0 );
    float v39_0= roundto( test39(0), 10000.0 );
    float v39_1= roundto( test39(1), 10000.0 );
    float v39_2= roundto( test39(2), 10000.0 );
    fprintf(m_pFile, "%9.4f %9.4f %9.4f\n",  v0_0, v0_1, v0_2);
    fprintf(m_pFile, "%9.4f %9.4f %9.4f\n",  v1_0, v1_1, v1_2);
    fprintf(m_pFile, "%9.4f %9.4f %9.4f\n", v10_0, v10_1, v10_2);
    fprintf(m_pFile, "%9.4f %9.4f %9.4f\n", v19_0, v19_1, v19_2);
    fprintf(m_pFile, "%9.4f %9.4f %9.4f\n", v39_0, v39_1, v39_2);
    //sleep(1);
  }
  
  if( getTime() >0.3 )
  {
      //fprintf( m_pFile, "closing file" );
      fclose (m_pFile);
      exit(0); 
  }
  
  return true;
}

float  ObjTranslatorWeta::roundto(float val, float prec) { return floor( prec * val + 0.5 ) / prec; }

void CollisionTestWeta::MovingSphereSetup()
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
  objparser.loadTriangularMesh( "assets/TriangulatedSphere.obj", *tri_mesh );
  for( TriangleMesh::vertex_iter itr = tri_mesh->vertices_begin(); itr != tri_mesh->vertices_end(); ++itr )
  {
    tri_mesh->getVertex(*itr) -= Vec3d(-20.0,5.0,0.0);  //was (0.0,5.0,0.0
  }
  m_world->addObject(tri_mesh);
  m_tri_objs.push_back(tri_mesh);
  FILE * pFile;

  pFile=fopen("./testing/testBASimulator/run7_t10_test.txt","w");
  if ( pFile == NULL )
  {
       std::cerr << "Couldn't open the file!" << std::endl;
   }

  m_scripting_controllers.push_back(new ObjTranslatorWeta(*m_rods[0],*tri_mesh,getTime(),GetScalarOpt("dt"), pFile));

  m_br_stepper = new BridsonStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);


}



void CollisionTestWeta::RodObjectSnaggingMovingObjectAtEachTimestep()
{
}

void CollisionTestWeta::RodObjectSnaggingVertFaceSetup()
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
  
  m_scripting_controllers.push_back(new ObjTranslator(*tri_mesh,getTime(),GetScalarOpt("dt")));
  
  m_br_stepper = new BridsonStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

void CollisionTestWeta::RodObjectSnaggingVertFaceTwoSetup()
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
  FILE * pFile;

  pFile=fopen("./testing/testBASimulator/run7_t3_test.txt","w");
   if ( pFile == NULL )
  {
       std::cerr << "Couldn't open the file!" << std::endl;
   }

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
  
  m_scripting_controllers.push_back(new RodObjSnaggingController(*m_rods[0],getTime(),GetScalarOpt("dt"),pFile));

  m_br_stepper = new BridsonStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}


/////////////////////////////////////////////////////////////////////////////
// ROD-ROD SNAGGING

RodRodSnagController::RodRodSnagController( ElasticRod& rod, double time, double dt )
: ScriptingController(time,dt)
, m_rod(rod)
{}

bool RodRodSnagController::execute()
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
  
  if( ScriptingController::getTime() > 10.5 ) exit(0);
  
  return true;
}

void CollisionTestWeta::RodRodSnaggingSetup()
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
  
  m_scripting_controllers.push_back(new RodRodSnagController(*m_rods[0],getTime(),GetScalarOpt("dt")));
  
  m_br_stepper = new BridsonStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

void CollisionTestWeta::RodRodSnaggingAtEachTimestep()
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

void CollisionTestWeta::RodRodSnaggingTwoSetup()
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

  m_scripting_controllers.push_back(new RodRodSnagController(*m_rods[0],getTime(),GetScalarOpt("dt")));

  m_br_stepper = new BridsonStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

void CollisionTestWeta::RodRodSnaggingTwoAtEachTimestep()
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

void CollisionTestWeta::RodRodSnaggingThreeSetup()
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
  
  m_scripting_controllers.push_back(new RodRodSnagController(*m_rods[0],getTime(),GetScalarOpt("dt")));

  m_br_stepper = new BridsonStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

void CollisionTestWeta::RodRodSnaggingThreeAtEachTimestep()
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

RodFixedRodSnagScriptingController::RodFixedRodSnagScriptingController( ElasticRod& rod, double time, double dt )
: ScriptingController(time,dt)
, m_rod(rod)
{}

bool RodFixedRodSnagScriptingController::execute()
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

void CollisionTestWeta::RodFixedRodSnaggingSetup()
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
  
  m_scripting_controllers.push_back(new RodFixedRodSnagScriptingController(*m_rods[0],getTime(),GetScalarOpt("dt")));
  
  m_br_stepper = new BridsonStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

void CollisionTestWeta::RodFixedRodSnaggingAtEachTimestep()
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














void CollisionTestWeta::RodRodSnaggingSmallSetup()
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
  
  m_scripting_controllers.push_back(new RodFixedRodSnagScriptingController(*m_rods[0],getTime(),GetScalarOpt("dt")));
  
  m_br_stepper = new BridsonStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}

void CollisionTestWeta::RodRodSnaggingSmallAtEachTimestep()
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

void CollisionTestWeta::RodRodSnaggingDifferentSizeSetup()
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
  
  m_scripting_controllers.push_back(new RodFixedRodSnagScriptingController(*m_rods[0],getTime(),GetScalarOpt("dt")));
  
  m_br_stepper = new BridsonStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
}








