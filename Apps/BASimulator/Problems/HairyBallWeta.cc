/**
 * \file HairyBallWeta.cc
 *
 * \author smith@cs.columbia.edu
 * \date 05/11/2010
 */

#include "HairyBallWeta.hh"

SphereRotatorW::SphereRotatorW( TriangleMesh& mesh, std::vector<ElasticRod*>& rods, double time, double dt )
: ScriptingController(time,dt)
, m_mesh(mesh)
, m_rods(rods)
{
  m_pFile=fopen("./testing/run8_t0_test.txt","w");
  if ( m_pFile == NULL )
  {
     std::cerr << "Couldn't open the file!" << std::endl;
  }
}

SphereRotatorW::~SphereRotatorW()
{
  fclose (m_pFile);
}


bool SphereRotatorW::execute()
{
  double thetax = 2.0*pi;
  double thetay = 2.0*pi;
  double thetaz = 2.0*pi;
  
  double txshake = 2.0;
  double thetaxrate = thetax/(txshake/2.0);
  
  double tyshake = 2.0;
  double thetayrate = thetay/(tyshake/2.0);
  
  double tzshake = 2.0;
  double thetazrate = thetaz/(tzshake/2.0);
  
  
  //std::cout << getTime() << "  -  " << getTime()*thetazrate << std::endl;
  //std::cout << thetazrate*GetScalarOpt("dt") << std::endl;
  
  /////////////////////////////////////////////////////////////////////////////
  // Allow the hair to fall for two seconds
  if( getTime() <= 1.0 )
  {
    // Do nothing
  }
  /////////////////////////////////////////////////////////////////////////////
  // "Shake" the sphere about the z axis
  else if( getTime() <= 2.0 )
  {
    //g_testthetaz += thetazrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetaz " << g_testthetaz << std::endl;
    
    Mat3d rotz;
    rotz << cos(thetazrate*getDt()),-sin(thetazrate*getDt()),0,
    sin(thetazrate*getDt()),cos(thetazrate*getDt()),0,
    0,0,1;
    assert( approxEq(rotz.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_mesh, rotz );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], rotz );
  }
  // Settle for a second
  else if( getTime() <= 2.5 )
  {
    // Do nothing
  }
  else if( getTime() <= 3.5)
  {
    //g_testthetaz -= thetazrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetaz " << g_testthetaz << std::endl;
    
    Mat3d rotz;
    rotz << cos(-thetazrate*getDt()),-sin(-thetazrate*getDt()),0,
    sin(-thetazrate*getDt()),cos(-thetazrate*getDt()),0,
    0,0,1;
    assert( approxEq(rotz.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_mesh, rotz );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], rotz );
  }
  // Settle for a second
  else if( getTime() <= 4.0 )
  {
    // Do nothing
  }
  /////////////////////////////////////////////////////////////////////////////
  // "Shake" the sphere about the x axis
  else if( getTime() <= 5.0 )
  {
    //g_testthetax += thetaxrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetax " << g_testthetax << std::endl;
    
    Mat3d rotx;
    rotx << 1,0,0,
    0,cos(thetaxrate*getDt()),-sin(thetaxrate*getDt()),
    0,sin(thetaxrate*getDt()),cos(thetaxrate*getDt());
    assert( approxEq(rotx.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_mesh, rotx );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], rotx );
  }
  // Settle for a second
  else if( getTime() <= 5.5 )
  {
    // Do nothing
  }
  else if( getTime() <= 6.5 )
  {
    //g_testthetax -= thetaxrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetax " << g_testthetax << std::endl;
    
    
    Mat3d rotx;
    rotx << 1,0,0,
    0,cos(-thetaxrate*getDt()),-sin(-thetaxrate*getDt()),
    0,sin(-thetaxrate*getDt()),cos(-thetaxrate*getDt());    
    assert( approxEq(rotx.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_mesh, rotx );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], rotx );
  }
  // Settle for a second
  else if( getTime() <= 7.0 )
  {
    // Do nothing
  }
  /////////////////////////////////////////////////////////////////////////////
  // "Shake" the sphere about the y axis
  else if( getTime() <= 8.0 )
  {
    //g_testthetay += thetayrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetay " << g_testthetay << std::endl;
    
    Mat3d roty;
    roty << cos(thetayrate*getDt()),0,sin(thetayrate*getDt()),
    0,1,0,
    -sin(thetayrate*getDt()),0,cos(thetayrate*getDt());
    assert( approxEq(roty.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_mesh, roty );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], roty );
  }
  // Settle for a second
  else if( getTime() <= 8.5 )
  {
    // Do nothing
  }
  else if( getTime() <= 9.5 )
  {
    //g_testthetay -= thetayrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetay " << g_testthetay << std::endl;
    
    Mat3d roty;
    roty << cos(-thetayrate*getDt()),0,sin(-thetayrate*getDt()),
    0,1,0,
    -sin(-thetayrate*getDt()),0,cos(-thetayrate*getDt());
    assert( approxEq(roty.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_mesh, roty );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], roty );
  }
  // Settle for a few seconds
  else if( getTime() <= 10.5 )
  {
    // Do nothing
  }
  else 
  {
    std::cout << "Simulation complete. Exiting." << std::endl;
    exit(0);
  }
  
  return true;
}

//void SphereRotator::setTime( double time )
//{
//  m_time = time;
//}

//void SphereRotator::setDt( double dt )
//{
//  m_dt = dt;
//}

//double SphereRotator::getTime()
//{
//  return m_time;
//}

//double SphereRotator::getDt()
//{
//  return m_dt;
//}

void SphereRotatorW::transformTriangleObject( TriangleMesh& tri_mesh, const Mat3d& transformation )
{
  for( TriangleMesh::vertex_iter itr = tri_mesh.vertices_begin(); itr != tri_mesh.vertices_end(); ++itr )
  {
    tri_mesh.getVertex(*itr) = transformation*tri_mesh.getVertex(*itr);
  } 
}

void SphereRotatorW::transformRodRoot( ElasticRod* rod, const Mat3d& transformation )
{
  RodBoundaryCondition* boundary = rod->getBoundaryCondition();
  rod->setVertex(0,transformation*rod->getVertex(0));
  rod->setVertex(1,transformation*rod->getVertex(1));
  boundary->setDesiredVertexPosition(0, rod->getVertex(0));
  boundary->setDesiredVertexPosition(1, rod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);

  //Vec3d test0=rod->getVertex(0);
  //Vec3d test1=rod->getVertex(1);
  Vec3d test10=rod->getVertex(10);
  Vec3d test19=rod->getVertex(19);
  Vec3d test39=rod->getVertex(39);
  //std::cout << "rod vertex 0  position "<<test0(0)<< std::endl;
  //std::cout << "rod vertex 0  position "<< rod->getVertex(0)<< std::endl;
  //std::cout << "rod vertex 1  position "<< rod->getVertex(1)<< std::endl;
  //std::cout << "rod vertex 10  position "<< rod->getVertex(10)<< std::endl;
  //sleep(1)
  //fprintf(m_pFile, "this is a test");
  float v10_0= roundto( test10(0), 10000.0 );
  float v10_1= roundto( test10(1), 10000.0 );
  float v10_2= roundto( test10(2), 10000.0 );
  float v19_0= roundto( test19(0), 10000.0 );
  float v19_1= roundto( test19(1), 10000.0 );
  float v19_2= roundto( test19(2), 10000.0 );
  float v39_0= roundto( test39(0), 10000.0 );
  float v39_1= roundto( test39(1), 10000.0 );
  float v39_2= roundto( test39(2), 10000.0 );
  fprintf(m_pFile, "%9.4f %9.4f %9.4f\n", v10_0, v10_1, v10_2);
  fprintf(m_pFile, "%9.4f %9.4f %9.4f\n", v19_0, v19_1, v19_2);
  fprintf(m_pFile, "%9.4f %9.4f %9.4f\n", v39_0, v39_1, v39_2);

}

float SphereRotatorW::roundto(float val, float prec) { return floor( prec * val + 0.5 ) / prec; }




HairyBallWeta::HairyBallWeta()
: Problem("Hairy Ball Weta", "An sphere with a number of elastic rods protruding from the surface.")
, m_tri_meshes()
, m_rods()
, m_scripting_controllers()
, m_steppers()
, m_pr_stepper(NULL)
, m_br_stepper(NULL)
{
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();
  
  // Set Defaults for built-in options
  GetStringOpt("integrator") = "implicit";
  // Gravity in CGS
  GetVecOpt("gravity") = Vec3d(0,-981,0);
  // Timestep
  GetScalarOpt("dt") = 0.001;
  // Maximum number of implicit solver itations
  GetIntOpt("iterations") = 100;
  // Number of vertices in each rod
  GetIntOpt("nv") = 40;
  // Assume a quasistatic material frame
  GetBoolOpt("quasistatic") = false;
  // Scale the radii of hairs for rendering/contact so we can see them
  GetScalarOpt("radius-scale") = 10.0;
  
  GetScalarOpt("mass-damping") = 0;
  GetScalarOpt("viscosity") = 0;
  
  AddOption("hairrootfile","file with locations (specified as point on sphere) to grow hair from","assets/hairyball/32normals.txt");
  //AddOption("hairrootfile","file with locations (specified as point on sphere) to grow hair from","assets/hairyball/100normals.txt");
  AddOption("curlyhair","curly if true, straight if false",false);
  AddOption("enablecollisions","true if collisions enabled, false otherwise",true);
  AddOption("adaptivestepper","true if adaptive timestepping should be used, failse otherwise",true);
}

HairyBallWeta::~HairyBallWeta()
{
  assert( m_rods.size() == m_steppers.size() );

  for( int i = 0; i < (int) m_tri_meshes.size(); ++i )
  {
    assert( m_tri_meshes[i] != NULL );
    delete m_tri_meshes[i];
    m_tri_meshes[i] = NULL;
  }

  for( int i = 0; i < (int) m_rods.size(); ++i )
  {
    assert( m_rods[i] != NULL );
    delete m_rods[i];
    m_rods[i] = NULL;
  }

  for( int i = 0; i < (int) m_scripting_controllers.size(); ++i )
  {
    assert( m_scripting_controllers[i] != NULL );
    delete m_scripting_controllers[i];
    m_scripting_controllers[i] = NULL;
  }  
  
  for( int i = 0; i < (int) m_steppers.size(); ++i )
  {
    assert( m_steppers[i] != NULL );
    delete m_steppers[i];
    m_steppers[i] = NULL;
  }  

  if( m_pr_stepper != NULL )
  {
    delete m_pr_stepper;
    m_pr_stepper = NULL;
  }
  
  if( m_br_stepper != NULL )
  {
    delete m_br_stepper;
    m_br_stepper = NULL;
  }
}

void HairyBallWeta::loadHandCodedVectors( std::vector<Vec3d>& normals )
{
  normals.clear();
  
  Vec3d nhat00(1,1,1);    nhat00.normalize(); normals.push_back(nhat00);
  Vec3d nhat01(1,1,-1);   nhat01.normalize(); normals.push_back(nhat01);
  Vec3d nhat02(1,-1,1);   nhat02.normalize(); normals.push_back(nhat02);
  Vec3d nhat03(-1,1,1);   nhat03.normalize(); normals.push_back(nhat03);
  Vec3d nhat04(1,-1,-1);  nhat04.normalize(); normals.push_back(nhat04);
  Vec3d nhat05(-1,1,-1);  nhat05.normalize(); normals.push_back(nhat05);
  Vec3d nhat06(-1,-1,1);  nhat06.normalize(); normals.push_back(nhat06);
  Vec3d nhat07(-1,-1,-1); nhat07.normalize(); normals.push_back(nhat07);
  
  
  Vec3d nhat08(1,0,0);    nhat08.normalize(); normals.push_back(nhat08);
  Vec3d nhat09(0,1,0);    nhat09.normalize(); normals.push_back(nhat09);
  Vec3d nhat10(0,0,1);    nhat10.normalize(); normals.push_back(nhat10);
  
  Vec3d nhat11(-1,0,0);   nhat11.normalize(); normals.push_back(nhat11);
  Vec3d nhat12(0,-1,0);   nhat12.normalize(); normals.push_back(nhat12);
  Vec3d nhat13(0,0,-1);   nhat13.normalize(); normals.push_back(nhat13);
  
  
  Vec3d nhat14(1,1,0);    nhat14.normalize(); normals.push_back(nhat14);
  Vec3d nhat15(-1,1,0);   nhat15.normalize(); normals.push_back(nhat15);
  Vec3d nhat16(1,-1,0);   nhat16.normalize(); normals.push_back(nhat16);
  Vec3d nhat17(-1,-1,0);  nhat17.normalize(); normals.push_back(nhat17);
  
  Vec3d nhat18(1,0,1);    nhat18.normalize(); normals.push_back(nhat18);
  Vec3d nhat19(-1,0,1);   nhat19.normalize(); normals.push_back(nhat19);
  Vec3d nhat20(1,0,-1);   nhat20.normalize(); normals.push_back(nhat20);
  Vec3d nhat21(-1,0,-1);  nhat21.normalize(); normals.push_back(nhat21);
  
  Vec3d nhat22(0,1,1);    nhat22.normalize(); normals.push_back(nhat22);
  Vec3d nhat23(0,-1,1);   nhat23.normalize(); normals.push_back(nhat23);
  Vec3d nhat24(0,1,-1);   nhat24.normalize(); normals.push_back(nhat24);
  Vec3d nhat25(0,-1,-1);  nhat25.normalize(); normals.push_back(nhat25);
}

void HairyBallWeta::loadNormalsFile( const std::string& filename, std::vector<Vec3d>& normals )
{
  normals.clear();
  
  std::ifstream infile(filename.c_str());
  if( !infile.is_open() )
  {
    std::cerr << "\033[31;1mERROR IN HAIRYBALL:\033[m Failed to open file with normals " << filename << std::endl;
    exit(0);
  }

  while( !infile.eof() )
  {
    Vec3d newnormal;
    if( (infile >> newnormal.x() >> newnormal.y() >> newnormal.z() >> std::ws) )
    {
      assert( newnormal.norm() != 0.0 );
      newnormal.normalize();
      normals.push_back(newnormal);
    }
    else
    {
      std::cerr << "\033[31;1mERROR IN HAIRYBALL:\033[m Failed to parse normal " << normals.size()+1 << " in file " << filename << std::endl;
      exit(0);
    }
  }

  #ifdef DEBUG
  for( int i = 0; i < (int) normals.size(); ++i ) assert(approxEq(normals[i].norm(),1.0,1.0e-6));
  for( int i = 0; i < (int) normals.size(); ++i ) for( int j = i+1; j < (int) normals.size(); ++j ) assert( !approxEq(normals[i],normals[j],1.0e-6) );
  #endif

  std::cout << "\033[35;1mHAIRYBALL MESSAGE:\033[m Successfully parsed: " << filename << std::endl;
}

void HairyBallWeta::generateStraightHair( const Vec3d& initnorm, const Vec3d& startpoint, const double& dL, const int& nv, std::vector<Vec3d>& vertices )
{
  vertices.push_back(startpoint);
  for( int j = 0; j < nv-1; ++j ) vertices.push_back(vertices.back()+dL*initnorm);
}

void HairyBallWeta::generateCurlyHair( const Vec3d& initnorm, const Vec3d& startpoint, const double& dL, const int& nv, std::vector<Vec3d>& vertices )
{
  // Generate an orthonormal frame
  Vec3d p1;
  findOrthogonal(p1,initnorm);
  Vec3d p2 = p1.cross(initnorm);
  
  double radius = 0.25;
  double curl_scale = 2.0;
  
  vertices.push_back(startpoint);
  vertices.push_back(startpoint+dL*initnorm);
  for( int j = 2; j < nv; ++j ) 
  {
    vertices.push_back(startpoint + j*dL*initnorm + radius*cos(j*curl_scale*dL)*p1 + radius*sin(j*curl_scale*dL)*p2 );
  }

}


void HairyBallWeta::Setup()
{
  assert( m_rods.size() == 0 );
  assert( m_steppers.size() == 0 );

  loadDynamicsProps();

  //RodTimeStepper::Method integrator;
  //if( GetStringOpt("integrator") == "implicit" ) integrator = RodTimeStepper::IMPL_EULER;
  //else assert( false );

  RodOptions opts;
  getRodOptions(opts);

  // cm
  double sphereradius = 4.0;

  // 12 inches == 60.96 cm
  double L = 60.96/4.0;

  //  double dx = 0.3;
  //  double dy = 0.3;
  double dL = L/((double)opts.numVertices-1);

  std::vector<Vec3d> initialnormals;
  loadNormalsFile( GetStringOpt("hairrootfile"), initialnormals );

  //std::vector<std::string> rod_labels;
  for( int i = 0; i < (int) initialnormals.size(); ++i )
  {
    //if( i > 30 ) continue;
        
    assert( approxEq(initialnormals[i].norm(),1.0,1.0e-6) );
    std::vector<Vec3d> vertices;
    if( !GetBoolOpt("curlyhair") ) generateStraightHair( initialnormals[i], sphereradius*initialnormals[i], dL, opts.numVertices, vertices );
    else generateCurlyHair( initialnormals[i], sphereradius*initialnormals[i], dL, opts.numVertices, vertices );
    //vertices.push_back(sphereradius*initialnormals[i]);
    //for( int j = 0; j < opts.numVertices-1; ++j ) vertices.push_back(vertices.back()+dL*initialnormals[i]);
    assert( (int) vertices.size() == opts.numVertices );


    //std::cout << initialnormals[i].transpose() << std::endl;
    //for( int j = 0; j < vertices.size(); ++j ) std::cout << vertices[j].transpose() << std::endl;


    ElasticRod* newrod = setupRod(opts, vertices, vertices);
    RodBoundaryCondition* boundary = newrod->getBoundaryCondition();
    boundary->setDesiredVertexPosition(0, newrod->getVertex(0));
    boundary->setDesiredVertexPosition(1, newrod->getVertex(1));
    boundary->setDesiredEdgeAngle(0,0.0);
    
    m_rods.push_back(newrod);
    
    //const ElasticRod::RodForces& frcs = newrod->getForces();
    //for( int q = 0; q < (int) frcs.size(); ++q ) std::cout << frcs[q]->getName() << std::endl;  
    
    RodTimeStepper* tstep = getRodTimeStepper(*newrod);
    if( !GetBoolOpt("adaptivestepper") ) 
    {
      m_steppers.push_back(tstep);
    }
    else
    {
      AdaptiveBinaryStepper* adpstep = new AdaptiveBinaryStepper( newrod, tstep );
      m_steppers.push_back(adpstep);
    }
    
    //rod_labels.push_back("rod_"+toString(i));
  }

  
  std::cout << "\033[35;1mHAIRYBALL MESSAGE:\033[m Simulating " << m_rods.size() << " rods." << std::endl;
  
  // Create a sphere in the center of the screen
  ObjParser objparser;
  TriangleMesh* tri_mesh = new TriangleMesh();
  objparser.loadTriangularMesh( "assets/UniformlyTriangulatedSphere.obj", *tri_mesh );
  // Compute the sphere's center
  Vec3d center(0.0,0.0,0.0);
  for( TriangleMesh::vertex_iter itr = tri_mesh->vertices_begin(); itr != tri_mesh->vertices_end(); ++itr )
  {
    center += tri_mesh->getVertex(*itr);
  }
  center /= ((double)tri_mesh->nv());
  // Translate the sphere so its center is at the origin
  for( TriangleMesh::vertex_iter itr = tri_mesh->vertices_begin(); itr != tri_mesh->vertices_end(); ++itr )
  {
    tri_mesh->getVertex(*itr) -= center;
  }
  // Make the sphere the desired radius
  double meshradius = (tri_mesh->getVertex(*tri_mesh->vertices_begin())).norm();
  double meshscale = sphereradius/meshradius;
  for( TriangleMesh::vertex_iter itr = tri_mesh->vertices_begin(); itr != tri_mesh->vertices_end(); ++itr )
  {
    tri_mesh->getVertex(*itr) *= meshscale;
  } 
  m_tri_meshes.push_back(tri_mesh);

  m_rotator = new SphereRotatorW(*tri_mesh,m_rods,getTime(),GetScalarOpt("dt"));
  m_scripting_controllers.push_back(m_rotator);

  for( int i = 0; i < (int) m_rods.size(); ++i ) m_world->addObject(m_rods[i]);
  for( int i = 0; i < (int) m_tri_meshes.size(); ++i ) m_world->addObject(m_tri_meshes[i]);

  if( !GetBoolOpt("enablecollisions") )
  {
    std::cout << "\033[35;1mHAIRYBALL MESSAGE:\033[m No collisions in use." << std::endl;
    m_pr_stepper = new ParallelStepper;
    for( int i = 0; i < (int) m_steppers.size(); ++i ) m_pr_stepper->addController(m_steppers[i]);
    for( int i = 0; i < (int) m_scripting_controllers.size(); ++i ) m_pr_stepper->addScriptingController(m_scripting_controllers[i]);
    m_world->addController(m_pr_stepper);
  }
  else
  {
    std::cout << "\033[35;1mHAIRYBALL MESSAGE:\033[m Collisions in use." << std::endl;
    m_br_stepper = new BridsonStepper( m_rods, m_tri_meshes, m_scripting_controllers, m_steppers, GetScalarOpt("dt") );
    //m_br_stepper->setRodLabels(rod_labels);
    m_world->addController(m_br_stepper);
  }

  assert( m_rods.size() == m_steppers.size() );

  
  // Add a few debugging "breakpoints"
  //insertBreakpoint(0.54);
  
  //insertBreakpoint(3.7);
  //insertBreakpoint(3.793);
  //insertBreakpoint(3.175);
}

void HairyBallWeta::AtEachTimestep()
{
  //if( approxEq(getTime(),0.318,1.0e-6) )
  //{
  //  GetScalarOpt("dt") -= 1.0e-4;
  //  std::cout << "CHANGING TIMESTEP TO: " << GetScalarOpt("dt") << std::endl;
  //  this->setDt(GetScalarOpt("dt"));
  //}

  //for( int i = 0; i < (int) m_steppers.size(); ++i ) m_steppers[i]->setTimeStep(GetScalarOpt("dt"));
  //if( m_br_stepper != NULL ) m_br_stepper->setDt(GetScalarOpt("dt"));
  //if( m_pr_stepper != NULL ) m_pr_stepper->setDt(GetScalarOpt("dt"));
  m_rotator->setDt(GetScalarOpt("dt"));
  m_rotator->setTime(getTime());
}

void HairyBallWeta::serialize( std::ofstream& of )
{
  // Serialize the parent class
  Problem::serializeProblem(of);
  
  TopologicalObjectSerializer objserializer;

  // Serialize the number of triangle meshes
  int numtrimeshes = (int) m_tri_meshes.size();
  assert( numtrimeshes == 1 );
  serializeVal(of,numtrimeshes);

  // Serialize each triangle mesh
  for( std::vector<TriangleMesh*>::size_type i = 0; i < m_tri_meshes.size(); ++i )
  {
    assert( m_tri_meshes[i] != NULL );
    objserializer.appendTopologicalObjectToFile(*m_tri_meshes[i],of);
  }

  // Serialize the number of rods
  int numrods = (int) m_rods.size();
  serializeVal(of,numrods);  

  // Serialize each rod
  for( std::vector<ElasticRod*>::size_type i = 0; i < m_rods.size(); ++i )
  {
    assert( m_rods[i] != NULL );
    objserializer.appendTopologicalObjectToFile(*m_rods[i],of);
  }
}

void HairyBallWeta::resumeFromfile( std::ifstream& ifs )
{
  Problem::resumeProblem(ifs);

  TopologicalObjectSerializer objserializer;

  // Load the number of triangle meshes
  int numtrimeshes = -1;
  loadVal(ifs,numtrimeshes);
  assert( numtrimeshes == 1 );

  std::cout << "About to load one mesh" << std::endl;
  
  // Load the serialized triangle meshes
  for( int i = 0; i < numtrimeshes; ++i )
  {
    TriangleMesh* loadedmesh = NULL;
    TopologicalObject* tempobj = NULL;
    objserializer.loadGenericTopologicalObject(&tempobj,ifs);
    loadedmesh = dynamic_cast<TriangleMesh*>(const_cast<TopologicalObject*>(tempobj));
    assert( loadedmesh != NULL );
    m_tri_meshes.push_back(loadedmesh);
  }
    
  // Load the number of rods
  int numrods = -1;
  loadVal(ifs,numrods);
  assert( numrods >= 0 );
  
  // Load the serialized rods
  for( int i = 0; i < numrods; ++i )
  {
    ElasticRod* erod = NULL;
    objserializer.loadTopologicalObjectFromFile(&erod,ifs);
    assert( erod != NULL );
    m_rods.push_back(erod);
  }
  
  // Create time steppers for each rod
  for( std::vector<ElasticRod*>::size_type i = 0; i < m_rods.size(); ++i )
  {
    RodTimeStepper* tstep = getRodTimeStepper(*m_rods[i]);
    if( !GetBoolOpt("adaptivestepper") ) 
    {
      m_steppers.push_back(tstep);
    }
    else
    {
      AdaptiveBinaryStepper* adpstep = new AdaptiveBinaryStepper( m_rods[i], tstep );
      m_steppers.push_back(adpstep);
    }
  }

  // Create a controller to rotate the sphere and the hair roots
  m_rotator = new SphereRotatorW(*m_tri_meshes[0],m_rods,getTime(),GetScalarOpt("dt"));
  m_scripting_controllers.push_back(m_rotator);

  for( int i = 0; i < (int) m_rods.size(); ++i ) m_world->addObject(m_rods[i]);
  for( int i = 0; i < (int) m_tri_meshes.size(); ++i ) m_world->addObject(m_tri_meshes[i]);
  
  if( !GetBoolOpt("enablecollisions") )
  {
    //std::cout << "\033[35;1mHAIRYBALL MESSAGE:\033[m No collisions in use." << std::endl;
    m_pr_stepper = new ParallelStepper;
    for( int i = 0; i < (int) m_steppers.size(); ++i ) m_pr_stepper->addController(m_steppers[i]);
    for( int i = 0; i < (int) m_scripting_controllers.size(); ++i ) m_pr_stepper->addScriptingController(m_scripting_controllers[i]);
    m_world->addController(m_pr_stepper);
  }
  else
  {
    //std::cout << "\033[35;1mHAIRYBALL MESSAGE:\033[m Collisions in use." << std::endl;
    m_br_stepper = new BridsonStepper( m_rods, m_tri_meshes, m_scripting_controllers, m_steppers, GetScalarOpt("dt") );
    m_world->addController(m_br_stepper);
  }
  
  assert( m_rods.size() == m_steppers.size() );
}



















