/**
 * \file HairyBall.cc
 *
 * \author smith@cs.columbia.edu
 * \date 05/11/2010
 */

#include "HairyBall.hh"

HairyBall::HairyBall()
: Problem("Hairy Ball", "An sphere with a number of elastic rods protruding from the surface.")
, m_tri_meshes()
, m_rods()
, m_steppers()
, m_pr_stepper(NULL)
{
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();
  
  // Set Defaults for built-in options
  GetStringOpt("integrator") = "implicit";
  // Gravity in CGS
  GetVecOpt("gravity") = Vec3d(0,-981,0);
  // Timestep
  GetScalarOpt("dt") = 0.01;
  // Maximum number of implicit solver itations
  GetIntOpt("iterations") = 1000;
  // Number of vertices in each rod
  GetIntOpt("nv") = 40;
  // Assume a quasistatic material frame
  GetBoolOpt("quasistatic") = true;
  // Scale the radii of hairs for rendering/contact so we can see them
  GetScalarOpt("radius-scale") = 10.0;
  
  AddOption("hairrootfile","file with locations (specified as point on sphere) to grow hair from","assets/hairyball/32normals.txt");
  AddOption("curlyhair","curly if true, straight if false",false);
}

HairyBall::~HairyBall()
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
  
  //if( m_br_stepper != NULL )
  //{
  //  delete m_br_stepper;
  //  m_br_stepper = NULL;
  //}
}

void HairyBall::loadHandCodedVectors( std::vector<Vec3d>& normals )
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

void HairyBall::loadNormalsFile( const std::string& filename, std::vector<Vec3d>& normals )
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
  
  std::cout << "Successfully parsed: " << filename << std::endl;
  std::cout << "   Loaded " << normals.size() << " normal directions" << std::endl;
}

void HairyBall::generateStraightHair( const Vec3d& initnorm, const Vec3d& startpoint, const double& dL, const int& nv, std::vector<Vec3d>& vertices )
{
  vertices.push_back(startpoint);
  for( int j = 0; j < nv-1; ++j ) vertices.push_back(vertices.back()+dL*initnorm);
}

void HairyBall::generateCurlyHair( const Vec3d& initnorm, const Vec3d& startpoint, const double& dL, const int& nv, std::vector<Vec3d>& vertices )
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


void HairyBall::Setup()
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

  m_pr_stepper = new ParallelStepper;
  
  std::vector<Vec3d> initialnormals;
  loadNormalsFile( GetStringOpt("hairrootfile"), initialnormals );

  for( int i = 0; i < (int) initialnormals.size(); ++i )
  {
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
    m_world->addObject(newrod);
    
    RodTimeStepper* tstep = getRodTimeStepper(*newrod);
    m_steppers.push_back(tstep);
    m_pr_stepper->addController(tstep);
  }
  
  
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
  m_world->addObject(tri_mesh);
  //m_br_stepper->addTriangleMesh(tri_mesh);
  
  //m_world->addController(m_br_stepper);
  //m_br_stepper->prepareForExecution();

  m_world->addController(m_pr_stepper);

  assert( m_rods.size() == m_steppers.size() );
}

void HairyBall::transformTriangleObject( TriangleMesh* tri_mesh, const Mat3d& transformation )
{
  for( TriangleMesh::vertex_iter itr = tri_mesh->vertices_begin(); itr != tri_mesh->vertices_end(); ++itr )
  {
    tri_mesh->getVertex(*itr) = transformation*tri_mesh->getVertex(*itr);
  } 
}

void HairyBall::transformRodRoot( ElasticRod* rod, const Mat3d& transformation )
{
  RodBoundaryCondition* boundary = rod->getBoundaryCondition();
  boundary->setDesiredVertexPosition(0, transformation*rod->getVertex(0));
  boundary->setDesiredVertexPosition(1, transformation*rod->getVertex(1));
  boundary->setDesiredEdgeAngle(0,0.0);
}

void HairyBall::AtEachTimestep()
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
  if( getTime() <= 2.0 )
  {
    // Do nothing
  }
  /////////////////////////////////////////////////////////////////////////////
  // "Shake" the sphere about the z axis
  else if( getTime() <= 3.0 )
  {
    //g_testthetaz += thetazrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetaz " << g_testthetaz << std::endl;
    
    Mat3d rotz;
    rotz << cos(thetazrate*GetScalarOpt("dt")),-sin(thetazrate*GetScalarOpt("dt")),0,
    sin(thetazrate*GetScalarOpt("dt")),cos(thetazrate*GetScalarOpt("dt")),0,
    0,0,1;
    assert( approxEq(rotz.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_tri_meshes[0], rotz );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], rotz );
  }
  // Settle for a second
  else if( getTime() <= 4.0 )
  {
    // Do nothing
  }
  else if( getTime() <= 5.0 )
  {
    //g_testthetaz -= thetazrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetaz " << g_testthetaz << std::endl;
    
    Mat3d rotz;
    rotz << cos(-thetazrate*GetScalarOpt("dt")),-sin(-thetazrate*GetScalarOpt("dt")),0,
    sin(-thetazrate*GetScalarOpt("dt")),cos(-thetazrate*GetScalarOpt("dt")),0,
    0,0,1;
    assert( approxEq(rotz.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_tri_meshes[0], rotz );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], rotz );
  }
  // Settle for a second
  else if( getTime() <= 6.0 )
  {
    // Do nothing
  }
  /////////////////////////////////////////////////////////////////////////////
  // "Shake" the sphere about the x axis
  else if( getTime() <= 7.0 )
  {
    //g_testthetax += thetaxrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetax " << g_testthetax << std::endl;
    
    Mat3d rotx;
    rotx << 1,0,0,
    0,cos(thetaxrate*GetScalarOpt("dt")),-sin(thetaxrate*GetScalarOpt("dt")),
    0,sin(thetaxrate*GetScalarOpt("dt")),cos(thetaxrate*GetScalarOpt("dt"));
    assert( approxEq(rotx.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_tri_meshes[0], rotx );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], rotx );
  }
  // Settle for a second
  else if( getTime() <= 8.0 )
  {
    // Do nothing
  }
  else if( getTime() <= 9.0 )
  {
    //g_testthetax -= thetaxrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetax " << g_testthetax << std::endl;
    
    
    Mat3d rotx;
    rotx << 1,0,0,
    0,cos(-thetaxrate*GetScalarOpt("dt")),-sin(-thetaxrate*GetScalarOpt("dt")),
    0,sin(-thetaxrate*GetScalarOpt("dt")),cos(-thetaxrate*GetScalarOpt("dt"));    
    assert( approxEq(rotx.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_tri_meshes[0], rotx );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], rotx );
  }
  // Settle for a second
  else if( getTime() <= 10.0 )
  {
    // Do nothing
  }
  /////////////////////////////////////////////////////////////////////////////
  // "Shake" the sphere about the y axis
  else if( getTime() <= 11.0 )
  {
    //g_testthetay += thetayrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetay " << g_testthetay << std::endl;
    
    Mat3d roty;
    roty << cos(thetayrate*GetScalarOpt("dt")),0,sin(thetayrate*GetScalarOpt("dt")),
    0,1,0,
    -sin(thetayrate*GetScalarOpt("dt")),0,cos(thetayrate*GetScalarOpt("dt"));
    assert( approxEq(roty.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_tri_meshes[0], roty );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], roty );
  }
  // Settle for a second
  else if( getTime() <= 12.0 )
  {
    // Do nothing
  }
  else if( getTime() <= 13.0 )
  {
    //g_testthetay -= thetayrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetay " << g_testthetay << std::endl;
    
    Mat3d roty;
    roty << cos(-thetayrate*GetScalarOpt("dt")),0,sin(-thetayrate*GetScalarOpt("dt")),
    0,1,0,
    -sin(-thetayrate*GetScalarOpt("dt")),0,cos(-thetayrate*GetScalarOpt("dt"));
    assert( approxEq(roty.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_tri_meshes[0], roty );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], roty );
  }
  // Settle for a few seconds
  else if( getTime() <= 20.0 )
  {
    // Do nothing
  }
  else 
  {
    std::cout << "Simulation complete. Exiting." << std::endl;
    exit(0);
  }
}

void HairyBall::AtEachTimestepOLD()
{
  double thetax = pi/4.0;
  double thetay = pi/4.0;
  double thetaz = pi/4.0;

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
  if( getTime() <= 2.0 )
  {
    // Do nothing
  }
  /////////////////////////////////////////////////////////////////////////////
  // "Shake" the sphere about the z axis
  else if( getTime() <= 3.0 )
  {
    //g_testthetaz += thetazrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetaz " << g_testthetaz << std::endl;

    Mat3d rotz;
    rotz << cos(thetazrate*GetScalarOpt("dt")),-sin(thetazrate*GetScalarOpt("dt")),0,
            sin(thetazrate*GetScalarOpt("dt")),cos(thetazrate*GetScalarOpt("dt")),0,
            0,0,1;
    assert( approxEq(rotz.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_tri_meshes[0], rotz );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], rotz );
  }
  // Settle for a second
  else if( getTime() <= 4.0 )
  {
    // Do nothing
  }
  else if( getTime() <= 5.0 )
  {
    //g_testthetaz -= thetazrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetaz " << g_testthetaz << std::endl;

    Mat3d rotz;
    rotz << cos(-thetazrate*GetScalarOpt("dt")),-sin(-thetazrate*GetScalarOpt("dt")),0,
            sin(-thetazrate*GetScalarOpt("dt")),cos(-thetazrate*GetScalarOpt("dt")),0,
            0,0,1;
    assert( approxEq(rotz.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_tri_meshes[0], rotz );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], rotz );
  }
  // Settle for a second
  else if( getTime() <= 6.0 )
  {
    // Do nothing
  }
  /////////////////////////////////////////////////////////////////////////////
  // "Shake" the sphere about the x axis
  else if( getTime() <= 7.0 )
  {
    //g_testthetax += thetaxrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetax " << g_testthetax << std::endl;
    
    Mat3d rotx;
    rotx << 1,0,0,
            0,cos(thetaxrate*GetScalarOpt("dt")),-sin(thetaxrate*GetScalarOpt("dt")),
            0,sin(thetaxrate*GetScalarOpt("dt")),cos(thetaxrate*GetScalarOpt("dt"));
    assert( approxEq(rotx.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_tri_meshes[0], rotx );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], rotx );
  }
  // Settle for a second
  else if( getTime() <= 8.0 )
  {
    // Do nothing
  }
  else if( getTime() <= 9.0 )
  {
    //g_testthetax -= thetaxrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetax " << g_testthetax << std::endl;
    
    
    Mat3d rotx;
    rotx << 1,0,0,
            0,cos(-thetaxrate*GetScalarOpt("dt")),-sin(-thetaxrate*GetScalarOpt("dt")),
            0,sin(-thetaxrate*GetScalarOpt("dt")),cos(-thetaxrate*GetScalarOpt("dt"));    
    assert( approxEq(rotx.determinant(),1.0,1.0e-6) );

    transformTriangleObject( m_tri_meshes[0], rotx );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], rotx );
  }
  // Settle for a second
  else if( getTime() <= 10.0 )
  {
    // Do nothing
  }
  /////////////////////////////////////////////////////////////////////////////
  // "Shake" the sphere about the y axis
  else if( getTime() <= 11.0 )
  {
    //g_testthetay += thetayrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetay " << g_testthetay << std::endl;
    
    Mat3d roty;
    roty << cos(thetayrate*GetScalarOpt("dt")),0,sin(thetayrate*GetScalarOpt("dt")),
            0,1,0,
            -sin(thetayrate*GetScalarOpt("dt")),0,cos(thetayrate*GetScalarOpt("dt"));
    assert( approxEq(roty.determinant(),1.0,1.0e-6) );
    
    transformTriangleObject( m_tri_meshes[0], roty );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], roty );
  }
  // Settle for a second
  else if( getTime() <= 12.0 )
  {
    // Do nothing
  }
  else if( getTime() <= 13.0 )
  {
    //g_testthetay -= thetayrate*GetScalarOpt("dt");
    //std::cout << "At time: " << getTime() << "     thetay " << g_testthetay << std::endl;
    
    Mat3d roty;
    roty << cos(-thetayrate*GetScalarOpt("dt")),0,sin(-thetayrate*GetScalarOpt("dt")),
            0,1,0,
            -sin(-thetayrate*GetScalarOpt("dt")),0,cos(-thetayrate*GetScalarOpt("dt"));
    assert( approxEq(roty.determinant(),1.0,1.0e-6) );

    transformTriangleObject( m_tri_meshes[0], roty );
    for( int i = 0; i < (int) m_rods.size(); ++i ) transformRodRoot( m_rods[i], roty );
  }
  // Settle for a few seconds
  else if( getTime() <= 16.0 )
  {
    // Do nothing
  }
  else 
  {
    std::cout << "Simulation complete. Exiting." << std::endl;
    exit(0);
  }
}
