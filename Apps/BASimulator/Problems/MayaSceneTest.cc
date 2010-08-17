/**
 * \file MayaSceneTest.hh
 *
 * \author 
 * \date 
 */

#include "MayaSceneTest.hh"

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

MayaSceneTest::MayaSceneTest()
  : Problem("Maya Scene Test", "Maya Scene Test.")
, m_rods()
, m_controllers()
{
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();
  
  collision_method = WETA;
  
  //AddOption("rodfile-path", "rodfile-path", "/usr/home/jjoo/scenes/mau160.");
  AddOption("rodfile-path", "rodfile-path", "/local1/jjoo/scenes/sac2/sac.");
  AddOption("rodfile-firstindex", "rodfile-firstindex", -39);
  AddOption("rod-input-rodnum", "rod-input-rodnum", 417);
  AddOption("maya-substeps", "maya-substeps", 20);
  
  AddOption("reverse", "reverse", false);

  AddOption("load-meshes", "load-meshes", false);
  AddOption("meshfile-path", "meshfile-path", "");
  AddOption("collisionsEnabled", "i_collisionsEnabled", false);
  
  AddOption("selected-running", "selected-running", false);
  AddOption("selected-rod", "selected-rod", 0);
  
  AddOption("impulse_enabled", "impulse_enabled", true);
  
  AddOption("figaro-friction", "figaro-friction", 0.0);
  AddOption("figaro-thickness", "figaro-thickness", 0.1);
  AddOption("figaro-full-collision", "figaro-full-collision", false);
  AddOption("figaro-separationStrength", "figaro-separationStrength", 40.0);
  AddOption("figaro-cor", "figaro-cor", 0.1);
  
  AddOption("collision_method", "0 - weta, 1 - columbia", 0);
  
  GetScalarOpt("viscosity") = 10;
  GetScalarOpt("mass-damping") = 10;
  
  
  AddOption("selected-running", "selected-running", false);
  AddOption("selected-rod", "selected-rod", 0);
  
  br_mesh = NULL;

  /*
  is_selected_running = true;
  selected_rods.push_back(130);
  selected_rods.push_back(167);
  selected_rods.push_back(184);
  selected_rods.push_back(383);
  selected_rods.push_back(397);
*/
  
}

MayaSceneTest::~MayaSceneTest()
{
  for( int i = 0; i < (int) m_rods.size(); ++i )
  {
    assert( m_rods[i] != NULL );
    delete m_rods[i];
    m_rods[i] = NULL;
  }
  
  for( int i = 0; i < (int) m_controllers.size(); ++i )
  {
    assert( m_controllers[i] != NULL );
    delete m_controllers[i];
    m_controllers[i] = NULL;
  }
}


FILE* MayaSceneTest::readNumberOfRodsFromFile( const std::string i_cacheFilename, size_t& o_numRodsInFile)
{
  FILE *fp = NULL;

  fp = fopen( i_cacheFilename.c_str(), "r" );
  if ( fp == NULL )
  {
    return NULL; 
  }

  std::cerr << "Reading file " << i_cacheFilename << std::endl;

  size_t ret = fread( &o_numRodsInFile, sizeof( size_t ), 1, fp );

  return fp;
}

void MayaSceneTest::LoadNextFrameFromCache() {
  // Copy data from NEXT to PREV
  prevBoundaryCondition.clear();
  prevBoundaryCondition = nextBoundaryCondition;
  nextBoundaryCondition.clear();
  
  prevMayaPosition.clear();
  prevMayaPosition = nextMayaPosition;
  nextMayaPosition.clear();

  // Load data for NEXT
  std::stringstream names;
  std::string name;
  names << GetStringOpt ("rodfile-path") << GetIntOpt("rodfile-firstindex") + current_frame + 1<< ".fig";
  name = names.str();

  size_t nr = 0;
  int max_r = GetIntOpt("rod-input-rodnum");

  FILE *fp = readNumberOfRodsFromFile(name, nr);
  
  if (!fp) {
    std::cout << "end of record\n";
    
    nextBoundaryCondition = prevBoundaryCondition;
    nextMayaPosition = prevMayaPosition;
     
    return;
  }
  
  for ( size_t r=0; r<nr && r<(size_t)max_r; r++ )
  {
    size_t numVertices;
    size_t ret = fread( &numVertices, sizeof( size_t ), 1, fp );

    std::vector<Vec3d> vertices;
    VecXd fixedVerts(6);
    VecXd allVerts(numVertices * 3);

    //std::cout << numVertices <<"\n";
    
    for ( size_t v=0; v<numVertices; v++ )
    {
      Vec3d posvec;
      double pos[3];

      ret = fread( &pos[0], sizeof( double ), 3, fp );
      posvec = Vec3d(pos[0], pos[1], pos[2]);
//      vertices.push_back(posvec);
      allVerts.segment<3> (3 * v) = posvec;

      //std::cout << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
      
//      posvec = Vec3d(pos[0], pos[1], pos[2]);
//      vertices.push_back(posvec);
      
      ret = fread( &pos[0], sizeof( double ), 3, fp );
      //std::cout << pos[0] << " " << pos[1] << " " << pos[2] << "\n";
      if (v < 2) {
        fixedVerts(v*3+0) = pos[0];
        fixedVerts(v*3+1) = pos[1];
        fixedVerts(v*3+2) = pos[2];
      }
    }
    nextBoundaryCondition.push_back(fixedVerts);
    nextMayaPosition.push_back(allVerts);
  }
  fclose(fp);
  
  
  if (GetBoolOpt("load-meshes")) {
    size_t nv = 0;
    size_t nf = 0;
    std::vector<Vec3d> verts;
    std::vector<uint> faces;
  
    ReadMeshFile(mesh_frame, nv, nf, verts, faces);
  
    UpdateCollisionMesh(verts, faces);
    mesh_frame++;

  }
}

void MayaSceneTest::Setup()
{
  switch (GetIntOpt("collision_method")) {
    case 0:
      collision_method = WETA;
      break;
    case 1:
      collision_method = COLUMBIA;
      break;
  }
  
  if (GetBoolOpt("reverse")) {
    reverse_test();
    return;
  }

  impulse_enabled = GetBoolOpt("impulse_enabled");
  
  std::cout << "impulse_enabled " << impulse_enabled << "\n";
  
  nSubsteps = GetIntOpt("maya-substeps");
  fps = 24;
 
  current_frame = 0;
  current_step = 0;

  i_collisionsEnabled = GetBoolOpt("collisionsEnabled");
      

  is_selected_running = GetBoolOpt("selected-running");
  if (is_selected_running) {
    selected_rods.push_back(GetIntOpt("selected-rod"));
  }
  
  GetScalarOpt("dt") = 1.0 / (double) fps / (double)nSubsteps;
  
  //  GetVecOpt("gravity") = Vec3d(0, GetScalarOpt("vertical-gravity") ,0); // 981 cm/s^2 = 9.81 m/s^2
  
  loadDynamicsProps();
  RodOptions opts;
  
  getRodOptions(opts);
  
//  opts.viscosity = GetScalarOpt("viscosity");
  //opts.radiusScale = 1;
  
  std::stringstream names;
  std::string name;
  names << GetStringOpt ("rodfile-path") << GetIntOpt("rodfile-firstindex") << ".fig";
  name = names.str();
  
//	std::ifstream myfile;
	
//	myfile.open(name.c_str() , std::ifstream::in);

	size_t nr = 0;
	int max_r = GetIntOpt("rod-input-rodnum");

  FILE *fp = readNumberOfRodsFromFile(name, nr);
  
  for ( size_t r=0; r<nr && r<(size_t)max_r; r++ )
  {
    size_t numVertices;
    size_t ret = fread( &numVertices, sizeof( size_t ), 1, fp );

    std::vector<Vec3d> vertices;
    VecXd fixedVerts(6);
    VecXd allVerts(numVertices * 3);

//    o_rodVertices[ r ].resize( numVertices );
//    o_unsimulatedRodVertices[ r ].resize( numVertices );

    for ( size_t v=0; v<numVertices; v++ )
    {
      Vec3d posvec;
      double pos[3];

      ret = fread( &pos[0], sizeof( double ), 3, fp );
      posvec = Vec3d(pos[0], pos[1], pos[2]);
      vertices.push_back(posvec);
      allVerts.segment<3> (3 * v) = posvec;
      
//      o_rodVertices[ r ][ v ] = Vec3d( pos[0], pos[1], pos[2] );
      
      ret = fread( &pos[0], sizeof( double ), 3, fp );
      if (v < 2) {
        fixedVerts(v*3+0) = pos[0];
        fixedVerts(v*3+1) = pos[1];
        fixedVerts(v*3+2) = pos[2];
      }
//      o_unsimulatedRodVertices[ r ][ v ] = Vec3d( pos[0], pos[1], pos[2] );
    }
    
    nextBoundaryCondition.push_back(fixedVerts);
    nextMayaPosition.push_back(allVerts);
    
    bool found = false;
    if (is_selected_running) {
      for(int i=0; i<(int)selected_rods.size(); i++) {
        if (selected_rods[i] == r) {
          found = true;
          break;
        }
      }
      
      if (!found) continue;
    }
    opts.numVertices = (int)numVertices;

    ElasticRod* newrod = setupRod(opts, vertices, vertices);
    newrod->global_rodID = r;
    
//    if (numVertices < 3) 
//    std::cout << r << " rod created " << allVerts << "\n";
      
    newrod->getBoundaryCondition()->setDesiredVertexPosition(0, newrod->getVertex(0));
    newrod->getBoundaryCondition()->setDesiredVertexPosition(1, newrod->getVertex(1));
    newrod->getBoundaryCondition()->setDesiredEdgeAngle(0, newrod->getTheta(0));
      
    newrod->fixVert(0);
    newrod->fixVert(1);
    newrod->fixEdge(0);
    
//    if (newrod->nv() > 2) {
//      newrod->getBoundaryCondition()->setDesiredVertexPosition(2, newrod->getVertex(2));
//      newrod->getBoundaryCondition()->setDesiredEdgeAngle(1, newrod->getTheta(1));
//    }

    RodTimeStepper* stepper = new RodTimeStepper( *newrod );
    stepper->setDiffEqSolver( RodTimeStepper::SYM_IMPL_EULER );
  
    stepper->setMaxIterations(GetIntOpt("iterations"));
    stepper->set_stol(GetScalarOpt("stol"));
    stepper->set_atol(GetScalarOpt("atol"));
    stepper->set_rtol(GetScalarOpt("rtol"));
    stepper->set_inftol(GetScalarOpt("inftol"));

    // These get deleted by RodTimeStepper
    stepper->addExternalForce( new RodMassDamping( GetScalarOpt("mass-damping") ) );
    
//    if ( i_gravity.norm() > 0)
    {
      stepper->addExternalForce( new RodGravity( GetVecOpt("gravity") ) );
      //m_gravity = i_gravity;
    }

    AdaptiveBinaryStepper* adpstep = new AdaptiveBinaryStepper( newrod, stepper );
    //m_steppers.push_back(adpstep);
  
    RodCollisionTimeStepper *cstepper = new RodCollisionTimeStepper( adpstep, newrod );
    
    cstepper->impulse_enabled = impulse_enabled;
    
//    RodTimeStepper* newstepper = getRodTimeStepper(*newrod);
    m_rods.push_back(newrod);
    m_controllers.push_back(adpstep);
    m_world->addObject(newrod);
    m_collision_controllers.push_back(cstepper);
  }
  fclose(fp);
  
  InitializeMeshes();
  
  LoadNextFrameFromCache();
  
}

void MayaSceneTest::ReadMeshFile(int frame, size_t &numberOfVerts, size_t &numberOfFaces, std::vector<Vec3d> &verts, std::vector<uint> &faces)
{
  std::cout << "reading meshes\n";
  
  std::stringstream ss;
  std::string filename;
  
  ss << GetStringOpt("meshfile-path");
  ss << frame;
  ss >> filename;
  
  std::cout << filename << "\n";
  
  verts.clear();
  faces.clear();
  
  FILE *fp;
  fp = fopen( filename.c_str(), "r" );
  if ( fp == NULL )
  {
    std::cout << "can't open file\n";
    return;
  }

//    size_t numberOfVerts = 0;
//    size_t numberOfFaces = 0;
  
  size_t ret = fread( &numberOfVerts, sizeof( size_t ), 1, fp );
  ret = fread( &numberOfFaces, sizeof( size_t ), 1, fp );

  for ( size_t r=0; r<numberOfVerts; r++ )
  {
    double pos[3];
        
    fread( &pos[0], sizeof( double ), 3, fp );
    
    Vec3d p(pos[0], pos[1], pos[2]);
    verts.push_back(p);

  }
  
  for ( size_t r=0; r<numberOfFaces; r++ )
  {
    uint idx[3];
  
    fread( &idx[0], sizeof( uint ), 3, fp );
    
    faces.push_back(idx[0]);
    faces.push_back(idx[1]);
    faces.push_back(idx[2]);
  }
  
  fclose ( fp );
}


void MayaSceneTest::initialiseCollisionMesh( BASim::CollisionMeshData *collisionMeshData, size_t id )
{
  m_collisionMeshMap[ id ] = collisionMeshData;
  m_collisionMeshMap[ id ]->initialize();
}

void MayaSceneTest::UpdateCollisionMesh( std::vector<Vec3d> &points, std::vector<uint> &faces )
{    
    
  if ( m_collisionMeshData->currPositions.size() != points.size() ||
       m_collisionMeshData->oldPositions.size() != points.size() ||
       m_collisionMeshData->newPositions.size() != points.size() ||
       m_collisionMeshData->velocities.size() != points.size() ) {
    
    cerr << "Initialising size and data for collision mesh ( " << points.size() << " points\n";
        
    m_collisionMeshData->currPositions.resize( points.size() );
    m_collisionMeshData->prevPositions.resize( points.size() );
    m_collisionMeshData->oldPositions.resize( points.size() );
    m_collisionMeshData->newPositions.resize( points.size() );
    m_collisionMeshData->velocities.resize( points.size() );
       }

  
  m_collisionMeshData->triangleIndices.resize( faces.size() );
  
  for ( size_t i=0; i<faces.size(); ++i ) {
    m_collisionMeshData->triangleIndices[i] = faces[i];
  }

  if ( mesh_frame == 1 )
  {
    m_collisionMeshData->reset( points );
  }
  else
  {
    m_collisionMeshData->update( points, "", 0 );
  }

}


void MayaSceneTest::InitializeMeshes() {
  
  if (!GetBoolOpt("load-meshes")) return;
  
  mesh_frame = 1;
  
  size_t nv = 0;
  size_t nf = 0;
  std::vector<Vec3d> verts;
  std::vector<uint> faces;
  
  ReadMeshFile(mesh_frame, nv, nf, verts, faces);
  
  m_collisionMeshData = new CollisionMeshData;
  
  UpdateCollisionMesh(verts, faces);
  
  initialiseCollisionMesh(m_collisionMeshData, 0);
  
 
  Scalar m_friction = GetScalarOpt("figaro-friction");
  m_collisionMeshData->setFriction( m_friction );

  Scalar m_thickness = GetScalarOpt("figaro-thickness");
  m_collisionMeshData->setThickness( m_thickness );
        
  bool m_fullCollisions = GetBoolOpt("figaro-full-collision");
  m_collisionMeshData->setFullCollisions( m_fullCollisions );
        
  Scalar separationStrength = GetScalarOpt("figaro-separationStrength");
  m_collisionMeshData->setSeparationStrength( separationStrength );
        
  Scalar coefficientOfRestitution = GetScalarOpt("figaro-cor");
  m_collisionMeshData->setCoefficientOfRestitution( coefficientOfRestitution );

  
  
  MakeTriangleMesh(verts, faces);
  
  mesh_frame++;
//  ReadMeshFile(mesh_frame, nv, nf, verts, faces);
//  UpdateCollisionMesh(verts, faces);
  
  //  void setCollisionMeshesMap(CollisionMeshDataHashMap* collisionMeshes)
//m_collisionMeshMap
  //triangleIndices -> [ 0 2 4 ] [ 1 3 4 ] ....  vector<uint>
  
}


void MayaSceneTest::AtEachTimestep()
{
  if (GetBoolOpt("reverse")) {
    
    reverse_timestep();
    return;
  }
  
  current_step++;
  
  if (current_step >= nSubsteps) {
    current_frame ++;
    current_step = 0;
    LoadNextFrameFromCache();
  }
  
  std::cout << current_step << " step " << current_frame << " frame \n";
  
//  std::cout << "MAYA time step \n";
  
// Set boundary conditions for fixed DOFs
  for(int i=0; i<(int)m_rods.size(); i++) {
    int idbc = m_rods[i]->global_rodID;
    VecXd prev = prevBoundaryCondition[idbc];
    VecXd next = nextBoundaryCondition[idbc];
    for(int v=0; v<2; v++) {
      Vec3d prev3 = prev.segment<3> (3 * v);
      Vec3d next3 = next.segment<3> (3 * v);
      Vec3d target_pos = prev3 + (double) current_step / (double) nSubsteps * (next3 - prev3);
      
      //std::cout << m_rods[i]->getBoundaryCondition()->getDesiredVertexPosition(v);
      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(v, target_pos);
      //std::cout << target_pos;
    }
    //std::cout << "\n";
  }
  
  Scalar interpolateFactor = ( (double)(current_step) / nSubsteps );
  if ( i_collisionsEnabled )
  {
    for ( CollisionMeshDataHashMapIterator cmItr = m_collisionMeshMap.begin();
          cmItr != m_collisionMeshMap.end(); ++cmItr )
    {
      cmItr->second->interpolate( interpolateFactor );
      UpdateTriangleMesh(cmItr->second->currPositions);
    }
  }  
  
//  if ( !i_selfCollisionPenaltyForcesEnabled && !i_fullSelfCollisionsEnabled )
  {
            // Fantastic let's just run it all threaded!
//    timeval threadFrameTimer;
//    startTimer( threadFrameTimer );
    double dt = GetScalarOpt("dt");

#pragma omp parallel for num_threads( 8 )
    for ( int i=0; i<m_collision_controllers.size(); ++i )
    {                
//      RodCollisionTimeStepper* rodCollisionTimeStepper = dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ]);
      
      RodCollisionTimeStepper* rodCollisionTimeStepper = m_collision_controllers[i];

      rodCollisionTimeStepper->shouldDoCollisions( i_collisionsEnabled );
      rodCollisionTimeStepper->setTimeStep( dt );
      rodCollisionTimeStepper->setCollisionMeshesMap( &m_collisionMeshMap );
      rodCollisionTimeStepper->setClumping( false, 0.0 );
      
      m_controllers[i]->setTimeStep(dt);

      m_collision_controllers[i]->initialiseCollisionsAndApplyObjectCollisionForces();
    }
    
            // Clumping is like self collisions, has to all be done at the same time
//    if ( m_isClumpingEnabled )
//    {
//      RodCollisionTimeStepper::getClumpingPairs( m_rods );
//    }

#pragma omp parallel for num_threads( 8 )
    for ( int i=0; i<m_collision_controllers.size(); ++i )
    {
      int idbc = m_rods[i]->global_rodID;
      //std::cout << "Now simulating rod # " << i << "\n";
      
      RodCollisionTimeStepper* rodCollisionTimeStepper = m_collision_controllers[i];

      ElasticRod* elasticRod = rodCollisionTimeStepper->getRod();

      rodCollisionTimeStepper->execute();
      
      rodCollisionTimeStepper->respondToObjectCollisions();
      rodCollisionTimeStepper->tidyUpCollisionStructuresForNextStep();
      
//      VecXd maya_pos = prevMayaPosition[idbc] + (double) current_step / (double) nSubsteps * (nextMayaPosition[idbc] - prevMayaPosition[idbc]);
//      VecXd real_pos (elasticRod->nv() * 3);
      
//      for(int j=0; j<elasticRod->nv(); j++) {
//        real_pos.segment<3> (3 * j) = elasticRod->getVertex(j);
//      }
      
//      VecXd diff = maya_pos - real_pos;
//      double dif = diff.norm();
//      std::cout << idbc << " " << dif / (double)elasticRod->nv()  << "\n";

    }
            
//    frameTime += stopTimer( threadFrameTimer );
  }  
  
}

void MayaSceneTest::UpdateTriangleMesh( std::vector<Vec3d> &points ) {
  
  if (!br_mesh) {
    return;
  }
  
  size_t i=0;
  for( TriangleMesh::vertex_iter itr = br_mesh->vertices_begin(); itr != br_mesh->vertices_end(); ++itr, ++i )
  {
    br_mesh->getVertex(*itr) = points[i];
  }
}

// class TriangleMesh
void MayaSceneTest::MakeTriangleMesh( std::vector<Vec3d> &points, std::vector<uint> &faces ) {
  
  if (br_mesh) {
    delete br_mesh;
  }
  
  br_mesh = new TriangleMesh();

  for(size_t i=0; i<points.size(); i++) {
    TriangleMesh::vertex_handle vh = br_mesh->addVertex();
    br_mesh->getVertex(vh) = points[i];
  }
  
  for(size_t i=0; i<faces.size(); i+=3) {
  
    TriangleMesh::vertex_handle vhx(faces[i+0]);
    TriangleMesh::vertex_handle vhy(faces[i+1]);
    TriangleMesh::vertex_handle vhz(faces[i+2]);
    
    // Generate three edges for the face
    br_mesh->addEdge(vhx,vhy);
    br_mesh->addEdge(vhy,vhz);
    br_mesh->addEdge(vhz,vhx);
    
    // Generate a triangular face
    br_mesh->addFace(vhx,vhy,vhz);  
  }
  
  m_tri_objs.push_back(br_mesh);
  m_world->addObject(br_mesh);
  
}


void MayaSceneTest::reverse_test() {
  int nv = 3;
  std::vector <Vec3d> vertices;
  
  Vec3d v ( 0,0,0);
  
  for(int i=0; i<nv; i++) {
    vertices.push_back(v);
    v += Vec3d(1,1,0);
  }
  
  loadDynamicsProps();
  RodOptions opts;
  getRodOptions(opts);
  
  opts.numVertices = nv;

  ElasticRod* newrod = setupRod(opts, vertices, vertices);
  newrod->global_rodID = 0;
  
  newrod->getBoundaryCondition()->setDesiredVertexPosition(0, newrod->getVertex(0));
  newrod->getBoundaryCondition()->setDesiredVertexPosition(1, newrod->getVertex(1));
  newrod->getBoundaryCondition()->setDesiredEdgeAngle(0, newrod->getTheta(0));
      
//    if (newrod->nv() > 2) {
//      newrod->getBoundaryCondition()->setDesiredVertexPosition(2, newrod->getVertex(2));
//      newrod->getBoundaryCondition()->setDesiredEdgeAngle(1, newrod->getTheta(1));
//    }

  RodTimeStepper* stepper = new RodTimeStepper( *newrod );
  stepper->setDiffEqSolver( RodTimeStepper::SYM_IMPL_EULER );
  
    // These get deleted by RodTimeStepper
//  stepper->addExternalForce( new RodMassDamping( GetScalarOpt("mass-damping") ) );
    
//    if ( i_gravity.norm() > 0)
  {
    stepper->addExternalForce( new RodGravity( GetVecOpt("gravity") ) );
  }

//  AdaptiveBinaryStepper* adpstep = new AdaptiveBinaryStepper( newrod, stepper );
    //m_steppers.push_back(adpstep);
  
  RodCollisionTimeStepper *cstepper = new RodCollisionTimeStepper( stepper, newrod );
    
//    RodTimeStepper* newstepper = getRodTimeStepper(*newrod);
  m_rods.push_back(newrod);
  m_controllers.push_back(stepper);
  m_world->addObject(newrod);
  m_collision_controllers.push_back(cstepper);
  
  
}

void MayaSceneTest::reverse_timestep()
{
  static int fr = 0;
  
  fr++;
  if (fr == 2) {
    
    reverse_solve();
    return;
  }
  
  reverse_solve();

    double dt = GetScalarOpt("dt");

   //std::cout <<"wef\n";
  for ( int i=0; i<m_collision_controllers.size(); ++i )
  {                
  //      RodCollisionTimeStepper* rodCollisionTimeStepper = dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ]);
        
    RodCollisionTimeStepper* rodCollisionTimeStepper = m_collision_controllers[i];
  
    rodCollisionTimeStepper->shouldDoCollisions( false );
    rodCollisionTimeStepper->setTimeStep( dt );
    rodCollisionTimeStepper->setCollisionMeshesMap( &m_collisionMeshMap );
    rodCollisionTimeStepper->setClumping( false, 0.0 );
        
    m_controllers[i]->setTimeStep(dt);
  
    m_collision_controllers[i]->initialiseCollisionsAndApplyObjectCollisionForces();
  }
      
  //std::cout <<"wef2\n";
  
  for ( int i=0; i<m_collision_controllers.size(); ++i )
  {
    int idbc = m_rods[i]->global_rodID;
    std::cout << "Now simulating rod # " << i << "\n";
        
    RodCollisionTimeStepper* rodCollisionTimeStepper = m_collision_controllers[i];
  
    ElasticRod* elasticRod = rodCollisionTimeStepper->getRod();
  
    rodCollisionTimeStepper->execute();
        
    rodCollisionTimeStepper->respondToObjectCollisions();
    rodCollisionTimeStepper->tidyUpCollisionStructuresForNextStep();
        
  
  }
}

void MayaSceneTest::reverse_solve() 
{
  std::cout << "reverse solve\n";
  ElasticRod *rod = m_rods[0];
  
  std::vector<RodForce*> forces = rod->getForces();
  
  for(int i=0; i<forces.size(); i++) {
    std::cout << i << " " << typeid(forces[i]).name() << "\n";
    RodBendingForceSym *force = dynamic_cast<RodBendingForceSym *> (forces[i]);
    
    if (force) {
      
      if (force->viscous()) continue;
      
      force->computeGradKappa();
      
      ElasticRod::vertex_handle vh(1);
      
      Mat2d B = force->getB(vh);
      Scalar len = force->getRefVertexLength(vh);

      const Vec2d& kappa = force->getKappa(vh);
//      const Vec2d& kappaBar = getKappaBar(vh);
      VecXd f(11);
      MatXd m(11,2);
      
      m = -1.0/len * force->getGradKappa(vh) * B;

      f = -1.0/len * force->getGradKappa(vh) * B * (kappa);
      
      std::cout << m << "\n" << f <<"\n" << kappa << "\n";
      
      Scalar m00 = m(8,0);
      Scalar m01 = m(8,1);
      Scalar m10 = m(9,0);
      Scalar m11 = m(9,1);
      Scalar m20 = m(10,0);
      Scalar m21 = m(10,1);
      
      Scalar bx = m00 * kappa(0) + m01 * kappa(1);
      Scalar by = m10 * kappa(0) + m11 * kappa(1) - 981.0;
      Scalar bz = m20 * kappa(0) + m21 * kappa(1);
      
      Scalar k0;
      Scalar k1;
      
      std::cout << bx << " " << by << " " << bz << "\n";
      
      if (fabs(m00 * m11 - m01 * m10) < 1e-4) {
        k0 = (by * m21 - bz * m11) / (m10 * m21 - m11 * m20);
        k1 = (bz * m10 - by * m20) / (m10 * m21 - m11 * m20);
      } else {
        k0 = (bx * m11 - by * m01) / (m00 * m11 - m01 * m10);
        k1 = (by * m00 - bx * m10) / (m00 * m11 - m01 * m10);
      }
      
      Vec2d kappaBar(k0, k1);
      
      std::cout << kappaBar << "\n";
    }
    
  }
  
  
}




