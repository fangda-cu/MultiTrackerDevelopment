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



#ifdef COLUMBIA_COLLISION

MeshController::MeshController( TriangleMesh* mesh ) : m_mesh(mesh)
{
}

bool MeshController::execute() 
{
  if (m_mesh) {
    size_t i=0;
    for( TriangleMesh::vertex_iter itr = m_mesh->vertices_begin(); itr != m_mesh->vertices_end(); ++itr, ++i )
    {
      m_mesh->getVertex(*itr) = m_verts[i];
    }  
  }
  
  return true;
}

void MeshController::updateMesh(  std::vector<Vec3d> & prevPosition, std::vector<Vec3d> & nextPosition, double s )
{
  m_verts.clear();
  
  for(int i=0; i<(int)prevPosition.size(); i++) {
    Vec3d p = prevPosition[i] * (1.0 - s) + nextPosition[i] * (s);
    m_verts.push_back(p);
  }
  
}


#endif


MayaSceneTest::MayaSceneTest()
  : Problem("Maya Scene Test", "Maya Scene Test.")
, m_rods()
, m_controllers()
{
  
#ifdef COLUMBIA_COLLISION
  std::cout << "COLUMBIA_COLLISION enbled\n";
#endif
  
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();
  
  AddOption("read-camera", "read camera parameters from config file", false);
  AddOption("eye", "eye", Vec3d(0,0,0));
  AddOption("up", "up", Vec3d(0,1,0));
  AddOption("center", "view center", Vec3d(0,0,-1));
  
  
  //AddOption("rodfile-path", "rodfile-path", "/usr/home/jjoo/scenes/mau160.");
  AddOption("rodfile-path", "rodfile-path", "/local1/jjoo/scenes/sac2/sac.");
  AddOption("rodfile-firstindex", "rodfile-firstindex", -39);
  AddOption("rod-input-rodnum", "rod-input-rodnum", 417);
  AddOption("maya-substeps", "maya-substeps", 20);
  
  AddOption("reverse", "reverse", false);
  AddOption("reverse-internal", "reverse", false);

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
  
  AddOption("group-spring", "group-spring", false);
  AddOption("group-num", "group-num", 10);
  AddOption("group-layer-gap", "group-layer-gap", 1.0);
  AddOption("bulk-stiffness", "bulk-stiffness", 10.0);
  
  AddOption("multi-stepper", "multi-stepper", false);
  
//  GetScalarOpt("viscosity") = 10;
//  GetScalarOpt("mass-damping") = 10;
  
  br_mesh = NULL;
//  rodGroupManager = NULL;
  
  multi_stepper = NULL;

//  group_spring = NULL;
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

// ###################SETUP
void MayaSceneTest::Setup()
{
  const bool outputForMathe = false;
  
  nvs.clear();
  
  if (GetBoolOpt("reverse")) {
    if (0) {
      reverse_test();
      return;
    } else {
      // from list
//      reverse_setup_list();
    }
  }
  
  is_multi_stepper = GetBoolOpt("multi-stepper");

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
  setDt(1.0 / (double) fps / (double)nSubsteps);
  
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

  std::cout << name << " reading\n";
  FILE *fp = readNumberOfRodsFromFile(name, nr);
  
  if (fp == NULL) {
    std::cout << "failed to open\n";
  }
  
  std::ofstream out_file;
  if (outputForMathe) {
    out_file.open("sakr_mathematica.txt");
    out_file.precision(15);
    
    out_file << max_r << " ";
  }
    
  
  
  std::cout << nr << " rods\n";
  
      
  for ( size_t r=0; r<nr && r<(size_t)max_r; r++ )
  {
    size_t numVertices;
    size_t ret = fread( &numVertices, sizeof( size_t ), 1, fp );

//    if (GetBoolOpt("reverse-internal")) {
//      numVertices = 3;
//    }
    
    nvs.push_back(numVertices);
    
    std::vector<Vec3d> vertices;
    VecXd fixedVerts(6);
    VecXd allVerts(numVertices * 3);

    if (outputForMathe) {
      out_file << numVertices << " ";
    }
    
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
      
      out_file << pos[0] << " ";
      out_file << pos[1] << " ";
      out_file << pos[2] << " ";
      
      
//      std::cout << " { "  << pos[0] << "," <<pos[1]<<" ," <<pos[2] <<" }, ";
      
//      o_rodVertices[ r ][ v ] = Vec3d( pos[0], pos[1], pos[2] );
      
      ret = fread( &pos[0], sizeof( double ), 3, fp );
      if (v < 2) {
        fixedVerts(v*3+0) = pos[0];
        fixedVerts(v*3+1) = pos[1];
        fixedVerts(v*3+2) = pos[2];
      }
//      o_unsimulatedRodVertices[ r ][ v ] = Vec3d( pos[0], pos[1], pos[2] );
    }
    
//    std::cout <<"\n";
    
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
      
//    for(int i=1; i<(numVertices-1); i++) {
//        newrod->getBoundaryCondition()->setDesiredEdgeAngle(i, newrod->getTheta(i));
//    }
    
    
#ifdef COLUMBIA_COLLISION
#else 
    newrod->fixVert(0);
    newrod->fixVert(1);
    newrod->fixEdge(0);
#endif

    if (is_multi_stepper) {
      RodCollisionTimeStepper *cstepper = new RodCollisionTimeStepper( NULL, newrod ); // multi_stepper will handle dynamic
      cstepper->impulse_enabled = impulse_enabled;
      m_collision_controllers.push_back(cstepper);
      
    } else {
      RodTimeStepper* stepper = new RodTimeStepper( *newrod );
      stepper->setDiffEqSolver( RodTimeStepper::SYM_IMPL_EULER );
      
//      stepper->setDiffEqSolver( RodTimeStepper::SYMPL_EULER );
    
      stepper->setMaxIterations(GetIntOpt("iterations"));
      stepper->set_stol(GetScalarOpt("stol"));
      stepper->set_atol(GetScalarOpt("atol"));
      stepper->set_rtol(GetScalarOpt("rtol"));
      stepper->set_inftol(GetScalarOpt("inftol"));
  
      // These get deleted by RodTimeStepper
      stepper->addExternalForce( new RodMassDamping( GetScalarOpt("mass-damping") ) );
      stepper->addExternalForce( new RodGravity( GetVecOpt("gravity") ) );
  
      AdaptiveBinaryStepper* adpstep = new AdaptiveBinaryStepper( newrod, stepper );
      adpstep->setTimeStep(GetScalarOpt("dt"));
    
#ifdef COLUMBIA_COLLISION
#else 
      RodCollisionTimeStepper *cstepper = new RodCollisionTimeStepper( adpstep, newrod );
      
      cstepper->impulse_enabled = impulse_enabled;
      m_collision_controllers.push_back(cstepper);
#endif
    
      if (GetBoolOpt("reverse-internal")) {
        newrod->doReverseHairdo(stepper);
      }
      
      m_controllers.push_back(adpstep);
    }
          
    
    
    m_rods.push_back(newrod);
    m_world->addObject(newrod);
    
  }
  fclose(fp);
  
  if (outputForMathe) {
    std::cout << " rods were recorded into a file for mathematica\n";
    out_file.close();
  }
  
  if (GetBoolOpt("reverse")) { // from list
    reverse_findudf_list();
  }
  
  
/*  if (GetBoolOpt("group-spring") ) {
//    rodGroupManager = new RodGroupManager (m_rods, GetIntOpt("group-num"), GetScalarOpt("group-layer-gap"), GetScalarOpt("bulk-stiffness"));
    
    group_spring = new RodGroupSpringForce(m_rods);
  } else {
    group_spring = NULL;
  }*/
  
  if (is_multi_stepper) {
    makeMultipleRodStepper();
    
    for( int i = 0; i < (int) m_rods.size(); ++i ) multi_stepper->addRod(m_rods[i]);
    
//    if (group_spring)
//      multi_stepper->addRodRodExternalForce(group_spring);
    //m_world->addController(multi_stepper);
    
  }  
  
  InitializeMeshes();
  
#ifdef COLUMBIA_COLLISION
  std::cout << "GetScalarOpt dt " << GetScalarOpt("dt") << "\n";
  
  m_br_stepper = new BridsonStepper( m_rods, m_tri_objs, m_obj_controllers, m_controllers, GetScalarOpt("dt") );
  m_br_stepper->setEdgeEdgePenalty(1.0);
  m_br_stepper->setVertexFacePenalty(1.0);
  m_br_stepper->disablePenaltyImpulses();
  m_br_stepper->enableImplicitPenaltyImpulses();
  //m_br_stepper->setRodLabels(rod_labels);
  m_world->addController(m_br_stepper);
#endif  
  
  LoadNextFrameFromCache();
  
}

void MayaSceneTest::makeMultipleRodStepper()
{
  multi_stepper = NULL;
  
  MultipleRodTimeStepper* stepper = new MultipleRodTimeStepper();

  string& integrator = GetStringOpt("integrator");

  if (integrator == "symplectic") {
    stepper->setDiffEqSolver(MultipleRodTimeStepper::SYMPL_EULER);
  } 
  if (integrator == "implicit")
  {
    stepper->setDiffEqSolver(MultipleRodTimeStepper::SYM_IMPL_EULER);
  } 
  else 
  {
    cerr << "Unknown integrator " << integrator
        << ". Using default instead." << endl;
  }
  
  stepper->setTimeStep(getDt());
  
  Scalar massDamping = GetScalarOpt("mass-damping");
  if (massDamping != 0) {
    stepper->addExternalForce(new RodMassDamping(massDamping));
  }
  
  if (getGravity().norm() > 0) {
    stepper->addExternalForce(new RodGravity(getGravity()));
  }
  
  stepper->setMaxIterations(GetIntOpt("iterations"));
  stepper->set_stol(GetScalarOpt("stol"));
  stepper->set_atol(GetScalarOpt("atol"));
  stepper->set_rtol(GetScalarOpt("rtol"));
  stepper->set_inftol(GetScalarOpt("inftol"));
  
  multi_stepper = stepper;
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


#ifdef COLUMBIA_COLLISION
#else COLUMBIA_COLLISION
void MayaSceneTest::initialiseCollisionMesh( BASim::CollisionMeshData *collisionMeshData, size_t id )
{
  m_collisionMeshMap[ id ] = collisionMeshData;
  m_collisionMeshMap[ id ]->initialize();
}
#endif COLUMBIA_COLLISION

void MayaSceneTest::UpdateCollisionMesh( std::vector<Vec3d> &points, std::vector<uint> &faces )
{    
#ifdef COLUMBIA_COLLISION
  
  prevMeshPosition = nextMeshPosition;
  nextMeshPosition = points;

#else COLUMBIA_COLLISION
    
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

#endif COLUMBIA_COLLISION
  
}


void MayaSceneTest::InitializeMeshes() {
  
  if (!GetBoolOpt("load-meshes")) return;
  
  mesh_frame = 1;
  
  size_t nv = 0;
  size_t nf = 0;
  std::vector<Vec3d> verts;
  std::vector<uint> faces;
  
  ReadMeshFile(mesh_frame, nv, nf, verts, faces);
  MakeTriangleMesh(verts, faces);
  
#ifdef COLUMBIA_COLLISION
  nextMeshPosition = verts;

  for(int i=0; i<(int)m_tri_objs.size(); i++) {
    MeshController *new_controller = new MeshController(m_tri_objs[i]);
    m_mesh_controllers.push_back(new_controller);
    m_obj_controllers.push_back(new_controller);
  }
#else COLUMBIA_COLLISION
  
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
  

  
#endif COLUMBIA_COLLISION  

  mesh_frame++;
//  ReadMeshFile(mesh_frame, nv, nf, verts, faces);
//  UpdateCollisionMesh(verts, faces);
  
  //  void setCollisionMeshesMap(CollisionMeshDataHashMap* collisionMeshes)
//m_collisionMeshMap
  //triangleIndices -> [ 0 2 4 ] [ 1 3 4 ] ....  vector<uint>
  
}


// ####################TIMESTEP
void MayaSceneTest::AtEachTimestep()
{ 
  
  if (0) {
    std::cout.precision(15);
  
    std::cout << "getVertex\n";
    for(int i=0; i<m_rods.size(); i++) {
      std::cout <<  "{ ";
      for(int j=0; j<m_rods[i]->nv(); j++) {
        std::cout << "{" << m_rods[i]->getVertex(j)[0] << ", " << m_rods[i]->getVertex(j)[1] << ", " << m_rods[i]->getVertex(j)[2]  << "}, ";
      }
      std::cout << "};\n";
    }
      
  }
    

  if (0) 
  {
    static int fr = 0;
    fr++;
    
    if (fr > 200) {
      for ( int i=0; i<m_collision_controllers.size(); ++i )
      {                
        RodCollisionTimeStepper* rodCollisionTimeStepper = m_collision_controllers[i];
        
        rodCollisionTimeStepper->clearVertexPositionPenalty(-1);
      }
    } else {
      for ( int i=0; i<m_collision_controllers.size(); ++i )
      {                
        RodCollisionTimeStepper* rodCollisionTimeStepper = m_collision_controllers[i];
        
        Vec3d xx =  m_rods[i]->getVertex(m_rods[i]->nv()-1) + Vec3d(0.1, 0.1, 0.0);
        rodCollisionTimeStepper->setVertexPositionPenalty(m_rods[i]->nv()-1,xx , 100);
        
        rodCollisionTimeStepper->setVertexPositionPenalty(3,xx , 1);
        rodCollisionTimeStepper->setVertexPositionPenalty(4,xx , 2);
        rodCollisionTimeStepper->setVertexPositionPenalty(5,xx , 3);
        rodCollisionTimeStepper->setVertexPositionPenalty(3,xx , 4);
        rodCollisionTimeStepper->setVertexPositionPenalty(5,xx , 5);
        rodCollisionTimeStepper->setVertexPositionPenalty(3,xx , 6);
        
        rodCollisionTimeStepper->clearVertexPositionPenalty(5);
        rodCollisionTimeStepper->setVertexPositionPenalty(3,xx , 7);
        rodCollisionTimeStepper->clearVertexPositionPenalty(3);
        rodCollisionTimeStepper->clearVertexPositionPenalty(-1);
        
        
      }
    }
    
  }
  
  
  if (0)
  {
    ElasticRod *rod = m_rods[0];
    std::vector<RodForce*>& forces = rod->getForces();
    std::vector<RodForce*>::iterator fIt;
    for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
      if (dynamic_cast<RodBendingForceSym*> (*fIt) != NULL ) {
        RodBendingForceSym *f = dynamic_cast<RodBendingForceSym*> (*fIt);
        
        if (f->viscous()) {
          std::cout << "viscous bending stf\n";
        
          ElasticRod::vertex_iter vit = rod->vertices_begin();
          uint i=0;
          
          ++vit;
    
          for (; vit != rod->vertices_end() && i < rod->nv() - 1; ++vit, ++i) {
            ElasticRod::vertex_handle vh = *vit;
            std::cout << f->getKappaBar(vh) << "\n";
          }
        }
      }
    }
  }
  
  if (0 && GetBoolOpt("reverse")) {
    reverse_timestep();
    return;
  }
  
  /*
  for(int i=0; i<m_rods.size(); i++) {
      
    std::cout << "\ntheta\n";
    for(int j=0; j<m_rods[i]->ne(); j++) {
      std::cout << m_rods[i]->getTheta(j) << " ";
    }
    std::cout << "\n";
    
    std::cout.precision(15);

    std::cout << "getVertex\n";
    std::cout <<  "{ ";
    for(int j=0; j<m_rods[i]->nv(); j++) {
      std::cout << "{" << m_rods[i]->getVertex(j)[0] << ", " << m_rods[i]->getVertex(j)[1] << ", " << m_rods[i]->getVertex(j)[2]  << "}, ";
    }
    std::cout << "};\n";
    
    std::cout << "\ngetEdge\n";
    for(int j=0; j<m_rods[i]->ne(); j++) {
      std::cout << m_rods[i]->getEdge(j) << ", ";
    }
    std::cout << "\n";
    
    
    
    std::cout << "twist\n";
    std::vector<RodForce*>& forces = m_rods[i]->getForces();
    std::vector<RodForce*>::iterator fIt;
    for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
      if (dynamic_cast<RodTwistingForceSym*> (*fIt) != NULL ) {
        RodTwistingForceSym *f = dynamic_cast<RodTwistingForceSym*> (*fIt);
        
        if (f->viscous()) continue;
        ElasticRod::vertex_iter vit = m_rods[i]->vertices_begin();
        uint j=0;
        ++vit;

        for (; vit != m_rods[i]->vertices_end() && j < m_rods[i]->nv() - 1; ++vit, ++j) {
          ElasticRod::vertex_handle vh = *vit;
          std::cout << f->getTwist(vh) << " ";
        }
        std::cout << "\n";
      
      }
    }
    
    std::cout << " reference\n";
    for(int j=0; j<m_rods[i]->ne(); j++) {
      std::cout <<  m_rods[i]->getTangent(j) << " " ;
      std::cout <<  m_rods[i]->getReferenceDirector1(j) << " " ;
      std::cout <<  m_rods[i]->getReferenceDirector2(j) << "\n " ;
    }
    std::cout << "\n";
    
    std::cout << " material\n";
    for(int j=0; j<m_rods[i]->ne(); j++) {
      std::cout <<  m_rods[i]->getTangent(j) << " " ;
      std::cout <<  m_rods[i]->getMaterial1(j) << " " ;
      std::cout <<  m_rods[i]->getMaterial2(j) << "\n" ;
    }
    std::cout << "\n";
    
  }
  */
  
  
/*  for(int i=0; i<m_rods.size(); i++) {
    bool print = false;
    for(int j=0; j<m_rods[i]->nv(); j++) {
      Vec3d vvv = m_rods[i]->getVelocity(j);
      if (vvv.norm() > 2.0) print = true;
    }
    
    if (print) {
      
      std::cout << i << "  {";
      for(int j=0; j<m_rods[i]->nv(); j++) {
        std::cout << "{" << m_rods[i]->getVertex(j)[0] << ", " << m_rods[i]->getVertex(j)[1] << ", " << m_rods[i]->getVertex(j)[2]  << "}, ";
      }
      std::cout << "};\n";
    }
  }*/
  
  current_step++;
  
  if (current_step >= nSubsteps) {
    current_frame ++;
    current_step = 0;
    LoadNextFrameFromCache();
  }
  
  std::cout << current_step << " step " << current_frame << " frame \n";
  
//  if (group_spring) {
//    group_spring->checkActivatingCondition();
//  }
  
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
    }
  }
  
  Scalar interpolateFactor = ( (double)(current_step) / nSubsteps );
  if ( 1 || i_collisionsEnabled )
  {
#ifdef COLUMBIA_COLLISION
    for(int i=0; i<(int)m_mesh_controllers.size(); i++) {
      m_mesh_controllers[i]->updateMesh( prevMeshPosition, nextMeshPosition, interpolateFactor);
    }
    
#else 
    for ( CollisionMeshDataHashMapIterator cmItr = m_collisionMeshMap.begin();
          cmItr != m_collisionMeshMap.end(); ++cmItr )
    {
      cmItr->second->interpolate( interpolateFactor );
      UpdateTriangleMesh(cmItr->second->currPositions);
    }
#endif
  }  
  
//  if ( !i_selfCollisionPenaltyForcesEnabled && !i_fullSelfCollisionsEnabled )
  {
    double dt = GetScalarOpt("dt");
    
//    if (rodGroupManager) {
//      rodGroupManager->ActivateSpring(dt);
//    }
    
    if (is_multi_stepper) {
      multi_stepper->setTimeStep(dt);
      
      for ( int i=0; i<m_collision_controllers.size(); ++i )
      {                
        RodCollisionTimeStepper* rodCollisionTimeStepper = m_collision_controllers[i];
  
        rodCollisionTimeStepper->shouldDoCollisions( i_collisionsEnabled );
        rodCollisionTimeStepper->setTimeStep( dt );
        rodCollisionTimeStepper->setCollisionMeshesMap( &m_collisionMeshMap );
        rodCollisionTimeStepper->setClumping( false, 0.0 );
        
        m_collision_controllers[i]->initialiseCollisionsAndApplyObjectCollisionForces();
      }
      
      multi_stepper->execute();
      
      for ( int i=0; i<m_collision_controllers.size(); ++i )
      {                
        RodCollisionTimeStepper* rodCollisionTimeStepper = m_collision_controllers[i];
        rodCollisionTimeStepper->respondToObjectCollisions();
        rodCollisionTimeStepper->tidyUpCollisionStructuresForNextStep();
      }
      
    } else {

    
#ifdef COLUMBIA_COLLISION
#else 

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

#pragma omp parallel for num_threads( 8 )
      for ( int i=0; i<m_collision_controllers.size(); ++i )
      {
        int idbc = m_rods[i]->global_rodID;
        //std::cout << "Now simulating rod # " << i << "\n";
        
        RodCollisionTimeStepper* rodCollisionTimeStepper = m_collision_controllers[i];
  
//        ElasticRod* elasticRod = rodCollisionTimeStepper->getRod();
  
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
#endif

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

void reverse_setup_list() {
  

}

void MayaSceneTest::reverse_test() {
  //int nv = 3;
  std::vector <Vec3d> vertices;
  
  Vec3d v ( 0,0,0);
  
  
  double v0[][3] = {{-331.055, 185.907, 1284.26}, {-331.103, 186.868, 1283.63}, {-331.469, 187.419, 1282.69}, {-331.703, 187.485, 1281.59}, {-331.991, 187.344, 1280.5}, {-332.339, 186.879, 1279.52}, {-332.463, 186.236, 1278.6}, {-332.533, 185.692, 1277.61}, {-332.797, 185.24, 1276.61}, {-333.545, 184.861, 1275.85}, {-334.015, 184.195, 1275.06}, {-334.384, 183.422, 1274.32}, {-334.552, 182.571, 1273.58}, {-334.606, 181.678, 1272.99}, {-334.301, 180.719, 1272.44} };
  double v1[][3] = {{-330.043, 183.617, 1294.71}, {-329.649, 184.663, 1295.07}, {-329.452, 185.794, 1295.28}, {-329.656, 186.937, 1295.23}, {-330.208, 187.796, 1294.68}, {-330.858, 188.244, 1293.82}, {-331.658, 188.503, 1293.01}, {-332.581, 188.896, 1292.41}, {-333.37, 189.385, 1291.69}, {-334.158, 189.809, 1290.93}};
  double v2[][3] = {{-346.48, 181.787, 1288.76}, {-347.081, 181.466, 1287.89}, {-348.081, 180.948, 1288.13}, {-348.894, 180.282, 1288.65}, {-349.595, 179.416, 1288.99}, {-350.331, 178.566, 1288.68}, {-351.013, 177.735, 1288.21}, {-351.561, 176.815, 1287.74}, {-352.083, 175.892, 1287.23}, {-352.666, 175.021, 1286.73}, {-353.535, 174.551, 1286.08}, {-354.519, 174.009, 1285.75}, };
  double v3[][3] = {{-329.725, 180.429, 1281.8}, {-329.148, 179.542, 1281.36}, {-328.882, 178.581, 1280.75}, {-328.916, 177.792, 1279.88}, {-329.097, 177.066, 1278.97}, {-329.387, 176.546, 1277.96}, {-329.785, 176.355, 1276.87}, {-330.224, 176.543, 1275.81}, {-331, 176.92, 1275.01}, {-331.822, 177.31, 1274.26}, };
  double v4[][3] = {{-335.405, 185.113, 1280.32}, {-335.853, 185.129, 1279.25}, {-336.293, 184.813, 1278.23}, {-336.122, 184.252, 1277.25}, {-335.704, 183.589, 1276.4}, {-335.372, 182.814, 1275.61}, {-335.726, 181.921, 1275}, {-336.48, 181.153, 1274.57}, {-337.265, 180.359, 1274.26}, {-338.077, 179.613, 1273.9}, {-338.849, 178.972, 1273.33}, {-339.331, 178.757, 1272.32}, };
  double v5[][3] = {{-331.986, 185.928, 1292.86}, {-332.324, 186.73, 1291.97}, {-332.691, 187.679, 1291.25}, {-333.207, 188.736, 1290.86}, {-334.006, 189.648, 1290.59}, {-334.896, 190.211, 1289.94}, {-335.67, 190.583, 1289.03}, {-336.571, 190.62, 1288.2}, {-337.614, 190.129, 1287.73}, };
  double v6[][3] = {{-338.7, 172.448, 1278.32}, {-338.509, 171.528, 1277.64}, {-338.333, 170.763, 1276.79}, {-338.356, 170.064, 1275.82}, {-338.535, 169.313, 1274.91}, {-338.663, 168.387, 1274.18}, {-338.923, 167.299, 1273.77}, {-339.486, 166.283, 1273.52}, };
  double v7[][3] = {{-330.916, 185.982, 1284.69}, {-330.956, 186.11, 1283.51}, {-330.834, 186.411, 1282.38}, {-330.829, 186.256, 1281.2}, {-331.101, 185.96, 1280.09}, {-331.512, 185.472, 1279.08}, {-331.519, 184.781, 1278.12}, {-331.062, 184.177, 1277.21}, {-331.273, 183.807, 1276.12}, {-331.996, 183.595, 1275.19}, };
  double v8[][3] = {{-338.969, 186.931, 1292.91}, {-340.101, 187.212, 1292.95}, {-341.167, 187.71, 1292.83}, {-342.255, 188.16, 1292.75}, {-343.416, 188.341, 1292.64}, {-344.52, 188.134, 1292.31}, {-345.306, 187.432, 1291.79}, {-345.993, 186.673, 1291.21}, {-346.827, 186.035, 1290.67}, {-347.854, 185.566, 1290.34}, {-348.978, 185.268, 1290.15}, };
  double v9[][3] = {{-327.945, 183.337, 1291.47}, {-327.711, 184.371, 1291.04}, {-327.731, 185.352, 1290.44}, {-327.697, 186.283, 1289.77}, {-327.659, 187.123, 1288.99}, {-327.782, 187.884, 1288.15}, {-328.13, 188.498, 1287.24}, {-328.721, 188.938, 1286.37}, {-329.551, 189.219, 1285.63}, {-330.343, 189.27, 1284.81}, {-330.891, 189.256, 1283.8}, {-330.999, 189.043, 1282.69}, };

  int nvs[] = {15, 10, 12, 10, 12, 9, 8, 10, 11, 12};
  
  //std::vector < double** > rvecs;
  //rvecs.push_back(v0);  rvecs.push_back(v1);  rvecs.push_back(v2);  rvecs.push_back(v3);  rvecs.push_back(v4);  rvecs.push_back(v5);  rvecs.push_back(v6);  rvecs.push_back(v7);  rvecs.push_back(v8);  rvecs.push_back(v9);
    
  
#define REV_ORG v7
#define REV_UDF w7
  
  int rid = 7;
  
  int nv = nvs[rid];
      
  for(int i=0; i<nv; i++) {
//    v = Vec3d(rvecs[rid][i][0], rvecs[rid][1], rvecs[rid][2]);
    v = Vec3d(REV_ORG[i][0], REV_ORG[i][1], REV_ORG[i][2]);
    vertices.push_back(v);
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
      
  newrod->fixVert(0);
  newrod->fixVert(1);
  newrod->fixEdge(0);
  
//    if (newrod->nv() > 2) {
//      newrod->getBoundaryCondition()->setDesiredVertexPosition(2, newrod->getVertex(2));
//      newrod->getBoundaryCondition()->setDesiredEdgeAngle(1, newrod->getTheta(1));
//    }

//  RodTimeStepper* stepper = new RodTimeStepper( *newrod );
//  stepper->setDiffEqSolver( RodTimeStepper::SYM_IMPL_EULER );
    
//    if ( i_gravity.norm() > 0)
  {
//    stepper->addExternalForce( new RodGravity( GetVecOpt("gravity") ) );
  }

//  RodCollisionTimeStepper *cstepper = new RodCollisionTimeStepper( stepper, newrod );
    
//    RodTimeStepper* newstepper = getRodTimeStepper(*newrod);
  m_rods.push_back(newrod);
//  m_controllers.push_back(stepper);
  m_world->addObject(newrod);
  //m_collision_controllers.push_back(cstepper);
  
  {  
    ElasticRod* newrod = setupRod(opts, vertices, vertices);
    newrod->global_rodID = 1;
    
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
      
  //    if ( i_gravity.norm() > 0)
    {
      stepper->addExternalForce( new RodGravity( GetVecOpt("gravity") ) );
    }
  
    RodCollisionTimeStepper *cstepper = new RodCollisionTimeStepper( stepper, newrod );
      
  //    RodTimeStepper* newstepper = getRodTimeStepper(*newrod);
    m_rods.push_back(newrod);
    m_controllers.push_back(stepper);
    m_world->addObject(newrod);
    m_collision_controllers.push_back(cstepper);
  }
  
  
  if(1)
  {
    double w6[][3] = {{-338.7, 172.448, 1278.32}, {-338.509, 171.528, 1277.64}, {-338.338, 
      170.867, 1276.71}, {-338.372, 170.376, 1275.62}, {-338.564, 169.879,
      1274.55}, {-338.721, 169.188, 1273.6}, {-339.026, 168.269, 
      1272.9}, {-339.63, 167.404, 1272.36}};
    
    //double w7[][3] =   {{-330.916, 185.982, 1284.69}, {-330.956, 186.11, 1283.51}, {-330.847,        186.69, 1282.49}, {-330.859, 187.059, 1281.36}, {-331.119, 187.445,        1280.28}, {-331.515, 187.75, 1279.19}, {-331.632, 187.875,         1278.02}, {-331.417, 188.09, 1276.88}, {-331.709, 188.604,         1275.87}, {-332.308, 189.287, 1275.09}};
        
    double w7[][3] =      {{-330.916, 185.982, 1284.69}, {-330.956, 186.11, 1283.51}, {-330.912,      187.126, 1282.92}, {-330.921, 188.265, 1282.57}, {-330.797,       189.439, 1282.59},{-330.538, 190.599, 1282.71}, {-330.556, 191.768,      1282.89}, {-330.949, 192.852, 1283.16}, {-330.831, 193.816,       1283.81}, {-330.261, 194.678, 1284.42}};
      
    
//      double w7[][3]= {{-330.916, 185.982, 1284.69}, {-330.956, 186.11, 1283.51}, {-330.844,        186.642, 1282.47}, {-330.853, 186.918, 1281.31}, {-331.125,         187.184, 1280.19}, {-331.541,187.344, 1279.08}, {-331.644, 187.312,        1277.91}, {-331.371, 187.372, 1276.76}, {-331.661, 187.733,         1275.68}, {-332.316, 188.269, 1274.83}};
        
      std::vector <Vec3d> vertices_udf;
//      vertices.clear();
      
      for(int i=0; i<nv; i++) {
        v = Vec3d(REV_UDF[i][0], REV_UDF[i][1], REV_UDF[i][2]);
        if (0) {
/*          v = Vec3d(q[i][0], q[i][1], q[i][2]);
          
          if (i > 1) {
            Vec3d v0 = Vec3d(q[i-1][0], q[i-1][1], q[i-1][2]);
            double d = (v-v0).norm();
            Vec3d w0 = vertices[i-1];
            Vec3d e0 = vertices[1] - vertices[0];
            
            e0.normalize();
            
            v = w0 + d * e0;
            
        }*/
        }
        vertices_udf.push_back(v);
      }
      
//      vertices = vertices_udf;
  
//      loadDynamicsProps();
      RodOptions opts;
      getRodOptions(opts);
  
      opts.numVertices = nv;
      opts.radiusScale = 20.0;

      ElasticRod* newrod = setupRod(opts, vertices, vertices);
      newrod->global_rodID = 2;
  
      newrod->getBoundaryCondition()->setDesiredVertexPosition(0, newrod->getVertex(0));
      newrod->getBoundaryCondition()->setDesiredVertexPosition(1, newrod->getVertex(1));
      newrod->getBoundaryCondition()->setDesiredEdgeAngle(0, newrod->getTheta(0));
      
      for(int i=1; i<nv-1; i++) {
//        newrod->getBoundaryCondition()->setDesiredEdgeAngle(i, newrod->getTheta(i));
        
      }
      
      newrod->fixVert(0);
      newrod->fixVert(1);
      newrod->fixEdge(0);

      RodTimeStepper* stepper = new RodTimeStepper( *newrod );
      stepper->setDiffEqSolver( RodTimeStepper::SYM_IMPL_EULER );
  
    // These get deleted by RodTimeStepper
  stepper->addExternalForce( new RodMassDamping( GetScalarOpt("mass-damping") ) );
    
//    if ( i_gravity.norm() > 0)
      {
        stepper->addExternalForce( new RodGravity( GetVecOpt("gravity") ) );
      }

//  AdaptiveBinaryStepper* adpstep = new AdaptiveBinaryStepper( newrod, stepper );
    //m_steppers.push_back(adpstep);
  

//      for(int i=2; i<nv; i++) {
//        newrod->setVertex(i, vertices_udf[i]);
//      }
//      newrod->updateProperties();
            
      RodCollisionTimeStepper *cstepper = new RodCollisionTimeStepper( stepper, newrod );
    /*
      double mm[] = {0.000867657, 0.000994337, 0.000233227, -0.000367572, -0.000286615, 0.000114513, 0.000163147, 0.0000268727};
      std::vector<RodForce*>& forces = newrod->getForces();
      std::vector<RodForce*>::iterator fIt;
      for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
        if (dynamic_cast<RodTwistingForceSym*> (*fIt) != NULL ) {
          RodTwistingForceSym *f = dynamic_cast<RodTwistingForceSym*> (*fIt);
      
          ElasticRod::vertex_iter vit = newrod->vertices_begin();
          uint i=0;
      
          ++vit;

          for (; vit != newrod->vertices_end() && i < newrod->nv() - 1; ++vit, ++i) {
            ElasticRod::vertex_handle vh = *vit;
            f->setUndeformedTwist(vh, mm[i]);
//        f->setRefVertexLength(vh, 0.5*(len[i]+len[i+1]));
            f->updateProperties();
          }
        }
      }
      
      newrod->updateProperties();*/
                    
//    RodTimeStepper* newstepper = getRodTimeStepper(*newrod);
      m_rods.push_back(newrod);
      m_controllers.push_back(stepper);
      m_world->addObject(newrod);
      m_collision_controllers.push_back(cstepper);    
    
      reverse_findudf(2);
    
  }
}


// #####UDFLIST
void MayaSceneTest::reverse_findudf_list() {
  int nr =0;
      
  std::ifstream myfile;
  myfile.open("att.txt");
  
  myfile>>nr;
  
  if (is_selected_running) {
    for (int i=0; i<selected_rods[0]; i++) {
      int nv = 0;
      myfile >> nv;
      std::cout << nv << " ";
      for(int j=0; j< (nv-2)*4; j++) {
        double k1;
        myfile >> k1;
      }
    }
  }
  
  for(int i=0; i<m_rods.size(); i++) {
    
    ElasticRod* rod = m_rods[i];
    
    int nv = 0;
    myfile >> nv;
      
    if (nv < 3) continue;
    
    std::vector <double> ldata;
    std::vector <double> kdata;
    std::vector <double> mdata;
    
    double k;
    for(int j=0; j<nv-2; j++) {
      myfile >> k;  ldata.push_back(k);
    }
    
    for(int j=0; j<nv*2-4; j++) {
      myfile >> k;  kdata.push_back(k);
    }
    
    for(int j=0; j<nv-2; j++) {
      myfile >> k;  mdata.push_back(k);
    }
    
    std::vector<RodForce*>& forces = rod->getForces();
    std::vector<RodForce*>::iterator fIt;
    for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
      if (dynamic_cast<RodStretchingForce*> (*fIt) != NULL ) {
        (*fIt)->updateUndeformedConfiguration(ldata);
      }
      if (dynamic_cast<RodBendingForceSym*> (*fIt) != NULL ) {
        (*fIt)->updateUndeformedConfiguration(kdata);
      }
      if (dynamic_cast<RodTwistingForceSym*> (*fIt) != NULL ) {
        (*fIt)->updateUndeformedConfiguration(mdata);
      }
    
    }
    
    rod->updateProperties();
    
  }
  myfile.close();
  
}


void MayaSceneTest::reverse_findudf(uint id) {
  ElasticRod* rod = m_rods[id];
  double kappa[][2] = 
  //{{-0.127306, 1.15698}, {0.102148, 0.341424}, {0.235131, 0.413327},  {0.107683, 0.187215}, {-0.372765, 0.0319196}, {-0.434292, 0.136436},   {0.558839, 0.392574}, {0.46502, 0.159427}};
  //{{-0.135107, 0.403972}, {0.0977328, -0.212217}, {0.23592, 0.00347672},   {0.119101, -0.0878528}, {-0.345602, -0.139429}, {-0.415669,       0.0364725}, {0.562019, 0.345449}, {0.464473, 0.146755}};
            
  {{-0.119968, 1.01617}, {0.164064, 0.232882}, {0.308329, 0.333767}, {0.157027, 0.123967}, {-0.345193, -0.0067079}, {-0.418527, 0.109367}, {0.563588, 0.383028}, {0.46528, 0.158578}};
  
    double theta[] = { 0.000867657, 0.000994337, 0.000233227, -0.000367572, -0.000286615, 
      0.000114513, 0.000163147, 0.0000268727 };
    
      
//    double len[] = {1.1501, 1.14942, 1.12655, 1.13619, 1.13918, 1.12926, 1.13178,      1.12872, 1.1317, 1.13515, 1.13194, 1.14019, 1.07166, 1.14682};
      
  std::vector<RodForce*>& forces = rod->getForces();
  std::vector<RodForce*>::iterator fIt;
  for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
    if (dynamic_cast<RodBendingForceSym*> (*fIt) != NULL ) {
      RodBendingForceSym *f = dynamic_cast<RodBendingForceSym*> (*fIt);
      
      ElasticRod::vertex_iter vit = rod->vertices_begin();
      uint i=0;
      
      ++vit;

      for (; vit != rod->vertices_end() && i < rod->nv() - 1; ++vit, ++i) {
        Vec2d kb (kappa[i][0], kappa[i][1]);
        ElasticRod::vertex_handle vh = *vit;
        f->setKappaBar(vh, kb);
//        f->setRefVertexLength(vh, 0.5*(len[i]+len[i+1]));
        f->updateProperties();
        
      }
      
    }
    if (dynamic_cast<RodTwistingForceSym*> (*fIt) != NULL ) {
      RodTwistingForceSym *f = dynamic_cast<RodTwistingForceSym*> (*fIt);
      
      ElasticRod::vertex_iter vit = rod->vertices_begin();
      uint i=0;
      
      ++vit;

      for (; vit != rod->vertices_end() && i < rod->nv() - 1; ++vit, ++i) {
        ElasticRod::vertex_handle vh = *vit;
//        f->setUndeformedTwist(vh, theta[i]);
//        f->setRefVertexLength(vh, 0.5*(len[i]+len[i+1]));
        f->updateProperties();
      }
      
    }
    if (dynamic_cast<RodStretchingForce*> (*fIt) != NULL ) {
      RodStretchingForce *f = dynamic_cast<RodStretchingForce*> (*fIt);
      
      ElasticRod::edge_iter eit = rod->edges_begin();
      uint i=0;
      
      ++eit;
      ++i;
      
      for (; eit != rod->edges_end() && i < rod->nv() - 1; ++eit, ++i) {
        ElasticRod::edge_handle eh = *eit;
//        f->setRefLength(eh, rod->getEdgeLength(eh)); //len[i]);
        
//        std::cout << rod->getEdgeLength(eh) << " " << len[i] << "\n";
        f->updateProperties();
      }
      
    }
    
    
    rod->updateProperties();
    
    //cout << (*fIt)->getName() << " = " << curr_force.norm() << endl;
  }
  
  
  
}

#include <iostream>
#include <iomanip>
             
void MayaSceneTest::reverse_timestep()
{
  static int fr = 0;
  
//  std::cout.precision(10);
  //std::cout.width(10);
  
  //std::cout << setprecision(20) << setw(10) << 2343.0 << "\n";
  
  
  for(int i=0; i<m_rods.size(); i++) {
    std::cout << "{";
    for(int j=0; j<m_rods[i]->nv(); j++) {
      std::cout << "{" << m_rods[i]->getVertex(j)[0] << ", " << m_rods[i]->getVertex(j)[1] << ", " << m_rods[i]->getVertex(j)[2]  << "}, ";
    }
    std::cout << "};\n";
    
    for(int j=0; j<m_rods[i]->ne(); j++) {
      std::cout <<  m_rods[i]->getTheta(j) << " " ;
    }
    std::cout << "};\n";
    
    for(int j=0; j<m_rods[i]->ne(); j++) {
      std::cout <<  "{" << m_rods[i]->getTangent(j) << " " ;
      std::cout <<  m_rods[i]->getReferenceDirector1(j) << " " ;
      std::cout <<  m_rods[i]->getReferenceDirector2(j) << "},  " ;
    }
    
    
    std::cout << "};\n";
    std::cout << "};\n";
    
    
  }

  
  fr++;
  if (fr == 200) {
    
    std::cout << "set udf\n";
    //reverse_findudf(1);
    std::cout << "set udf done\n";
//    reverse_solve();
//    return;
  }
  
//  reverse_solve();

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
  std::cout.precision(10);

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







