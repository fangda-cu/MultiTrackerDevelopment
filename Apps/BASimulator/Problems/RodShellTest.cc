//
//  RodShellTest.cc
//  BASim
//
//  Created by Fang Da (fang@cs.columbia.edu) on 5/12/12.
//  Copyright (c) 2012 Columbia University. All rights reserved.
//

#include "RodShellTest.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/CSTMembraneForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/DSBendingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/MNBendingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellGravityForce.hh"
#include "BASim/src/Render/ShellRenderer.hh"
#include "BASim/src/Render/RodModelRenderer.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"
#include "BASim/src/Physics/DeformableObjects/Rods/RodModelForce.hh"

RodShellTest::RodShellTest() : 
  Problem("Rod Shell Test", "Rod shell integration tests"), 
  obj(NULL), 
  rod(NULL), 
  shell(NULL), 
  stepper(NULL),
  m_active_scene(1)
{
  addDynamicsProps();
  
  //Choice of scene
  AddOption("rodshell-scene", "the scene to test", 1);

  //Basic rod options
  AddOption("rod-radius-a", "major radius of the rod", 0.01);
  AddOption("rod-radius-b", "minor radius of the rod", 0.01);
  AddOption("rod-density", "volumetric density of the rod ", 1.0);

  //Properties for the thickness-dependent linear elasticity & viscosity
  AddOption("rod-Youngs", "the Young's modulus of the rodl material", 0.0f);
  AddOption("rod-Shear", "the shear modulus of the rod material", 0.0f);
  AddOption("rod-Youngs-damping", "the damping coefficient associated with the rod's Young's modulus", 0.0f);
  AddOption("rod-Shear-damping", "the damping coefficient associated to the rod's shear modulus", 0.0f);

  //Basic shell options
  AddOption("shell-thickness", "the (initial) thickness of the shell", 0.01);
  AddOption("shell-density", "volumetric density of the shell ", 1.0);
  
  //Scene specific stuff (geometry, scene-specific forces, etc.)
  AddOption("shell-width", "the horizontal side length of the shell", 1.0);
  AddOption("shell-height", "the vertical side length of the shell", 1.0);
  AddOption("shell-x-resolution", "the number of segments along first dimension", 30);
  AddOption("shell-y-resolution", "the number of segments along second dimension", 30);
  
  AddOption("shell-initial-velocity", "starting velocity in the y direction", -0.1);
  AddOption("shell-inflate-sphere-coeff", "coefficient for inflating sphere", 0.0);
  AddOption("shell-inflate-sphere-const-pressure", "whether to use constant pressure version", false);
  
  //sheared wrinkling parameters
  AddOption("shell-rotation-rate", "the rate at which inner cylinder rotates for sheared wrinkling test", 0.0);
  AddOption("shell-bath-density", "the density of water fluid bath for sheared wrinkling test", 1.0);
  
  //Remeshing options
  AddOption("shell-remeshing", "whether to perform remeshing", 0);
  AddOption("shell-remeshing-resolution", "target edge-length", 0.1);
  AddOption("shell-remeshing-iterations", "number of remeshing iterations to run", 2);
  
  //Area-based surface tension force
  AddOption("shell-surface-tension", "surface tension coefficient of the shell", 0.0);
  
  //Properties for thickness-dependent linear elasticity & viscosity
  AddOption("shell-Poisson", "the Poisson ratio of the shell material", 0.0f);
  AddOption("shell-Youngs", "the Young's modulus of the shell material", 0.0f);
  AddOption("shell-Poisson-damping", "the damping coefficient associated to the shell's Poisson ratio", 0.0f);
  AddOption("shell-Youngs-damping", "the damping coefficient associated with the shell's Young's modulus", 0.0f);
  
  //Shell forces
  AddOption("shell-CST-stretching", "whether to apply constant-strain-triangle in-plane stretching", true);
  AddOption("shell-DS-bending", "whether to apply \"Discrete Shells\" hinge-based bending", false);
  AddOption("shell-MN-bending", "whether to apply Mid-edge Normal based bending (Morley elt)", false);
  AddOption("shell-stretching-factor", "extra scale factor to multiply stretching coefficient by", 1.0);
  AddOption("shell-bending-factor", "extra scale factor to multiple bending coefficient by", 1.0);
  
  //Collision options
  AddOption("shell-self-collision", "whether to add self-collision springs", false);
  AddOption("shell-ground-plane", "whether to add ground plane collision springs", false);
  AddOption("shell-ground-plane-height", "height of the ground plane", 0.0);
  AddOption("shell-ground-plane-velocity", "the rate at which particles colliding with the ground (bath) are pulled into it.", 0.0);
  
  AddOption("shell-collision-spring-stiffness", "stiffness coefficient of the collision springs", 0.0);
  AddOption("shell-collision-spring-damping", "damping coefficient of the collision springs", 0.0);
  AddOption("shell-collision-proximity", "the collision spring rest length and distance at which to add springs", 0.0);
  
  AddOption("shell-collision-object-file", "source SDF for object collision", "");
  AddOption("shell-collision-object-offset", "translation of the object", Vec3d(0,0,0));
  
  //Tearing options
  AddOption("shell-tearing", "whether to add tearing to the model", false);
  AddOption("shell-tearing-threshold", "the thickness threshold to use for tearing", 0.0 );
  AddOption("shell-tearing-randomness", "percent of fracture edges that will actually tear apart", 1.0 );
  AddOption("shell-ring-velocity", "velocity in the x direction for the rings", 0.25);
  
  //Timestepper options
  AddOption("integrator", "type of integrator to use for the shell", "implicit");
  
  //Solver options
  AddOption("iterations", "maximum number of iterations for the implicit method", (int) 100);
  AddOption("atol", "absolute convergence tolerance", 1e-8);
  AddOption("rtol", "relative convergence tolerance", 1e-8);
  AddOption("stol", "convergence tolerance in terms of the norm of the change in the solution between steps", 1e-8);
  AddOption("inftol", "infinity norm convergence tolerance", 1e-8);
    
}

RodShellTest::~RodShellTest()
{
  if (obj) delete obj;
  if (rod) delete rod;
  if (shell) delete shell;
  if (stepper) delete stepper;
}

typedef void (RodShellTest::*sceneFunc)();

sceneFunc rodshell_scenes[] = 
{
  0,
  &RodShellTest::setupScene1,  // shell test: vertical flat sheet
  &RodShellTest::setupScene2,  // rod test: bent twisting
  &RodShellTest::setupScene3,  // simple rod shell test: scene 1 with a rod rib
  &RodShellTest::setupScene4   // umbrella opening
};

void RodShellTest::Setup()
{
  loadDynamicsProps();
  
  //////////////////////////////////////////////////////////////////////////
  //
  // load options
  //
  //////////////////////////////////////////////////////////////////////////
  //General shell forces and properties
  Scalar shell_density = GetScalarOpt("shell-density");
  Scalar shell_thickness = GetScalarOpt("shell-thickness");
  Scalar initial_thickness = shell_thickness;
  Scalar rod_density = GetScalarOpt("rod-density");
  Scalar rod_radius_a = GetScalarOpt("rod-radius-a");
  Scalar rod_radius_b = GetScalarOpt("rod-radius-b");
  
  Vec3d gravity = GetVecOpt("gravity");
  
  Scalar surface_tension = GetScalarOpt("shell-surface-tension");
  
  Scalar shell_Youngs_modulus = GetScalarOpt("shell-Youngs");
  Scalar shell_Poisson_ratio = GetScalarOpt("shell-Poisson");
  Scalar shell_Youngs_damping = GetScalarOpt("shell-Youngs-damping");
  Scalar shell_Poisson_damping = GetScalarOpt("shell-Poisson-damping");
  
  bool cst_stretch = GetBoolOpt("shell-CST-stretching");
  bool ds_bend = GetBoolOpt("shell-DS-bending");
  bool mn_bend = GetBoolOpt("shell-MN-bending");
  
  //fudge factors to modify the elastic-viscous coefficients (so as to manipulate them separately)
  Scalar cst_scale = GetScalarOpt("shell-stretching-factor");
  Scalar ds_scale = GetScalarOpt("shell-bending-factor");
  
  std::string integrator = GetStringOpt("integrator");
  
  Scalar timestep = getDt(); //Our Rayleigh damping model relies on knowing the timestep (folds it into the damping stiffness, as in Viscous Threads)
  m_timestep = timestep;
  
  //////////////////////////////////////////////////////////////////////////
  //
  // scene setup
  //
  //////////////////////////////////////////////////////////////////////////
  //Geometry/scene specific
  int sceneChoice = GetIntOpt("rodshell-scene");
  m_active_scene = sceneChoice;
  
  //Create the base deformable object (mesh)
  obj = new DeformableObject();
  
  //Call the appropriate scene setup function.
  (*this.*rodshell_scenes[sceneChoice])();
  
  //compute the dof indexing for use in the diff_eq solver
  obj->computeDofIndexing();
  
  //////////////////////////////////////////////////////////////////////////
  //
  // rod forces
  //
  //////////////////////////////////////////////////////////////////////////
  Scalar rod_Youngs_modulus = GetScalarOpt("rod-Youngs");
  Scalar rod_Shear_modulus = GetScalarOpt("rod-Shear");
  Scalar rod_Youngs_damping = GetScalarOpt("rod-Youngs-damping");
  Scalar rod_Shear_damping = GetScalarOpt("rod-Shear-damping");
  rod->setup(rod_Youngs_modulus, rod_Youngs_damping, rod_Shear_modulus, rod_Shear_damping, m_timestep);

  //////////////////////////////////////////////////////////////////////////
  //
  // shell forces
  //
  //////////////////////////////////////////////////////////////////////////
  //Stretching and bending forces
  if(shell_Youngs_modulus != 0 || shell_Youngs_damping != 0) {
    
    //Stretching force (Constant Strain Triangle, i.e. linear FEM)
    if(cst_stretch)
      shell->addForce(new CSTMembraneForce(*shell, "CSTMembrane", cst_scale*shell_Youngs_modulus, shell_Poisson_ratio, cst_scale*shell_Youngs_damping, shell_Poisson_damping, timestep));
    
    //Bending force (Hinge-based Bending, a la Discrete Shells)
    if(ds_bend)
      shell->addForce(new DSBendingForce(*shell, "DSBending", ds_scale*shell_Youngs_modulus, shell_Poisson_ratio, ds_scale*shell_Youngs_damping, shell_Poisson_damping, timestep));
    
    //Better bending model, not currently functional.
    if(mn_bend)
      shell->addForce(new MNBendingForce(*shell, "MNBending", shell_Youngs_modulus, shell_Poisson_ratio, shell_Youngs_damping, shell_Poisson_damping, timestep));
  }
  
  //Gravity force (handled by shell only, not rod because we only need one copy of the force)
  shell->addForce(new ShellGravityForce(*shell, "Gravity", gravity)); // TODO: move gravity force to the PositionDofsModel?
  
  //////////////////////////////////////////////////////////////////////////
  //
  // rod and shell mass
  //
  //////////////////////////////////////////////////////////////////////////
  shell->setThickness(shell_thickness);
  shell->setDensity(shell_density);
  rod->setRadii(rod_radius_a, rod_radius_b);
  rod->setDensity(rod_density);
  
  shell->computeMasses();
  rod->computeMasses();
  
  for (size_t i = 0; i < rod->getForces().size(); i++)
    rod->getForces()[i]->updateStiffness();
  
  //////////////////////////////////////////////////////////////////////////
  //
  // shell specific (remeshing, collision, tearing etc)
  //
  //////////////////////////////////////////////////////////////////////////
  bool remeshing = GetIntOpt("shell-remeshing") == 1?true:false;
  Scalar remeshing_res = GetScalarOpt("shell-remeshing-resolution");
  int remeshing_its = GetIntOpt("shell-remeshing-iterations");
  shell->setRemeshing(remeshing, remeshing_res, remeshing_its);
  
  //shell->remesh(remeshing_res);
  
  Scalar stiffness = GetScalarOpt("shell-collision-spring-stiffness");
  Scalar damping = GetScalarOpt("shell-collision-spring-damping");
  Scalar proximity = GetScalarOpt("shell-collision-proximity");
  shell->setCollisionParams(proximity, stiffness, damping);
  
  bool groundPlane = GetBoolOpt("shell-ground-plane");
  Scalar gpHeight = GetScalarOpt("shell-ground-plane-height");
  Scalar gpSpeed = GetScalarOpt("shell-ground-plane-velocity");
  shell->setGroundPlane(groundPlane, gpHeight, gpSpeed);
  bool selfCollide = GetBoolOpt("shell-self-collision");
  shell->setSelfCollision(selfCollide);
  
  bool tearing = GetBoolOpt("shell-tearing");
  Scalar tearingThres = GetScalarOpt("shell-tearing-threshold");
  Scalar tearingRand = GetScalarOpt( "shell-tearing-randomness");
  //  tearingRand = clamp(tearingRand, 0.0, 1.0);
  shell->setTearing(tearing, tearingThres, tearingRand);
  
  //for(int i = 0; i < 3; ++i)
  //  shell->remesh(remeshing_res);
  
  //////////////////////////////////////////////////////////////////////////
  //
  // create time stepper
  //
  //////////////////////////////////////////////////////////////////////////
  stepper = new DefoObjTimeStepper(*obj);
  if(integrator == "symplectic")
    stepper->setDiffEqSolver(DefoObjTimeStepper::SYMPL_EULER);
  else if(integrator == "implicit")
    stepper->setDiffEqSolver(DefoObjTimeStepper::IMPL_EULER);
  else {
    std::cout << "No valid integrator specified. Existing options are 'symplectic' and 'implicit'.\n";
    assert(false);
  }
  stepper->setTimeStep(getDt());
  stepper->setMaxIterations(GetIntOpt("iterations"));
  stepper->set_stol(GetScalarOpt("stol"));
  stepper->set_atol(GetScalarOpt("atol"));
  stepper->set_rtol(GetScalarOpt("rtol"));
  stepper->set_inftol(GetScalarOpt("inftol"));
  
  //////////////////////////////////////////////////////////////////////////
  //
  // add renderers and controllers to the simulator
  //
  //////////////////////////////////////////////////////////////////////////
  m_world->addObject(obj);
  m_world->addController(stepper);
  RenderBase* shellRender = new ShellRenderer(*shell, initial_thickness);
  m_world->addRenderer(shellRender);
  RenderBase* rodRender = new RodModelRenderer(*rod);
  m_world->addRenderer(rodRender);
  
}

void RodShellTest::AtEachTimestep()
{

}

void RodShellTest::setupScene1() 
{
  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");
  
  Scalar dx = (Scalar)width / (Scalar)xresolution;
  Scalar dy = (Scalar)height / (Scalar)yresolution;
  
  //build a rectangular grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(obj);
  VertexProperty<Vec3d> positions(obj);
  VertexProperty<Vec3d> velocities(obj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(obj);
  EdgeProperty<Scalar> edgeAngle(obj);
  EdgeProperty<Scalar> edgeVel(obj);
  
  Vec3d start_vel(0,0,0);
  
  for(int j = 0; j <= yresolution; ++j) {
    for(int i = 0; i <= xresolution; ++i) {
      Vec3d vert(i*dx, j*dy, 0);//0.01*dx*sin(100*j*dy + 17*i*dx));
      if(j < 0.5*yresolution) {
        int k = j;
        int j_mod = (int)(0.5*yresolution);
        vert(1) = j_mod*dx;
        vert(2) = (k-j_mod)*dx;
      }
      Vec3d undef = vert;
      
      VertexHandle h = obj->addVertex();
      
      positions[h] = vert;
      velocities[h] = start_vel;
      undeformed[h] = undef;
      vertHandles.push_back(h);
    }
  }
  
  
  //build the faces in a 4-8 pattern
  std::vector<Vec3i> tris;
  for(int i = 0; i < xresolution; ++i) {
    for(int j = 0; j < yresolution; ++j) {
      int tl = i + (xresolution+1)*j;
      int tr = i+1 + (xresolution+1)*j;
      int bl = i + (xresolution+1)*(j+1);
      int br = i+1 + (xresolution+1)*(j+1);
      
      if((i+j)%2 == 0) {
        obj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[br]);
        obj->addFace(vertHandles[tl], vertHandles[br], vertHandles[bl]);
      }
      else {
        obj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[bl]);
        obj->addFace(vertHandles[bl], vertHandles[tr], vertHandles[br]);
      }
    }
  }
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(obj); 
  DeformableObject::face_iter fIt;
  for(fIt = obj->faces_begin(); fIt != obj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(obj, shellFaces, m_timestep);
  obj->addModel(shell);
  
  //positions
  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  //mid-edge normal variables
  shell->setEdgeUndeformed(undefAngle);
  shell->setEdgeXis(edgeAngle);
  shell->setEdgeVelocities(edgeVel);
  
  //Find highest vertex
  VertexIterator vit = obj->vertices_begin();
  Scalar highest = -10000;
  for(;vit!= obj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest) {
      highest = pos[1];
    }
  }

  //Pin all verts at or near that height
  for(vit = obj->vertices_begin();vit!= obj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest - 1e-4)
      obj->constrainVertex(*vit, pos);
  }
  
  // create an empty rod model
  rod = new ElasticRodModel(obj, std::vector<EdgeHandle>(), m_timestep);
  obj->addModel(rod);
  
  // set init dofs for edges
  EdgeProperty<Scalar> zeros(obj);
  zeros.assign(0);
  rod->setEdgeThetas(zeros);
  rod->setEdgeThetaVelocities(zeros);
  rod->setEdgeUndeformedThetas(zeros);

}

void RodShellTest::setupScene2() 
{  
  Scalar radius = 20.0;
  int nv = 50;

  // allocate vertices
  std::vector<VertexHandle> vertex_handles;
  for (int i = 0; i < nv; i++)
    vertex_handles.push_back(obj->addVertex());

  // compute init configuration
  VertexProperty<Vec3d> positions(obj);
  VertexProperty<Vec3d> velocities(obj);
  VertexProperty<Vec3d> undeformed(obj);
  for (int i = 0; i < nv; i++)
  {
    VertexHandle & h = vertex_handles[i];
    positions[h] = Vec3d(radius * cos(i * M_PI / (nv - 1)), radius * sin(i * M_PI / (nv - 1)), 0);
    velocities[h] = Vec3d::Zero();
    undeformed[h] = positions[h];
  }

  obj->setVertexPositions(positions);
  obj->setVertexVelocities(velocities);
  obj->setVertexUndeformedPositions(undeformed);
  
  // allocate edges, and set the undeformed reference frame
  std::vector<EdgeHandle> rodEdges;
  EdgeProperty<Vec3d> ref_dir(obj);
  for (int i = 0; i < nv - 1; i++)
  {
    EdgeHandle h = obj->addEdge(vertex_handles[i], vertex_handles[i + 1]);
    rodEdges.push_back(h);
    ref_dir[h] = Vec3d(0, 0, 1);
  }

  // create rod model
  rod = new ElasticRodModel(obj, rodEdges, m_timestep);
  obj->addModel(rod);

  // set init rod dofs for edges
  EdgeProperty<Scalar> zeros(obj);
  zeros.assign(0);
  rod->setEdgeThetas(zeros);
  rod->setEdgeThetaVelocities(zeros);
  rod->setEdgeUndeformedThetas(zeros);

  // rod forces, and init ref frame
  rod->setUndeformedReferenceDirector1(ref_dir);
  
  obj->constrainVertex(vertex_handles[0], positions[vertex_handles[0]]);  // fix head position
  obj->constrainVertex(vertex_handles[1], positions[vertex_handles[1]]);  // fix head tangent
  obj->constrainVertex(vertex_handles[nv - 1], positions[vertex_handles[nv - 1]]);  // fix tail position
  obj->constrainVertex(vertex_handles[nv - 2], positions[vertex_handles[nv - 2]]);  // fix tail tangent
  rod->constrainEdge(rodEdges[0], 0);      // no twist at the head
  rod->constrainEdgeVel(rodEdges[nv - 2], 0, -1, 0);  // twist rate = 1 at the tail
  
  // create an empty shell model
  FaceProperty<char> shellFaces(obj); 
  shellFaces.assign(false); // no face anyway  
  shell = new ElasticShell(obj, shellFaces, m_timestep);
  obj->addModel(shell);

  // set init shell dofs for edges
  shell->setEdgeXis(zeros);
  shell->setEdgeVelocities(zeros);
  shell->setEdgeUndeformed(zeros);
  
}

void RodShellTest::setupScene3() 
{
  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");
  
  Scalar dx = (Scalar)width / (Scalar)xresolution;
  Scalar dy = (Scalar)height / (Scalar)yresolution;
  
  //build a rectangular grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(obj);
  VertexProperty<Vec3d> positions(obj);
  VertexProperty<Vec3d> velocities(obj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(obj);
  EdgeProperty<Scalar> edgeAngle(obj);
  EdgeProperty<Scalar> edgeVel(obj);
  
  Vec3d start_vel(0,0,0);
  
  for(int j = 0; j <= yresolution; ++j) {
    for(int i = 0; i <= xresolution; ++i) {
      Vec3d vert(i*dx, j*dy, 0);//0.01*dx*sin(100*j*dy + 17*i*dx));
      if(j < 0.5*yresolution) {
        int k = j;
        int j_mod = (int)(0.5*yresolution);
        vert(1) = j_mod*dx;
        vert(2) = (k-j_mod)*dx;
      }
      Vec3d undef = vert;
      
      VertexHandle h = obj->addVertex();
      
      positions[h] = vert;
      velocities[h] = start_vel;
      undeformed[h] = undef;
      vertHandles.push_back(h);
    }
  }
  
  //build the faces in a 4-8 pattern
  std::vector<Vec3i> tris;
  for(int i = 0; i < xresolution; ++i) {
    for(int j = 0; j < yresolution; ++j) {
      int tl = i + (xresolution+1)*j;
      int tr = i+1 + (xresolution+1)*j;
      int bl = i + (xresolution+1)*(j+1);
      int br = i+1 + (xresolution+1)*(j+1);
      
      if((i+j)%2 == 0) {
        obj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[br]);
        obj->addFace(vertHandles[tl], vertHandles[br], vertHandles[bl]);
      }
      else {
        obj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[bl]);
        obj->addFace(vertHandles[bl], vertHandles[tr], vertHandles[br]);
      }
    }
  }
  
  std::cout << "resolution = " << xresolution << "x" << yresolution << std::endl;
  std::cout << "mesh nv = " << obj->nv() << " ne = " << obj->ne() << " nf = " << obj->nf() << std::endl;
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(obj); 
  DeformableObject::face_iter fIt;
  for(fIt = obj->faces_begin(); fIt != obj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(obj, shellFaces, m_timestep);
  obj->addModel(shell);
  
  //positions
  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  //mid-edge normal variables
  shell->setEdgeUndeformed(undefAngle);
  shell->setEdgeXis(edgeAngle);
  shell->setEdgeVelocities(edgeVel);
  
  //Find highest vertex
  VertexIterator vit = obj->vertices_begin();
  Scalar highest = -10000;
  for(;vit!= obj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest) {
      highest = pos[1];
    }
  }

  //Pin all verts at or near that height
  for(vit = obj->vertices_begin();vit!= obj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest - 1e-4)
      obj->constrainVertex(*vit, pos);
  }
  
  // find vertical edges in the center
  std::vector<EdgeHandle> rodEdges;
  EdgeHandle highestedge;
  Scalar bw = 0.001;
  for (EdgeIterator eit = obj->edges_begin(); eit != obj->edges_end(); ++eit)
  {
    EdgeVertexIterator evit = obj->ev_iter(*eit);
    VertexHandle v1 = *evit; ++evit;
    VertexHandle v2 = *evit; ++evit;
    Vec3d pos1 = obj->getVertexPosition(v1);
    Vec3d pos2 = obj->getVertexPosition(v2);
    
    if (pos1[0] >= width * (0.5 - bw) && pos1[0] <= width * (0.5 + bw) && pos2[0] >= width * (0.5 - bw) && pos2[0] <= width * (0.5 + bw))
    {
      rodEdges.push_back(*eit);
      if (pos1[1] >= highest - 1e-4 || pos2[1] >= highest - 1e-4)
      {
        highestedge = *eit;
      }
    }
    
  }
  
  std::cout << "rod edge count = " << rodEdges.size() << std::endl;
  
  // create an empty rod model
  rod = new ElasticRodModel(obj, rodEdges, m_timestep);
  obj->addModel(rod);
  
  // set init dofs for edges
  EdgeProperty<Scalar> zeros(obj);
  zeros.assign(0);
  rod->setEdgeThetas(zeros);
  rod->setEdgeThetaVelocities(zeros);
  rod->setEdgeUndeformedThetas(zeros);
  
  if (highestedge.isValid())
  {
    std::cout << highestedge.idx() << std::endl;
//    rod->constrainEdgeVel(highestedge, 0, 0.1, 0);
    rod->constrainEdge(highestedge, 0);
  }
  
}

void RodShellTest::setupScene4() 
{
  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");
  
  Scalar dx = (Scalar)width / (Scalar)xresolution;
  Scalar dy = (Scalar)height / (Scalar)yresolution;
  
  //build a rectangular grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(obj);
  VertexProperty<Vec3d> positions(obj);
  VertexProperty<Vec3d> velocities(obj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(obj);
  EdgeProperty<Scalar> edgeAngle(obj);
  EdgeProperty<Scalar> edgeVel(obj);
  
  Vec3d start_vel(0,0,0);
  
  for(int j = 0; j <= yresolution; ++j) {
    for(int i = 0; i <= xresolution; ++i) {
      Vec3d vert(i*dx, j*dy, 0);//0.01*dx*sin(100*j*dy + 17*i*dx));
      if(j < 0.5*yresolution) {
        int k = j;
        int j_mod = (int)(0.5*yresolution);
        vert(1) = j_mod*dx;
        vert(2) = (k-j_mod)*dx;
      }
      Vec3d undef = vert;
      
      VertexHandle h = obj->addVertex();
      
      positions[h] = vert;
      velocities[h] = start_vel;
      undeformed[h] = undef;
      vertHandles.push_back(h);
    }
  }
  
  //build the faces in a 4-8 pattern
  std::vector<Vec3i> tris;
  for(int i = 0; i < xresolution; ++i) {
    for(int j = 0; j < yresolution; ++j) {
      int tl = i + (xresolution+1)*j;
      int tr = i+1 + (xresolution+1)*j;
      int bl = i + (xresolution+1)*(j+1);
      int br = i+1 + (xresolution+1)*(j+1);
      
      if((i+j)%2 == 0) {
        obj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[br]);
        obj->addFace(vertHandles[tl], vertHandles[br], vertHandles[bl]);
      }
      else {
        obj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[bl]);
        obj->addFace(vertHandles[bl], vertHandles[tr], vertHandles[br]);
      }
    }
  }
  
  std::cout << "resolution = " << xresolution << "x" << yresolution << std::endl;
  std::cout << "mesh nv = " << obj->nv() << " ne = " << obj->ne() << " nf = " << obj->nf() << std::endl;
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(obj); 
  DeformableObject::face_iter fIt;
  for(fIt = obj->faces_begin(); fIt != obj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(obj, shellFaces, m_timestep);
  obj->addModel(shell);
  
  //positions
  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  //mid-edge normal variables
  shell->setEdgeUndeformed(undefAngle);
  shell->setEdgeXis(edgeAngle);
  shell->setEdgeVelocities(edgeVel);
  
  //Find highest vertex
  VertexIterator vit = obj->vertices_begin();
  Scalar highest = -10000;
  for(;vit!= obj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest) {
      highest = pos[1];
    }
  }
  
  //Pin all verts at or near that height
  for(vit = obj->vertices_begin();vit!= obj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest - 1e-4)
      obj->constrainVertex(*vit, pos);
  }
  
  // find vertical edges in the center
  std::vector<EdgeHandle> rodEdges;
  EdgeHandle highestedge;
  Scalar bw = 0.001;
  for (EdgeIterator eit = obj->edges_begin(); eit != obj->edges_end(); ++eit)
  {
    EdgeVertexIterator evit = obj->ev_iter(*eit);
    VertexHandle v1 = *evit; ++evit;
    VertexHandle v2 = *evit; ++evit;
    Vec3d pos1 = obj->getVertexPosition(v1);
    Vec3d pos2 = obj->getVertexPosition(v2);
    
    if (pos1[0] >= width * (0.5 - bw) && pos1[0] <= width * (0.5 + bw) && pos2[0] >= width * (0.5 - bw) && pos2[0] <= width * (0.5 + bw))
    {
      rodEdges.push_back(*eit);
      if (pos1[1] >= highest - 1e-4 || pos2[1] >= highest - 1e-4)
      {
        highestedge = *eit;
      }
    }
    
  }
  
  std::cout << "rod edge count = " << rodEdges.size() << std::endl;
  
  // create an empty rod model
  rod = new ElasticRodModel(obj, rodEdges, m_timestep);
  obj->addModel(rod);
  
  // set init dofs for edges
  EdgeProperty<Scalar> zeros(obj);
  zeros.assign(0);
  rod->setEdgeThetas(zeros);
  rod->setEdgeThetaVelocities(zeros);
  rod->setEdgeUndeformedThetas(zeros);
  
  if (highestedge.isValid())
  {
    std::cout << highestedge.idx() << std::endl;
    //    rod->constrainEdgeVel(highestedge, 0, 0.1, 0);
    rod->constrainEdge(highestedge, 0);
  }
  
}