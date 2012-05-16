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
  m_active_scene(1),
  m_time(0)
{
  addDynamicsProps();
  
  //Choice of scene
  AddOption("rodshell-scene", "the scene to test", 1);

  //Basic rod options
  AddOption("rod-radius-a", "major radius of the rod", 0.01);
  AddOption("rod-radius-b", "minor radius of the rod", 0.01);
  AddOption("rod-density", "volumetric density of the rod ", 1.0);

  //Properties for the thickness-dependent linear elasticity & viscosity
  AddOption("rod-Youngs", "the Young's modulus of the rod material", 0.0f);
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
  AddOption("shell-collision-epsilon", "the distance tolerance for El Topo to flag a collision", 1e-5);

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
  &RodShellTest::setupScene4,  // umbrella opening
  &RodShellTest::setupScene5   // car sunshade folding
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
  Scalar epsilon = GetScalarOpt("shell-collision-epsilon");
  shell->setCollisionParams(proximity, epsilon, stiffness, damping);
  
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
  if (obj->nf() > 0)
  {
    // no shell render if no shell at all (for pure rod tests)
    RenderBase* shellRender = new ShellRenderer(*shell, initial_thickness);
    m_world->addRenderer(shellRender); 
  }
  RenderBase* rodRender = new RodModelRenderer(*rod);
  m_world->addRenderer(rodRender);
  
}

void RodShellTest::AtEachTimestep()
{
  m_time += m_timestep;
  
  if (m_active_scene == 5)
  {
    const Vec3d l_vel(0.2, 0, 0);
    const Vec3d r_vel(-0.2, 0, 0);
    const Scalar l_rot_rate = 0.0;
    const Scalar r_rot_rate = 1.0;
    const Scalar l_rot_max = M_PI;
    const Scalar r_rot_max = M_PI;
    const Scalar l_twist_rate = 0.5;
    const Scalar r_twist_rate = -0.5;
    const Scalar l_twist_max = M_PI * 0.5;
    const Scalar r_twist_max = -M_PI * 0.5;
    
    Vec3d l1 = obj->getVertexUndeformedPosition(m_s5_l1);
    Vec3d l2 = obj->getVertexUndeformedPosition(m_s5_l2);
    Vec3d r1 = obj->getVertexUndeformedPosition(m_s5_r1);
    Vec3d r2 = obj->getVertexUndeformedPosition(m_s5_r2);

    Vec3d lc = (l1 + l2) / 2 + l_vel * std::min(m_time, 15.0);
    Vec3d rc = (r1 + r2) / 2 + r_vel * std::min(m_time, 15.0);
    
    Scalar l_angle = (l_rot_rate >= 0 ? std::min(l_rot_max, l_rot_rate * m_time) : std::max(l_rot_max, l_rot_rate * m_time));
    Scalar r_angle = (r_rot_rate >= 0 ? std::min(r_rot_max, r_rot_rate * m_time) : std::max(r_rot_max, r_rot_rate * m_time));
    Mat3d l_rot, r_rot;
    l_rot << cos(l_angle), 0, sin(l_angle), 0, 1, 0, -sin(l_angle), 0, cos(l_angle);
//    r_rot << cos(r_angle), 0, sin(r_angle), 0, 1, 0, -sin(r_angle), 0, cos(r_angle);
    r_rot << 1, 0, 0, 0, cos(r_angle), sin(r_angle), 0, -sin(r_angle), cos(r_angle);
    
    obj->constrainVertex(m_s5_l1, lc + l_rot * (l1 - l2) / 2);
    obj->constrainVertex(m_s5_l2, lc + l_rot * (l2 - l1) / 2);
    obj->constrainVertex(m_s5_r1, rc + r_rot * (r1 - r2) / 2);
    obj->constrainVertex(m_s5_r2, rc + r_rot * (r2 - r1) / 2);
    
    Scalar l_twist = (l_twist_rate > 0 ? std::min(l_twist_max, l_twist_rate * m_time) : std::max(l_twist_max, l_twist_rate * m_time));
    Scalar r_twist = (r_twist_rate > 0 ? std::min(r_twist_max, r_twist_rate * m_time) : std::max(r_twist_max, r_twist_rate * m_time));
    
    rod->constrainEdge(m_s5_le, l_twist);
    rod->constrainEdge(m_s5_re, r_twist);
  }

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
  
  // create a rod model
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
  int resolution = GetIntOpt("shell-x-resolution"); // only one resolution
  
  //build a hexagonal grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> positions(obj);
  VertexProperty<Vec3d> velocities(obj);
  VertexProperty<Vec3d> undeformed(obj);
  VertexProperty<Vec3d> rodundeformed(obj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(obj);
  EdgeProperty<Scalar> edgeAngle(obj);
  EdgeProperty<Scalar> edgeVel(obj);
  
  VertexHandle h = obj->addVertex();
  positions[h] = Vec3d(0, 0, 0);
  velocities[h] = Vec3d(0, 0, 0);
  undeformed[h] = Vec3d(0, 0, 0);
  vertHandles.push_back(h);
  
  int nblock = 6;
  int nv_block = resolution * (resolution + 1) / 2;
  int nv_total = nv_block * nblock + 1;
  
  // hexagonal pyramid umbrella, hard-coded for now
  for (int k = 0; k < nblock; k++)
  {
    Mat2d R;
    R(0, 0) = cos(k * 2 * M_PI / nblock);
    R(0, 1) = sin(k * 2 * M_PI / nblock);
    R(1, 0) = -sin(k * 2 * M_PI / nblock);
    R(1, 1) = cos(k * 2 * M_PI / nblock);
    for (int j = 1; j <= resolution; j++) 
    {
      for (int i = 0; i < j; i++) 
      {
        Vec2d v2d;
        v2d.x() = width / resolution * ((i - j * 0.5) * sqrt(3.0) * 2 / 3);
        v2d.y() = width / resolution * j;
        v2d = R * v2d;
        
        Vec3d vert(v2d.x(), -height / resolution * j, -v2d.y());
        Vec3d undef = vert;
        
        v2d *= (vert.norm() / v2d.norm());
        Vec3d rodundef(v2d.x(), 0, -v2d.y()); // fully opened state as rod's undeformed configuration
        
        VertexHandle h = obj->addVertex();
        
        positions[h] = vert;
        velocities[h] = Vec3d(0, 0, 0);
        undeformed[h] = undef;
        rodundeformed[h] = rodundef;
        vertHandles.push_back(h);
      }
    }
  }
  
  // connect the vertices into hexagonal mesh
  std::vector<Vec3i> tris;
  for (int k = 0; k < nblock; k++)
  {
    for (int j = 0; j < resolution; j++) 
    {
      for (int i = 0; i < j; i++) 
      {
        int v1 = 1 + k * nv_block + j * (j - 1) / 2 + i;
        int v2 = 1 + k * nv_block + (j + 1) * j / 2 + i;
        int v3 = 1 + k * nv_block + (j + 1) * j / 2 + i + 1;
        obj->addFace(vertHandles[v1], vertHandles[v2], vertHandles[v3]);
      }

      for (int i = 0; i < j - 1; i++) 
      {
        int v1 = 1 + k * nv_block + j * (j - 1) / 2 + i;
        int v2 = 1 + k * nv_block + (j + 1) * j / 2 + i + 1;
        int v3 = 1 + k * nv_block + j * (j - 1) / 2 + i + 1;
        obj->addFace(vertHandles[v1], vertHandles[v2], vertHandles[v3]);
      }
      
      int v1 = (j == 0 ? 0 : (1 + ((k + 1) % nblock) * nv_block + j * (j - 1) / 2 + 0));
      int v2 = 1 + k * nv_block + (j + 1) * j / 2 + j;
      int v3 = 1 + ((k + 1) % nblock) * nv_block + (j + 1) * j / 2 + 0;
      obj->addFace(vertHandles[v1], vertHandles[v2], vertHandles[v3]);
      
      if (j > 0)
      {
        int v1 = 1 + k * nv_block + j * (j - 1) / 2 + j - 1;
        int v2 = 1 + k * nv_block + (j + 1) * j / 2 + j;
        int v3 = 1 + ((k + 1) % nblock) * nv_block + j * (j - 1) / 2 + 0;
        obj->addFace(vertHandles[v1], vertHandles[v2], vertHandles[v3]);
      }
      
    }
  }
  
  std::cout << "resolution = " << resolution << std::endl;
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
  
  //Pin the center vertex
  obj->constrainVertex(vertHandles[0], positions[vertHandles[0]]);
  
  // collect rod edges
  std::vector<EdgeHandle> rodEdges;  
  for (int k = 0; k < nblock; k++)
  {
    for (int j = 0; j < resolution; j++)
    {
      VertexHandle v1 = (j == 0 ? vertHandles[0] : (vertHandles[1 + k * nv_block + j * (j - 1) / 2 + 0]));
      VertexHandle v2 = vertHandles[1 + k * nv_block + (j + 1) * j / 2 + 0];
      for (VertexEdgeIterator veit = obj->ve_iter(v1); veit; ++veit)
      {
        VertexHandle vother = (obj->toVertex(*veit) == v1 ? obj->fromVertex(*veit) : obj->toVertex(*veit));
        if (vother == v2)
          rodEdges.push_back(*veit);
      }
    }
  }
  
  // create an empty rod model
  rod = new ElasticRodModel(obj, rodEdges, m_timestep);
  obj->addModel(rod);
  
  // set init dofs for edges
  EdgeProperty<Scalar> zeros(obj);
  zeros.assign(0);
  rod->setEdgeThetas(zeros);
  rod->setEdgeThetaVelocities(zeros);
  rod->setEdgeUndeformedThetas(zeros);
  
  rod->setUndeformedPositions(rodundeformed);
  
}

void RodShellTest::setupScene5()
{
	std::vector<Vec3d> vertices;
	std::vector<Vec3i> faces;
  std::ifstream objfile("assets/rodshelltest/tescircle400.obj");
  
	char c;
	while (objfile >> c)
	{
		switch (c)
		{
      case 'v':
			{
				Vec3d pos;
				objfile >> pos.x() >> pos.y() >> pos.z();
				vertices.push_back(pos);
			}
        break;
      case 'f':
			{
				Vec3i face;
				objfile >> face.x() >> face.y() >> face.z();
        face.x()--;
        face.y()--;
        face.z()--;
				faces.push_back(face);
			}
        break;
		}
	}
  
  objfile.close();
  
  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  
  //build a hexagonal grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> positions(obj);
  VertexProperty<Vec3d> velocities(obj);
  VertexProperty<Vec3d> undeformed(obj);
  VertexProperty<Vec3d> rodundeformed(obj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(obj);
  EdgeProperty<Scalar> edgeAngle(obj);
  EdgeProperty<Scalar> edgeVel(obj);
  
  // hexagonal pyramid umbrella, hard-coded for now
  for (size_t i = 0; i < vertices.size(); i++)
  {
    VertexHandle h = obj->addVertex();
      
    Vec3d vert = vertices[i] * width;
    Vec3d rodundef = vert;
    
    positions[h] = vert;
    velocities[h] = Vec3d(0, 0, 0);
    undeformed[h] = vert;
    rodundeformed[h] = rodundef;
    vertHandles.push_back(h);
  }
  
  // connect the vertices into hexagonal mesh
  for (size_t i = 0; i < faces.size(); i++)
  {
    obj->addFace(vertHandles[faces[i][0]], vertHandles[faces[i][1]], vertHandles[faces[i][2]]);
  }
  
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
  
  // find the two end points in the x direction
  m_s5_l1 = *obj->vertices_begin();
  m_s5_r1 = m_s5_l1;
  for (VertexIterator i = obj->vertices_begin(); i != obj->vertices_end(); ++i)
  {
    if (obj->getVertexPosition(*i).x() < obj->getVertexPosition(m_s5_l1).x())
      m_s5_l1 = *i;
    if (obj->getVertexPosition(*i).x() > obj->getVertexPosition(m_s5_r1).x())
      m_s5_r1 = *i;
  }
  
  for (VertexEdgeIterator veit = obj->ve_iter(m_s5_l1); veit; ++veit)
    if (obj->isBoundary(*veit))
      if (obj->fromVertex(*veit) == m_s5_l1)
        m_s5_l2 = obj->toVertex(*veit), m_s5_le = *veit;
  
  for (VertexEdgeIterator veit = obj->ve_iter(m_s5_r1); veit; ++veit)
    if (obj->isBoundary(*veit))
      if (obj->fromVertex(*veit) == m_s5_r1)
        m_s5_r2 = obj->toVertex(*veit), m_s5_re = *veit;
  
  std::cout << "left 1 = " << obj->getVertexPosition(m_s5_l1) << std::endl;
  std::cout << "left 2 = " << obj->getVertexPosition(m_s5_l2) << std::endl;
  std::cout << "right 1 = " << obj->getVertexPosition(m_s5_r1) << std::endl;
  std::cout << "right 2 = " << obj->getVertexPosition(m_s5_r2) << std::endl;
  
  //Pin the two end points
//  obj->constrainVertex(m_s5_l1, positions[m_s5_l1]);
//  obj->constrainVertex(m_s5_l2, positions[m_s5_l2]);
//  obj->constrainVertex(m_s5_r1, positions[m_s5_r1]);
//  obj->constrainVertex(m_s5_r2, positions[m_s5_r2]);
  
  // collect rod edges
  std::vector<EdgeHandle> rodEdges;  
  for (EdgeIterator i = obj->edges_begin(); i != obj->edges_end(); ++i)
    if (obj->isBoundary(*i))
      rodEdges.push_back(*i);
  
  // create an empty rod model
  rod = new ElasticRodModel(obj, rodEdges, m_timestep);
  obj->addModel(rod);
  
  // set init dofs for edges
  EdgeProperty<Scalar> zeros(obj);
  zeros.assign(0);
  rod->setEdgeThetas(zeros);
  rod->setEdgeThetaVelocities(zeros);
  rod->setEdgeUndeformedThetas(zeros);
  
  rod->setUndeformedPositions(rodundeformed);
  
//  rod->constrainEdge(m_s5_le, 0);
//  rod->constrainEdge(m_s5_re, 0);
  
}

