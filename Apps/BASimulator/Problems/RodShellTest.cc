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
#include "BASim/src/Physics/DeformableObjects/Coupling/RodTwistShellFaceCouplingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Coupling/RodTwistSolidTetCouplingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Coupling/RodShellVertexJointCouplingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Coupling/RodSolidVertexJointCouplingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Coupling/ShellShellVertexJointCouplingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Coupling/ShellSolidVertexJointCouplingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Coupling/ShellSolidEdgeJointCouplingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Coupling/SolidSolidVertexJointCouplingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Coupling/SolidSolidEdgeJointCouplingForce.hh"
#include "BASim/src/Render/ShellRenderer.hh"
#include "BASim/src/Render/RodModelRenderer.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"
#include "BASim/src/Physics/DeformableObjects/Rods/RodModelForce.hh"
#include "BASim/src/Physics/DeformableObjects/Solids/SolidElasticityForce.hh"

using namespace BASim;

void load_obj(std::string file_location, std::vector<Vec3d>& vertices, std::vector<Vec3i>& faces) {

  std::ifstream objfile(file_location.c_str());

  char c;
  while (objfile >> c)
  {
    switch (c)
    {
    case 'v':
      {
        Vec3d pos;
        objfile >> pos.x() >> pos.y() >> pos.z();
        vertices.push_back(pos / 50);
      }
      break;
    case 'f':
      {
        Vec3i face;
        Vec3i normals;
        Vec3i texcoord;
        objfile >> face.x(); objfile.ignore(1); objfile >> normals.x(); objfile.ignore(1); objfile >> texcoord.x();
        objfile >> face.y(); objfile.ignore(1); objfile >> normals.y(); objfile.ignore(1); objfile >> texcoord.y();
        objfile >> face.z(); objfile.ignore(1); objfile >> normals.z(); objfile.ignore(1); objfile >> texcoord.z();
        face.x()--;
        face.y()--;
        face.z()--;
        faces.push_back(face);
      }
      break;
    case ' ':
      break;
    }
  }
  objfile.close();
  std::cout << "obj load: nv = " << vertices.size() << " nf = " << faces.size() << std::endl;
}


void load_obj2(std::string file_location, std::vector<Vec3d>& vertices, std::vector<Vec3i>& faces) {

  std::ifstream objfile(file_location.c_str());

  char c;
  while (objfile >> c)
  {
    switch (c)
    {
    case 'v':
      {
        Vec3d pos;
        objfile >> pos.x() >> pos.y() >> pos.z();
        vertices.push_back(pos / 50);
      }
      break;
    case 'f':
      {
        Vec3i face;
        Vec3i normals;
        Vec3i texcoord;
        objfile >> face.x();
        objfile >> face.y();
        objfile >> face.z();
        face.x()--;
        face.y()--;
        face.z()--;
        faces.push_back(face);
      }
      break;
    case ' ':
      break;
    }
  }
  objfile.close();
  std::cout << "obj load: nv = " << vertices.size() << " nf = " << faces.size() << std::endl;
}

RodShellTest::RodShellTest() : 
  Problem("Rod Shell Test", "Rod shell integration tests"), 
  obj(NULL), 
  rod(NULL), 
  shell(NULL), 
  solid(NULL),
  stepper(NULL),
  m_active_scene(1),
  m_time(0)
{
  addDynamicsProps();
  
  //Choice of scene
  AddOption("rodshell-scene", "the scene to test", 1);

  //Scene specific options
  AddOption("s9-rod-border", "whether or nor to put rod edges on the border of the bag model (scene 9)", true);
  AddOption("s9-ball-thickness", "thickness of the ball shell (scene 9)", 0.1f);
  
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
  
  //Coupling forces
  AddOption("rod-twist-shell-face-coupling-stiffness", "rod twist-shell face coupling force stiffness", 1.0);
  AddOption("rod-twist-shell-face-coupling-stiffness-damping", "rod twist-shell face coupling force viscous stiffness", 0.0);
  AddOption("rod-twist-solid-tet-coupling-stiffness", "rod twist-solid tet coupling force stiffness", 1.0);
  AddOption("rod-twist-solid-tet-coupling-stiffness-damping", "rod twist-solid tet coupling force viscous stiffness", 0.0);
  AddOption("rod-shell-vertex-joint-coupling-stiffness", "rod shell vertex joint coupling force stiffness", 1.0);
  AddOption("rod-shell-vertex-joint-coupling-stiffness-damping", "rod shell vertex joint coupling force viscous stiffness", 0.0);
  AddOption("rod-solid-vertex-joint-coupling-stiffness", "rod solid vertex joint coupling force stiffness", 1.0);
  AddOption("rod-solid-vertex-joint-coupling-stiffness-damping", "rod solid vertex joint coupling force viscous stiffness", 0.0);
  AddOption("shell-shell-vertex-joint-coupling-stiffness", "shell shell vertex joint coupling force stiffness", 1.0);
  AddOption("shell-shell-vertex-joint-coupling-stiffness-damping", "shell shell vertex joint coupling force viscous stiffness", 0.0);
  AddOption("shell-solid-vertex-joint-coupling-stiffness", "shell solid vertex joint coupling force stiffness", 1.0);
  AddOption("shell-solid-vertex-joint-coupling-stiffness-damping", "shell solid vertex joint coupling force viscous stiffness", 0.0);
  AddOption("shell-solid-edge-joint-coupling-stiffness", "shell solid edge joint coupling force stiffness", 1.0);
  AddOption("shell-solid-edge-joint-coupling-stiffness-damping", "shell solid edge joint coupling force viscous stiffness", 0.0);
  AddOption("solid-solid-vertex-joint-coupling-stiffness", "solid solid vertex joint coupling force stiffness", 1.0);
  AddOption("solid-solid-vertex-joint-coupling-stiffness-damping", "solid solid vertex joint coupling force viscous stiffness", 0.0);
  AddOption("solid-solid-edge-joint-coupling-stiffness", "solid solid edge joint coupling force stiffness", 1.0);
  AddOption("solid-solid-edge-joint-coupling-stiffness-damping", "solid solid edge joint coupling force viscous stiffness", 0.0);
  
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
  
  //Basic solid options
  AddOption("solid-density", "volumetric density of the solid", 1.0);

  //Properties for three-dimensional solid linear elasticity & viscosity
  AddOption("solid-Poisson", "the Poisson ratio of the solid material", 0.0f);
  AddOption("solid-Youngs", "the Young's modulus of the solid material", 0.0f);
  AddOption("solid-Poisson-damping", "the damping coefficient associated to the solid's Poisson ratio", 0.0f);
  AddOption("solid-Youngs-damping", "the damping coefficient associated with the solid's Young's modulus", 0.0f);


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
  if (solid) delete solid;
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
  &RodShellTest::setupScene5,  // car sunshade folding
  &RodShellTest::setupScene6,  // shell contraction
  &RodShellTest::setupScene7,  // collapsible tunnel folding
  &RodShellTest::setupScene8,  // collapsible tunnel opening
  &RodShellTest::setupScene9,  // balls in bag
  &RodShellTest::setupScene10,  // Christopher's first Houdini export test
  &RodShellTest::setupScene11,  // solid elasticity test
  &RodShellTest::setupScene12,  // rod twist shell face coupling test
  &RodShellTest::setupScene13,  // rod twist solid tet coupling test
  &RodShellTest::setupScene14   // rod shell joint coupling test

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
  
  Scalar solid_Youngs_modulus = GetScalarOpt("solid-Youngs");
  Scalar solid_Poisson_ratio = GetScalarOpt("solid-Poisson");
  Scalar solid_Youngs_damping = GetScalarOpt("solid-Youngs-damping");
  Scalar solid_Poisson_damping = GetScalarOpt("solid-Poisson-damping");

  Scalar solid_density = GetScalarOpt("solid-density");

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
  // solid forces
  //
  //////////////////////////////////////////////////////////////////////////

  if(solid_Youngs_modulus != 0 || solid_Youngs_damping != 0) {

    //Linear elastic forces (Constant Strain tetrahedra, i.e. linear FEM)
    solid->addForce(new SolidElasticityForce(*solid, "SolidElasticity", 
      solid_Youngs_modulus, solid_Poisson_ratio, 
      solid_Youngs_damping, solid_Poisson_damping, timestep));
  }

  //////////////////////////////////////////////////////////////////////////
  //
  // rod and shell mass
  //
  //////////////////////////////////////////////////////////////////////////
  shell->setThickness(shell_thickness);
  shell->setDensity(shell_density);
  
  rod->setRadii(rod_radius_a, rod_radius_b);
  rod->setDensity(rod_density);
  
  if(solid != 0) {
    solid->setDensity(solid_density);
  
    solid->computeMasses();
  }

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
    
  } else if (m_active_scene == 7)
  {
    const Scalar rot_rate = m_s7_rot_rate;
    const Scalar rot_tmax = m_s7_tmax;
    const Vec3d vel(0.0, 0.0, m_s7_vel);
    const Scalar vel_tmax = m_s7_tmax;
    
    Scalar angle = std::min(rot_tmax, m_time) * rot_rate;
    Mat3d lrot;
    lrot << cos(angle), -sin(angle), 0, sin(angle), cos(angle), 0, 0, 0, 1;
    Mat3d rrot;
    rrot << cos(angle), sin(angle), 0, -sin(angle), cos(angle), 0, 0, 0, 1;
    
    Vec3d lcenter = Vec3d::Zero();
    Vec3d rcenter = Vec3d::Zero();
    for (size_t i = 0; i < m_s7_leftring.size(); i++)
      lcenter += obj->getVertexUndeformedPosition(m_s7_leftring[i]) / m_s7_leftring.size();
    for (size_t i = 0; i < m_s7_rightring.size(); i++)
      rcenter += obj->getVertexUndeformedPosition(m_s7_rightring[i]) / m_s7_rightring.size();
    
    for (size_t i = 0; i < m_s7_leftring.size(); i++)
      obj->constrainVertex(m_s7_leftring[i], lcenter + std::min(vel_tmax, m_time) * vel + lrot * (obj->getVertexUndeformedPosition(m_s7_leftring[i]) - lcenter));
    for (size_t i = 0; i < m_s7_rightring.size(); i++)
      obj->constrainVertex(m_s7_rightring[i], rcenter - std::min(vel_tmax, m_time) * vel + rrot * (obj->getVertexUndeformedPosition(m_s7_rightring[i]) - rcenter));
  } else if (m_active_scene == 9)
  {
    static bool init = true;
    if (init)
    {
      init = false;
      Scalar thickness = GetScalarOpt("s9-ball-thickness");
      for (size_t i = 0; i < m_s9_ball_faces.size(); i++)
      {
        shell->setThickness(m_s9_ball_faces[i], thickness);
      }
    }
  } else if (m_active_scene == 12)
  {
    static bool init = true;
    if (init)
    {
      init = false;
      std::vector<RodTwistShellFaceCouplingForce::Stencil> stencils;
      for (size_t i = 0; i < m_s12_rod_edges.size(); i++)
      {
        RodTwistShellFaceCouplingForce::Stencil s;
        s.e = m_s12_rod_edges[i];
        VertexHandle ev1 = obj->fromVertex(s.e);
        VertexHandle ev2 = obj->toVertex(s.e);
        EdgeFaceIterator efit = obj->ef_iter(s.e);
        FaceHandle f = *efit; ++efit; assert(!efit);
        FaceVertexIterator fvit = obj->fv_iter(f);
        VertexHandle v1 = *fvit; ++fvit; assert(fvit);
        VertexHandle v2 = *fvit; ++fvit; assert(fvit);
        VertexHandle v3 = *fvit; ++fvit; assert(!fvit);
        s.v = (ev1 == v1 ? (ev2 == v2 ? v3 : v2) : (ev1 == v2 ? (ev2 == v3 ? v1 : v3) : (ev2 == v1 ? v2 : v1)));
        stencils.push_back(s);
      }
      
      Scalar stiffness =          GetScalarOpt("rod-twist-shell-face-coupling-stiffness");
      Scalar stiffness_damping =  GetScalarOpt("rod-twist-shell-face-coupling-stiffness-damping");
      obj->addForce(new RodTwistShellFaceCouplingForce(*rod, *shell, stencils, stiffness, stiffness_damping, m_timestep));
    }
  } else if (m_active_scene == 13)
  {
    static bool init = true;
    if (init)
    {
      init = false;
      std::vector<RodTwistSolidTetCouplingForce::Stencil> stencils;
      for (size_t i = 0; i < m_s13_rod_edges.size(); i++)
      {
        RodTwistSolidTetCouplingForce::Stencil s;
        s.e = m_s13_rod_edges[i];
        VertexHandle ev1 = obj->fromVertex(s.e);
        VertexHandle ev2 = obj->toVertex(s.e);
        EdgeFaceIterator efit = obj->ef_iter(s.e);
        if (efit)
        {
          FaceHandle f = *efit; ++efit; assert(efit);
          FaceVertexIterator fvit = obj->fv_iter(f);
          VertexHandle v1 = *fvit; ++fvit; assert(fvit);
          VertexHandle v2 = *fvit; ++fvit; assert(fvit);
          VertexHandle v3 = *fvit; ++fvit; assert(!fvit);
          s.v1 = (ev1 == v1 ? (ev2 == v2 ? v3 : v2) : (ev1 == v2 ? (ev2 == v3 ? v1 : v3) : (ev2 == v1 ? v2 : v1)));
          
          f = *efit; ++efit; assert(!efit);
          FaceVertexIterator fvit2 = obj->fv_iter(f);
          VertexHandle w1 = *fvit2; ++fvit2; assert(fvit2);
          VertexHandle w2 = *fvit2; ++fvit2; assert(fvit2);
          VertexHandle w3 = *fvit2; ++fvit2; assert(!fvit2);
          s.v2 = (ev1 == w1 ? (ev2 == w2 ? w3 : w2) : (ev1 == w2 ? (ev2 == w3 ? w1 : w3) : (ev2 == w1 ? w2 : w1)));
          
          std::cout << "stencil " << ev1.idx() << " " << ev2.idx() << " " << s.v1.idx() << " " << s.v2.idx() << std::endl;
          stencils.push_back(s);
        } else
        {
          std::cout << "dedge " << ev1.idx() << " " << ev2.idx() << std::endl;
        }
      }
      
      Scalar stiffness =          GetScalarOpt("rod-twist-solid-tet-coupling-stiffness");
      Scalar stiffness_damping =  GetScalarOpt("rod-twist-solid-tet-coupling-stiffness-damping");
      obj->addForce(new RodTwistSolidTetCouplingForce(*rod, *solid, stencils, stiffness, stiffness_damping, m_timestep));
    }
  } else if (m_active_scene == 14)
  {
    static bool init = true;
    if (init)
    {
      init = false;
      
      Scalar rsv_stiffness =          GetScalarOpt("rod-shell-vertex-joint-coupling-stiffness");
      Scalar rsv_stiffness_damping =  GetScalarOpt("rod-shell-vertex-joint-coupling-stiffness-damping");
      Scalar rdv_stiffness =          GetScalarOpt("rod-solid-vertex-joint-coupling-stiffness");
      Scalar rdv_stiffness_damping =  GetScalarOpt("rod-solid-vertex-joint-coupling-stiffness-damping");
      Scalar ssv_stiffness =          GetScalarOpt("shell-shell-vertex-joint-coupling-stiffness");
      Scalar ssv_stiffness_damping =  GetScalarOpt("shell-shell-vertex-joint-coupling-stiffness-damping");
      Scalar sdv_stiffness =          GetScalarOpt("shell-solid-vertex-joint-coupling-stiffness");
      Scalar sdv_stiffness_damping =  GetScalarOpt("shell-solid-vertex-joint-coupling-stiffness-damping");
      Scalar sde_stiffness =          GetScalarOpt("shell-solid-edge-joint-coupling-stiffness");
      Scalar sde_stiffness_damping =  GetScalarOpt("shell-solid-edge-joint-coupling-stiffness-damping");
      Scalar ddv_stiffness =          GetScalarOpt("solid-solid-vertex-joint-coupling-stiffness");
      Scalar ddv_stiffness_damping =  GetScalarOpt("solid-solid-vertex-joint-coupling-stiffness-damping");
      Scalar dde_stiffness =          GetScalarOpt("solid-solid-edge-joint-coupling-stiffness");
      Scalar dde_stiffness_damping =  GetScalarOpt("solid-solid-edge-joint-coupling-stiffness-damping");

      std::vector<RodShellVertexJointCouplingForce::Stencil> rsv_stencils;
      RodShellVertexJointCouplingForce::Stencil rsvs;
      rsvs.e = m_s14_edges[1];
      rsvs.f = m_s14_faces[0];
      rsv_stencils.push_back(rsvs);
      obj->addForce(new RodShellVertexJointCouplingForce(*rod, *shell, rsv_stencils, rsv_stiffness, rsv_stiffness_damping, m_timestep));
      
      std::vector<ShellShellVertexJointCouplingForce::Stencil> ssv_stencils;
      ShellShellVertexJointCouplingForce::Stencil ssvs;
      ssvs.f1 = m_s14_faces[0];
      ssvs.f2 = m_s14_faces[1];
      ssv_stencils.push_back(ssvs);
      obj->addForce(new ShellShellVertexJointCouplingForce(*shell, ssv_stencils, ssv_stiffness, ssv_stiffness_damping, m_timestep));
      
      std::vector<ShellSolidVertexJointCouplingForce::Stencil> sdv_stencils;
      ShellSolidVertexJointCouplingForce::Stencil sdvs;
      sdvs.f = m_s14_faces[1];
      sdvs.t = m_s14_tets[0];
      sdv_stencils.push_back(sdvs);
      obj->addForce(new ShellSolidVertexJointCouplingForce(*shell, *solid, sdv_stencils, sdv_stiffness, sdv_stiffness_damping, m_timestep));
      
      std::vector<SolidSolidVertexJointCouplingForce::Stencil> ddv_stencils;
      SolidSolidVertexJointCouplingForce::Stencil ddvs;
      ddvs.t1 = m_s14_tets[0];
      ddvs.t2 = m_s14_tets[1];
      ddv_stencils.push_back(ddvs);
      obj->addForce(new SolidSolidVertexJointCouplingForce(*solid, ddv_stencils, ddv_stiffness, ddv_stiffness_damping, m_timestep));
      
    }
  }
  
  std::cout << "====================================================" << std::endl;
  /*
  for (VertexIterator i = obj->vertices_begin(); i != obj->vertices_end(); ++i)
  {
    Vec3d v = obj->getVertexPosition(*i);
    std::cout << v.x() << " " << v.y() << " " << v.z() << std::endl;
  }
  for (EdgeIterator i = obj->edges_begin(); i != obj->edges_end(); ++i)
  {
    if (rod->isEdgeActive(*i))
    {
      Vec3d r1 = rod->getReferenceDirector1(*i);
      Vec3d r2 = rod->getReferenceDirector2(*i);
      Scalar t = rod->getEdgeTheta(*i);
      std::cout << r1.x() << " " << r1.y() << " " << r1.z() << " " << r2.x() << " " << r2.y() << " " << r2.z() << " " << t << std::endl;
    }
  }
  */
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
  // load an obj, hard coded path for now
	std::vector<Vec3d> vertices;
	std::vector<Vec3i> faces;
  load_obj("assets/rodshelltest/tescircle400.obj", vertices, faces);
  
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
  
  // vertex positions from obj file
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

void RodShellTest::setupScene6()
{
  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");

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
  
  for (int j = 0; j <= yresolution; j++)
  {
    for (int i = 0; i <= xresolution; i++)
    {
      VertexHandle h = obj->addVertex();
      
      Vec3d vert(width * (i - xresolution / 2) / xresolution, -height * j / yresolution, 0.01 * sin(i * M_PI / 2));
      Vec3d rodundef = vert;
      rodundef.x() *= 0.25;
      
      positions[h] = vert;
      velocities[h] = Vec3d(0, 0, 0);
      undeformed[h] = vert;
      rodundeformed[h] = rodundef;
      vertHandles.push_back(h);
    }
  }
  
  //build the faces in a 4-8 pattern
  for(int i = 0; i < xresolution; ++i) 
  {
    for(int j = 0; j < yresolution; ++j) 
    {
      int tl = i + (xresolution+1)*j;
      int tr = i+1 + (xresolution+1)*j;
      int bl = i + (xresolution+1)*(j+1);
      int br = i+1 + (xresolution+1)*(j+1);
      
      if((i+j)%2 == 0) 
      {
        obj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[br]);
        obj->addFace(vertHandles[tl], vertHandles[br], vertHandles[bl]);
      } else 
      {
        obj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[bl]);
        obj->addFace(vertHandles[bl], vertHandles[tr], vertHandles[br]);
      }
    }
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
  
  //Find highest vertex
  VertexIterator vit = obj->vertices_begin();
  Scalar highest = -10000;
  for(;vit!= obj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest) {
      highest = pos[1];
    }
  }
  
  //pin the vertices on the top border
  VertexHandle left = *(obj->vertices_begin());
  VertexHandle right = left;
  for (vit = obj->vertices_begin();vit!= obj->vertices_end(); ++vit) 
  {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest - 1e-4)
    {
      static int a = 0;
      
      if (obj->getVertexPosition(*vit).x() < obj->getVertexPosition(left).x())
        left = *vit;
      if (obj->getVertexPosition(*vit).x() > obj->getVertexPosition(right).x())
        right = *vit;
      
      if (a % 2 == 0)
        obj->constrainVertex(*vit, new FixedVelocityConstraint(obj->getVertexPosition(*vit), -obj->getVertexPosition(*vit) * 0.05, 0));
      a++;
    }
  }
//  obj->constrainVertex(left, new FixedVelocityConstraint(obj->getVertexPosition(left), Vec3d(width * 0.02, 0, 0), 0));
//  obj->constrainVertex(right, new FixedVelocityConstraint(obj->getVertexPosition(right), Vec3d(-width * 0.02, 0, 0), 0));
  
  // collect rod edges
  std::vector<EdgeHandle> rodEdges;  
  for (EdgeIterator i = obj->edges_begin(); i != obj->edges_end(); ++i)
    if (obj->isBoundary(*i))
      rodEdges.push_back(*i);
  rodEdges.clear();
  
  // create a rod model
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

EdgeHandle edgeBetweenVertices(TopologicalObject & obj, VertexHandle & v1, VertexHandle v2)
{
  for (VertexEdgeIterator veit = obj.ve_iter(v1); veit; ++veit)
  {
    assert (obj.fromVertex(*veit) == v1 || obj.toVertex(*veit) == v1);
    VertexHandle v_other = (obj.fromVertex(*veit) == v1 ? obj.toVertex(*veit) : obj.fromVertex(*veit));
    assert (v_other != v1);
    if (v_other == v2)
      return *veit;
  }
  return EdgeHandle();
}

void RodShellTest::setupScene7()
{
  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");
  
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> positions(obj);
  VertexProperty<Vec3d> velocities(obj);
  VertexProperty<Vec3d> undeformed(obj);
  VertexProperty<Vec3d> rodundeformed(obj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(obj);
  EdgeProperty<Scalar> edgeAngle(obj);
  EdgeProperty<Scalar> edgeVel(obj);
  
  for (int j = 0; j <= yresolution; j++)
  {
    for (int i = 0; i < xresolution; i++)
    {
      VertexHandle h = obj->addVertex();
      
      Scalar noise = 0.01;
      Vec3d vert(width * cos(i * 2 * M_PI / xresolution) + width * noise * ((double)rand() / RAND_MAX - 0.5), 
                 width * sin(i * 2 * M_PI / xresolution) + width * noise * ((double)rand() / RAND_MAX - 0.5), 
                 height * (j - yresolution / 2) / yresolution);
      Vec3d rodundef = vert;
      
      positions[h] = vert;
      velocities[h] = Vec3d(0, 0, 0);
      undeformed[h] = vert;
      rodundeformed[h] = rodundef;
      vertHandles.push_back(h);
    }
  }
  
  //build the faces in a 4-8 pattern
  for(int i = 0; i < xresolution; ++i) 
  {
    for(int j = 0; j < yresolution; ++j) 
    {
      int tl = i + (xresolution)*j;
      int tr = (i+1)%xresolution + (xresolution)*j;
      int bl = i + (xresolution)*(j+1);
      int br = (i+1)%xresolution + (xresolution)*(j+1);
      
      if((i+j)%2 == 0) 
      {
        obj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[br]);
        obj->addFace(vertHandles[tl], vertHandles[br], vertHandles[bl]);
      } else 
      {
        obj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[bl]);
        obj->addFace(vertHandles[bl], vertHandles[tr], vertHandles[br]);
      }
    }
  }
  
  std::cout << "mesh nv = " << obj->nv() << " ne = " << obj->ne() << " nf = " << obj->nf() << std::endl;
  
  // compute the compressing parameters
  m_s7_tmax = 20;
  m_s7_vel = height / m_s7_tmax;
  Scalar dydx = (height / yresolution) / (2 * M_PI * width / xresolution);
  Scalar init_rot = yresolution * 2 * M_PI / xresolution;
  m_s7_rot_rate = -init_rot * (sqrt(dydx * dydx + 1) - 1) / m_s7_tmax;
  
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
  
  // pin (and compress) the two rings on both ends
  for (int i = 0; i < xresolution; i++)
  {
    int idx = i;
//    obj->constrainVertex(vertHandles[idx], new FixedVelocityConstraint(obj->getVertexPosition(vertHandles[idx]), Vec3d(0, 0, height * 0.05), 0));
//    obj->constrainVertex(vertHandles[idx], obj->getVertexPosition(vertHandles[idx]));
    m_s7_leftring.push_back(vertHandles[idx]);
    idx = i + yresolution * xresolution;
//    obj->constrainVertex(vertHandles[idx], new FixedVelocityConstraint(obj->getVertexPosition(vertHandles[idx]), Vec3d(0, 0, -height * 0.05), 0));
    m_s7_rightring.push_back(vertHandles[idx]);
  }
  
  // collect rod edges
  std::vector<EdgeHandle> rodEdges;  
  for (int i = 0; i < xresolution; i++)
  {
    EdgeHandle e;
    e = edgeBetweenVertices(*obj, vertHandles[i], vertHandles[(i + 1) % xresolution]);
    assert (e.isValid());
//    if (i == (xresolution - 1) % xresolution)
      rodEdges.push_back(e);
    
    e = edgeBetweenVertices(*obj, vertHandles[i + yresolution * xresolution], vertHandles[(i + 1) % xresolution + yresolution * xresolution]);
    assert (e.isValid());
//    if (i == (yresolution) % xresolution)
      rodEdges.push_back(e);
  }
  
  for (int i = 0; i < yresolution; i++)
  {
    EdgeHandle e;
    e = edgeBetweenVertices(*obj, vertHandles[i * xresolution + (i + 0) % xresolution], vertHandles[(i + 1) * xresolution + (i + 1) % xresolution]);
    assert (e.isValid());
    rodEdges.push_back(e);
  }
  
  // create a rod model
  rod = new ElasticRodModel(obj, rodEdges, m_timestep);
  obj->addModel(rod);
  
  // set init dofs for edges
  EdgeProperty<Scalar> zeros(obj);
  zeros.assign(0);
  rod->setEdgeThetas(zeros);
  rod->setEdgeThetaVelocities(zeros);
  rod->setEdgeUndeformedThetas(zeros);
  
  rod->setUndeformedPositions(rodundeformed);

  EdgeProperty<Vec3d> undefref(obj);
  for (size_t i = 0; i < rodEdges.size(); i++)
  {
    Vec3d v1 = obj->getVertexPosition(obj->fromVertex(rodEdges[i]));
    Vec3d v2 = obj->getVertexPosition(obj->toVertex(rodEdges[i]));
    Vec3d t = v2 - v1;
    undefref[rodEdges[i]] = Vec3d(-t.y(), t.x(), 0);
    undefref[rodEdges[i]].normalize();
  }
  rod->setUndeformedReferenceDirector1(undefref);
  
}

void RodShellTest::setupScene8()
{
  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");
  
  // compute the compressing parameters
  Scalar dy = height / yresolution;
  Scalar dx = 2 * M_PI * width / xresolution;
  Scalar undef_rot = yresolution * 2 * M_PI / xresolution;
  Scalar init_len = 0.1;
  Scalar extra_rot = -undef_rot * (sqrt(square(dx) + square(dy) - square(dy * init_len)) / dx - 1);
  
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> positions(obj);
  VertexProperty<Vec3d> velocities(obj);
  VertexProperty<Vec3d> undeformed(obj);
  VertexProperty<Vec3d> rodundeformed(obj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(obj);
  EdgeProperty<Scalar> edgeAngle(obj);
  EdgeProperty<Scalar> edgeVel(obj);

  // hard-coded init configuration file for now
  std::ifstream initconfig_file("assets/rodshelltest/scene7_s_t7.65.txt");
  
  for (int j = 0; j <= yresolution; j++)
  {
    for (int i = 0; i < xresolution; i++)
    {
      VertexHandle h = obj->addVertex();
      
//      Vec3d vert(width * cos(i * 2 * M_PI / xresolution - j * extra_rot / yresolution), 
//                 width * sin(i * 2 * M_PI / xresolution - j * extra_rot / yresolution), 
//                 height * init_len * (j - yresolution / 2) / yresolution);
      Vec3d vert;
      initconfig_file >> vert.x() >> vert.y() >> vert.z();
      
      Vec3d undef(width * cos(i * 2 * M_PI / xresolution), 
                 width * sin(i * 2 * M_PI / xresolution), 
                 height * (j - yresolution / 2) / yresolution);
      
      Vec3d rodundef = undef;
      
      positions[h] = vert;
      velocities[h] = Vec3d(0, 0, 0);
      undeformed[h] = undef;
      rodundeformed[h] = rodundef;
      vertHandles.push_back(h);
    }
  }
    
  //build the faces in a 4-8 pattern
  for(int i = 0; i < xresolution; ++i) 
  {
    for(int j = 0; j < yresolution; ++j) 
    {
      int tl = i + (xresolution)*j;
      int tr = (i+1)%xresolution + (xresolution)*j;
      int bl = i + (xresolution)*(j+1);
      int br = (i+1)%xresolution + (xresolution)*(j+1);
      
      if((i+j)%2 == 0) 
      {
        obj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[br]);
        obj->addFace(vertHandles[tl], vertHandles[br], vertHandles[bl]);
      } else 
      {
        obj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[bl]);
        obj->addFace(vertHandles[bl], vertHandles[tr], vertHandles[br]);
      }
    }
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
  
  // collect rod edges
  std::vector<EdgeHandle> rodEdges;  
  for (int i = 0; i < xresolution; i++)
  {
    EdgeHandle e;
    e = edgeBetweenVertices(*obj, vertHandles[i], vertHandles[(i + 1) % xresolution]);
    assert (e.isValid());
    rodEdges.push_back(e);
    
    e = edgeBetweenVertices(*obj, vertHandles[i + yresolution * xresolution], vertHandles[(i + 1) % xresolution + yresolution * xresolution]);
    assert (e.isValid());
    rodEdges.push_back(e);
  }

  for (int i = 0; i < yresolution; i++)
  {
    EdgeHandle e;
    e = edgeBetweenVertices(*obj, vertHandles[i * xresolution + (i + 0) % xresolution], vertHandles[(i + 1) * xresolution + (i + 1) % xresolution]);
    assert (e.isValid());
    rodEdges.push_back(e);
  }
  
  // create a rod model
  rod = new ElasticRodModel(obj, rodEdges, m_timestep);
  obj->addModel(rod);
  
  // set init dofs for edges
  EdgeProperty<Scalar> zeros(obj);
  zeros.assign(0);
  rod->setEdgeThetas(zeros);
  rod->setEdgeThetaVelocities(zeros);
  rod->setEdgeUndeformedThetas(zeros);
  
  rod->setUndeformedPositions(rodundeformed);
  
  EdgeProperty<Vec3d> undefref(obj);
  for (size_t i = 0; i < rodEdges.size(); i++)
  {
    Vec3d v1 = rodundeformed[obj->fromVertex(rodEdges[i])];
    Vec3d v2 = rodundeformed[obj->toVertex(rodEdges[i])];
    Vec3d t = v2 - v1;
    undefref[rodEdges[i]] = Vec3d(-t.y(), t.x(), 0);
    undefref[rodEdges[i]].normalize();
  }
  rod->setUndeformedReferenceDirector1(undefref);
 
  obj->computeDofIndexing();

  Scalar rod_Youngs_modulus = GetScalarOpt("rod-Youngs");
  Scalar rod_Shear_modulus = GetScalarOpt("rod-Shear");
  Scalar rod_Youngs_damping = GetScalarOpt("rod-Youngs-damping");
  Scalar rod_Shear_damping = GetScalarOpt("rod-Shear-damping");
  rod->setup(rod_Youngs_modulus, rod_Youngs_damping, rod_Shear_modulus, rod_Shear_damping, m_timestep);

  EdgeProperty<Scalar> thetas(obj);
  for (EdgeIterator i = obj->edges_begin(); i != obj->edges_end(); ++i)
  {
    if (rod->isEdgeActive(*i))
    {
      thetas[*i] = rod->getEdgeTheta(*i);
    }
  }
  
  for (EdgeIterator i = obj->edges_begin(); i != obj->edges_end(); ++i)
  {
    if (rod->isEdgeActive(*i))
    {
      std::cout << "active edge: " << i->idx() << std::endl;
      
      Vec3d t = obj->getVertexPosition(obj->toVertex(*i)) - obj->getVertexPosition(obj->fromVertex(*i));
      Vec3d init_r1, init_r2;
      Scalar init_theta;
      initconfig_file >> init_r1.x() >> init_r1.y() >> init_r1.z() >> init_r2.x() >> init_r2.y() >> init_r2.z() >> init_theta;
      if (!approxEq(t.dot(init_r1), 0, 1e-4) || !approxEq(t.dot(init_r2), 0, 1e-4))
        std::cout << "t = " << t << " init r1 = " << init_r1 << " init r2 = " << init_r2 << " dot1 = " << t.dot(init_r1) << " dot2 = " << t.dot(init_r2) << std::endl;
      Vec3d r1 = rod->getReferenceDirector1(*i);
      Vec3d r2 = rod->getReferenceDirector2(*i);
      if (!approxEq(t.dot(r1), 0, 1e-4) || !approxEq(t.dot(r2), 0, 1e-4))
        std::cout << "t = " << t << " r1 = " << r1 << " r2 = " << r2 << " dot1 = " << t.dot(r1) << " dot2 = " << t.dot(r2) << std::endl;

      Scalar ca = cos(init_theta);
      Scalar sa = sin(init_theta);
      Vec3d m1 = ca * init_r1 + sa * init_r2;
      Vec3d m2 = -sa * init_r1 + ca * init_r2;
      
      Scalar theta = signedAngle(r1, m1, t);
      thetas[*i] = theta;
      rod->setEdgeThetas(thetas);
      rod->updateProperties();
      Vec3d rodm1 = rod->getMaterialDirector1(*i);
      Vec3d rodm2 = rod->getMaterialDirector2(*i);
      if (!approxEq(rodm1.cross(m1).norm(), 0, 1e-6) || !approxEq(rodm2.cross(m2).norm(), 0, 1e-6))
        std::cout << "m1 = " << m1 << " m2 = " << m2 << " rod m1 = " << rodm1 << " rod m2 = " << rodm2 << std::endl;
      std::cout << "theta = " << theta << "/" << init_theta << std::endl;
    }
  }
}

void RodShellTest::setupScene9()
{
  // load an obj, hard coded path for now
  std::vector<Vec3d> vertices;
  std::vector<Vec3i> faces;
  load_obj("assets/rodshelltest/bag4.obj", vertices, faces);

  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");
  bool include_rod_border = GetBoolOpt("s9-rod-border");
  
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
  
  // vertex positions from obj file
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
  int nv = vertices.size();

  // create K UV spheres of radius R with N tessellation (thin shell spheres, i.e. christmas balls)
  int K = 2;
  int N = 6;
  Scalar R = 1.0;
  Scalar D = 2.0;
  for (int k = 0; k < K; k++)
    for (int i = 0; i <= N; i++)
      for (int j = 0; j < ((i == 0 || i == N) ? 1 : 2 * N); j++)
      {
        VertexHandle h = obj->addVertex();
        Scalar delta = k * 2 * M_PI / K + 0.2 * M_PI;
        Vec3d ballcenter = Vec3d(0, 15, -2) + Vec3d(1, 0, 0) * D * cos(delta) + Vec3d(0, 1, 0) * D * sin(delta);
        Scalar theta = j * 2 * M_PI / 2 / N;
        Scalar alpha = i * M_PI / N;
        Vec3d vert = ballcenter + Vec3d(R * cos(theta) * sin(alpha), R * sin(theta) * sin(alpha), -R * cos(alpha));
        
        positions[h] = vert;
        velocities[h] = Vec3d(0, 0, 0);
        undeformed[h] = vert;
        rodundeformed[h] = vert;
        vertHandles.push_back(h);
      }

  // connect the vertices into hexagonal mesh
  for (size_t i = 0; i < faces.size(); i++)
  {
    Vec3d v1 = positions[vertHandles[faces[i][0]]];
    Vec3d v2 = positions[vertHandles[faces[i][1]]];
    Vec3d v3 = positions[vertHandles[faces[i][2]]];
    if (v1.y() > 11.699 && v2.y() > 11.699 && v3.y() > 11.699)
    {
      if (v1.z() == v2.z()) obj->addEdge(vertHandles[faces[i][0]], vertHandles[faces[i][1]]);
      if (v1.z() == v3.z()) obj->addEdge(vertHandles[faces[i][0]], vertHandles[faces[i][2]]);
      if (v2.z() == v3.z()) obj->addEdge(vertHandles[faces[i][1]], vertHandles[faces[i][2]]);
    } else
    {
      obj->addFace(vertHandles[faces[i][0]], vertHandles[faces[i][1]], vertHandles[faces[i][2]]); 
    }
  }
  
  for (int k = 0; k < K; k++)
    for (int i = 0; i < N; i++)
      for (int j = 0; j < 2 * N; j++)
      {
        if (i == 0)
        {
          int v1 = nv + (2 * (N - 1) * N + 2) * k + 0;
          int v2 = nv + (2 * (N - 1) * N + 2) * k + 1 + j;
          int v3 = nv + (2 * (N - 1) * N + 2) * k + 1 + (j + 1) % (2 * N);
          FaceHandle h = obj->addFace(vertHandles[v1], vertHandles[v2], vertHandles[v3]); 
          m_s9_ball_faces.push_back(h);
        } else if (i == N - 1)
        {
          int v1 = nv + (2 * (N - 1) * N + 2) * k + 1 + 2 * (N - 2) * N + j;
          int v2 = nv + (2 * (N - 1) * N + 2) * k + 1 + 2 * (N - 1) * N;
          int v3 = nv + (2 * (N - 1) * N + 2) * k + 1 + 2 * (N - 2) * N + (j + 1) % (2 * N);
          FaceHandle h = obj->addFace(vertHandles[v1], vertHandles[v2], vertHandles[v3]); 
          m_s9_ball_faces.push_back(h);
        } else
        {
          int v1 = nv + (2 * (N - 1) * N + 2) * k + 1 + 2 * (i - 1) * N + j;
          int v2 = nv + (2 * (N - 1) * N + 2) * k + 1 + 2 * i * N + j;
          int v3 = nv + (2 * (N - 1) * N + 2) * k + 1 + 2 * (i - 1) * N + (j + 1) % (2 * N);
          FaceHandle h = obj->addFace(vertHandles[v1], vertHandles[v2], vertHandles[v3]); 
          m_s9_ball_faces.push_back(h);
          v1 = nv + (2 * (N - 1) * N + 2) * k + 1 + 2 * (i - 1) * N + (j + 1) % (2 * N);
          v2 = nv + (2 * (N - 1) * N + 2) * k + 1 + 2 * i * N + j;
          v3 = nv + (2 * (N - 1) * N + 2) * k + 1 + 2 * i * N + (j + 1) % (2 * N);
          h = obj->addFace(vertHandles[v1], vertHandles[v2], vertHandles[v3]); 
          m_s9_ball_faces.push_back(h);
        }
      }
  
//  std::map<EdgeHandle, int> edgevalence;
//  for (int i = 0; i < faces.size(); i++)
//  {
//    EdgeHandle e1 = edgeBetweenVertices(*obj, vertHandles[faces[i].x()], vertHandles[faces[i].y()]);
//    EdgeHandle e2 = edgeBetweenVertices(*obj, vertHandles[faces[i].y()], vertHandles[faces[i].z()]);
//    EdgeHandle e3 = edgeBetweenVertices(*obj, vertHandles[faces[i].z()], vertHandles[faces[i].x()]);
//    int a;
//    a = edgevalence[e1];
//    edgevalence[e1] = a + 1;
//    a = edgevalence[e2];
//    edgevalence[e2] = a + 1;
//    a = edgevalence[e3];
//    edgevalence[e3] = a + 1;
//  }
//  std::vector<EdgeHandle> hv;
//  for (std::map<EdgeHandle, int>::iterator i = edgevalence.begin(); i != edgevalence.end(); i++)
//    if (i->second > 2)
//      hv.push_back(i->first);
//  assert(hv.size() == 2);
//  std::cout << positions[obj->fromVertex(hv[0])] << " " << positions[obj->toVertex(hv[0])] << std::endl;
//  std::cout << positions[obj->fromVertex(hv[1])] << " " << positions[obj->toVertex(hv[1])] << std::endl;
  
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
  
  // collect rod edges
  std::vector<EdgeHandle> rodEdges;  
  for (EdgeIterator i = obj->edges_begin(); i != obj->edges_end(); ++i)
  {
    if (obj->isBoundary(*i))
      continue;
        
    EdgeFaceIterator efit = obj->ef_iter(*i);
    if (efit)
    {
      FaceHandle f1 = *efit; ++efit;
      FaceHandle f2 = *efit; ++efit;
      assert(!efit);
      VertexHandle v1 = obj->fromVertex(*i);
      VertexHandle v2 = obj->toVertex(*i);
      VertexHandle v1o, v2o;
      FaceVertexIterator fvit;
      fvit = obj->fv_iter(f1);
      if (*fvit != v1 && *fvit != v2) v1o = *fvit; ++fvit;
      if (*fvit != v1 && *fvit != v2) v1o = *fvit; ++fvit;
      if (*fvit != v1 && *fvit != v2) v1o = *fvit; ++fvit;
      assert(!fvit);
      assert(v1o.isValid());
      fvit = obj->fv_iter(f2);
      if (*fvit != v1 && *fvit != v2) v2o = *fvit; ++fvit;
      if (*fvit != v1 && *fvit != v2) v2o = *fvit; ++fvit;
      if (*fvit != v1 && *fvit != v2) v2o = *fvit; ++fvit;
      assert(!fvit);
      assert(v2o.isValid());
      if (v1.idx() < nv && obj->getVertexPosition(v1).z() == obj->getVertexPosition(v2).z() && obj->getVertexPosition(v1o).z() != obj->getVertexPosition(v2o).z())
      {
        if (include_rod_border)
          rodEdges.push_back(*i);
      }
    } else
    {
      rodEdges.push_back(*i);
    }
    
  }
    
  // create a rod model
  rod = new ElasticRodModel(obj, rodEdges, m_timestep);
  obj->addModel(rod);
  
  // set init dofs for edges
  EdgeProperty<Scalar> zeros(obj);
  zeros.assign(0);
  rod->setEdgeThetas(zeros);
  rod->setEdgeThetaVelocities(zeros);
  rod->setEdgeUndeformedThetas(zeros);
  
  rod->setUndeformedPositions(rodundeformed);
  
  EdgeProperty<Vec3d> ref_dir(obj);
  for (unsigned int i = 0; i < rodEdges.size(); i++)
  {
    ref_dir[rodEdges[i]] = Vec3d(0, 0, 1);
  }
  rod->setUndeformedReferenceDirector1(ref_dir);
}


void RodShellTest::setupScene10()
{
  // load an obj, hard coded path for now
  std::vector<Vec3d> vertices;
  std::vector<Vec3i> faces;
  load_obj2("../../assets/rodshelltest/scene10_data/shellgeom.obj", vertices, faces);

  //load the rod vertices list
  std::ifstream infile("../../assets/rodshelltest/scene10_data/rodvertices.txt");
  std::vector<int> rod_vert_indices;
  do 
  {
    int vertex_number;
    infile >> vertex_number;
    if(infile.eof()) break;
    rod_vert_indices.push_back(vertex_number);
  } while (!infile.eof());
  
  std::cout << "Rod verts:";
  for(unsigned int i = 0; i < rod_vert_indices.size(); ++i)
    std::cout << rod_vert_indices[i] << " ";
  std::cout << std::endl;
  
  
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> positions(obj);
  VertexProperty<Vec3d> velocities(obj);
  VertexProperty<Vec3d> undeformed(obj);
  VertexProperty<Vec3d> rodundeformed(obj);

  //edge properties
  EdgeProperty<Scalar> undefAngle(obj);
  EdgeProperty<Scalar> edgeAngle(obj);
  EdgeProperty<Scalar> edgeVel(obj);

  // vertex positions from obj file
  for (size_t i = 0; i < vertices.size(); i++)
  {
    VertexHandle h = obj->addVertex();

    Vec3d vert = vertices[i];
    Vec3d rodundef = vert;

    positions[h] = vert;
    velocities[h] = Vec3d(0, 0, 0);
    undeformed[h] = vert;
    rodundeformed[h] = rodundef;
    vertHandles.push_back(h);
  }
  int nv = vertices.size();


  std::cout << "mesh nv = " << obj->nv() << " ne = " << obj->ne() << " nf = " << obj->nf() << std::endl;

  for(unsigned int i = 0; i < faces.size(); ++i) {
    
    FaceHandle newFace = obj->addFace(vertHandles[faces[i][0]], vertHandles[faces[i][1]], vertHandles[faces[i][2]]);

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

  
  // collect rod edges
  std::vector<EdgeHandle> rodEdges;  
  for (EdgeIterator i = obj->edges_begin(); i != obj->edges_end(); ++i)
  {
    VertexHandle v0 = obj->fromVertex(*i);
    VertexHandle v1 = obj->toVertex(*i);
    if(std::find(rod_vert_indices.begin(), rod_vert_indices.end(), v0.idx()) != rod_vert_indices.end() &&
       std::find(rod_vert_indices.begin(), rod_vert_indices.end(), v1.idx()) != rod_vert_indices.end())
       rodEdges.push_back(*i);

  }
  std::cout << "Rod edges? " << rodEdges.size() << std::endl;

  // create a rod model
  rod = new ElasticRodModel(obj, rodEdges, m_timestep);
  obj->addModel(rod);

  // set init dofs for edges
  EdgeProperty<Scalar> zeros(obj);
  zeros.assign(0);
  rod->setEdgeThetas(zeros);
  rod->setEdgeThetaVelocities(zeros);
  rod->setEdgeUndeformedThetas(zeros);

  rod->setUndeformedPositions(rodundeformed);

  EdgeProperty<Vec3d> ref_dir(obj);
  for (unsigned int i = 0; i < rodEdges.size(); i++)
  {
    ref_dir[rodEdges[i]] = Vec3d(0, 0, 1);
  }
  rod->setUndeformedReferenceDirector1(ref_dir);
  
}

void RodShellTest::setupScene11()
{
  // load an obj, hard coded path for now
  std::vector<Vec3d> vertices;

  // Create 4 vertices
  vertices.push_back(Vec3d(0,0,0));
  vertices.push_back(Vec3d(1,0,0));
  vertices.push_back(Vec3d(0,1,0));
  vertices.push_back(Vec3d(0,0,1));
  
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> positions(obj);
  VertexProperty<Vec3d> velocities(obj);
  VertexProperty<Vec3d> undeformed(obj);
  
  // vertex positions from obj file
  for (size_t i = 0; i < vertices.size(); i++)
  {
    VertexHandle h = obj->addVertex();

    Vec3d vert = vertices[i];
    positions[h] = vert;
    if(i == 0)
      positions[h] = Vec3d(-0.2, -0.2, -0.2);
    velocities[h] = Vec3d(0, 0, 0);
    undeformed[h] = vert;
    vertHandles.push_back(h);
  }
  int nv = vertices.size();

  FaceHandle face0 = obj->addFace(vertHandles[0], vertHandles[1], vertHandles[2]);
  FaceHandle face1 = obj->addFace(vertHandles[0], vertHandles[3], vertHandles[1]);
  FaceHandle face2 = obj->addFace(vertHandles[0], vertHandles[2], vertHandles[3]);
  FaceHandle face3 = obj->addFace(vertHandles[1], vertHandles[3], vertHandles[2]);
  
  TetHandle tet = obj->addTet(face0, face1, face2, face3);

  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  TetProperty<char> solidTets(obj); 
  DeformableObject::tet_iter tIt;
  for(tIt = obj->tets_begin(); tIt != obj->tets_end(); ++tIt)
    solidTets[*tIt] = true;

  //now create the physical model to hang on the mesh
  solid = new ElasticSolid(obj, solidTets, m_timestep);
  obj->addModel(solid);

  //positions
  solid->setVertexUndeformed(undeformed);
  solid->setVertexPositions(positions);
  solid->setVertexVelocities(velocities);

  //Find highest vertex
  VertexIterator vit = obj->vertices_begin();
  Scalar highest = -10000;
  for(;vit!= obj->vertices_end(); ++vit) {
    Vec3d pos = solid->getVertexPosition(*vit);
    if(pos[1] >= highest) {
      highest = pos[1];
    }
  }

  //Pin all verts at or near that height
  for(vit = obj->vertices_begin();vit!= obj->vertices_end(); ++vit) {
    Vec3d pos = solid->getVertexPosition(*vit);
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

void RodShellTest::setupScene12() 
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
  Scalar leftmost = 10000;
  for(;vit!= obj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest)
      highest = pos[1];
    if(pos[0] <= leftmost)
      leftmost = pos[0];
  }
  
  //Pin all verts at or near that height, on the left
//  for(vit = obj->vertices_begin();vit!= obj->vertices_end(); ++vit) {
//    Vec3d pos = shell->getVertexPosition(*vit);
//    if(pos[1] >= highest - 1e-4)
//      obj->constrainVertex(*vit, pos);
//  }
  
  // find top edges
  std::vector<EdgeHandle> rodEdges;
  EdgeHandle leftmostedge;
  for (EdgeIterator eit = obj->edges_begin(); eit != obj->edges_end(); ++eit)
  {
    EdgeVertexIterator evit = obj->ev_iter(*eit);
    VertexHandle v1 = *evit; ++evit;
    VertexHandle v2 = *evit; ++evit;
    Vec3d pos1 = obj->getVertexPosition(v1);
    Vec3d pos2 = obj->getVertexPosition(v2);
    
    if (pos2[1] >= highest - 1e-4 && pos1[1] >= highest - 1e-4)
    {
      rodEdges.push_back(*eit);
      if (pos1[0] <= leftmost + 1e-4 || pos2[0] <= leftmost + 1e-4)
        leftmostedge = *eit;
    }
  }
  
  std::cout << "rod edge count = " << rodEdges.size() << std::endl;
  
  if (leftmostedge.isValid())
  {
    obj->constrainVertex(obj->fromVertex(leftmostedge), obj->getVertexPosition(obj->fromVertex(leftmostedge)));
    obj->constrainVertex(obj->toVertex(leftmostedge),   obj->getVertexPosition(obj->toVertex(leftmostedge)));
  }
  
  // create a rod model
  rod = new ElasticRodModel(obj, rodEdges, m_timestep);
  obj->addModel(rod);
  
  // set init dofs for edges
  EdgeProperty<Scalar> zeros(obj);
  zeros.assign(0);
  rod->setEdgeThetas(zeros);
  rod->setEdgeThetaVelocities(zeros);
  rod->setEdgeUndeformedThetas(zeros);
  
  if (leftmostedge.isValid())
  {
    std::cout << leftmostedge.idx() << std::endl;
    rod->constrainEdgeVel(leftmostedge, 0, 0.1, 0);
  }
  
  m_s12_rod_edges = rodEdges;
}

void RodShellTest::setupScene13() 
{
  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");
  
  Scalar dx = (Scalar)width / (Scalar)xresolution;
  Scalar dy = (Scalar)height / (Scalar)yresolution;
  
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(obj);
  VertexProperty<Vec3d> positions(obj);
  VertexProperty<Vec3d> velocities(obj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(obj);
  EdgeProperty<Scalar> edgeAngle(obj);
  EdgeProperty<Scalar> edgeVel(obj);
  
  Vec3d start_vel(0,0,0);
  
  // pyramids
  int n = xresolution + 1;
  Scalar d = dx * 1.5;
  VertexHandle basep1;
  VertexHandle basep2;
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n - i; ++j)
    {
      for (int k = 0; k < j + 1; ++k)
      {
        Vec3d vert(j * dx - k * dx * 0.5 + i * dx * 0.5 + d, i * dy, -(k * dx * 0.8660254 + i * dx * 0.2886751));
        Vec3d undef = vert;
        
        VertexHandle h = obj->addVertex();
        if (i == 0 && j == 0 && k == 0)
          basep1 = h;
        
        positions[h] = vert;
        velocities[h] = start_vel;
        undeformed[h] = undef;
        vertHandles.push_back(h);
      }
    }
  }
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n - i; ++j)
    {
      for (int k = 0; k < j + 1; ++k)
      {
        Vec3d vert(-(j * dx - k * dx * 0.5 + i * dx * 0.5 + d), i * dy, k * dx * 0.8660254 + i * dx * 0.2886751);
        Vec3d undef = vert;
        
        VertexHandle h = obj->addVertex();
        if (i == 0 && j == 0 && k == 0)
          basep2 = h;
        
        positions[h] = vert;
        velocities[h] = start_vel;
        undeformed[h] = undef;
        vertHandles.push_back(h);
      }
    }
  }
  Vec3d o1(dx * 0.5, 0, 0);
  Vec3d o2(-dx * 0.5, 0, 0);
  VertexHandle origin1 = obj->addVertex();
  VertexHandle origin2 = obj->addVertex();
  positions[origin1] = o1;
  positions[origin2] = o2;
  velocities[origin1] = start_vel;
  velocities[origin2] = start_vel;
  undeformed[origin1] = o1;
  undeformed[origin2] = o2;
  vertHandles.push_back(origin1);
  vertHandles.push_back(origin2);
  obj->addEdge(basep1, origin1);
  obj->addEdge(basep2, origin2);
  obj->addEdge(origin1, origin2);
  
  FaceProperty<char> shellFaces(obj);
  int base_i = 0;
  for (int i = 0; i < n - 1; ++i)
  {
    int stride_i = (n - i) * (n - i + 1) / 2;
    int base_j = 0;
    for (int j = 0; j < n - i - 1; ++j)
    {
      int stride_j = j + 1;
      int base_k = 0;
      for (int k = 0; k < j + 1; ++k)
      {
        int stride_k = 1;
        FaceHandle f1 = obj->addFace(vertHandles[base_i + base_j + base_k], vertHandles[base_i + base_j + base_k + stride_j + stride_k], vertHandles[base_i + base_j + base_k + stride_j]);
        FaceHandle f2 = obj->addFace(vertHandles[base_i + base_j + base_k], vertHandles[base_i + base_j + base_k + stride_j], vertHandles[base_i + base_j + base_k + stride_i]);
        FaceHandle f3 = obj->addFace(vertHandles[base_i + base_j + base_k], vertHandles[base_i + base_j + base_k + stride_i], vertHandles[base_i + base_j + base_k + stride_j + stride_k]);
        FaceHandle f4 = obj->addFace(vertHandles[base_i + base_j + base_k + stride_j], vertHandles[base_i + base_j + base_k + stride_j + stride_k], vertHandles[base_i + base_j + base_k + stride_i]);
        obj->addTet(f1, f2, f3, f4);
        if (i == 0)
        {
          shellFaces[f1] = true;
        }
        base_k += stride_k;
      }
      base_j += stride_j;
    }
    base_i += stride_i;
  }
  base_i++;
  for (int i = 0; i < n - 1; ++i)
  {
    int stride_i = (n - i) * (n - i + 1) / 2;
    int base_j = 0;
    for (int j = 0; j < n - i - 1; ++j)
    {
      int stride_j = j + 1;
      int base_k = 0;
      for (int k = 0; k < j + 1; ++k)
      {
        int stride_k = 1;
        int base = base_i + base_j + base_k;
        FaceHandle f1 = obj->addFace(vertHandles[base], vertHandles[base + stride_j + stride_k], vertHandles[base + stride_j]);
        FaceHandle f2 = obj->addFace(vertHandles[base], vertHandles[base + stride_j], vertHandles[base + stride_i]);
        FaceHandle f3 = obj->addFace(vertHandles[base], vertHandles[base + stride_i], vertHandles[base + stride_j + stride_k]);
        FaceHandle f4 = obj->addFace(vertHandles[base + stride_j], vertHandles[base + stride_j + stride_k], vertHandles[base + stride_i]);
        obj->addTet(f1, f2, f3, f4);
        base_k += stride_k;
      }
      base_j += stride_j;
    }
    base_i += stride_i;
  }

  std::cout << "resolution = " << xresolution << "x" << yresolution << std::endl;
  std::cout << "mesh nv = " << obj->nv() << " ne = " << obj->ne() << " nf = " << obj->nf() << std::endl;
  
  //create a tet property to flag which of the tets are part of the object. (All of them, in this case.)
  TetProperty<char> solidTets(obj); 
  DeformableObject::tet_iter tIt;
  for (tIt = obj->tets_begin(); tIt != obj->tets_end(); ++tIt)
    solidTets[*tIt] = true;
  
  //now create the physical model to hang on the mesh
  solid = new ElasticSolid(obj, solidTets, m_timestep);
  obj->addModel(solid);
  shell = new ElasticShell(obj, shellFaces, m_timestep);
  obj->addModel(shell);
  
  //positions
  solid->setVertexUndeformed(undeformed);
  solid->setVertexPositions(positions);
  solid->setVertexVelocities(velocities);
  
  //Find lowest vertex
  VertexIterator vit = obj->vertices_begin();
  Scalar lowest = 10000;
  for(;vit!= obj->vertices_end(); ++vit) {
    Vec3d pos = solid->getVertexPosition(*vit);
    if(pos[1] <= lowest)
      lowest = pos[1];
  }
  
  //Pin all verts at or near that height, on the left
  //  for(vit = obj->vertices_begin();vit!= obj->vertices_end(); ++vit) {
  //    Vec3d pos = shell->getVertexPosition(*vit);
  //    if(pos[1] >= highest - 1e-4)
  //      obj->constrainVertex(*vit, pos);
  //  }
  
  // find top edges
  std::vector<EdgeHandle> rodEdges;
  EdgeHandle originedge;
  for (EdgeIterator eit = obj->edges_begin(); eit != obj->edges_end(); ++eit)
  {
    EdgeVertexIterator evit = obj->ev_iter(*eit);
    VertexHandle v1 = *evit; ++evit;
    VertexHandle v2 = *evit; ++evit;
    Vec3d pos1 = obj->getVertexPosition(v1);
    Vec3d pos2 = obj->getVertexPosition(v2);
    
    if (pos2[1] <= lowest + 1e-4 && pos1[1] <= lowest + 1e-4 && pos2[2] <= 1e-4 && pos1[2] <= 1e-4 && pos2[2] >= -1e-4 && pos1[2] >= -1e-4)
    {
      rodEdges.push_back(*eit);
      if ((v1 == origin1 && v2 == origin2) || (v1 == origin2 && v2 == origin1))
      {
        originedge = *eit;
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

  obj->constrainVertex(origin1, obj->getVertexPosition(origin1));
  obj->constrainVertex(origin2, obj->getVertexPosition(origin2));
//  obj->constrainVertex(basep1, obj->getVertexPosition(basep1));
//  obj->constrainVertex(basep2, obj->getVertexPosition(basep2));
  
  if (originedge.isValid())
  {
    rod->constrainEdge(originedge, 0);
  }
  
  m_s13_rod_edges = rodEdges;
}

void RodShellTest::setupScene14() 
{
  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar dx = width;
  
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
  VertexHandle h;
  Vec3d vert;
  
  vert = Vec3d(0, 0, 0);                      // vertex 0
  h = obj->addVertex(); positions[h] = vert;  velocities[h] = start_vel;  undeformed[h] = vert;   vertHandles.push_back(h);
  vert = Vec3d(dx, 0, 0);                     // vertex 1
  h = obj->addVertex(); positions[h] = vert;  velocities[h] = start_vel;  undeformed[h] = vert;   vertHandles.push_back(h);
  vert = Vec3d(2 * dx, 0, 0);                 // vertex 2
  h = obj->addVertex(); positions[h] = vert;  velocities[h] = start_vel;  undeformed[h] = vert;   vertHandles.push_back(h);
  vert = Vec3d(3 * dx, dx * 0.3, -dx * 0.5);  // vertex 3
  h = obj->addVertex(); positions[h] = vert;  velocities[h] = start_vel;  undeformed[h] = vert;   vertHandles.push_back(h);
  vert = Vec3d(3 * dx, dx * 0.3, dx * 0.5);   // vertex 4
  h = obj->addVertex(); positions[h] = vert;  velocities[h] = start_vel;  undeformed[h] = vert;   vertHandles.push_back(h);
  vert = Vec3d(4 * dx, dx * 0.3, dx * 0.5);   // vertex 5
  h = obj->addVertex(); positions[h] = vert;  velocities[h] = start_vel;  undeformed[h] = vert;   vertHandles.push_back(h);
  vert = Vec3d(3.5 * dx, dx * 0.3, -dx * 0.5);// vertex 6
  h = obj->addVertex(); positions[h] = vert;  velocities[h] = start_vel;  undeformed[h] = vert;   vertHandles.push_back(h);
  vert = Vec3d(4 * dx, dx * 0.3, -dx * 1.5);  // vertex 7
  h = obj->addVertex(); positions[h] = vert;  velocities[h] = start_vel;  undeformed[h] = vert;   vertHandles.push_back(h);
  vert = Vec3d(3 * dx, dx * 0.3, -dx * 1.5);  // vertex 8
  h = obj->addVertex(); positions[h] = vert;  velocities[h] = start_vel;  undeformed[h] = vert;   vertHandles.push_back(h);
  vert = Vec3d(3.5 * dx, dx * 1.0, -dx * 1.0);// vertex 9
  h = obj->addVertex(); positions[h] = vert;  velocities[h] = start_vel;  undeformed[h] = vert;   vertHandles.push_back(h);
  vert = Vec3d(4 * dx, dx * 1.7, -dx * 1.5);  // vertex 10
  h = obj->addVertex(); positions[h] = vert;  velocities[h] = start_vel;  undeformed[h] = vert;   vertHandles.push_back(h);
  vert = Vec3d(3 * dx, dx * 1.7, -dx * 1.5);  // vertex 11
  h = obj->addVertex(); positions[h] = vert;  velocities[h] = start_vel;  undeformed[h] = vert;   vertHandles.push_back(h);
  vert = Vec3d(3.5 * dx, dx * 1.7, -dx * 0.5);// vertex 12
  h = obj->addVertex(); positions[h] = vert;  velocities[h] = start_vel;  undeformed[h] = vert;   vertHandles.push_back(h);

  std::vector<EdgeHandle> edges;
  std::vector<FaceHandle> faces;
  std::vector<TetHandle> tets;

  edges.push_back(obj->addEdge(vertHandles[0], vertHandles[1]));
  edges.push_back(obj->addEdge(vertHandles[1], vertHandles[2]));
  
  faces.push_back(obj->addFace(vertHandles[2], vertHandles[3], vertHandles[4]));
  faces.push_back(obj->addFace(vertHandles[4], vertHandles[5], vertHandles[6]));
  
  faces.push_back(obj->addFace(vertHandles[6], vertHandles[8], vertHandles[7]));
  faces.push_back(obj->addFace(vertHandles[6], vertHandles[9], vertHandles[8]));
  faces.push_back(obj->addFace(vertHandles[6], vertHandles[7], vertHandles[9]));
  faces.push_back(obj->addFace(vertHandles[7], vertHandles[8], vertHandles[9]));
  tets.push_back(obj->addTet(faces[2], faces[3], faces[4], faces[5]));
  faces.push_back(obj->addFace(vertHandles[9], vertHandles[10], vertHandles[11]));
  faces.push_back(obj->addFace(vertHandles[9], vertHandles[11], vertHandles[12]));
  faces.push_back(obj->addFace(vertHandles[9], vertHandles[12], vertHandles[10]));
  faces.push_back(obj->addFace(vertHandles[10], vertHandles[12], vertHandles[11]));
  tets.push_back(obj->addTet(faces[6], faces[7], faces[8], faces[9]));

  FaceProperty<char> shellFaces(obj);
  shellFaces[faces[0]] = true;
  shellFaces[faces[1]] = true;
  TetProperty<char> solidTets(obj);
  solidTets[tets[0]] = true;
  solidTets[tets[1]] = true;
  
  // create a solid model
  solid = new ElasticSolid(obj, solidTets, m_timestep);
  obj->addModel(solid);
  
  //positions
  solid->setVertexUndeformed(undeformed);
  solid->setVertexPositions(positions);
  solid->setVertexVelocities(velocities);

  // create a shell model
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
    
  // create a rod model
  rod = new ElasticRodModel(obj, edges, m_timestep);
  obj->addModel(rod);
  
  // set init dofs for edges
  EdgeProperty<Scalar> zeros(obj);
  zeros.assign(0);
  rod->setEdgeThetas(zeros);
  rod->setEdgeThetaVelocities(zeros);
  rod->setEdgeUndeformedThetas(zeros);

  obj->constrainVertex(vertHandles[0], obj->getVertexPosition(vertHandles[0]));
  obj->constrainVertex(vertHandles[1], obj->getVertexPosition(vertHandles[1]));
  rod->constrainEdge(edges[0], 0);

  m_s14_vertices = vertHandles;
  m_s14_edges = edges;
  m_s14_faces = faces;
  m_s14_tets = tets;
}

