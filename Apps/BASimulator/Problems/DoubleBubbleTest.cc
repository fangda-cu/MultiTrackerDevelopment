/**
 * \file DoubleBubbleTest.cc
 *
 * \author fang@cs.columbia.edu
 * \date Nov 20, 2012
 */

#include "DoubleBubbleTest.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShell.hh"
#include "BASim/src/Physics/DeformableObjects/DefoObjTimeStepper.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/CSTMembraneForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/DSBendingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/MNBendingForce.hh"
#include "BASim/src/Render/ShellRenderer.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellRadialForce.hh"
#include "BASim/src/Physics/DeformableObjects/GravityForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellSurfaceTensionForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellLinearSurfaceTensionForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellLinearSurfaceTensionForce2.hh"

#include "BASim/src/Physics/DeformableObjects/Shells/ShellVolumeForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/DrainingBubblePressureForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellBathForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellVerticalForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellPointForce.hh"
#include "BASim/src/Collisions/ElTopo/util.hh"

//An ElTopo global variable.
#include "runstats.h"

#include <set>
#include <fstream>

//For OBJDUMP
#include <sstream>
#include "BASim/src/IO/ObjWriter.hh"
#include "BASim/src/IO/PlyWriter.hh"
#include <iomanip>
#ifdef _MSC_VER
#include <direct.h>
#endif
#include <sys/stat.h>
#include <sys/types.h>

// Toggle switch to output obj files
extern bool g_obj_dump;
int db_current_obj_frame = 0;
extern bool g_ply_dump;
int db_current_ply_frame = 0;

extern std::string outputdirectory;

DoubleBubbleTest::DoubleBubbleTest() : 
  Problem("Double Bubble Test", "Various bubble dynamics tests"), 
  shell(NULL), 
  shellObj(NULL), 
  stepper(NULL)
{
  addDynamicsProps();
  
  //Choice of scene
  AddOption("shell-scene", "the shell scene to test", 1);

  //Basic shell options
  AddOption("shell-thickness", "the (initial) thickness of the shell", 0.01);
  AddOption("shell-density", "volumetric density of the shell ", 1.0);
  
  AddOption("shell-update-thickness", "whether the shell thickness should be dynamically updated to maintain volume", true);
  

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
  AddOption("shell-remeshing", "whether to perform remeshing", false);
  AddOption("shell-remeshing-resolution", "target edge-length", 0.0); //for backwards compatibility
  AddOption("shell-remeshing-max-length", "upper bound on edge-length", 0.5);
  AddOption("shell-remeshing-min-length", "lower bound on edge-length", 0.1);
  AddOption("shell-remeshing-iterations", "number of remeshing iterations to run", 2);

  //Area-based surface tension force
  AddOption("shell-surface-tension", "surface tension coefficient of the shell", 0.0);
  
  //Properties for thickness-dependent linear elasticity & viscosity
  AddOption("shell-Poisson", "the Poisson ratio of the shell material", 0.0f);
  AddOption("shell-Youngs", "the Young's modulus of the shell material", 0.0f);
  AddOption("shell-Poisson-damping", "the damping coefficient associated to the shell's Poisson ratio", 0.0f);
  AddOption("shell-Youngs-damping", "the damping coefficient associated with the shell's Young's modulus", 0.0f);
  AddOption("shell-uniform-load", "a vertical force per unit area applied to the shell surface", 0.0);

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

  AddOption("shell-eltopo-collisions", "whether to apply bridson/harmon-style CCD collision handling, via the El Topo library", true);

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

  //Time stepper options
  AddOption("integrator", "type of integrator to use for the shell", "implicit");

  //Solver options
  AddOption("iterations", "maximum number of iterations for the implicit method", (int) 100);
  AddOption("atol", "absolute convergence tolerance", 1e-8);
  AddOption("rtol", "relative convergence tolerance", 1e-8);
  AddOption("stol", "convergence tolerance in terms of the norm of the change in the solution between steps", 1e-8);
  AddOption("inftol", "infinity norm convergence tolerance", 1e-8);
  
  //OBJ file dump
  AddOption("generate-OBJ", "Generate an OBJ file at each timestep", g_obj_dump);
  AddOption("generate-PLY", "Generate a PLY file at each timestep", g_ply_dump);


}

DoubleBubbleTest::~DoubleBubbleTest()
{
  if (shellObj != NULL) delete shellObj;
  //the shell is deleted by the obj destructor
  if (stepper != NULL) delete stepper;
}

typedef void (DoubleBubbleTest::*sceneFunc)();


sceneFunc db_scenes[] = 
{
  0,
  &DoubleBubbleTest::setupScene1,
  &DoubleBubbleTest::setupScene2,
  &DoubleBubbleTest::setupScene3,
};

Scalar db_bubbleThicknessFunction(Vec3d pos) {
  //NOTE: Assumes radius of 0.01;
  Scalar rad = 0.01f;
  Scalar height_frac = pos[1] / 0.01f;
  Scalar transition_point = 1.0f;
  if(height_frac > transition_point) {
    return 0.0001;
  }
  else {
    Scalar s = height_frac/transition_point;
    return s*0.0001 + (1-s)*0.00015;
  }

}

void DoubleBubbleTest::Setup()
{

  loadDynamicsProps();

  //General shell forces and properties
  Scalar density = GetScalarOpt("shell-density");
  Scalar thickness = GetScalarOpt("shell-thickness");
  m_initial_thickness = thickness;

  Vec3d gravity = GetVecOpt("gravity");
  
  Scalar surface_tension = GetScalarOpt("shell-surface-tension");

  Scalar Youngs_modulus = GetScalarOpt("shell-Youngs");
  Scalar Poisson_ratio = GetScalarOpt("shell-Poisson");
  Scalar Youngs_damping = GetScalarOpt("shell-Youngs-damping");
  Scalar Poisson_damping = GetScalarOpt("shell-Poisson-damping");
  
  bool cst_stretch = GetBoolOpt("shell-CST-stretching");
  bool ds_bend = GetBoolOpt("shell-DS-bending");
  bool mn_bend = GetBoolOpt("shell-MN-bending");

  //fudge factors to modify the elastic-viscous coefficients (so as to manipulate them separately)
  Scalar cst_scale = GetScalarOpt("shell-stretching-factor");
  Scalar ds_scale = GetScalarOpt("shell-bending-factor");
  
  std::string integrator = GetStringOpt("integrator");

  Scalar timestep = getDt(); //Our Rayleigh damping model relies on knowing the timestep (folds it into the damping stiffness, as in Viscous Threads)
  m_timestep = timestep;

  //Geometry/scene specific
  int sceneChoice = GetIntOpt("shell-scene");
  m_active_scene = sceneChoice;

  //Create the base deformable object (mesh)
  shellObj = new DeformableObject();

  //Call the appropriate scene setup function.
  (*this.*db_scenes[sceneChoice])();

  
  shell->setThickness(thickness);

  //now add forces to the model
  MNBendingForce* bender;

  //Stretching and bending forces
  if(Youngs_modulus != 0 || Youngs_damping != 0) {
    
    //Stretching force (Constant Strain Triangle, i.e. linear FEM)
    if(cst_stretch)
      shell->addForce(new CSTMembraneForce(*shell, "CSTMembrane", cst_scale*Youngs_modulus, Poisson_ratio, cst_scale*Youngs_damping, Poisson_damping, timestep));
    
    //Bending force (Hinge-based Bending, a la Discrete Shells)
    if(ds_bend)
      shell->addForce(new DSBendingForce(*shell, "DSBending", ds_scale*Youngs_modulus, Poisson_ratio, ds_scale*Youngs_damping, Poisson_damping, timestep));

    //Better bending model (Mid-Edge normals, a la Computing discrete shape operators on general meshes)
    if(mn_bend) {
      bender = new MNBendingForce(*shell, "MNBending", ds_scale*Youngs_modulus, Poisson_ratio, ds_scale*Youngs_damping, Poisson_damping, timestep);
      shell->addForce(bender);
    }

  }

  //Gravity force
  shellObj->addForce(new GravityForce(*shellObj, timestep, "Gravity", gravity)); 

  Scalar loadForce = GetScalarOpt("shell-uniform-load");
  if(loadForce != 0)
    shell->addForce(new ShellVerticalForce(*shell, "Load", Vec3d(0,-1,0), loadForce));

  
  //Surface tension force
  if(surface_tension != 0) {
    //Viscous sheets-style surface tension
    shell->addForce(new ShellSurfaceTensionForce(*shell, "Surface Tension", surface_tension));
    
    //Experimental piecewise linear surface tension
    //shell->addForce(new ShellLinearSurfaceTensionForce(*shell, "Surface Tension", surface_tension));
    //shell->addForce(new ShellLinearSurfaceTensionForce2(*shell, "Surface Tension", surface_tension));
  }

  shell->setDensity(density);
  
  bool remeshing = GetBoolOpt("shell-remeshing");
  bool self = GetBoolOpt("shell-self-collision");
  std::cout << "Remeshing: " << remeshing << std::endl;
  Scalar remeshing_rez = GetScalarOpt("shell-remeshing-resolution");
  Scalar remeshing_min, remeshing_max;
  if(remeshing_rez == 0.0) {
    remeshing_min = GetScalarOpt("shell-remeshing-min-length");
    remeshing_max = GetScalarOpt("shell-remeshing-max-length");
  }
  else {
      remeshing_min = 0.5*remeshing_rez;
      remeshing_max = 1.5*remeshing_rez;
  }
  int remeshing_its = GetIntOpt("shell-remeshing-iterations");
  shell->setRemeshing(remeshing, remeshing_min, remeshing_max, remeshing_its);
  
  bool eltopo_collisions = GetBoolOpt("shell-eltopo-collisions");
  shell->setElTopoCollisions(eltopo_collisions);

  bool thickness_evolution = GetBoolOpt("shell-update-thickness");
  shell->setThicknessUpdating(thickness_evolution);

  //shell->remesh(remeshing_res);

  shell->computeMasses();
 

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
    
  //compute the dof indexing for use in the diff_eq solver
  shellObj->computeDofIndexing();

  stepper = new DefoObjTimeStepper(*shellObj);
  if(integrator == "symplectic")
    stepper->setDiffEqSolver(DefoObjTimeStepper::SYMPL_EULER);
  else if(integrator == "implicit")
    stepper->setDiffEqSolver(DefoObjTimeStepper::IMPL_EULER);
  else if(integrator == "statics")
    stepper->setDiffEqSolver(DefoObjTimeStepper::STATICS);
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
  
  m_world->addObject(shellObj);
  m_world->addController(stepper);
  RenderBase* shellRender = new ShellRenderer(*shell, m_initial_thickness);
  m_world->addRenderer(shellRender);
  
  g_obj_dump = GetBoolOpt("generate-OBJ");
  g_ply_dump = GetBoolOpt("generate-PLY");

  if (g_obj_dump || g_ply_dump){
#ifdef _MSC_VER
    _mkdir(outputdirectory.c_str());
#else
    mkdir(outputdirectory.c_str(), 0755);
#endif
  }

}

void DoubleBubbleTest::AtEachTimestep()
{
  
    
  //dump PLY files if needed
    if ( g_ply_dump ){
        std::stringstream name;
        int file_width = 20;

        name << std::setfill('0');
        name << outputdirectory << "/frame" << std::setw(file_width) << db_current_ply_frame << ".PLY";

        PlyWriter::write(name.str(), *shell);
        std::cout << "Frame: " << db_current_ply_frame << "   Time: " << getTime() << "   PLYDump: "
                            << name.str() << std::endl;

        ++db_current_ply_frame;
    }
    //Dump OBJ files if needed
    //Start OBJ file stuff
    if ( g_obj_dump ){

        std::stringstream name;
        int file_width = 20;

        name << std::setfill('0');
        name << outputdirectory << "/frame" << std::setw(file_width) << db_current_obj_frame << ".OBJ";

        ObjWriter::write(name.str(), *shell);
        std::cout << "Frame: " << db_current_obj_frame << "   Time: " << getTime() << "   OBJDump: "
                            << name.str() << std::endl;

        ++db_current_obj_frame;
    }

    std::cout << "Time: " << this->getTime() << std::endl; 

}

void DoubleBubbleTest::setupScene1() {
  
  Vec3d start_vel(0,0,0);
  
  //vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);
  
  //create a cube
  std::vector<VertexHandle> vertList;
  
  for(int i = 0; i < 12; ++i) {
    vertList.push_back(shellObj->addVertex());
    velocities[vertList[i]] = start_vel;
  }
  
  //create positions
  positions[vertList[0]] = Vec3d(0,0,0);
  positions[vertList[1]] = Vec3d(0,0,1);
  positions[vertList[2]] = Vec3d(0,1,0);
  positions[vertList[3]] = Vec3d(0,1,1);
  
  positions[vertList[4]] = Vec3d(1,0,0);
  positions[vertList[5]] = Vec3d(1,0,1);
  positions[vertList[6]] = Vec3d(1,1,0);
  positions[vertList[7]] = Vec3d(1,1,1);
  
  positions[vertList[8]] = Vec3d(2,0,0);
  positions[vertList[9]] = Vec3d(2,0,1);
  positions[vertList[10]] = Vec3d(2,1,0);
  positions[vertList[11]] = Vec3d(2,1,1);
  
  
  for(int i = 0; i < 8; ++i) {
    undeformed[vertList[i]] = positions[vertList[i]];
  }
  
  std::vector<FaceHandle> faceList;
  
  //first cube
  faceList.push_back(shellObj->addFace(vertList[0], vertList[2], vertList[4]));
  faceList.push_back(shellObj->addFace(vertList[6], vertList[4], vertList[2]));
  
  faceList.push_back(shellObj->addFace(vertList[7], vertList[6], vertList[3]));
  faceList.push_back(shellObj->addFace(vertList[3], vertList[6], vertList[2]));
  
  faceList.push_back(shellObj->addFace(vertList[3], vertList[2], vertList[1]));
  faceList.push_back(shellObj->addFace(vertList[1], vertList[2], vertList[0]));
  
  faceList.push_back(shellObj->addFace(vertList[0], vertList[4], vertList[5]));
  faceList.push_back(shellObj->addFace(vertList[0], vertList[5], vertList[1]));
  
  faceList.push_back(shellObj->addFace(vertList[7], vertList[3], vertList[1]));
  faceList.push_back(shellObj->addFace(vertList[5], vertList[7], vertList[1]));
  
  //shared
  faceList.push_back(shellObj->addFace(vertList[4], vertList[6], vertList[5]));
  faceList.push_back(shellObj->addFace(vertList[5], vertList[6], vertList[7]));
  
  //second cube
  faceList.push_back(shellObj->addFace(vertList[4], vertList[6], vertList[10]));
  faceList.push_back(shellObj->addFace(vertList[4], vertList[10], vertList[8]));
  
  faceList.push_back(shellObj->addFace(vertList[6], vertList[7], vertList[11]));
  faceList.push_back(shellObj->addFace(vertList[6], vertList[11], vertList[10]));
  
  faceList.push_back(shellObj->addFace(vertList[8], vertList[10], vertList[9]));
  faceList.push_back(shellObj->addFace(vertList[10], vertList[11], vertList[9]));
  
  faceList.push_back(shellObj->addFace(vertList[4], vertList[8], vertList[5]));
  faceList.push_back(shellObj->addFace(vertList[5], vertList[8], vertList[9]));
  
  faceList.push_back(shellObj->addFace(vertList[9], vertList[11], vertList[7]));
  faceList.push_back(shellObj->addFace(vertList[9], vertList[7], vertList[5]));
  
  FaceProperty<Vec2i> faceLabels(shellObj); //label face regions to do volume constrained bubbles
  
  //cube 1
  for(int i = 0; i < 10; ++i)
    faceLabels[faceList[i]] = Vec2i(0,-1);
  
  //cube 2
  for(int i = 12; i < 22; ++i)
    faceLabels[faceList[i]] = Vec2i(1, -1);
  
  //connector faces
  faceLabels[faceList[10]] = Vec2i(0, 1);
  faceLabels[faceList[11]] = Vec2i(0, 1);
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(shellObj); 
  DeformableObject::face_iter fIt;
  for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(shellObj, shellFaces, m_timestep);
  shellObj->addModel(shell);
  
  //positions
  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  //mid-edge normal variables
  shell->setEdgeUndeformed(undefAngle);
  shell->setEdgeXis(edgeAngle);
  shell->setEdgeVelocities(edgeVel);
  
  shell->setFaceLabels(faceLabels);
  
  shell->addForce(new ShellVolumeForce(*shell, "Volume", 1000));
  
}

void DoubleBubbleTest::setupScene2() 
{
  Vec3d start_vel(0,0,0);
  
  //vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);
  
  //create a cube
  std::vector<VertexHandle> vertList;
  
  for(int i = 0; i < 8; ++i) 
  {
    vertList.push_back(shellObj->addVertex());
    velocities[vertList[i]] = start_vel;
  }
  
  //create positions
  positions[vertList[0]] = Vec3d(0,0,0);
  positions[vertList[1]] = Vec3d(0,0,1);
  positions[vertList[2]] = Vec3d(0,1,0);
  positions[vertList[3]] = Vec3d(0,1,1);
  
  positions[vertList[4]] = Vec3d(1,0,0);
  positions[vertList[5]] = Vec3d(1,0,1);
  positions[vertList[6]] = Vec3d(1,1,0);
  positions[vertList[7]] = Vec3d(1,1,1);
  
  for(int i = 0; i < 8; ++i)
    undeformed[vertList[i]] = positions[vertList[i]];
  
  std::vector<FaceHandle> faceList;
  
  faceList.push_back(shellObj->addFace(vertList[0], vertList[2], vertList[4]));
  faceList.push_back(shellObj->addFace(vertList[6], vertList[4], vertList[2]));
  
  faceList.push_back(shellObj->addFace(vertList[7], vertList[6], vertList[3]));
  faceList.push_back(shellObj->addFace(vertList[3], vertList[6], vertList[2]));
  
  faceList.push_back(shellObj->addFace(vertList[3], vertList[2], vertList[1]));
  faceList.push_back(shellObj->addFace(vertList[1], vertList[2], vertList[0]));
  
  faceList.push_back(shellObj->addFace(vertList[0], vertList[4], vertList[5]));
  faceList.push_back(shellObj->addFace(vertList[0], vertList[5], vertList[1]));
  
  faceList.push_back(shellObj->addFace(vertList[7], vertList[3], vertList[1]));
  faceList.push_back(shellObj->addFace(vertList[5], vertList[7], vertList[1]));
  
  faceList.push_back(shellObj->addFace(vertList[4], vertList[6], vertList[5]));
  faceList.push_back(shellObj->addFace(vertList[5], vertList[6], vertList[7]));
  
  FaceProperty<Vec2i> faceLabels(shellObj); //label face regions to do volume constrained bubbles  
  for(int i = 0; i < 12; ++i)
    faceLabels[faceList[i]] = Vec2i(0,-1);
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(shellObj); 
  DeformableObject::face_iter fIt;
  for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(shellObj, shellFaces, m_timestep);
  shellObj->addModel(shell);
  
  //positions
  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  //mid-edge normal variables
  shell->setEdgeUndeformed(undefAngle);
  shell->setEdgeXis(edgeAngle);
  shell->setEdgeVelocities(edgeVel);
  
  shell->setFaceLabels(faceLabels);
  
  shell->addForce(new ShellVolumeForce(*shell, "Volume", 1000));
  
}

void DoubleBubbleTest::setupScene3() 
{
  Vec3d start_vel(0,0,0);
  
  //vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);
  
  //create a cube
  std::vector<VertexHandle> vertList;
  
  for(int i = 0; i < 8; ++i) 
  {
    vertList.push_back(shellObj->addVertex());
    velocities[vertList[i]] = start_vel;
  }
  
  //create positions
  positions[vertList[0]] = Vec3d(0,0,0);
  positions[vertList[1]] = Vec3d(0,0,1);
  positions[vertList[2]] = Vec3d(0,1,0);
  positions[vertList[3]] = Vec3d(0,1,1);
  
  positions[vertList[4]] = Vec3d(1,0,0);
  positions[vertList[5]] = Vec3d(1,0,1);
  positions[vertList[6]] = Vec3d(1,1,0);
  positions[vertList[7]] = Vec3d(1,1,1);
  
  for(int i = 0; i < 8; ++i)
    undeformed[vertList[i]] = positions[vertList[i]];
  
  std::vector<FaceHandle> faceList;
  
  faceList.push_back(shellObj->addFace(vertList[0], vertList[2], vertList[4]));
  faceList.push_back(shellObj->addFace(vertList[6], vertList[4], vertList[2]));
  
  faceList.push_back(shellObj->addFace(vertList[7], vertList[6], vertList[3]));
  faceList.push_back(shellObj->addFace(vertList[3], vertList[6], vertList[2]));
  
  faceList.push_back(shellObj->addFace(vertList[3], vertList[2], vertList[1]));
  faceList.push_back(shellObj->addFace(vertList[1], vertList[2], vertList[0]));
  
  faceList.push_back(shellObj->addFace(vertList[0], vertList[4], vertList[5]));
  faceList.push_back(shellObj->addFace(vertList[0], vertList[5], vertList[1]));
  
  faceList.push_back(shellObj->addFace(vertList[7], vertList[3], vertList[1]));
  faceList.push_back(shellObj->addFace(vertList[5], vertList[7], vertList[1]));
  
  faceList.push_back(shellObj->addFace(vertList[4], vertList[6], vertList[5]));
  faceList.push_back(shellObj->addFace(vertList[5], vertList[6], vertList[7]));
  
  FaceProperty<Vec2i> faceLabels(shellObj); //label face regions to do volume constrained bubbles  
  for(int i = 0; i < 12; ++i)
    faceLabels[faceList[i]] = Vec2i(0,-1);
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(shellObj); 
  DeformableObject::face_iter fIt;
  for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(shellObj, shellFaces, m_timestep);
  shellObj->addModel(shell);
  
  //positions
  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  //mid-edge normal variables
  shell->setEdgeUndeformed(undefAngle);
  shell->setEdgeXis(edgeAngle);
  shell->setEdgeVelocities(edgeVel);
  
  shell->setFaceLabels(faceLabels);
  
  shell->addForce(new ShellVolumeForce(*shell, "Volume", 1000));
  
}
