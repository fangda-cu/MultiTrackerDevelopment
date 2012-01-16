/**
 * \file ShellTest.cc
 *
 * \author batty@cs.columbia.edu
 * \date April 18, 2011
 */

#include "ShellTest.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShell.hh"
#include "BASim/src/Physics/DeformableObjects/DefoObjTimeStepper.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/CSTMembraneForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/DSBendingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/MNBendingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellGravityForce.hh"
#include "BASim/src/Render/ShellRenderer.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellRadialForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellSurfaceTensionForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellVolumeForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellVertexTriSpringForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellVertexPointSpringForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/DrainingBubblePressureForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellBathForce.hh"
#include "BASim/src/Collisions/ElTopo/util.hh"

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
int current_obj_frame = 0;
extern bool g_ply_dump;
int current_ply_frame = 0;

extern std::string outputdirectory;

ShellTest::ShellTest()
: Problem("Shell Test", "Various viscous and elastic sheet/shell tests"), 
  shell(NULL), shellObj(NULL), stepper(NULL)
{
  addDynamicsProps();
  
  //Choice of scene
  AddOption("shell-scene", "the shell scene to test", 1);

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
  
  //Tearing options
  AddOption("shell-tearing", "wheter to add tearing to the model", false);
  AddOption("shell-tearing-threshold", "the thickness threshold to use for tearing", 0.0 );

    //Timestepper options
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

ShellTest::~ShellTest()
{
  if (shellObj != NULL) delete shellObj;
  if (shell != NULL) delete shell;
  if (stepper != NULL) delete stepper;
}

typedef void (ShellTest::*sceneFunc)();


sceneFunc scenes[] = {0,
                      &ShellTest::setupScene1, //vertical flat sheet
                      &ShellTest::setupScene2, //vertical cylindrical sheet
                      &ShellTest::setupScene3, //spherical sheet
                      &ShellTest::setupScene4, //two-triangle bending test
                      &ShellTest::setupScene5, //catenary 
                      &ShellTest::setupScene6, //hemispherical bubble
                      &ShellTest::setupScene7, //sheet sheared between two circles
                      &ShellTest::setupScene8, //torus
                      &ShellTest::setupScene9, //non-manifold edge / bending
                      &ShellTest::setupScene10, //pouring inflow with deletion 
                      &ShellTest::setupScene11, //a cube with surface tension collapsing to a sphere  
                      &ShellTest::setupScene12, //hemispherical bubble popping with low viscosity
                      &ShellTest::setupScene13 };  //an constant inflow hitting a solid floor and buckling

Scalar bubbleThicknessFunction(Vec3d pos) {
  //NOTE: Assumes radius of 0.01;
  float rad = 0.01f;
  float height_frac = pos[1] / 0.01f;
  float transition_point = 1.0f;
  if(height_frac > transition_point) {
    return 0.0001;
  }
  else {
    float s = height_frac/transition_point;
    return s*0.0001 + (1-s)*0.00015;
  }

}

void ShellTest::Setup()
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

  //fudge factors to modify the elastic-viscous coefficients (so as to manipulate them separately)
  Scalar cst_scale = GetScalarOpt("shell-stretching-factor");
  Scalar ds_scale = GetScalarOpt("shell-bending-factor");
  
  std::string integrator = GetStringOpt("integrator");

  Scalar timestep = getDt(); //Our Rayleigh damping model relies on knowing the timestep (folds it into the damping stiffness, as in Viscous Threads)
  m_timestep = timestep;

  //Geometry/scene specific
  int sceneChoice = GetIntOpt("shell-scene");
  m_active_scene = sceneChoice;

  if(sceneChoice == 3) {
    Scalar pressureStrength = GetScalarOpt("shell-inflate-sphere-coeff");

    sphere_data.open("sphere_data.txt");
    sphere_data << "Input parameter data..." << std::endl;
    sphere_data << "Viscosity: " << Youngs_damping/3 << std::endl;
    sphere_data << "Initial thickness: " << thickness << std::endl;
    sphere_data << "Pressure coefficient, k: " << pressureStrength << std::endl;
    sphere_data << "Density: " << density << std::endl;
    sphere_data << std::endl;
    sphere_data << "Simulation data, in comma-separated format... radius,velocity and thickness are area-weighted averages over all triangles" << std::endl;
    sphere_data << "time,radius,computed pressure,surface area,current thickness,velocity" << std::endl;

  }

  //Create the base deformable object (mesh)
  shellObj = new DeformableObject();

  //Call the appropriate scene setup function.
  (*this.*scenes[sceneChoice])();


  //now add forces to the model

  //Stretching and bending forces
  if(Youngs_modulus != 0 || Youngs_damping != 0) {
    
    //Stretching force (Constant Strain Triangle, i.e. linear FEM)
    if(cst_stretch)
      shell->addForce(new CSTMembraneForce(*shell, "CSTMembrane", cst_scale*Youngs_modulus, Poisson_ratio, cst_scale*Youngs_damping, Poisson_damping, timestep));
    
    //Bending force (Hinge-based Bending, a la Discrete Shells)
    if(ds_bend)
      shell->addForce(new DSBendingForce(*shell, "DSBending", ds_scale*Youngs_modulus, Poisson_ratio, ds_scale*Youngs_damping, Poisson_damping, timestep));

    //Better bending model, not currently functional.
    //shell->addForce(new MNBendingForce(*shell, "MNBending", Youngs_modulus, Poisson_ratio, Youngs_damping, Poisson_damping, timestep));
  }


  //Gravity force
  shell->addForce(new ShellGravityForce(*shell, "Gravity", gravity));

  //Surface tension force
  if(surface_tension != 0)
    shell->addForce(new ShellSurfaceTensionForce(*shell, "Surface Tension", surface_tension));

  //and set its standard properties
  if(sceneChoice == 6) {
    //make the thickness vary for the bubble example
    for(FaceIterator fit = shellObj->faces_begin(); fit != shellObj->faces_end(); ++fit) {
      FaceHandle fh = *fit;
      Vec3d barycentre(0,0,0);
      for(FaceVertexIterator fvit = shellObj->fv_iter(fh); fvit; ++fvit) {
        VertexHandle vh = *fvit;
        Vec3d position = shell->getVertexPosition(vh);
        barycentre += position;
      }
      barycentre/=3.0;
      shell->setThickness(fh, bubbleThicknessFunction(barycentre));
    }
  }
  else {
    shell->setThickness(thickness);
  }
  shell->setDensity(density);
  
  bool remeshing = GetIntOpt("shell-remeshing") == 1?true:false;
  Scalar remeshing_res = GetScalarOpt("shell-remeshing-resolution");
  int remeshing_its = GetIntOpt("shell-remeshing-iterations");
  shell->setRemeshing(remeshing, remeshing_res, remeshing_its);
  
  //shell->remesh(remeshing_res);

  shell->computeMasses();

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
  shell->setTearing(tearing, tearingThres);

  //for(int i = 0; i < 3; ++i)
  //  shell->remesh(remeshing_res);
    
  //compute the dof indexing for use in the diff_eq solver
  shellObj->computeDofIndexing();

  stepper = new DefoObjTimeStepper(*shellObj);
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

void ShellTest::AtEachTimestep()
{
  
  std::cout << "Time: " << std::endl;
    
  //dump PLY files if needed
    if ( g_ply_dump ){
        std::stringstream name;
        int file_width = 20;

        name << std::setfill('0');
        name << outputdirectory << "/frame" << std::setw(file_width) << current_ply_frame << ".PLY";

        PlyWriter::write(name.str(), *shell);
        std::cout << "Frame: " << current_ply_frame << "   Time: " << getTime() << "   PLYDump: "
                            << name.str() << std::endl;

        ++current_ply_frame;
    }
    //Dump OBJ files if needed
    //Start OBJ file stuff
    if ( g_obj_dump ){

        std::stringstream name;
        int file_width = 20;

        name << std::setfill('0');
        name << outputdirectory << "/frame" << std::setw(file_width) << current_obj_frame << ".OBJ";

        ObjWriter::write(name.str(), *shell);
        std::cout << "Frame: " << current_obj_frame << "   Time: " << getTime() << "   OBJDump: "
                            << name.str() << std::endl;

        ++current_obj_frame;
    }

    if(m_active_scene) {
      //Write out the data for the analytical example

      Scalar velocity = 0;
      int count = 0;
      Scalar rad = 0;
      Scalar area = 0;
      Scalar thickness = 0;
      for(FaceIterator fit = shellObj->faces_begin(); fit!= shellObj->faces_end(); ++fit) {
        FaceHandle fh = *fit;
        Scalar faceArea = shell->getArea(fh);
        Scalar faceRad = 0;
        Scalar faceVel = 0;
        for(FaceVertexIterator fvit = shellObj->fv_iter(fh); fvit; ++fvit) {
          VertexHandle vh = *fvit;
          faceVel += shell->getVertexVelocity(vh).norm();
          faceRad += shell->getVertexPosition(vh).norm();
        }
        thickness += shell->getThickness(fh)*faceArea;
        velocity += faceVel/3*faceArea;
        rad += faceRad/3*faceArea;
        area += faceArea;
      }
      rad /= area;
      velocity /= area;
      thickness /= area;
      std::cout << "Average velocity: " << velocity << std::endl;
      std::cout << "Average radius: " << rad << std::endl;
      std::cout << "Average thickness: " << thickness << std::endl;
      Scalar viscosity = GetScalarOpt("shell-Youngs-damping") / 3.0;
      Scalar thickness0 = GetScalarOpt("shell-thickness");
      Scalar pressureStrength = GetScalarOpt("shell-inflate-sphere-coeff");
      bool constPressure = GetBoolOpt("shell-inflate-sphere-const-pressure");

      Scalar radMul = constPressure? rad*rad : 1;
      double expected_vel = pressureStrength * radMul /(12.0*viscosity*thickness0);
      Scalar pressure_val = pressureStrength / (constPressure? 1 : rad*rad);
      //time,radius,computed pressure,surface area,current thickness,velocity
      sphere_data << getTime() << "," << rad << "," << pressure_val << "," << area << "," << thickness << "," << velocity << std::endl;
      sphere_data.flush();

      std::cout << "Velocity: " << velocity << std::endl;
      std::cout << "Expected velocity: " << expected_vel << std::endl;
    }
    std::cout << "Time: " << this->getTime() << std::endl; 

}

//vertical flat sheet
void ShellTest::setupScene1() {

  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");

  Scalar dx = (Scalar)width / (Scalar)xresolution;
  Scalar dy = (Scalar)height / (Scalar)yresolution;

  //build a rectangular grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);

  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);

  Vec3d start_vel(0,0,0);
  
  for(int j = 0; j <= yresolution; ++j) {
    for(int i = 0; i <= xresolution; ++i) {
      Vec3d vert(i*dx, j*dy, 0.01*dx*sin(100*j*dy + 17*i*dx));
     /* if(j < 0.5*yresolution) {
        int k = j;
        int j_mod = (int)(0.5*yresolution);
        vert(1) = j_mod*dx;
        vert(2) = (k-j_mod)*dx;
      }*/
      Vec3d undef = vert;

      VertexHandle h = shellObj->addVertex();

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
        shellObj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[br]);
        shellObj->addFace(vertHandles[tl], vertHandles[br], vertHandles[bl]);
      }
      else {
        shellObj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[bl]);
        shellObj->addFace(vertHandles[bl], vertHandles[tr], vertHandles[br]);
      }
    }
  }

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



  ////find the top left and right corners for adding a constraint
  //VertexIterator vit = shellObj->vertices_begin();
  //Vec3d maxPos(-100,-100,-100), minPos(100,100,100);
  //VertexHandle minH, maxH;
  //for(;vit!= shellObj->vertices_end(); ++vit) {
  //  Vec3d pos = shell->getVertexPosition(*vit);
  //  if(pos[0] <= minPos[0]) {
  //    minH = *vit;
  //    minPos = pos;
  //  }
  //  if(pos[0] >= maxPos[0]) {
  //    maxH = *vit;
  //    maxPos = pos;
  //  }
  //}
  //shell->constrainVertex(minH, minPos);
  //shell->constrainVertex(maxH, maxPos);


  //Find highest vertex
  VertexIterator vit = shellObj->vertices_begin();
  Scalar highest = -10000;
  for(;vit!= shellObj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest) {
      highest = pos[1];
    }
  }
  //Pin all verts at or near that height
  for(vit = shellObj->vertices_begin();vit!= shellObj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest - 1e-4)
      shell->constrainVertex(*vit, pos);
  }

}

//vertical cylinder, pinned or flowing at top
void ShellTest::setupScene2() {

  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");

  //Scalar dx = (Scalar)width / (Scalar)xresolution;
  Scalar dy = (Scalar)height / (Scalar)yresolution;

  //build a rectangular grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);

  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);
  
  Vec3d start_vel(0,0,0);
  for(int j = 0; j <= yresolution; ++j) {
    for(int i = 0; i <= xresolution; ++i) {
      Scalar circumference = width;
      Scalar radius = circumference / 2 / pi;
      Scalar angle = ((Scalar)i / (Scalar)(xresolution+1)) * 2 * pi;
      Scalar xpos = radius * cos(angle);
      Scalar zpos = radius * sin(angle);
      Vec3d vert(xpos, -1 + j*dy, zpos);

      Vec3d undef = vert;

      VertexHandle h = shellObj->addVertex();

      positions[h] = vert;
      velocities[h] = start_vel;
      undeformed[h] = undef;
      vertHandles.push_back(h);
    }
  }


  //build the faces
  std::vector<Vec3i> tris;

  for(int i = 0; i < xresolution+1; ++i) {
    for(int j = 0; j < yresolution; ++j) {
      int tl = i + (xresolution+1)*j;
      int tr = (i+1)%(xresolution+1) + (xresolution+1)*j;
      int bl = i + (xresolution+1)*(j+1);
      int br = (i+1)%(xresolution+1) + (xresolution+1)*(j+1);

      shellObj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[br]);
      shellObj->addFace(vertHandles[tl], vertHandles[br], vertHandles[bl]);
    }
    //close the circle
  }

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


  //Find highest vertex
  VertexIterator vit = shellObj->vertices_begin();
  Scalar highest = -10000;
  for(;vit!= shellObj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest) {
      highest = pos[1];
    }
  }
  //Pin all verts at or near that height
  for(vit = shellObj->vertices_begin();vit!= shellObj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest - 1e-4)
      shell->constrainVertex(*vit, pos);
  }
}

//spherical shell
void ShellTest::setupScene3() {

  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");

  //build a rectangular grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);

  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);

  //create a sphere
  int layers = yresolution;
  int slices = xresolution;
  Vec3d centre(0,0,0);
  Scalar radius = 1.0;
  Vec3d start_vel(0,0,0);

  
  /*
  //create top pole
  VertexHandle topV = shellObj->addVertex();
  positions[topV] = centre + Vec3d(0,0,radius);
  velocities[topV] = start_vel;
  undeformed[topV] = positions[topV];

  //create bottom pole
  VertexHandle botV = shellObj->addVertex();
  positions[botV] = centre + Vec3d(0,0,-radius);
  velocities[botV] = start_vel;
  undeformed[botV] = positions[botV];

  //fill in the interior
  vertList.resize(layers-1);
  for(int j = 0; j < layers-1; ++j) {
    Scalar heightAngle = -pi/2 + (j+1) * pi/(Scalar)layers;
    for(int i = 0; i < slices; ++i) {
      Scalar rotAngle = 2*pi * (Scalar)i / (Scalar)slices;
      Scalar zVal = radius*sin(heightAngle);
      Scalar newRad = radius*cos(heightAngle);
      Scalar xVal = newRad*cos(rotAngle);
      Scalar yVal = newRad*sin(rotAngle);

      VertexHandle vNew = shellObj->addVertex();
      positions[vNew] = centre + Vec3d(xVal,yVal,zVal);
      velocities[vNew] = start_vel;
      undeformed[vNew] = positions[vNew];
      vertList[j].push_back(vNew);
    }
  }

  //construct faces
  for(int j = 0; j < layers; ++j) {
    for(int i = 0; i < slices; ++i) {
      if(j == 0) {
        shellObj->addFace(botV, vertList[j][i], vertList[j][(i+1)%slices]);
      }
      else if(j == layers-1) {
        shellObj->addFace(topV, vertList[j-1][(i+1)%slices], vertList[j-1][i]);
      }
      else {
        shellObj->addFace(vertList[j-1][i], vertList[j][i], vertList[j][(i+1)%slices]);
        shellObj->addFace(vertList[j-1][i], vertList[j][(i+1)%slices], vertList[j-1][(i+1)%slices]);
      }

    }
  }
  */

  //load a nicer sphere from disk
  ifstream infile("sphere1.obj");
  if(!infile)
    std::cout << "Error loading file\n";
  std::string line;
  std::vector<VertexHandle> vertList;
  while(!infile.eof()) {
    std::getline(infile, line);
    if(line.substr(0,1) == std::string("v")) {
      std::stringstream data(line);
      char c;
      Vec3d point;
      data >> c >> point[0] >> point[1] >> point[2];

      VertexHandle vNew = shellObj->addVertex();
      positions[vNew] = point;
      velocities[vNew] = start_vel;
      undeformed[vNew] = point;
      vertList.push_back(vNew);
    }
    else {
      std::stringstream data(line);
      char c;
      int v0,v1,v2;
      data >> c >> v0 >> v1 >> v2;
      shellObj->addFace(vertList[v0-1],vertList[v1-1],vertList[v2-1]);
    }
  }

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

  //inflate the sphere by some fixed amount
  /*Scalar inflateDist = 0.1;
  for(VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit) {  
    Vec3d oldPos = positions[*vit];
    Vec3d normVec = (oldPos-centre);
    normVec.normalize();
    shell->setVertexPosition(*vit, oldPos + inflateDist*normVec);
  }*/
  
  //Add an outward pressure force to inflate the sphere
  Scalar pressureStrength = GetScalarOpt("shell-inflate-sphere-coeff");
  bool constPressure = GetBoolOpt("shell-inflate-sphere-const-pressure");
  //Scalar pressureStrength = 0.1;
  shell->addForce(new ShellRadialForce(*shell, "Radial", Vec3d(0,0,0), pressureStrength, constPressure));

  //shell->addForce(new ShellVolumeForce(*shell, "Volume", 0.5));
}

//simple square with two triangles, one pinned in place, to test bending
void ShellTest::setupScene4() {
  
  //build a rectangular grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);

  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);

  //set up a small chunk of shell for testing bending
  VertexHandle v0 = shellObj->addVertex();
  VertexHandle v1 = shellObj->addVertex();
  VertexHandle v2 = shellObj->addVertex();
  VertexHandle v3 = shellObj->addVertex();
  FaceHandle f0 = shellObj->addFace(v0, v1, v2);
  FaceHandle f1 = shellObj->addFace(v2, v1, v3);

  //set up a square
  positions[v0] = undeformed[v0] = Vec3d(0,0,0);
  positions[v1] = undeformed[v1] = Vec3d(0,0,-1);
  positions[v2] = undeformed[v2] = Vec3d(1,0,0);
  positions[v3] = undeformed[v3] = Vec3d(1,0,-1);
  velocities[v0] = velocities[v1] = velocities[v2] = velocities[v3] = Vec3d(0,0,0);

  //initialize all edges to zero angle for now
  for(EdgeIterator eit = shellObj->edges_begin(); eit!= shellObj->edges_end(); ++eit) {
    EdgeHandle eh = *eit;
    undefAngle[*eit] = edgeAngle[*eit] = edgeVel[*eit] = 0;
  }
  
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


  //CONSTRAINTS

  //Just pin the first triangle right where it is.
  shell->constrainVertex(v0, shell->getVertexPosition(v0));
  shell->constrainVertex(v1, shell->getVertexPosition(v1));
  shell->constrainVertex(v2, shell->getVertexPosition(v2));

}

//a catenary sheet
void ShellTest::setupScene5() {

  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");

  Scalar dx = (Scalar)width / (Scalar)xresolution;
  Scalar dy = (Scalar)height / (Scalar)yresolution;

  //build a rectangular grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);

  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);

  Vec3d start_vel(0,0,0);
  for(int j = 0; j <= yresolution; ++j) {
    for(int i = 0; i <= xresolution; ++i) {
      Vec3d vert(i*dx, 0, j*dy);
    /*  if(j < 0.5*yresolution) {
        int k = j;
        int j_mod = (int)(0.5*yresolution);
        vert(1) = j_mod*dx;
        vert(2) = (k-j_mod)*dx;
      }*/
      Vec3d undef = vert;

      VertexHandle h = shellObj->addVertex();

      positions[h] = vert;
      velocities[h] = start_vel;
      undeformed[h] = undef;
      vertHandles.push_back(h);
    }
  }


  //build the faces
  std::vector<Vec3i> tris;
  for(int i = 0; i < xresolution; ++i) {
    for(int j = 0; j < yresolution; ++j) {
      int tl = i + (xresolution+1)*j;
      int tr = i+1 + (xresolution+1)*j;
      int bl = i + (xresolution+1)*(j+1);
      int br = i+1 + (xresolution+1)*(j+1);

      shellObj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[br]);
      shellObj->addFace(vertHandles[tl], vertHandles[br], vertHandles[bl]);
    }
  }

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

  //Find the leftest vertex
  VertexIterator vit = shellObj->vertices_begin();
  Scalar lowest = 10000;
  for(;vit!= shellObj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[0] <= lowest) {
      lowest = pos[0];
    }
  }
  //Find the rightest vertex
  vit = shellObj->vertices_begin();
  Scalar highest = -10000;
  for(;vit!= shellObj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[0] >= highest) {
      highest = pos[0];
    }
  }

  //Pin all verts at or near that height
  for(vit = shellObj->vertices_begin();vit!= shellObj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[0] >= highest - 1e-4) {
      shell->constrainVertex(*vit, pos);
    }
    if(pos[0] <= lowest + 1e-4) {
      shell->constrainVertex(*vit, pos);
    }
  }
 

}

//a hemispherical bubble with a hole in the top and pinned at the bottom ring
void ShellTest::setupScene6() {
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");

  //vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);

  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);

  //create a sphere
  int layers = yresolution;
  int slices = xresolution;
  Vec3d centre(0,0,0);
  Scalar radius = 0.01;
  Vec3d start_vel(0,0,0);

  std::vector<std::vector<VertexHandle> > vertList;

  
   int seed = 0;
  //fill in the interior
  vertList.resize(layers-1);
  for(int j = 0; j < layers-1; ++j) {
    Scalar heightAngle = (j+1) * 0.95 * pi / 2 /(Scalar)layers;
    for(int i = 0; i < slices; ++i) {
      Scalar rotAngle = 2*pi * (Scalar)i / (Scalar)slices;
      Scalar zVal = radius*sin(heightAngle);
      Scalar newRad = radius*cos(heightAngle);
      Scalar xVal = newRad*cos(rotAngle);
      Scalar yVal = newRad*sin(rotAngle);

      float x_n,y_n,z_n;
      x_n = 0.00001*ElTopo::randhashd(seed++);
      y_n = 0.00001*ElTopo::randhashd(seed++);
      z_n = 0.00001*ElTopo::randhashd(seed++);
      VertexHandle vNew = shellObj->addVertex();
      positions[vNew] = centre + Vec3d(xVal,zVal,yVal) + Vec3d(x_n, y_n, z_n);
      velocities[vNew] = start_vel;
      undeformed[vNew] = positions[vNew];
      vertList[j].push_back(vNew);
    }
  }

  //construct faces
  for(int j = 0; j < layers-2; ++j) {
    for(int i = 0; i < slices; ++i) {
      shellObj->addFace(vertList[j][i], vertList[j+1][i], vertList[j+1][(i+1)%slices]);
      shellObj->addFace(vertList[j][i], vertList[j+1][(i+1)%slices], vertList[j][(i+1)%slices]);
    }
  }

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

 
  float freeze_height = 0.003;
  for(VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit) {
    VertexHandle vh = *vit; 
    Vec3d position = shell->getVertexPosition(vh);
    if(position[1] < freeze_height)
      shell->constrainVertex(vh, position);
  }


}

//Planar rotating points
class XZPlaneRotatingConstraint : public PositionConstraint {
  Scalar m_rate;
  Vec3d m_position, m_centre;
  Scalar m_rad;
  Scalar m_baseangle;
public:
  XZPlaneRotatingConstraint(Vec3d& position, Vec3d& centre, Scalar rate):m_position(position), m_centre(centre), m_rate(rate) {
    Vec3d offset = m_position - m_centre;
    offset[1] = 0;
    m_rad = offset.norm();
    m_baseangle = atan2(m_position[2]-centre[2], m_position[0]-centre[0]);
  }

  Vec3d operator()(Scalar time) {
    Scalar angle = m_baseangle + m_rate*time;
    return m_centre + m_rad*Vec3d(cos(angle), 0, sin(angle));
  }

};

//Planar rotating points
class XZPlaneVariableRotationConstraint : public PositionConstraint {
  
  std::vector<Scalar> m_times, m_start_rates, m_accelerations;

  Vec3d m_position, m_centre;
  Scalar m_rad;
  Scalar m_baseangle;

public:
  XZPlaneVariableRotationConstraint(Vec3d& position, Vec3d& centre, 
                            std::vector<Scalar>& times, std::vector<Scalar>& startRates, std::vector<Scalar>& accels)
    :m_position(position), m_centre(centre), m_times(times), m_start_rates(startRates), m_accelerations(accels) {
    Vec3d offset = m_position - m_centre;
    offset[1] = 0;
    m_rad = offset.norm();
    m_baseangle = atan2(m_position[2]-centre[2], m_position[0]-centre[0]);
  }

  Vec3d operator()(Scalar time) {
    
    //integrate up from t=0 to t=time to get the current angle
    int segment = 0;
    Scalar accum_time = 0;
    Scalar accum_angle = m_baseangle;
    while(m_times[segment] < time) {
      Scalar dt = m_times[segment+1] < time? m_times[segment+1]-m_times[segment] : time - m_times[segment];
      accum_angle += m_start_rates[segment]*dt + 0.5*m_accelerations[segment]*square(dt);
      ++segment;
    }

    //work out the current target position
    return m_centre + m_rad*Vec3d(cos(accum_angle), 0, sin(accum_angle));
  }

};

//a horizontal sheet pinned between two circles
void ShellTest::setupScene7() {
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");

  //vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);

  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);

  //create a sphere
  int layers = yresolution;
  int slices = xresolution;
  Vec3d centre(0,0,0);
  Scalar out_radius = 0.25;
  Scalar in_radius = 0.1;
  Vec3d start_vel(0,0,0);
  Scalar rotation_rate = GetScalarOpt("shell-rotation-rate");

  std::vector<std::vector<VertexHandle> > vertList;

  Scalar dr = (out_radius - in_radius) / (Scalar) layers;
  
  //fill in the interior
  vertList.resize(layers+1);
  for(int j = 0; j < layers+1; ++j) {
    
    for(int i = 0; i < slices; ++i) {
      Scalar rotAngle = 2 * pi * (Scalar)i / (Scalar)slices;
      Scalar newRad = in_radius + j*dr;
      Scalar xVal = newRad*cos(rotAngle);
      Scalar yVal = newRad*sin(rotAngle);

      VertexHandle vNew = shellObj->addVertex();

      //perturbed vertically to generate some initial bending
      positions[vNew] = centre + Vec3d(xVal, 0.0001*sin(4*M_PI*(newRad-in_radius)/(out_radius-in_radius)), yVal);

      velocities[vNew] = start_vel;
      undeformed[vNew] = positions[vNew];
      vertList[j].push_back(vNew);
    }
  }

  //construct faces
  for(int j = 0; j < layers; ++j) {
    for(int i = 0; i < slices; ++i) {
      shellObj->addFace(vertList[j][i], vertList[j+1][i], vertList[j+1][(i+1)%slices]);
      shellObj->addFace(vertList[j][i], vertList[j+1][(i+1)%slices], vertList[j][(i+1)%slices]);
    }
  }

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

  //constrain inner and outer loops
  for(unsigned int i = 0; i < vertList[0].size(); ++i) {
    int outside = vertList.size()-1;
    int inside = 0;

    shell->constrainVertex(vertList[outside][i], shell->getVertexPosition(vertList[outside][i]));

    Vec3d pos = shell->getVertexPosition(vertList[inside][i]);
    
    //Constant rotation
    //XZPlaneRotatingConstraint*p = new XZPlaneRotatingConstraint(pos, centre, rotation_rate);

    //Use the rotation rate already inside
    std::vector<Scalar> times, rates, accels;
    times.push_back(0); rates.push_back(rotation_rate); accels.push_back(0);          //fixed velocity of 0.03 over [0, 400]
    times.push_back(400); rates.push_back(rotation_rate); accels.push_back(-rotation_rate/800.0); //decelerate linearly from 0.03 to 0 over [400,1200]
    times.push_back(1200); rates.push_back(0.0); accels.push_back(0);        //stop rotating
    times.push_back(2000); rates.push_back(0.0); accels.push_back(0);
    XZPlaneVariableRotationConstraint*p = new XZPlaneVariableRotationConstraint(pos, centre, times, rates, accels);
    
    shell->constrainVertex(vertList[inside][i], p);

  }
  
  Scalar bath_density = GetScalarOpt("shell-bath-density");
  //we assume gravity of -9.81 and bath height of 0.0
  shell->addForce(new ShellBathForce(*shell, "BathForce", Vec3d(0,-9.81,0), bath_density, 0.0));

}

//a torus
void ShellTest::setupScene8() {

  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");
  float height = (float)GetScalarOpt("shell-height");
  float width = (float)GetScalarOpt("shell-width");
  
  //width maps to the tube radius
  //height maps to the large radius
  float r0 = width;
  float r1 = height;

  //vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);

  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);

  //create a torus
  int resolution0 = xresolution;
  int resolution1 = yresolution;
  Vec3d centre(0,0,0);
  Vec3d start_vel(0,0,0);

  std::vector<std::vector<VertexHandle> > vertList;

  Scalar d_phi0 = 2*M_PI / (Scalar) (resolution0-1);
  Scalar d_phi1 = 2*M_PI / (Scalar) (resolution1-1);

  //fill in the interior
  vertList.resize(resolution1-1);
  for(int j = 0; j < resolution1-1; ++j) {
    Scalar phi1 = j*d_phi1;
    for(int i = 0; i < resolution0-1; ++i) {
      Scalar phi0 = i*d_phi0;  

      Scalar x = centre[0] + cos(phi1)*(r1 + r0*cos(phi0));
      Scalar y = centre[2] + sin(phi1)*(r1 + r0*cos(phi0));
      Scalar z = centre[1] + r0 * sin(phi0) + 0.3*r0*sin(10*phi1);
      Vec3d pos(x,z,y);
      
      VertexHandle vNew = shellObj->addVertex();
      positions[vNew] = pos;
      velocities[vNew] = start_vel;
      undeformed[vNew] = positions[vNew];
      vertList[j].push_back(vNew);
    }
  }

  
  //construct faces
  std::cout << "Testing\n";
  for(unsigned int j = 0; j < vertList.size(); ++j) {
    for(unsigned int i = 0; i < vertList[j].size(); ++i) {
      int j_next = (j+1)%vertList.size();
      int i_next = (i+1)%vertList[j].size();

      shellObj->addFace(vertList[j][i], vertList[j_next][i_next], vertList[j_next][i]);
      shellObj->addFace(vertList[j][i], vertList[j][i_next], vertList[j_next][i_next]);
    }
  }
  
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

  ////constrain inner and outer loops
  //for(unsigned int i = 0; i < vertList[0].size(); ++i) {
  //  int outside = vertList.size()-1;
  //  int inside = 0;

  //  shell->constrainVertex(vertList[outside][i], shell->getVertexPosition(vertList[outside][i]));

  //  Vec3d pos = shell->getVertexPosition(vertList[inside][i]);
  //  XZPlaneRotatingConstraint*p = new XZPlaneRotatingConstraint(pos, centre, 20);
  //  shell->constrainVertex(vertList[inside][i], p);

  //}
  shell->addForce(new ShellVolumeForce(*shell, "Volume", 10));

}


//the triangles connected about a single edge, in non-manifold way
void ShellTest::setupScene9() {

  //build a rectangular grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);

  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);

  //set up a small chunk of shell for testing bending
  VertexHandle v0 = shellObj->addVertex();
  VertexHandle v1 = shellObj->addVertex();
  VertexHandle v2 = shellObj->addVertex();
  VertexHandle v3 = shellObj->addVertex();
  VertexHandle v4 = shellObj->addVertex();
  VertexHandle v5 = shellObj->addVertex();

  FaceHandle f0 = shellObj->addFace(v0, v1, v2);
  FaceHandle f1 = shellObj->addFace(v2, v1, v3);
  FaceHandle f2 = shellObj->addFace(v2, v3, v4);
  FaceHandle f3 = shellObj->addFace(v3, v4, v5);

  /*VertexHandle v5 = shellObj->addVertex();
  VertexHandle v6 = shellObj->addVertex();
  VertexHandle v7 = shellObj->addVertex();

  FaceHandle f3 = shellObj->addFace(v5,v6,v7);*/

  //set up a square
  positions[v0] = undeformed[v0] = Vec3d(0,0,0);
  positions[v1] = undeformed[v1] = Vec3d(0,0,-1);
  positions[v2] = undeformed[v2] = Vec3d(1,0,0);
  positions[v3] = undeformed[v3] = Vec3d(1,0,-1);
  
  positions[v4] = undeformed[v4] = Vec3d(1.5,0,-0.5);
  positions[v5] = undeformed[v5] = Vec3d(1.5,0,-1);

  velocities[v0] = velocities[v1] = velocities[v2] = velocities[v3] = velocities[v4] = velocities[v5] = Vec3d(0,0,0);

  //initialize all edges to zero offsets
  for(EdgeIterator eit = shellObj->edges_begin(); eit!= shellObj->edges_end(); ++eit) {
    EdgeHandle eh = *eit;
    undefAngle[*eit] = edgeAngle[*eit] = edgeVel[*eit] = 0;
  }

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
  
  /*ShellVertexTriSpringForce* spring = new ShellVertexTriSpringForce(*shell, "SpringTest", m_timestep);
  spring->addSpring(f3, v0, Vec3d(0.4, 0.4, 0.2), 0.0, 0.001, 0.5);
  shell->addForce(spring);*/

  /*ShellVertexPointSpringForce* spring = new ShellVertexPointSpringForce(*shell, "Spring", m_timestep);
  spring->addSpring(v7, Vec3d(0.1, 0.3, 0.3), 0.1, 0.0, 0.0);
  shell->addForce(spring);*/


  //CONSTRAINTS

  //Just pin the first triangle right where it is.
  shell->constrainVertex(v0, shell->getVertexPosition(v0));
  shell->constrainVertex(v1, shell->getVertexPosition(v1));
  shell->constrainVertex(v2, shell->getVertexPosition(v2));

}

//vertical flat sheet with inflow at the top
void ShellTest::setupScene10() {

  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");

  Scalar dx = (Scalar)width / (Scalar)xresolution;
  Scalar dy = (Scalar)height / (Scalar)yresolution;

  //build a rectangular grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);

  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);

  Vec3d start_vel(0,-2,0);
  std::set<VertexHandle> topVerts;
  for(int j = 0; j <= yresolution; ++j) {
    for(int i = 0; i <= xresolution; ++i) {
      Vec3d vert(i*dx, j*dy, 0);
      Vec3d undef = vert;

      VertexHandle h = shellObj->addVertex();

      positions[h] = vert;
      velocities[h] = start_vel;
      undeformed[h] = undef;
      vertHandles.push_back(h);
      if(j == yresolution)
        topVerts.insert(h);
    }
  }

  //build the faces
  std::vector<Vec3i> tris;
  for(int i = 0; i < xresolution; ++i) {
    for(int j = 0; j < yresolution; ++j) {
      int tl = i + (xresolution+1)*j;
      int tr = i+1 + (xresolution+1)*j;
      int bl = i + (xresolution+1)*(j+1);
      int br = i+1 + (xresolution+1)*(j+1);

      shellObj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[br]);
      shellObj->addFace(vertHandles[tl], vertHandles[br], vertHandles[bl]);
    }
  }

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

  //add an inflow at the top of the sheet
  
  std::vector<EdgeHandle> extendEdgeList;
  EdgeIterator eit = shellObj->edges_begin();
  for(;eit != shellObj->edges_end(); ++eit) {
    EdgeHandle eh = *eit;
    VertexHandle vh0 = shellObj->fromVertex(eh);
    VertexHandle vh1 = shellObj->toVertex(eh);
    Vec3d pos0 = shell->getVertexPosition(vh0);
    Vec3d pos1 = shell->getVertexPosition(vh1);
    if(topVerts.find(vh0) != topVerts.end() && topVerts.find(vh1) != topVerts.end()) {
      extendEdgeList.push_back(eh);
    }
  }
  Vec3d inflow_vel = start_vel;
  
  shell->setInflowSection(extendEdgeList, inflow_vel, m_initial_thickness);
  shell->setDeletionBox(Vec3d(-1, -10, -1), Vec3d(2, -4.0, 1));

}

//a cube
void ShellTest::setupScene11() {

  float height = (float)GetScalarOpt("shell-height");
  float width = (float)GetScalarOpt("shell-width");
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

  Scalar offset = height;
  for(int i = 0; i < 8; ++i) {
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

  for(int i = 0; i < 8; ++i) {
    undeformed[vertList[i]] = positions[vertList[i]];
  }
  
  shellObj->addFace(vertList[0], vertList[2], vertList[4]);
  shellObj->addFace(vertList[6], vertList[4], vertList[2]);
  shellObj->addFace(vertList[4], vertList[6], vertList[5]);
  shellObj->addFace(vertList[5], vertList[6], vertList[7]);
  shellObj->addFace(vertList[7], vertList[6], vertList[3]);
  shellObj->addFace(vertList[3], vertList[6], vertList[2]);
  shellObj->addFace(vertList[3], vertList[2], vertList[1]);
  shellObj->addFace(vertList[1], vertList[2], vertList[0]);
  shellObj->addFace(vertList[0], vertList[4], vertList[5]);
  shellObj->addFace(vertList[0], vertList[5], vertList[1]);
  shellObj->addFace(vertList[7], vertList[3], vertList[5]);
  shellObj->addFace(vertList[5], vertList[3], vertList[1]);

  
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

  shell->addForce(new ShellVolumeForce(*shell, "Volume", 10));

}


void ShellTest::setupScene12() {
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");

  //vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);

  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);

  //create a sphere
  int layers = yresolution;
  int slices = xresolution;
  Vec3d centre(0,0,0);
  Scalar radius = 5.3;
  Vec3d start_vel(0,0,0);

  std::vector<std::vector<VertexHandle> > vertList;

  //fill in the interior
  vertList.resize(layers);
  for(int j = 0; j < layers; ++j) {
    Scalar heightAngle = j * 0.98 * pi / 2 /(Scalar)layers;
    for(int i = 0; i < slices; ++i) {
      Scalar rotAngle = 2*pi * (Scalar)i / (Scalar)slices;
      Scalar zVal = radius*sin(heightAngle);
      Scalar newRad = radius*cos(heightAngle);
      Scalar xVal = newRad*cos(rotAngle);
      Scalar yVal = newRad*sin(rotAngle);

      VertexHandle vNew = shellObj->addVertex();
      positions[vNew] = centre + Vec3d(xVal,zVal,yVal);
      velocities[vNew] = start_vel;
      undeformed[vNew] = positions[vNew];
      vertList[j].push_back(vNew);
    }
  }

  //construct faces
  for(int j = 0; j < layers-1; ++j) {
    for(int i = 0; i < slices; ++i) {
      if((i+j)%2 == 0) {
        shellObj->addFace(vertList[j][i], vertList[j+1][i], vertList[j+1][(i+1)%slices]);
        shellObj->addFace(vertList[j][i], vertList[j+1][(i+1)%slices], vertList[j][(i+1)%slices]);
      }
      else {
        shellObj->addFace(vertList[j][i], vertList[j+1][i], vertList[j][(i+1)%slices]);
        shellObj->addFace(vertList[j+1][i], vertList[j+1][(i+1)%slices], vertList[j][(i+1)%slices]);
      }
    }
  }

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

  //pin just the bottom layer
  for(unsigned int i = 0; i < vertList[0].size(); ++i)
    shell->constrainVertex(vertList[0][i], shell->getVertexPosition(vertList[0][i]));
  
  //construct list of hole and base edges
  std::vector<EdgeHandle> holeEdges, baseEdges;
  
  //Find the boundary loops
  int last = vertList.size()-1;
  for(EdgeIterator eit = shell->getDefoObj().edges_begin(); eit != shell->getDefoObj().edges_end(); ++eit) {
    EdgeHandle eh = *eit;
    VertexHandle v0 = shell->getDefoObj().fromVertex(eh);
    VertexHandle v1 = shell->getDefoObj().toVertex(eh);
    if(std::find(vertList[0].begin(), vertList[0].end(), v0) != vertList[0].end() && 
       std::find(vertList[0].begin(), vertList[0].end(), v1) != vertList[0].end()) {
      baseEdges.push_back(eh);
    }
    else if(std::find(vertList[last].begin(), vertList[last].end(), v0) != vertList[last].end() && 
            std::find(vertList[last].begin(), vertList[last].end(), v1) != vertList[last].end()) {
      holeEdges.push_back(eh);
    }
  }

  //Scalar air_density = 1.225e-9; 
  //shell->addForce(new DrainingBubblePressureForce(*shell, "DrainingBubblePressure", holeEdges, baseEdges, air_density, m_timestep));

}

//vertical flat sheet with inflow at the top
void ShellTest::setupScene13() {

  //get params
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("shell-x-resolution");
  int yresolution = GetIntOpt("shell-y-resolution");

  Scalar dx = (Scalar)width / (Scalar)xresolution;
  Scalar dy = (Scalar)height / (Scalar)yresolution;

  //build a rectangular grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);

  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);

  Scalar drop_vel = GetScalarOpt("shell-initial-velocity");

  Vec3d start_vel(0,drop_vel,0);
  std::set<VertexHandle> topVerts;
  for(int j = 0; j <= yresolution; ++j) {
    for(int i = 0; i <= xresolution; ++i) {
      Vec3d vert(i*dx, j*dy, 0.00001*sin(10000*(j*dy)));
      Vec3d undef = vert;

      VertexHandle h = shellObj->addVertex();

      positions[h] = vert;
      velocities[h] = start_vel;
      undeformed[h] = undef;
      vertHandles.push_back(h);
      if(j == yresolution)
        topVerts.insert(h);
    }
  }

  //build the faces
  std::vector<Vec3i> tris;
  for(int i = 0; i < xresolution; ++i) {
    for(int j = 0; j < yresolution; ++j) {
      int tl = i + (xresolution+1)*j;
      int tr = i+1 + (xresolution+1)*j;
      int bl = i + (xresolution+1)*(j+1);
      int br = i+1 + (xresolution+1)*(j+1);

      if((i+j)%2 == 0) {
        shellObj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[br]);
        shellObj->addFace(vertHandles[tl], vertHandles[br], vertHandles[bl]);
      }
      else {
        shellObj->addFace(vertHandles[tl], vertHandles[tr], vertHandles[bl]);
        shellObj->addFace(vertHandles[bl], vertHandles[tr], vertHandles[br]);
      }
    }
  }

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

  //add an inflow at the top of the sheet

  std::vector<EdgeHandle> extendEdgeList;
  EdgeIterator eit = shellObj->edges_begin();
  for(;eit != shellObj->edges_end(); ++eit) {
    EdgeHandle eh = *eit;
    VertexHandle vh0 = shellObj->fromVertex(eh);
    VertexHandle vh1 = shellObj->toVertex(eh);
    Vec3d pos0 = shell->getVertexPosition(vh0);
    Vec3d pos1 = shell->getVertexPosition(vh1);
    if(topVerts.find(vh0) != topVerts.end() && topVerts.find(vh1) != topVerts.end()) {
      extendEdgeList.push_back(eh);
    }
  }
  Vec3d inflow_vel = start_vel;

  shell->setInflowSection(extendEdgeList, inflow_vel, m_initial_thickness);
  Scalar ground = GetScalarOpt("shell-ground-plane-height");
  //shell->setDeletionBox(Vec3d(-2, -5, -2), Vec3d(2, ground-0.02, 2));

}
