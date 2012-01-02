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

#include <fstream>

ShellTest::ShellTest()
: Problem("Shell Test", "Various viscous and elastic sheet/shell tests"), 
  shell(NULL), shellObj(NULL), stepper(NULL)
{
  addDynamicsProps();
  
  //Choice of scene
  AddOption("shell-scene", "the shell scene to test", 1);

  //Basic shell options
  AddOption("shell-thickness", "the thickness of the shell", 0.01);
  AddOption("shell-density", "volumetric density of the shell ", 1.0);

  //Shell geometry (x/y/width/height may also refer to resolutions/sizes in non-cartesian scenarios)
  AddOption("shell-width", "the horizontal side length of the shell", 1.0);
  AddOption("shell-height", "the vertical side length of the shell", 1.0);
  AddOption("shell-x-resolution", "the number of segments along first dimension", 30);
  AddOption("shell-y-resolution", "the number of segments along second dimension", 30);
  
  //Remeshing options
  AddOption("shell-remeshing", "whether to perform remeshing", 0);
  AddOption("shell-remeshing-resolution", "target edge-length", 0.1);
  AddOption("shell-remeshing-iterations", "number of remeshing iterations to run", 2);

  //Area-based surface tension force
  AddOption("shell-surface-tension", "surface tension coefficient of the shell", 0.0);
  
  //Properties for proper thickness dependent elasticity & viscosity (just CSTMembrane so far)
  AddOption("shell-Poisson", "the Poisson ratio of the shell material", 0.0f);
  AddOption("shell-Youngs", "the Young's modulus of the shell material", 0.0f);
  AddOption("shell-Poisson-damping", "the damping coefficient associated to the shell's Poisson ratio", 0.0f);
  AddOption("shell-Youngs-damping", "the damping coefficient associated with the shell's Young's modulus", 0.0f);

  //Hinge bending (discrete shells) stiffness and damping
  AddOption("shell-bending-stiffness", "Hinge (Discrete shells) bending stiffness of the shell", 0.0);
  AddOption("shell-bending-damping", "Hinge (Discrete shells) bending damping coefficient of the shell ", 0.0);

  //Timestepper options
  AddOption("integrator", "type of integrator to use for the shell", "implicit");

  //Solver options
  AddOption("iterations", "maximum number of iterations for the implicit method", (int) 100);
  AddOption("atol", "absolute convergence tolerance", 1e-8);
  AddOption("rtol", "relative convergence tolerance", 1e-8);
  AddOption("stol", "convergence tolerance in terms of the norm of the change in the solution between steps", 1e-8);
  AddOption("inftol", "infinity norm convergence tolerance", 1e-8);
  

}

ShellTest::~ShellTest()
{
  if (shellObj != NULL) delete shellObj;
  if (shell != NULL) delete shell;
  if (stepper != NULL) delete stepper;
}

typedef void (ShellTest::*sceneFunc)();


sceneFunc scenes[] = {0,
                      &ShellTest::setupScene1,  //vertical flat sheet
                      &ShellTest::setupScene2, //vertical cylindrical sheet
                      &ShellTest::setupScene3, //spherical sheet
                      &ShellTest::setupScene4, //two-triangle bending test
                      &ShellTest::setupScene5, //catenary 
                      &ShellTest::setupScene6, //hemispherical bubble
                      &ShellTest::setupScene7, //sheet between two circles
                      &ShellTest::setupScene8, //torus
                      &ShellTest::setupScene9};  //non-manifold edge

void ShellTest::Setup()
{

  loadDynamicsProps();

  //General shell forces and properties
  Scalar density = GetScalarOpt("shell-density");
  Scalar thickness = GetScalarOpt("shell-thickness");
  
  Vec3d gravity = GetVecOpt("gravity");
  
  Scalar surface_tension = GetScalarOpt("shell-surface-tension");

  Scalar Youngs_modulus = GetScalarOpt("shell-Youngs");
  Scalar Poisson_ratio = GetScalarOpt("shell-Poisson");
  Scalar Youngs_damping = GetScalarOpt("shell-Youngs-damping");
  Scalar Poisson_damping = GetScalarOpt("shell-Poisson-damping");
  
  Scalar DSbendstiffness = GetScalarOpt("shell-bending-stiffness");
  Scalar DSbenddamping = GetScalarOpt("shell-bending-damping");

  std::string integrator = GetStringOpt("integrator");

  Scalar timestep = getDt(); //Our Rayleigh damping model relies on knowing the timestep (folds it into the damping stiffness, as in Viscous Threads)
  m_timestep = timestep;

  //Geometry/scene specific
  int sceneChoice = GetIntOpt("shell-scene");


  

  //Create the base deformable object (mesh)
  shellObj = new DeformableObject();

  (*this.*scenes[sceneChoice])();
  /*
  switch(sceneChoice) {
    case 1: setupScene1();
      break;
    case 2: setupScene2();
      break;
    case 3: setupScene3();
      break;
    case 4: setupScene4();
      break;
    default:
      break;
  }*/
  
  
  bool circular = false;

  
  


  /*
  Vec3d start_vel(0,0,0);
  if(!circular) {
     
     
  }
  else {
    
  }
  */
  
  
  /*
  
  */

  ////choose some initial velocities
  //for(VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit) {
  //  Vec3d radius = (positions[*vit] - centre);
  //  radius.normalize();
  //Scalar viscosity = Youngs_damping / 2 / (1 + Poisson_damping);
  //  //velocities[*vit] = pressureStrength/12.0/viscosity/thickness*radius;
  //  //velocities[*vit] = radius;
  //}

 

  //now add forces to the model

  //Stretching and bending forces
  if(Youngs_modulus != 0 || Youngs_damping != 0) {
    shell->addForce(new CSTMembraneForce(*shell, "CSTMembrane", Youngs_modulus, Poisson_ratio, Youngs_damping, Poisson_damping, timestep));
    //shell->addForce(new MNBendingForce(*shell, "MNBending", Youngs_modulus, Poisson_ratio, Youngs_damping, Poisson_damping, timestep));
  }

  if(DSbendstiffness != 0 || DSbenddamping !=0)
    shell->addForce(new DSBendingForce(*shell, "DSBending", DSbendstiffness, DSbenddamping, timestep));
  
  

  //Gravity force
  shell->addForce(new ShellGravityForce(*shell, "Gravity", gravity));

  //Surface tension force
  if(surface_tension != 0)
    shell->addForce(new ShellSurfaceTensionForce(*shell, "Surface Tension", surface_tension));
  

  //and set its standard properties
  shell->setThickness(thickness);
  shell->setDensity(density);
  
  bool remeshing = GetIntOpt("shell-remeshing") == 1?true:false;
  Scalar remeshing_res = GetScalarOpt("shell-remeshing-resolution");
  int remeshing_its = GetIntOpt("shell-remeshing-iterations");
  

  shell->setRemeshing(remeshing, remeshing_res, remeshing_its);

 /* std::vector<EdgeHandle> extendEdgeList;
  EdgeIterator eit = shellObj->edges_begin();
  for(;eit != shellObj->edges_end(); ++eit) {
    EdgeHandle eh = *eit;
    VertexHandle vh0 = shellObj->fromVertex(eh);
    VertexHandle vh1 = shellObj->toVertex(eh);
    Vec3d pos0 = shell->getVertexPosition(vh0);
    Vec3d pos1 = shell->getVertexPosition(vh1);
    if(pos0[1] >= -0.01 && pos1[1] >= -0.01) {
      extendEdgeList.push_back(eh);
      std::cout << "Edge: " << eh.idx() << std::endl;
    }
  }
  Vec3d inflow_vel = start_vel;
  shell->setInflowSection(extendEdgeList, inflow_vel);*/


  //for(int q = 0; q < 0; ++q) {
  //  int c = 0; 
  //  std::cout << "Setting up adjusted positions\n";
  //  VertexIterator vit2 = shellObj->vertices_begin();
  //  for(;vit2 != shellObj->vertices_end(); ++vit2) {
  //    //if(c++==0) {
  //    //  continue;
  //    //}
  //    
  //    Vec3d pos = shell->getVertexPosition(*vit2);
  //    pos += Vec3d(0,-dy,0);
  //    shell->setVertexPosition(*vit2, pos);
  //    //std::cout << "After: " << pos << std::endl;
  //  }

  //  shell->extendMesh();
  //}

  //shell->remesh(0.2);
 

  shell->computeMasses();

  shell->setCollisionParams(true, 0.005);

  //compute the dof indexing for use in the diff_eq solver
  shellObj->computeDofIndexing();

  std::cout << "Set up\n";
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
  RenderBase* shellRender = new ShellRenderer(*shell);
  m_world->addRenderer(shellRender);
  
}

void ShellTest::AtEachTimestep()
{

}

//vertical flat sheet, pinned or flowing at top
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
 /* for(vit = shellObj->vertices_begin();vit!= shellObj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= highest - 1e-4)
      shell->constrainVertex(*vit, pos);
  }*/

}

//vertical cylinder, pinned or flowing at top
void ShellTest::setupScene2() {

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
  Scalar radius = 1.234;
  Vec3d start_vel(0,0,0);

  std::vector<std::vector<VertexHandle> > vertList;
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
  //Scalar pressureStrength = 0.1;
  //shell->addForce(new ShellRadialForce(*shell, "Radial", Vec3d(0,0,0), pressureStrength));

  shell->addForce(new ShellVolumeForce(*shell, "Volume", 0.5));
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
  Scalar radius = 1.0;
  Vec3d start_vel(0,0,0);

  std::vector<std::vector<VertexHandle> > vertList;

  //fill in the interior
  vertList.resize(layers-1);
  for(int j = 0; j < layers-1; ++j) {
    Scalar heightAngle = (j+1) * 0.9* pi / 2 /(Scalar)layers;
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

 
  for(unsigned int i = 0; i < vertList[0].size(); ++i)
    shell->constrainVertex(vertList[0][i], shell->getVertexPosition(vertList[0][i]));
 


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
  Scalar out_radius = 1.0;
  Scalar in_radius = 0.2;
  Vec3d start_vel(0,0,0);

  std::vector<std::vector<VertexHandle> > vertList;

  Scalar dr = (out_radius - in_radius) / (Scalar) layers;
  
  //fill in the interior
  vertList.resize(layers-1);
  for(int j = 0; j < layers-1; ++j) {
    
    for(int i = 0; i < slices; ++i) {
      Scalar rotAngle = 2 * pi * (Scalar)i / (Scalar)slices;
      Scalar newRad = in_radius + j*dr;
      Scalar xVal = newRad*cos(rotAngle);
      Scalar yVal = newRad*sin(rotAngle);

      VertexHandle vNew = shellObj->addVertex();
      positions[vNew] = centre + Vec3d(xVal, 0, yVal);
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

  //constrain inner and outer loops
  for(unsigned int i = 0; i < vertList[0].size(); ++i) {
    int outside = vertList.size()-1;
    int inside = 0;

    shell->constrainVertex(vertList[outside][i], shell->getVertexPosition(vertList[outside][i]));

    Vec3d pos = shell->getVertexPosition(vertList[inside][i]);
    XZPlaneRotatingConstraint*p = new XZPlaneRotatingConstraint(pos, centre, 20);
    shell->constrainVertex(vertList[inside][i], p);

  }

  /*shell->remesh(0.05);
  shell->remesh(0.05);
  shell->remesh(0.05);*/

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


