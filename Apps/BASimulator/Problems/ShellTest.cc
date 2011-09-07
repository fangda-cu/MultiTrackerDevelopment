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
#include "BASim/src/Physics/DeformableObjects/Shells/ShellGravityForce.hh"
#include "BASim/src/Render/ShellRenderer.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"

#include <fstream>

ShellTest::ShellTest()
: Problem("Shell Test", "A rectangular shell suspended by its corner points"), 
  shell(NULL), shellObj(NULL), stepper(NULL)
{
  addDynamicsProps();
  
  //shell options
  AddOption("shell-thickness", "the thickness of the shell", 0.01);
  AddOption("shell-width", "the horizontal side length of the shell", 1.0);
  AddOption("shell-height", "the vertical side length of the shell", 1.0);
  AddOption("x-resolution", "the number of segments along the horizontal edge", 30);
  AddOption("y-resolution", "the number of segments along the vertical edge", 30);
  AddOption("shell-density", "volumetric density of the shell ", 1.0);
  
  //discrete shell bending stiffnesses
  AddOption("shell-bending-stiffness", "the bending stiffness of the shell", 0.0f);
  AddOption("shell-bending-damping", "the bending viscosity of the shell", 0.0f);

  //thickness-dependent elasticity with proper physical-parameters, only for CST membrane so far.
  AddOption("shell-Poisson", "the Poisson ratio of the shell material", 0.0f);
  AddOption("shell-Youngs", "the Young's modulus of the shell material", 0.0f);
  AddOption("shell-Poisson-damping", "the Poisson ratio of the shell material", 0.0f);
  AddOption("shell-Youngs-damping", "the Young's modulus of the shell material", 0.0f);

  //timestepper options
  AddOption("integrator", "type of integrator to use for the shell", "symplectic");
  AddOption("iterations", "maximum number of iterations for the implicit method", (int) 100);
  AddOption("atol", "absolute convergence tolerance", 1e-8);
  AddOption("rtol", "relative convergence tolerance", 1e-8);
  AddOption("stol", "convergence tolerance in terms of the norm of the change in the solution between steps", 1e-8);
  AddOption("inftol", "infinity norm convergence tolerance", 1e-8);
 
  // default to no gravity
  GetVecOpt("gravity") = Vec3d::Zero();
}

ShellTest::~ShellTest()
{
  if (shellObj != NULL) delete shellObj;
  if (shell != NULL) delete shell;
  if (stepper != NULL) delete stepper;
}

void ShellTest::Setup()
{

  loadDynamicsProps();

  Scalar density = GetScalarOpt("shell-density");
  Scalar thickness = GetScalarOpt("shell-thickness");
  Scalar width = GetScalarOpt("shell-width");
  Scalar height = GetScalarOpt("shell-height");
  int xresolution = GetIntOpt("x-resolution");
  int yresolution = GetIntOpt("y-resolution");
  Vec3d gravity = GetVecOpt("gravity");
  
  Scalar bend_stiffness = GetScalarOpt("shell-bending-stiffness");
  Scalar bend_damping = GetScalarOpt("shell-bending-damping");
  
  Scalar Youngs_modulus = GetScalarOpt("shell-Youngs");
  Scalar Poisson_ratio = GetScalarOpt("shell-Poisson");
  Scalar Youngs_damping = GetScalarOpt("shell-Youngs-damping");
  Scalar Poisson_damping = GetScalarOpt("shell-Poisson-damping");

  std::string integrator = GetStringOpt("integrator");

  Scalar dx = (Scalar)width / (Scalar)xresolution;
  Scalar dy = (Scalar)height / (Scalar)yresolution;

  //Create the base deformable object (mesh)
  shellObj = new DeformableObject();

  //build a rectangular grid of vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);
  VertexProperty<Vec3d> positions(shellObj);

  /*for(int i = 0; i <= xresolution; ++i) {
     for(int j = 0; j <= yresolution; ++j) {
        Vec3d vert(i*dx, 0, j*dy);
        Vec3d undef = vert;

        VertexHandle h = shellObj->addVertex();
        positions[h] = vert;
        velocities[h] = Vec3d(0,0,0);
        undeformed[h] = undef;
        vertHandles.push_back(h);
     }
  }*/

  
  for(int i = 0; i <= xresolution; ++i) {
    for(int j = 0; j <= yresolution; ++j) {
      Vec3d vert(i*dx, j*dy, 0);
     /* if(j < 0.5*yresolution) {
        int k = j;
        int j_mod = (int)(0.5*yresolution);
        vert(1) = j_mod*dx;
        vert(2) = (k-j_mod)*dx;
      }*/
      Vec3d undef = vert;
        
      VertexHandle h = shellObj->addVertex();

      positions[h] = vert;
      velocities[h] = Vec3d(0,0,0);
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
  shell = new ElasticShell(shellObj, shellFaces);
  shellObj->addModel(shell);

  //now add forces to the model
  Scalar timestep = getDt(); //Our Rayleigh damping model relies on knowing the timestep (folds it into the damping stiffness, as in Viscous Threads)
 
  if(bend_stiffness != 0 || bend_damping != 0)
    shell->addForce(new DSBendingForce(*shell, "DSBending", bend_stiffness, bend_damping, timestep));

  if(Youngs_modulus != 0 || Youngs_damping != 0) {
    shell->addForce(new CSTMembraneForce(*shell, "CSTMembrane", Youngs_modulus, Poisson_ratio, Youngs_damping, Poisson_damping, timestep));
  }
  

  shell->addForce(new ShellGravityForce(*shell, "Gravity", gravity));

  //and set its properties, including geometry
  shell->setThickness(thickness);
  shell->setDensity(density);
  shell->setUndeformedConfig(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);

  for(int i  = 0; i < 5; ++i)
    shell->remesh(0.05);

  std::cout << "Finding top vertices\n";
  //find the top left and right corners for adding a constraint
  VertexIterator vit = shellObj->vertices_begin();
  Vec3d minPos(0,0,0), maxPos(0,0,0);
  VertexHandle minH, maxH;
  for(;vit!= shellObj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[0] <= minPos[0] && pos[1] >= minPos[1]) {
      minH = *vit;
      minPos = pos;
    }
    if(pos[0] >= maxPos[0] && pos[1] >= maxPos[1]) {
      maxH = *vit;
      maxPos = pos;
    }
  }
  
  
  //Pin just the left and right top corners
 /* shell->constrainVertex(minH, minPos);
  shell->constrainVertex(maxH, maxPos);*/
 
  /*
  int count = 0;
  for(FaceIterator fit = shellObj->faces_begin(); fit!=shellObj->faces_end(); ++fit) {
    ++count;
    if(count % 3 ==0)
      shellObj->deleteFace(*fit, false);
  }*/

  std::cout << "Pinning vertices\n";
  //Pin all vertices in the top row.
  vit = shellObj->vertices_begin();
  for(;vit!= shellObj->vertices_end(); ++vit) {
    Vec3d pos = shell->getVertexPosition(*vit);
    if(pos[1] >= 0.99)
      shell->constrainVertex(*vit, pos);
  }

 
  shell->computeMasses();

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

  /*
  //Dump out data for a meshopt testcas
  std::ofstream output("crease.obj");
  
  //(undeformed vertices)
  vit = shellObj->vertices_begin();
  for(;vit!=shellObj->vertices_end(); ++vit) {
    VertexHandle vh = *vit;
    Vec3d vertPos = shell->getUndeformedPosition(vh);
    output << "v " << vertPos[0]  << " " << vertPos[1] << " " << vertPos[2] << std::endl;
  }

  //faces
  FaceIterator fit = shellObj->faces_begin();
  for(;fit!=shellObj->faces_end(); ++fit) {
    FaceHandle fh = *fit;
    FaceVertexIterator fvit = shellObj->fv_iter(fh);
    output << "f "; 
    for(;fvit; ++fvit) {
      VertexHandle vh = *fvit; //shift by one for OBJ numbering
      output << (vh.idx()+1) << " ";
    }
    output << std::endl;
  }

  //initial positions
  vit = shellObj->vertices_begin();
  for(;vit!=shellObj->vertices_end(); ++vit) {
    VertexHandle vh = *vit;
    Vec3d vertPos = shell->getVertexPosition(vh);
    output << "i " << (vh.idx()+1) << " " << vertPos[0]  << " " << vertPos[1] << " " << vertPos[2] << std::endl;
  }

  //pinned vertices
  output << "b " << (minH.idx()+1) << " 1 1 1\n";
  output << "b " << (maxH.idx()+1) << " 1 1 1\n";
  output.close();
  */
  
  m_world->addObject(shellObj);
  m_world->addController(stepper);
  RenderBase* shellRender = new ShellRenderer(*shell);
  m_world->addRenderer(shellRender);
  std::cout << "Entering main\n";
}

void ShellTest::AtEachTimestep()
{
  // Test of object serialization
  //RodState origrodstate;
  //origrodstate.copyState(*rod);
  //rodstate.print(*rod);
  
  //TopologicalObjectSerializer testserializer;
  //testserializer.saveTopologicalObject( *rod, "testoutput.bin" );

  //ElasticRod* testrod = NULL;
  //testserializer.loadTopologicalObject( &testrod, "testoutput.bin" );
  //assert( testrod != NULL );

  //origrodstate.compareProperties(*testrod);
  //delete testrod;
  
 /* if (m_maxTwist > 0 && m_currentTwist >= m_maxTwist) return;

  Scalar twistIncrement = getDt() * m_twistRate;
  int edge = rod->ne() - 1;
  Scalar theta = rod->getTheta(edge) - twistIncrement;
  stepper->getBoundaryCondition()->setDesiredEdgeAngle(edge, theta);
  m_currentTwist += twistIncrement;*/
}








