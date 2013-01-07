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
//#include "BASim/src/Physics/DeformableObjects/Shells/CSTMembraneForce.hh"
//#include "BASim/src/Physics/DeformableObjects/Shells/DSBendingForce.hh"
//#include "BASim/src/Physics/DeformableObjects/Shells/MNBendingForce.hh"
#include "BASim/src/Render/ShellRenderer.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellRadialForce.hh"
#include "BASim/src/Physics/DeformableObjects/GravityForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellSurfaceTensionForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellLinearSurfaceTensionForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellLinearSurfaceTensionForce2.hh"

#include "BASim/src/Physics/DeformableObjects/Shells/DrainingBubblePressureForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellBathForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellVerticalForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellPointForce.hh"
#include "BASim/src/Collisions/ElTopo/util.hh"

#include "DelaunayTriangulator.hh"

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

  AddOption("volume-force-stiffness", "stiffness of the penalty force maintaining volume of each region", 1000.0);
  
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
  0, //  &DoubleBubbleTest::setupScene3,
  0, //  &DoubleBubbleTest::setupScene4,
  &DoubleBubbleTest::setupScene5,
  &DoubleBubbleTest::setupScene6,
};

//Scalar db_bubbleThicknessFunction(Vec3d pos) {
//  //NOTE: Assumes radius of 0.01;
//  Scalar rad = 0.01f;
//  Scalar height_frac = pos[1] / 0.01f;
//  Scalar transition_point = 1.0f;
//  if(height_frac > transition_point) {
//    return 0.0001;
//  }
//  else {
//    Scalar s = height_frac/transition_point;
//    return s*0.0001 + (1-s)*0.00015;
//  }
//
//}

int DoubleBubbleTest::onBBWall(const Vec3d & pos) const
{
  int walls = 0;
  if (pos.x() < 0 + 1e-6)
    walls |= (1 << 0);
  if (pos.y() < 0 + 1e-6)
    walls |= (1 << 1);
  if (pos.z() < 0 + 1e-6)
    walls |= (1 << 2);
  if (pos.x() > 1 - 1e-6)
    walls |= (1 << 3);
  if (pos.y() > 1 - 1e-6)
    walls |= (1 << 4);
  if (pos.z() > 1 - 1e-6)
    walls |= (1 << 5);
  
  return walls;
}

void DoubleBubbleTest::Setup()
{

  loadDynamicsProps();

  //General shell forces and properties
//  Scalar density = GetScalarOpt("shell-density");
//  Scalar thickness = GetScalarOpt("shell-thickness");
//  m_initial_thickness = thickness;

  Vec3d gravity = GetVecOpt("gravity");
  
  Scalar surface_tension = GetScalarOpt("shell-surface-tension");

//  Scalar Youngs_modulus = GetScalarOpt("shell-Youngs");
//  Scalar Poisson_ratio = GetScalarOpt("shell-Poisson");
//  Scalar Youngs_damping = GetScalarOpt("shell-Youngs-damping");
//  Scalar Poisson_damping = GetScalarOpt("shell-Poisson-damping");
  
//  bool cst_stretch = GetBoolOpt("shell-CST-stretching");
//  bool ds_bend = GetBoolOpt("shell-DS-bending");
//  bool mn_bend = GetBoolOpt("shell-MN-bending");

//  //fudge factors to modify the elastic-viscous coefficients (so as to manipulate them separately)
//  Scalar cst_scale = GetScalarOpt("shell-stretching-factor");
//  Scalar ds_scale = GetScalarOpt("shell-bending-factor");
  
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

  shell->getVertexConstraintLabels().assign(0);
  
//  shell->setThickness(thickness); /////////////////////////

  //now add forces to the model
//  MNBendingForce* bender;
//
//  //Stretching and bending forces
//  if(Youngs_modulus != 0 || Youngs_damping != 0) {
//    
//    //Stretching force (Constant Strain Triangle, i.e. linear FEM)
//    if(cst_stretch)
//      shell->addForce(new CSTMembraneForce(*shell, "CSTMembrane", cst_scale*Youngs_modulus, Poisson_ratio, cst_scale*Youngs_damping, Poisson_damping, timestep));
//    
//    //Bending force (Hinge-based Bending, a la Discrete Shells)
//    if(ds_bend)
//      shell->addForce(new DSBendingForce(*shell, "DSBending", ds_scale*Youngs_modulus, Poisson_ratio, ds_scale*Youngs_damping, Poisson_damping, timestep));
//
//    //Better bending model (Mid-Edge normals, a la Computing discrete shape operators on general meshes)
//    if(mn_bend) {
//      bender = new MNBendingForce(*shell, "MNBending", ds_scale*Youngs_modulus, Poisson_ratio, ds_scale*Youngs_damping, Poisson_damping, timestep);
//      shell->addForce(bender);
//    }
//
//  }

//  //Gravity force
//  shellObj->addForce(new GravityForce(*shellObj, timestep, "Gravity", gravity)); 
//
//  Scalar loadForce = GetScalarOpt("shell-uniform-load");
//  if(loadForce != 0)
//    shell->addForce(new ShellVerticalForce(*shell, "Load", Vec3d(0,-1,0), loadForce));
  
  //Surface tension force
  if(surface_tension != 0) {
    //Viscous sheets-style surface tension
    shell->addForce(new ShellSurfaceTensionForce(*shell, "Surface Tension", surface_tension));
    
    //Experimental piecewise linear surface tension
    //shell->addForce(new ShellLinearSurfaceTensionForce(*shell, "Surface Tension", surface_tension));
    //shell->addForce(new ShellLinearSurfaceTensionForce2(*shell, "Surface Tension", surface_tension));
  }

//  shell->setDensity(density); ////////////////////////////
  
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

//  bool thickness_evolution = GetBoolOpt("shell-update-thickness");
//  shell->setThicknessUpdating(thickness_evolution);

  Scalar vf_stiffness = GetScalarOpt("volume-force-stiffness");
  if (vf_stiffness > 0)
  {
    svf = new ShellVolumeForce(*shell, "Volume", vf_stiffness);
    shell->addForce(svf);
  }
  
//  shell->computeMasses(); /////////////////////////////
 

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

//  shell->remesh();
    
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
  RenderBase* shellRender = new ShellRenderer(*shell);
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

class VertexHandleComp
{
public:
  bool operator() (const VertexHandle & v1, const VertexHandle & v2) const 
  { 
    return v1.idx() < v2.idx(); 
  }
};

class Vec3dComp
{
public:
  bool operator() (const Vec3d & v1, const Vec3d & v2) const 
  { 
    return v1.x() < v2.x() || 
          (v1.x() == v2.x() && v1.y() < v2.y()) || 
          (v1.x() == v2.x() && v1.y() == v2.y() && v1.z() < v2.z());
  }
};

void DoubleBubbleTest::AtEachTimestep()
{
//  for (size_t i = 0; i < triangulation_added_faces.size(); i++)
//    shellObj->deleteFace(triangulation_added_faces[i], false);
//  for (size_t i = 0; i < triangulation_added_edges.size(); i++)
//    shellObj->deleteEdge(triangulation_added_edges[i], false);
//  for (size_t i = 0; i < triangulation_added_vertices.size(); i++)
//    shellObj->deleteVertex(triangulation_added_vertices[i]);

  
    
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

//  if (m_active_scene == 4)
//  {
//    std::vector<std::set<Vec3d, Vec3dComp> > pos;               // bubble i vertex positions
//    std::vector<std::set<VertexHandle, VertexHandleComp> > vts; // bubble i vertices
//    pos.resize(m_s4_nbubble);
//    vts.resize(m_s4_nbubble);
//
//    for (FaceIterator fit = shellObj->faces_begin(); fit != shellObj->faces_end(); ++fit)
//    {
//      FaceHandle f = *fit;
//      FaceVertexIterator fvit = shellObj->fv_iter(f); assert(fvit);
//      Vec3d pos0 = shell->getVertexPosition(*fvit); VertexHandle v0 = *fvit; ++fvit; assert(fvit);
//      Vec3d pos1 = shell->getVertexPosition(*fvit); VertexHandle v1 = *fvit; ++fvit; assert(fvit);
//      Vec3d pos2 = shell->getVertexPosition(*fvit); VertexHandle v2 = *fvit; ++fvit; assert(!fvit);
//      
//      Vec2i label = shell->getFaceLabel(f);
//      for (int i = 0; i < m_s4_nbubble; i++)
//      {
//        if ((label.x() == i && label.y() == -1) || (label.y() == i && label.x() == -1))
//        {
//          pos[i].insert(pos0);
//          pos[i].insert(pos1);
//          pos[i].insert(pos2);
//          vts[i].insert(v0);
//          vts[i].insert(v1);
//          vts[i].insert(v2);
//        }
//      }
//    }
//    
//    for (int j = 0; j < m_s4_nbubble; j++)
//    {
//      // least square fitting of the sphere centers
//      Mat4d A;
//      VecXd b(4);
//      A.setZero();
//      b.setZero();
//      for (std::set<Vec3d>::iterator i = pos[j].begin(); i != pos[j].end(); i++)
//      {
//        A(0, 0) +=  8 * (*i).x() * (*i).x();
//        A(0, 1) +=  8 * (*i).y() * (*i).x();
//        A(0, 2) +=  8 * (*i).z() * (*i).x();
//        A(0, 3) += -4 * (*i).x();
//        A(1, 0) +=  8 * (*i).x() * (*i).y();
//        A(1, 1) +=  8 * (*i).y() * (*i).y();
//        A(1, 2) +=  8 * (*i).z() * (*i).y();
//        A(1, 3) += -4 * (*i).y();
//        A(2, 0) +=  8 * (*i).x() * (*i).z();
//        A(2, 1) +=  8 * (*i).y() * (*i).z();
//        A(2, 2) +=  8 * (*i).z() * (*i).z();
//        A(2, 3) += -4 * (*i).z();
//        A(3, 0) += -4 * (*i).x();
//        A(3, 1) += -4 * (*i).y();
//        A(3, 2) += -4 * (*i).z();
//        A(3, 3) +=  2;
//        
//        b(0) +=  4 * (*i).squaredNorm() * (*i).x();
//        b(1) +=  4 * (*i).squaredNorm() * (*i).y();
//        b(2) +=  4 * (*i).squaredNorm() * (*i).z();
//        b(3) += -2 * (*i).squaredNorm();
//      }
//      
//      VecXd x = A.fullPivLu().solve(b);
//      Scalar r = sqrt(x.segment<3>(0).dot(x.segment<3>(0)) - x.w());
//      
//      Vec3d v;
//      v.setZero();
//      for (std::set<VertexHandle>::iterator i = vts[j].begin(); i != vts[j].end(); i++)
//        v += shell->getVertexVelocity(*i);
//      v /= vts[j].size();
//      
//      std::cout << "Bubble " << j << " (" << pos[j].size() << ") : center = " << x.segment<3>(0).transpose() << " radius = " << r << " average velocity = " << v << std::endl;
//    }
//  } else if (m_active_scene == 5 || m_active_scene == 6)
  {
//    Scalar minz = 1;
//    VertexHandle minzv;
//    for (VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit)
//    {
//      bool foi = false;
//      for (VertexFaceIterator vfit = shellObj->vf_iter(*vit); vfit; ++vfit)
//      {
//        Vec2i label = shell->getFaceLabel(*vfit);
//        if (label.x() == 6 || label.y() == 6)
//          foi = true;
//      }
//      if (foi)
//        if (shell->getVertexPosition(*vit).z() < minz)
//        {
//          minz = shell->getVertexPosition(*vit).z();
//          minzv = *vit;
//        }
//    }
//    
//    std::cout << "min z = " << minz << " vertex = " << minzv.idx() << std::endl;
//    
//    shellObj->releaseAllVertices();
//    
//    for (VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit)
//    {
//      VertexHandle v = *vit;
//      
//      bool boundary = false;
//      for (VertexFaceIterator vfit = shellObj->vf_iter(v); vfit; ++vfit)
//      {
//        FaceHandle f = *vfit;
//        Vec2i labels = shell->getFaceLabel(f);
//        if (labels.x() < 0 || labels.y() < 0)
//          boundary = true;
//      }
//      
//      Vec3d pos = shell->getVertexPosition(v);
//      if (boundary)
//      {
//        bool x = false;
//        bool y = false;
//        bool z = false;
//        if (pos.x() < 1e-4)
//        {
//          pos.x() = 0;
//          x = true;
//        }
//        if (pos.x() > 1 - 1e-4)
//        {
//          pos.x() = 1;
//          x = true;
//        }
//        if (pos.y() < 1e-4)
//        {
//          pos.y() = 0;
//          y = true;
//        }
//        if (pos.y() > 1 - 1e-4)
//        {
//          pos.y() = 1;
//          y = true;
//        }
//        if (pos.z() < 1e-4)
//        {
//          pos.z() = 0;
//          z = true;
//        }
//        if (pos.z() > 1 - 1e-4)
//        {
//          pos.z() = 1;
//          z = true;
//        }
//        
//        assert(x || y || z);
//        PositionConstraint * pc = new PartialPositionConstraint(pos, x, y, z);
//        shellObj->constrainVertex(v, pc);
//      }
//    }

    shellObj->releaseAllVertices();
    shell->getVertexConstraintLabels().assign(0);
    
    for (VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit)
    {
      VertexHandle v = *vit;

      Vec3d pos = shell->getVertexPosition(v);
      int constraint = 0;
      bool x = false;
      bool y = false;
      bool z = false;
      if (pos.x() < 1e-4)
      {
        pos.x() = 0;
        x = true;
        constraint |= (1 << 0);
      }
      if (pos.x() > 1 - 1e-4)
      {
        pos.x() = 1;
        x = true;
        constraint |= (1 << 3);
      }
      if (pos.y() < 1e-4)
      {
        pos.y() = 0;
        y = true;
        constraint |= (1 << 1);
      }
      if (pos.y() > 1 - 1e-4)
      {
        pos.y() = 1;
        y = true;
        constraint |= (1 << 4);
      }
      if (pos.z() < 1e-4)
      {
        pos.z() = 0;
        z = true;
        constraint |= (1 << 2);
      }
      if (pos.z() > 1 - 1e-4)
      {
        pos.z() = 1;
        z = true;
        constraint |= (1 << 5);
      }

      if (x || y || z)
      {
        PositionConstraint * pc = new PartialPositionConstraint(pos, x, y, z);
        shellObj->constrainVertex(v, pc);
        shell->getVertexConstraintLabel(v) = constraint;
      }
    }
    
  }
}

void DoubleBubbleTest::AfterStep()
{
//  triangulation_added_vertices.clear();
//  triangulation_added_edges.clear();
//  triangulation_added_faces.clear();
//  svf->triangulateBBWalls(triangulation_added_vertices, triangulation_added_edges, triangulation_added_faces);
}

void DoubleBubbleTest::setupScene1() 
{
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
  
  for(int i = 0; i < 4; ++i) 
  {
    vertList.push_back(shellObj->addVertex());
    velocities[vertList[i]] = Vec3d(0,0,0);
  }
  
  //create positions
  positions[vertList[ 0]] = Vec3d(0.3,1,0);
  positions[vertList[ 1]] = Vec3d(0.3,1,1);
  positions[vertList[ 2]] = Vec3d(1,0.3,0);
  positions[vertList[ 3]] = Vec3d(1,0.3,1);
//  positions[vertList[ 0]] = Vec3d(0.5,1,0);
//  positions[vertList[ 1]] = Vec3d(0.5,1,1);
//  positions[vertList[ 2]] = Vec3d(0.5,0,0);
//  positions[vertList[ 3]] = Vec3d(0.5,0,1);
  
  for(int i = 0; i < 4; ++i)
    undeformed[vertList[i]] = positions[vertList[i]];
  
  std::vector<FaceHandle> faceList;
  FaceProperty<Vec2i> faceLabels(shellObj); //label face regions to do volume constrained bubbles  
  
  faceList.push_back(shellObj->addFace(vertList[ 0], vertList[ 1], vertList[ 3]));  faceLabels[faceList.back()] = Vec2i(0,1);
  faceList.push_back(shellObj->addFace(vertList[ 3], vertList[ 2], vertList[ 0]));  faceLabels[faceList.back()] = Vec2i(0,1);
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(shellObj); 
  DeformableObject::face_iter fIt;
  for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(shellObj, shellFaces, m_timestep);
  shellObj->addModel(shell);
  
  //positions
  //  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  //mid-edge normal variables
  //  shell->setEdgeUndeformed(undefAngle);
  //  shell->setEdgeXis(edgeAngle);
  //  shell->setEdgeVelocities(edgeVel);
  
  shell->setFaceLabels(faceLabels);
  
}

void DoubleBubbleTest::setupScene2()
{
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
  
  for(int i = 0; i < 20; ++i) 
  {
    vertList.push_back(shellObj->addVertex());
    velocities[vertList[i]] = Vec3d(0,0,0);
  }
  
  //create positions
  positions[vertList[ 0]] = Vec3d(0,0,0);
  positions[vertList[ 1]] = Vec3d(0,0,1);
  positions[vertList[ 2]] = Vec3d(0,1,0);
  positions[vertList[ 3]] = Vec3d(0,1,1);
  
  positions[vertList[ 4]] = Vec3d(1,0,0);
  positions[vertList[ 5]] = Vec3d(1,0,1);
  positions[vertList[ 6]] = Vec3d(1,1,0);
  positions[vertList[ 7]] = Vec3d(1,1,1);
  
  positions[vertList[ 8]] = Vec3d(0.3,0,0);
  positions[vertList[ 9]] = Vec3d(0.3,0,1);
  positions[vertList[10]] = Vec3d(0.3,1,0);
  positions[vertList[11]] = Vec3d(0.3,1,1);
  
  positions[vertList[12]] = Vec3d(0.7,0,0);
  positions[vertList[13]] = Vec3d(0.7,0,1);
  positions[vertList[14]] = Vec3d(0.7,1,0);
  positions[vertList[15]] = Vec3d(0.7,1,1);
  
  positions[vertList[16]] = Vec3d(0.3,0.5,0);
  positions[vertList[17]] = Vec3d(0.3,0.5,1);
  positions[vertList[18]] = Vec3d(0.7,0.5,0);
  positions[vertList[19]] = Vec3d(0.7,0.5,1);
  
  for(int i = 0; i < 20; ++i)
    undeformed[vertList[i]] = positions[vertList[i]];
  
  std::vector<FaceHandle> faceList;
  FaceProperty<Vec2i> faceLabels(shellObj); //label face regions to do volume constrained bubbles  
    
  // face with x = 0.3
  faceList.push_back(shellObj->addFace(vertList[ 8], vertList[16], vertList[17]));  faceLabels[faceList.back()] = Vec2i(0,1);
  faceList.push_back(shellObj->addFace(vertList[ 8], vertList[17], vertList[ 9]));  faceLabels[faceList.back()] = Vec2i(0,1);
  faceList.push_back(shellObj->addFace(vertList[16], vertList[10], vertList[11]));  faceLabels[faceList.back()] = Vec2i(0,3);
  faceList.push_back(shellObj->addFace(vertList[16], vertList[11], vertList[17]));  faceLabels[faceList.back()] = Vec2i(0,3);
  
  // face with x = 0.7
  faceList.push_back(shellObj->addFace(vertList[18], vertList[14], vertList[15]));  faceLabels[faceList.back()] = Vec2i(3,2);
  faceList.push_back(shellObj->addFace(vertList[18], vertList[15], vertList[19]));  faceLabels[faceList.back()] = Vec2i(3,2);
  faceList.push_back(shellObj->addFace(vertList[12], vertList[18], vertList[19]));  faceLabels[faceList.back()] = Vec2i(1,2);
  faceList.push_back(shellObj->addFace(vertList[12], vertList[19], vertList[13]));  faceLabels[faceList.back()] = Vec2i(1,2);
  
  // face with y = 0.5
  faceList.push_back(shellObj->addFace(vertList[16], vertList[17], vertList[19]));  faceLabels[faceList.back()] = Vec2i(1,3);
  faceList.push_back(shellObj->addFace(vertList[19], vertList[18], vertList[16]));  faceLabels[faceList.back()] = Vec2i(1,3);
  
  // remove unused vertices
  shellObj->deleteVertex(VertexHandle(vertList[0]));
  shellObj->deleteVertex(VertexHandle(vertList[1]));
  shellObj->deleteVertex(VertexHandle(vertList[2]));
  shellObj->deleteVertex(VertexHandle(vertList[3]));
  shellObj->deleteVertex(VertexHandle(vertList[4]));
  shellObj->deleteVertex(VertexHandle(vertList[5]));
  shellObj->deleteVertex(VertexHandle(vertList[6]));
  shellObj->deleteVertex(VertexHandle(vertList[7]));
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(shellObj); 
  DeformableObject::face_iter fIt;
  for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(shellObj, shellFaces, m_timestep);
  shellObj->addModel(shell);
  
  //positions
  //  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  //mid-edge normal variables
  //  shell->setEdgeUndeformed(undefAngle);
  //  shell->setEdgeXis(edgeAngle);
  //  shell->setEdgeVelocities(edgeVel);
  
  shell->setFaceLabels(faceLabels);
  
}

void DoubleBubbleTest::setupScene5()
{
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
  
  for(int i = 0; i < 20; ++i) 
  {
    vertList.push_back(shellObj->addVertex());
    velocities[vertList[i]] = Vec3d(0,0,0);
  }
  
  //create positions
  positions[vertList[ 0]] = Vec3d(0,0,0);
  positions[vertList[ 1]] = Vec3d(0,0,1);
  positions[vertList[ 2]] = Vec3d(0,1,0);
  positions[vertList[ 3]] = Vec3d(0,1,1);
  
  positions[vertList[ 4]] = Vec3d(1,0,0);
  positions[vertList[ 5]] = Vec3d(1,0,1);
  positions[vertList[ 6]] = Vec3d(1,1,0);
  positions[vertList[ 7]] = Vec3d(1,1,1);
  
  positions[vertList[ 8]] = Vec3d(0.3,0,0);
  positions[vertList[ 9]] = Vec3d(0.3,0,1);
  positions[vertList[10]] = Vec3d(0.3,1,0);
  positions[vertList[11]] = Vec3d(0.3,1,1);
  
  positions[vertList[12]] = Vec3d(0.3,0.3,0);
  positions[vertList[13]] = Vec3d(0.3,0.3,1);
  positions[vertList[14]] = Vec3d(1,0.3,0);
  positions[vertList[15]] = Vec3d(1,0.3,1);
  
  positions[vertList[16]] = Vec3d(0.3,0.3,0.5);
  positions[vertList[17]] = Vec3d(0.3,1,0.5);
  positions[vertList[18]] = Vec3d(1,0.3,0.5);
  positions[vertList[19]] = Vec3d(1,1,0.5);
  
  for(int i = 0; i < 20; ++i)
    undeformed[vertList[i]] = positions[vertList[i]];
  
  std::vector<FaceHandle> faceList;
  FaceProperty<Vec2i> faceLabels(shellObj); //label face regions to do volume constrained bubbles  

//  // face with z = 0
//  faceList.push_back(shellObj->addFace(vertList[ 0], vertList[12], vertList[ 8]));  faceLabels[faceList.back()] = Vec2i(0,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 0], vertList[ 2], vertList[12]));  faceLabels[faceList.back()] = Vec2i(0,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 2], vertList[10], vertList[12]));  faceLabels[faceList.back()] = Vec2i(0,-1);
//  faceList.push_back(shellObj->addFace(vertList[12], vertList[10], vertList[14]));  faceLabels[faceList.back()] = Vec2i(3,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 6], vertList[14], vertList[10]));  faceLabels[faceList.back()] = Vec2i(3,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 8], vertList[12], vertList[ 4]));  faceLabels[faceList.back()] = Vec2i(1,-1);
//  faceList.push_back(shellObj->addFace(vertList[14], vertList[ 4], vertList[12]));  faceLabels[faceList.back()] = Vec2i(1,-1);
  
//  // face with y = 1
//  faceList.push_back(shellObj->addFace(vertList[11], vertList[17], vertList[ 3]));  faceLabels[faceList.back()] = Vec2i(0,-1);
//  faceList.push_back(shellObj->addFace(vertList[17], vertList[ 2], vertList[ 3]));  faceLabels[faceList.back()] = Vec2i(0,-1);
//  faceList.push_back(shellObj->addFace(vertList[17], vertList[10], vertList[ 2]));  faceLabels[faceList.back()] = Vec2i(0,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 7], vertList[19], vertList[11]));  faceLabels[faceList.back()] = Vec2i(2,-1);
//  faceList.push_back(shellObj->addFace(vertList[17], vertList[11], vertList[19]));  faceLabels[faceList.back()] = Vec2i(2,-1);
//  faceList.push_back(shellObj->addFace(vertList[19], vertList[ 6], vertList[17]));  faceLabels[faceList.back()] = Vec2i(3,-1);
//  faceList.push_back(shellObj->addFace(vertList[10], vertList[17], vertList[ 6]));  faceLabels[faceList.back()] = Vec2i(3,-1);
  
//  // face with z = 1
//  faceList.push_back(shellObj->addFace(vertList[11], vertList[ 3], vertList[13]));  faceLabels[faceList.back()] = Vec2i(0,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 3], vertList[ 1], vertList[13]));  faceLabels[faceList.back()] = Vec2i(0,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 9], vertList[13], vertList[ 1]));  faceLabels[faceList.back()] = Vec2i(0,-1);
//  faceList.push_back(shellObj->addFace(vertList[11], vertList[13], vertList[ 7]));  faceLabels[faceList.back()] = Vec2i(2,-1);
//  faceList.push_back(shellObj->addFace(vertList[15], vertList[ 7], vertList[13]));  faceLabels[faceList.back()] = Vec2i(2,-1);
//  faceList.push_back(shellObj->addFace(vertList[13], vertList[ 9], vertList[15]));  faceLabels[faceList.back()] = Vec2i(1,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 5], vertList[15], vertList[ 9]));  faceLabels[faceList.back()] = Vec2i(1,-1);
  
//  // face with y = 0
//  faceList.push_back(shellObj->addFace(vertList[ 0], vertList[ 8], vertList[ 1]));  faceLabels[faceList.back()] = Vec2i(0,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 9], vertList[ 1], vertList[ 8]));  faceLabels[faceList.back()] = Vec2i(0,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 8], vertList[ 4], vertList[ 9]));  faceLabels[faceList.back()] = Vec2i(1,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 5], vertList[ 9], vertList[ 4]));  faceLabels[faceList.back()] = Vec2i(1,-1);
  
//  // face with x = 0
//  faceList.push_back(shellObj->addFace(vertList[ 3], vertList[ 2], vertList[ 1]));  faceLabels[faceList.back()] = Vec2i(0,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 0], vertList[ 1], vertList[ 2]));  faceLabels[faceList.back()] = Vec2i(0,-1);
  
//  // face with x = 1
//  faceList.push_back(shellObj->addFace(vertList[ 4], vertList[14], vertList[18]));  faceLabels[faceList.back()] = Vec2i(1,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 4], vertList[18], vertList[ 5]));  faceLabels[faceList.back()] = Vec2i(1,-1);
//  faceList.push_back(shellObj->addFace(vertList[18], vertList[15], vertList[ 5]));  faceLabels[faceList.back()] = Vec2i(1,-1);
//  faceList.push_back(shellObj->addFace(vertList[14], vertList[ 6], vertList[18]));  faceLabels[faceList.back()] = Vec2i(3,-1);
//  faceList.push_back(shellObj->addFace(vertList[19], vertList[18], vertList[ 6]));  faceLabels[faceList.back()] = Vec2i(3,-1);
//  faceList.push_back(shellObj->addFace(vertList[18], vertList[19], vertList[15]));  faceLabels[faceList.back()] = Vec2i(2,-1);
//  faceList.push_back(shellObj->addFace(vertList[ 7], vertList[15], vertList[19]));  faceLabels[faceList.back()] = Vec2i(2,-1);
  
  // face with x = 0.3
  faceList.push_back(shellObj->addFace(vertList[ 8], vertList[12], vertList[16]));  faceLabels[faceList.back()] = Vec2i(0,1);
  faceList.push_back(shellObj->addFace(vertList[ 8], vertList[16], vertList[ 9]));  faceLabels[faceList.back()] = Vec2i(0,1);
  faceList.push_back(shellObj->addFace(vertList[13], vertList[ 9], vertList[16]));  faceLabels[faceList.back()] = Vec2i(0,1);
  faceList.push_back(shellObj->addFace(vertList[12], vertList[10], vertList[16]));  faceLabels[faceList.back()] = Vec2i(0,3);
  faceList.push_back(shellObj->addFace(vertList[17], vertList[16], vertList[10]));  faceLabels[faceList.back()] = Vec2i(0,3);
  faceList.push_back(shellObj->addFace(vertList[16], vertList[17], vertList[13]));  faceLabels[faceList.back()] = Vec2i(0,2);
  faceList.push_back(shellObj->addFace(vertList[11], vertList[13], vertList[17]));  faceLabels[faceList.back()] = Vec2i(0,2);
  
  // face with y = 0.3
  faceList.push_back(shellObj->addFace(vertList[18], vertList[14], vertList[16]));  faceLabels[faceList.back()] = Vec2i(1,3);
  faceList.push_back(shellObj->addFace(vertList[12], vertList[16], vertList[14]));  faceLabels[faceList.back()] = Vec2i(1,3);
  faceList.push_back(shellObj->addFace(vertList[15], vertList[18], vertList[13]));  faceLabels[faceList.back()] = Vec2i(1,2);
  faceList.push_back(shellObj->addFace(vertList[16], vertList[13], vertList[18]));  faceLabels[faceList.back()] = Vec2i(1,2);
  
  // face with z = 0.5
  faceList.push_back(shellObj->addFace(vertList[17], vertList[16], vertList[19]));  faceLabels[faceList.back()] = Vec2i(3,2);
  faceList.push_back(shellObj->addFace(vertList[18], vertList[19], vertList[16]));  faceLabels[faceList.back()] = Vec2i(3,2);
  
  // remove unused vertices
  shellObj->deleteVertex(VertexHandle(vertList[0]));
  shellObj->deleteVertex(VertexHandle(vertList[1]));
  shellObj->deleteVertex(VertexHandle(vertList[2]));
  shellObj->deleteVertex(VertexHandle(vertList[3]));
  shellObj->deleteVertex(VertexHandle(vertList[4]));
  shellObj->deleteVertex(VertexHandle(vertList[5]));
  shellObj->deleteVertex(VertexHandle(vertList[6]));
  shellObj->deleteVertex(VertexHandle(vertList[7]));
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(shellObj); 
  DeformableObject::face_iter fIt;
  for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(shellObj, shellFaces, m_timestep);
  shellObj->addModel(shell);
  
  //positions
//  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  //mid-edge normal variables
//  shell->setEdgeUndeformed(undefAngle);
//  shell->setEdgeXis(edgeAngle);
//  shell->setEdgeVelocities(edgeVel);
  
  shell->setFaceLabels(faceLabels);
  
}

void DoubleBubbleTest::setupScene6()
{
  //vertices
  std::vector<VertexHandle> vertHandles;
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);
  
  //generate voronoi sites
  int nsite = 10;
  srand(100000);
  std::vector<Vec3d> sites;
  for (int i = 0; i < nsite; i++)
    sites.push_back(Vec3d((Scalar)rand() / RAND_MAX, (Scalar)rand() / RAND_MAX, (Scalar)rand() / RAND_MAX));
 
  //create an initial tet that contains the entire range [0, 1]^3
  TopologicalObject dtmesh;
  VertexProperty<Vec3d> dtpositions(&dtmesh);
  std::vector<VertexHandle> dtvertList;
  
  for (int i = 0; i < 4; i++)
    dtvertList.push_back(dtmesh.addVertex());
  
  dtpositions[dtvertList[0]] = Vec3d(-1,-1,-1);
  dtpositions[dtvertList[1]] = Vec3d(4,0,0);
  dtpositions[dtvertList[2]] = Vec3d(0,4,0);
  dtpositions[dtvertList[3]] = Vec3d(0,0,4);
  
  std::vector<FaceHandle> dtfaceList;
  dtfaceList.push_back(dtmesh.addFace(dtvertList[0], dtvertList[1], dtvertList[2])); // 0  // z = 0
  dtfaceList.push_back(dtmesh.addFace(dtvertList[0], dtvertList[1], dtvertList[3])); // 1  // y = 0
  dtfaceList.push_back(dtmesh.addFace(dtvertList[0], dtvertList[2], dtvertList[3])); // 2  // x = 0
  dtfaceList.push_back(dtmesh.addFace(dtvertList[1], dtvertList[2], dtvertList[3])); // 3  // x+y+z
  
  std::vector<TetHandle> dttetList;
  dttetList.push_back(dtmesh.addTet(dtfaceList[0], dtfaceList[1], dtfaceList[2], dtfaceList[3]));  
  
  DelaunayTriangulator::DelaunayTriangulator dt(&dtmesh, dtpositions);
  
  for (int i = 0; i < nsite; i++)
    dt.insertVertex(sites[i]);
  dtpositions = dt.positions();

  FaceProperty<Vec2i> faceLabels(shellObj);
  dt.extractVoronoiDiagram(shellObj, positions, faceLabels, Vec3d(0, 0, 0), Vec3d(1, 1, 1));
  
  for (VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit)
  {
    undeformed[*vit] = positions[*vit];
    velocities[*vit] = Vec3d(0, 0, 0);
  }
  
  // remove the bounding box wall faces
  std::vector<FaceHandle> faces_to_remove;
  for (FaceIterator fit = shellObj->faces_begin(); fit != shellObj->faces_end(); ++fit)
  {
    FaceVertexIterator fvit = shellObj->fv_iter(*fit); assert(fvit);
    VertexHandle v0 = *fvit; ++fvit; assert(fvit);
    VertexHandle v1 = *fvit; ++fvit; assert(fvit);
    VertexHandle v2 = *fvit; ++fvit; assert(!fvit);
    
    int wall0 = onBBWall(positions[v0]);
    int wall1 = onBBWall(positions[v1]);
    int wall2 = onBBWall(positions[v2]);
    
    if (wall0 & wall1 & wall2)
      faces_to_remove.push_back(*fit);
  }
  
  for (size_t i = 0; i < faces_to_remove.size(); i++)
    shellObj->deleteFace(faces_to_remove[i], true);
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(shellObj); 
  DeformableObject::face_iter fIt;
  for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(shellObj, shellFaces, m_timestep);
  shellObj->addModel(shell);
  
  //positions
//  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  //mid-edge normal variables
//  shell->setEdgeUndeformed(undefAngle);
//  shell->setEdgeXis(edgeAngle);
//  shell->setEdgeVelocities(edgeVel);
  
  shell->setFaceLabels(faceLabels);
  
 
}