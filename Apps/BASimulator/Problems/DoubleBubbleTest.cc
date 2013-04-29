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
#include "ElTopo/eltopo3d/iomesh.h"

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
Recording g_recording;

void Recording::writeSurfTrack(std::ostream & os, const ElTopo::SurfTrack & st)
{
  size_t n;
  n = st.m_mesh.nv();
  os.write((char *)&n, sizeof (size_t));
  for (size_t i = 0; i < n; i++)
  {
    ElTopo::Vec3d x = st.get_position(i);
    os.write((char *)&(x[0]), sizeof (x[0]));
    os.write((char *)&(x[1]), sizeof (x[1]));
    os.write((char *)&(x[2]), sizeof (x[2]));
    
//    x = st.get_newposition(i);
//    os.write((char *)&(x[0]), sizeof (x[0]));
//    os.write((char *)&(x[1]), sizeof (x[1]));
//    os.write((char *)&(x[2]), sizeof (x[2]));
//    
//    x = st.get_remesh_velocity(i);
//    os.write((char *)&(x[0]), sizeof (x[0]));
//    os.write((char *)&(x[1]), sizeof (x[1]));
//    os.write((char *)&(x[2]), sizeof (x[2]));
//    
//    double mass = st.m_masses[i];
//    os.write((char *)&mass, sizeof (mass));
//    
//    bool cl = st.m_mesh.m_vertex_constraint_labels[i];
//    os.write((char *)&cl, sizeof (cl));
  }
  
  n = st.m_mesh.nt();
  os.write((char *)&n, sizeof (size_t));
  for (size_t i = 0; i < n; i++)
  {
    ElTopo::Vec3st t = st.m_mesh.get_triangle(i);
    os.write((char *)&(t[0]), sizeof (t[0]));
    os.write((char *)&(t[1]), sizeof (t[1]));
    os.write((char *)&(t[2]), sizeof (t[2]));
    
    ElTopo::Vec2i l = st.m_mesh.get_triangle_label(i);
    os.write((char *)&(l[0]), sizeof (l[0]));
    os.write((char *)&(l[1]), sizeof (l[1]));
  }
  
}

void Recording::readSurfTrack(std::istream & is, ElTopo::SurfTrack & st)
{
  // clear the mesh
  for (size_t i = 0; i < st.m_mesh.nt(); i++)
  {
    if (st.m_mesh.get_triangle(i)[0] == st.m_mesh.get_triangle(i)[1])
      continue;
    st.remove_triangle(i);
  }
  
  for (size_t i = 0; i < st.m_mesh.nv(); i++)
    st.remove_vertex(i);
  
  size_t n;
  n = st.m_mesh.nv();
  is.read((char *)&n, sizeof (size_t));
//  std::cout << " nv = " << n << std::endl;
  st.m_mesh.set_num_vertices(n);
  std::vector<ElTopo::Vec3d> pos(n);
  for (size_t i = 0; i < n; i++)
  {
    ElTopo::Vec3d x;
    is.read((char *)&(x[0]), sizeof (x[0]));
    is.read((char *)&(x[1]), sizeof (x[1]));
    is.read((char *)&(x[2]), sizeof (x[2]));
    pos[i] = x;
    
//    is.read((char *)&(x[0]), sizeof (x[0]));
//    is.read((char *)&(x[1]), sizeof (x[1]));
//    is.read((char *)&(x[2]), sizeof (x[2]));
//    st.set_newposition(i, x);
//    
//    is.read((char *)&(x[0]), sizeof (x[0]));
//    is.read((char *)&(x[1]), sizeof (x[1]));
//    is.read((char *)&(x[2]), sizeof (x[2]));
//    st.set_remesh_velocity(i, x);
//    
//    double mass;
//    is.read((char *)&mass, sizeof (mass));
//    st.m_masses[i] = mass;
//    
//    bool cl;
//    is.read((char *)&cl, sizeof (cl));
//    st.m_mesh.m_vertex_constraint_labels[i] = cl;
  }
  
  st.m_masses.resize(n);
  for (size_t i = 0; i < n; i++)
    st.m_masses[i] = ElTopo::Vec3d(1, 1, 1);
  
  st.set_all_positions(pos);
  st.set_all_newpositions(pos);  
  st.set_all_remesh_velocities(std::vector<ElTopo::Vec3d>(n, ElTopo::Vec3d(0)));
  
//  std::cout << "nv: " << pos.size() << " " << st.pm_velocities.size() << " " << st.m_mesh.m_vertex_constraint_labels.size() << " " << st.m_mesh.m_vertex_to_triangle_map.size() << std::endl;
//  std::cout << "nv: " << st.m_mesh.nv() << std::endl;
  
  n = st.m_mesh.nt();
  is.read((char *)&n, sizeof (size_t));
//  std::cout << "nt = " << n << std::endl;
  std::vector<ElTopo::Vec3st> tris;
  std::vector<ElTopo::Vec2i> labels;
  for (size_t i = 0; i < n; i++)
  {
    ElTopo::Vec3st t;
    is.read((char *)&(t[0]), sizeof (t[0]));
    is.read((char *)&(t[1]), sizeof (t[1]));
    is.read((char *)&(t[2]), sizeof (t[2]));
    tris.push_back(t);
    
    ElTopo::Vec2i l;
    is.read((char *)&(l[0]), sizeof (l[0]));
    is.read((char *)&(l[1]), sizeof (l[1]));
    labels.push_back(l);
  }
  
  st.m_mesh.replace_all_triangles(tris, labels);
//  for (size_t i = 0; i < n; i++)
//    st.m_mesh.set_triangle_label(i, labels[i]);
  
  size_t nv = st.m_mesh.m_vertex_to_triangle_map.size();
  st.pm_positions.resize(nv);
  st.pm_newpositions.resize(nv);
  st.pm_velocities.resize(nv);
  st.m_velocities.resize(nv);

//  std::cout << "nv: " << pos.size() << " " << st.pm_velocities.size() << " " << st.m_mesh.m_vertex_constraint_labels.size() << " " << st.m_mesh.m_vertex_to_triangle_map.size() << std::endl;
//  std::cout << "nv: " << st.m_mesh.nv() << std::endl;
  
}

void Recording::recordSurfTrack(const ElTopo::SurfTrack & st)
{
  if (!isRecording())
    return;
  
  if (!m_of.is_open())
  {
    std::stringstream filename;
    filename << m_recording_name << "_" << m_current_frame << ".rec";
    m_of.open(filename.str().c_str());
      
    if (!m_of.is_open())
    {
      std::cout << "Cannot open recording frame file " << filename.str() << std::endl;
      return;
    }
  }
  
  m_of.write((char *)&m_current_step, sizeof (m_current_step));
  
  size_t loglen = m_log.str().size();
  m_of.write((char *)&loglen, sizeof (size_t));
  m_of.write((char *)m_log.str().c_str(), loglen);
  m_log.str("");
  
  writeSurfTrack(m_of, st);
  
  m_current_step++;
}

void Recording::loadRecording(ElTopo::SurfTrack & st, int next)
{
  if (!isPlaybackOn())
    return;
  
  if (!m_if.is_open())
  {
    std::stringstream filename;
    filename << m_recording_name << "_" << m_current_frame << ".rec";
    m_if.open(filename.str().c_str());
    
    if (!m_if.is_open())
    {
      std::cout << "Requested recording frame not found!" << std::endl;
      return;
    }
    
    m_step_pos.clear();
    m_step_pos.push_back(m_if.tellg());
    while (!m_if.eof())
    {
      int step = 0;
      m_if.read((char *)&step, sizeof (step));
      /////////
//      readSurfTrack(m_if, st);
      /////////
      size_t n;
      m_if.read((char *)&n, sizeof (n));  // skip log
      m_if.seekg(n, ios_base::cur);
      m_if.read((char *)&n, sizeof (n));  // skip vertices
      m_if.seekg(sizeof (double) * 3 * n, ios_base::cur);
      m_if.read((char *)&n, sizeof (n));  // skip faces
      m_if.seekg((sizeof (size_t) * 3 + sizeof (int) * 2) * n, ios_base::cur);
      /////////
      m_if.peek();
      if (!m_if.eof())
        m_step_pos.push_back(m_if.tellg());
    }
    std::cout << "Recording file " << filename << " contains " << m_step_pos.size() << " steps." << std::endl;
    
    m_if.seekg(m_step_pos.front());
    m_if.clear();
  }
  
  assert(m_current_step < m_step_pos.size());
//  std::cout << "Loading step " << (m_current_step + next) % m_step_pos.size() << std::endl;
  m_if.seekg(m_step_pos[(m_current_step + m_step_pos.size() + next) % m_step_pos.size()]);
  m_if.read((char *)&m_current_step, sizeof (m_current_step));
  
  size_t loglen;
  m_if.read((char *)&loglen, sizeof (size_t));
  char * log = new char[loglen + 1];
  m_if.read((char *)log, loglen);
  *(log + loglen) = NULL;
  
  readSurfTrack(m_if, st);

  m_if.peek();
  if (m_if.eof())
    m_if.close();
  
  std::cout << "Loaded recording: step " << m_current_step << "/" << m_step_pos.size() << " of frame " << m_current_frame << std::endl;
  std::cout << log << std::endl;
}

DoubleBubbleTest::DoubleBubbleTest() : 
  Problem("Double Bubble Test", "Various bubble dynamics tests"), 
  shell(NULL), 
  shellObj(NULL), 
  stepper(NULL),
  svf(NULL)
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
  
  AddOption("volume-targeting", "whether to set the desired volume of each region to equal", false);
    
  //sheared wrinkling parameters
  AddOption("shell-rotation-rate", "the rate at which inner cylinder rotates for sheared wrinkling test", 0.0);
  AddOption("shell-bath-density", "the density of water fluid bath for sheared wrinkling test", 1.0);

  //Remeshing options
  AddOption("shell-remeshing", "whether to perform remeshing", false);
  AddOption("shell-remeshing-resolution", "target edge-length", 0.0); //for backwards compatibility
  AddOption("shell-remeshing-max-length", "upper bound on edge-length", 0.5);
  AddOption("shell-remeshing-min-length", "lower bound on edge-length", 0.1);
  AddOption("shell-remeshing-iterations", "number of remeshing iterations to run", 2);
  
  AddOption("shell-init-remesh", "whether or not to run a remeshing pass before simulation starts", false);

  AddOption("t1-transition-enabled", "whether T1 transition operations are enabled", false);
  AddOption("smooth-subdivision-scheme", "whether or not to use the Modified Butterfly subdivision scheme (default is false, i.e. midpoint subdivision)", false);
    
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
  AddOption("obj-mode", "0 = one OBJ for entire mesh; 1 = one OBJ per region; 2 = one OBJ per region pair; 3 = regions specified in parameter 'obj-regions', plus an OBJ for the entire mesh with the interfaces involving those specified regions excluded", 0);
  AddOption("obj-regions", "a string listing region indices of the desired regions, separated by spaces", "");
  AddOption("record", "Generate a recording", false);
  AddOption("playback", "Playback a recording", false);
  AddOption("playback-path", "The path to the recording to playback", "");

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
  0, //  &DoubleBubbleTest::setupScene4,
  &DoubleBubbleTest::setupScene5,
  &DoubleBubbleTest::setupScene6,
  &DoubleBubbleTest::setupScene7,
  &DoubleBubbleTest::setupScene8,
  &DoubleBubbleTest::setupScene9,
  &DoubleBubbleTest::setupScene10,
  &DoubleBubbleTest::setupScene11,
  &DoubleBubbleTest::setupScene12,
  &DoubleBubbleTest::setupScene13,
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
    if (GetBoolOpt("record"))
    {
        assert(!GetBoolOpt("playback"));
        
        g_recording.turnOnRecording();
        g_recording.setRecordingName(outputdirectory + "/rec");
    }
    
    if (GetBoolOpt("playback"))
    {
        assert(!g_obj_dump);
        assert(!g_ply_dump);
        assert(!g_recording.isRecording());
        
        g_recording.turnOnPlayback();
        
        std::string playback_path = GetStringOpt("playback-path");
        assert(playback_path != "");
        
        g_recording.setRecordingName(playback_path + "/rec");
    }
    
    g_obj_dump = GetBoolOpt("generate-OBJ");
    g_ply_dump = GetBoolOpt("generate-PLY");
    
    if (g_obj_dump || g_ply_dump || g_recording.isRecording()){
#ifdef _MSC_VER
        _mkdir(outputdirectory.c_str());
#else
        mkdir(outputdirectory.c_str(), 0755);
#endif
    }
    
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

  shell->setMeshEventCallback(this);
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
  svf = NULL;
  if (vf_stiffness > 0)
  {
    svf = new ShellVolumeForce(*shell, "Volume", vf_stiffness);
    shell->addForce(svf);
      
    if (GetBoolOpt("volume-targeting"))
    {
      svf->m_target_volumes.assign(svf->m_target_volumes.size(), 1.0 / svf->m_target_volumes.size());
    }
  }
  
//  shell->computeMasses(); /////////////////////////////
    
  shell->m_remesh_smooth_subdivision = GetBoolOpt("smooth-subdivision-scheme");
  shell->m_remesh_t1transition = GetBoolOpt("t1-transition-enabled");
 

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

  updateBBWallConstraints();
  if (GetBoolOpt("shell-init-remesh"))
  {
    shell->remesh(1.0, true);
    updateBBWallConstraints();
  }
    
  // count total number of regions
  int max_region = 0;
  for (FaceIterator fit = shellObj->faces_begin(); fit != shellObj->faces_end(); ++fit)
  {
    Vec2i label = shell->getFaceLabel(*fit);
    if (label.x() > max_region)
      max_region = label.x();
    if (label.y() > max_region)
      max_region = label.y();
  }
    
  m_nregion = max_region + 1;
    
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

ElTopo::SurfTrack * DoubleBubbleTest::mesh2surftrack()
{
  // convert BASim mesh to ElTopo mesh (code copied from ElasticShell::remesh())
  std::vector<ElTopo::Vec3d> vert_data;
  std::vector<ElTopo::Vec3d> vert_vel;
  std::vector<ElTopo::Vec3st> tri_data;
  std::vector<ElTopo::Vec2i> tri_labels;
  std::vector<ElTopo::Vec3d> masses;
  
  DeformableObject& mesh = *shellObj;
  
  //Index mappings between us and El Topo, used in remesh()
  VertexProperty<int> vert_numbers(&mesh);
  FaceProperty<int> face_numbers(&mesh);
  std::vector<VertexHandle> reverse_vertmap;
  std::vector<FaceHandle> reverse_trimap;
  
  reverse_vertmap.reserve(shellObj->nv());
  reverse_trimap.reserve(shellObj->nt());  
  
  vert_data.reserve(shellObj->nv());
  tri_data.reserve(shellObj->nt());
  tri_labels.reserve(shellObj->nt());
  masses.reserve(shellObj->nv());
  
  //walk through vertices, create linear list, store numbering
  int id = 0;
  for(VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit) {
    VertexHandle vh = *vit;
    Vec3d vert = shell->getVertexPosition(vh);
    vert_data.push_back(ElTopo::Vec3d(vert[0], vert[1], vert[2]));
    Vec3d vel = shell->getVertexVelocity(vh);
    vert_vel.push_back(ElTopo::Vec3d(vel[0], vel[1], vel[2]));
    assert(shellObj->isConstrained(vh) == (shell->getVertexConstraintLabel(vh) != 0));
    ElTopo::Vec3d mass(1, 1, 1);
    for (int i = 0; i < 3; i++)
      if (shellObj->isConstrainedInDirection(vh, i))
        mass[i] = numeric_limits<Scalar>::infinity();
    masses.push_back(mass);
    vert_numbers[vh] = id;
    reverse_vertmap.push_back(vh);
    
    ++id;
  }
  
  //walk through tris, creating linear list, using the vertex numbering assigned above
  id = 0;
  for(FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit) {
    FaceHandle fh = *fit;
    ElTopo::Vec3st tri;
    int i = 0;
    for(FaceVertexIterator fvit = mesh.fv_iter(fh); fvit; ++fvit) {
      VertexHandle vh = *fvit;
      tri[i] = vert_numbers[vh];
      ++i;
    }
    tri_data.push_back(tri);
    tri_labels.push_back(ElTopo::Vec2i(shell->getFaceLabel(fh).x(), shell->getFaceLabel(fh).y()));
    face_numbers[fh] = id;
    reverse_trimap.push_back(fh);
    ++id;
  }
  
  ElTopo::SurfTrackInitializationParameters construction_parameters;
  construction_parameters.m_proximity_epsilon = 1e-10;
  construction_parameters.m_merge_proximity_epsilon = 1e-10;
  construction_parameters.m_allow_vertex_movement_during_collapse = true;
  construction_parameters.m_perform_smoothing = false;
  construction_parameters.m_min_edge_length = 0.00001;
  construction_parameters.m_max_edge_length = 1000;
  construction_parameters.m_max_volume_change = numeric_limits<double>::max();   
  construction_parameters.m_min_triangle_angle = 3;
  construction_parameters.m_max_triangle_angle = 177;
  construction_parameters.m_large_triangle_angle_to_split = 160;
  construction_parameters.m_verbose = false;
  construction_parameters.m_allow_non_manifold = true;
  construction_parameters.m_allow_topology_changes = true;
  construction_parameters.m_collision_safety = true;
  construction_parameters.m_remesh_boundaries = true;
  construction_parameters.m_t1_transition_enabled = false;
  
  ElTopo::SurfTrack * st = new ElTopo::SurfTrack( vert_data, tri_data, tri_labels, masses, construction_parameters ); 
  st->m_solid_vertices_callback = shell;
  st->set_all_remesh_velocities(vert_vel);
  
  return st;
}

void DoubleBubbleTest::surftrack2mesh(const ElTopo::SurfTrack & surface_tracker)
{
//  std::cout << "nv 3: " << surface_tracker.pm_positions.size() << " " << surface_tracker.pm_velocities.size() << " " << surface_tracker.m_mesh.m_vertex_constraint_labels.size() << " " << surface_tracker.m_mesh.m_vertex_to_triangle_map.size() << std::endl;
//  std::cout << "nv 3: " << surface_tracker.m_mesh.nv() << std::endl;
  
  for (FaceIterator fit = shellObj->faces_begin(); fit != shellObj->faces_end(); ++fit)
    shellObj->deleteFace(*fit, true);
  
  assert(shellObj->nv() == 0);
  assert(shellObj->ne() == 0);
  assert(shellObj->nf() == 0);
  
  VertexProperty<int> vert_numbers(shellObj);
  FaceProperty<int> face_numbers(shellObj);
  std::vector<VertexHandle> reverse_vertmap;
  std::vector<FaceHandle> reverse_trimap;
  
  reverse_vertmap.resize(surface_tracker.m_mesh.nv());
  reverse_trimap.resize(surface_tracker.m_mesh.nt());
  
  for (size_t i = 0; i < surface_tracker.m_mesh.nv(); i++)
  {
    if (surface_tracker.m_mesh.m_vertex_to_edge_map[i].size() == 0 ||
        surface_tracker.m_mesh.m_vertex_to_triangle_map[i].size() == 0)
    {
      // dead vertex
      reverse_vertmap[i] = VertexHandle();
      continue;
    }
    VertexHandle v = shellObj->addVertex();
    ElTopo::Vec3d x = surface_tracker.get_newposition(i);
    shell->setVertexPosition(v, Vec3d(x[0], x[1], x[2]));
    shell->setVertexVelocity(v, Vec3d(0, 0, 0));
    
    vert_numbers[v] = i;
    reverse_vertmap[i] = v;
  }
  
  for (size_t i = 0; i < surface_tracker.m_mesh.nt(); i++)
  {
    ElTopo::Vec3st tri = surface_tracker.m_mesh.get_triangle(i);
    if (tri[0] == tri[1] || tri[0] == tri[2] || tri[1] == tri[2])
    {
      // dead face
      reverse_trimap[i] = FaceHandle();
      continue; 
    }
    FaceHandle f = shellObj->addFace(reverse_vertmap[tri[0]], reverse_vertmap[tri[1]], reverse_vertmap[tri[2]]);
    shell->setFaceLabel(f, Vec2i(surface_tracker.m_mesh.get_triangle_label(i)[0], surface_tracker.m_mesh.get_triangle_label(i)[1]));
    shell->setFaceActive(f);
    
    face_numbers[f] = i;
    reverse_trimap[i] = f;
  }
  
  shellObj->computeDofIndexing();
  
  for (VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit)
    shell->getVertexConstraintLabel(*vit) = onBBWall(shell->getVertexPosition(*vit));
  
  updateBBWallConstraints();

}


void DoubleBubbleTest::AtEachTimestep()
{
//  for (size_t i = 0; i < triangulation_added_faces.size(); i++)
//    shellObj->deleteFace(triangulation_added_faces[i], false);
//  for (size_t i = 0; i < triangulation_added_edges.size(); i++)
//    shellObj->deleteEdge(triangulation_added_edges[i], false);
//  for (size_t i = 0; i < triangulation_added_vertices.size(); i++)
//    shellObj->deleteVertex(triangulation_added_vertices[i]);

  if (svf)
    svf->globalEnergy();
  
    
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
    if ( g_obj_dump )
    {
        ElTopo::SurfTrack * st = mesh2surftrack();
        ElTopo::SurfTrack & surface_tracker = *st;
        
#ifdef _MSC_VER
        _mkdir(outputdirectory.c_str());
#else
        mkdir(outputdirectory.c_str(), 0755);
#endif

        int RENDER_METHOD = GetIntOpt("obj-mode");
        if (RENDER_METHOD == 0)
        {
            // dump one OBJ for the entire mesh
            std::stringstream name;
            name << std::setfill('0');
            name << outputdirectory << "/" << "mesh_frame" << std::setw(6) << db_current_obj_frame << ".OBJ";
            
            write_objfile(surface_tracker.m_mesh, surface_tracker.get_positions(), name.str().c_str());
            std::cout << "Frame: " << db_current_obj_frame << "   Time: " << getTime() << "   OBJDump: " << name.str() << std::endl;
            
        } else if (RENDER_METHOD == 1)
        {
            // dump one OBJ per region
            std::set<int> excluding_regions;
            for (int i = 0; i < m_nregion; i++)
            {
                std::stringstream name;
                name << std::setfill('0');
                name << outputdirectory << "/" << "region" << std::setw(4) << i << "_frame" << std::setw(6) << db_current_obj_frame << ".OBJ";

                write_objfile_per_region(surface_tracker.m_mesh, surface_tracker.get_positions(), i, excluding_regions, name.str().c_str());
                std::cout << "Frame: " << db_current_obj_frame << "   Time: " << getTime() << "   OBJDump: " << name.str() << std::endl;
                
                excluding_regions.insert(i);
            }
        } else if (RENDER_METHOD == 2)
        {
            // dump one OBJ per region pair
            for (int i = 0; i < m_nregion; i++)
            {
                for (int j = i + 1; j < m_nregion; j++)
                {
                    std::stringstream name;
                    name << std::setfill('0');
                    name << outputdirectory << "/" << "label_" << std::setw(4) << i << "_" << std::setw(4) << j << "_frame" << std::setw(6) << db_current_obj_frame << ".OBJ";
                    
                    write_objfile_per_region_pair(surface_tracker.m_mesh, surface_tracker.get_positions(), ElTopo::Vec2i(i, j), name.str().c_str());
                    std::cout << "Frame: " << db_current_obj_frame << "   Time: " << getTime() << "   OBJDump: " << name.str() << std::endl;
                }
                
            }
            
        } else
        {
            // dump one OBJ per region listed in "obj-regions", plus one OBJ for the rest of the mesh
            std::stringstream regions_ss(GetStringOpt("obj-regions"));
            std::set<int> regions;
            while (!regions_ss.eof())
            {
                int r = -1;
                regions_ss >> r;
                if (r >= 0)
                    regions.insert(r);
            }
            
            std::set<int> excluding_regions;
            for (std::set<int>::iterator i = regions.begin(); i != regions.end(); i++)
            {
                std::stringstream name;
                name << std::setfill('0');
                name << outputdirectory << "/" << "region" << std::setw(4) << *i << "_frame" << std::setw(6) << db_current_obj_frame << ".OBJ";
                
                write_objfile_per_region(surface_tracker.m_mesh, surface_tracker.get_positions(), *i, excluding_regions, name.str().c_str());
                std::cout << "Frame: " << db_current_obj_frame << "   Time: " << getTime() << "   OBJDump: " << name.str() << std::endl;
                excluding_regions.insert(*i);
            }
            
            std::stringstream name;
            name << std::setfill('0');
            name << outputdirectory << "/" << "rest_of_the_mesh_frame" << std::setw(6) << db_current_obj_frame << ".OBJ";
            
            write_objfile_excluding_regions(surface_tracker.m_mesh, surface_tracker.get_positions(), regions, name.str().c_str());
            std::cout << "Frame: " << db_current_obj_frame << "   Time: " << getTime() << "   OBJDump: " << name.str() << std::endl;
          
        }
      
        delete(st);

        ++db_current_obj_frame;
    }

    std::cout << "Time: " << this->getTime() << std::endl; 

  
    g_recording.setCurrentFrame((int)(this->getTime() / this->getDt() + 0.5));

//    if (g_recording.isPlaybackOn())
//    {
//      ElTopo::SurfTrack * st = mesh2surftrack();
//      g_recording.loadRecording(*st);
//      surftrack2mesh(*st);
//      delete st;
//    }
  
}

void DoubleBubbleTest::updateBBWallConstraints()
{
    if (m_active_scene == 9)
      return;
    
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
            shell->setVertexPosition(v, pos);
        }
    }
    
}

void DoubleBubbleTest::afterStartStep()
{
  updateBBWallConstraints();
  
  if (g_recording.isRecording())
  {
    ElTopo::SurfTrack * st = mesh2surftrack();
    g_recording.log() << "Begin time step" << std::endl;
    g_recording.recordSurfTrack(*st);
    delete st;
  }
}

void DoubleBubbleTest::beforeEndStep()
{
  updateBBWallConstraints();
    
  Scalar dt = getDt();
  Scalar current_t = getTime();
  
  if (m_active_scene == 7 || m_active_scene == 13)
  {
    //
    //  RK4 integration of the Enright velocity field
    //  code adapted from El Topo's Enright driver:
    //    https://github.com/tysonbrochu/eltopo/blob/master/talpa/drivers/enrightdriver.h
    //
    
    //Enright test
    for (VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit)
    {
      Vec3d v;
      Vec3d x = shell->getVertexPosition(*vit);
      
      // RK4
      // -----------
      // k1 = dt * f( t, x );
      s7_enright_velocity(current_t, x, v);
      Vec3d k1 = v;
      
      // k2 = dt * f( t + 0.5*dt, x + 0.5*k1 );
      s7_enright_velocity(current_t + 0.5 * dt, x + 0.5 * dt * k1, v);
      Vec3d k2 = v;
      
      // k3 = dt * f( t + 0.5*dt, x + 0.5*k2 );
      s7_enright_velocity(current_t + 0.5 * dt, x + 0.5 * dt * k2, v);
      Vec3d k3 = v;
      
      // k4 = dt * f( t + dt, x + k3 );
      s7_enright_velocity(current_t + dt, x + dt * k3, v);
      Vec3d k4 = v;
      
      v = (1./6. * (k1 + k4) + 1./3. * (k2 + k3));
      shell->setVertexVelocity(*vit, v);
      shell->setVertexPosition(*vit, x + v * dt);
    }  
  }
  else if (m_active_scene == 12)
  {
    //zalesak disk test
    for (VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit)
    {
      Vec3d v;
      Vec3d x = shell->getVertexPosition(*vit);
      
      // RK4
      // -----------
      // k1 = dt * f( t, x );
      s12_zalesak_velocity(current_t, x, v);
      Vec3d k1 = v;
      
      // k2 = dt * f( t + 0.5*dt, x + 0.5*k1 );
      s12_zalesak_velocity(current_t + 0.5 * dt, x + 0.5 * dt * k1, v);
      Vec3d k2 = v;
      
      // k3 = dt * f( t + 0.5*dt, x + 0.5*k2 );
      s12_zalesak_velocity(current_t + 0.5 * dt, x + 0.5 * dt * k2, v);
      Vec3d k3 = v;
      
      // k4 = dt * f( t + dt, x + k3 );
      s12_zalesak_velocity(current_t + dt, x + dt * k3, v);
      Vec3d k4 = v;
      
      v = (1./6. * (k1 + k4) + 1./3. * (k2 + k3));
      shell->setVertexVelocity(*vit, v);
      shell->setVertexPosition(*vit, x + v * dt);
    }  
  }
  else if (m_active_scene == 9 || m_active_scene == 10)
  {
    //normal flow examples
    static Scalar speeds_scene9[3][3] = 
    {
      { 0, -1, -1 },
      { 1, 0, 0 },
      { 1, 0, 0 }
    };
    
    static Scalar speeds_scene10[3][3] = 
    {
      { 0, -1, 1 },
      { 1, 0, -2 },
      { -1, 2, 0 }
    };
    
    Scalar speeds[3][3];
    memcpy(speeds, m_active_scene == 9 || m_active_scene == 11 ? speeds_scene9 : speeds_scene10, sizeof (Scalar) * 9);
    
    VertexProperty<Vec3d> velocities(shellObj);
    for (VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit)
    {
      if (shellObj->vertexIncidentEdges(*vit) == 0) 
      { 
        velocities[*vit] = Vec3d(0, 0, 0);
        continue;
      }
      
      //      std::cout << "vertex (" << shell->getVertexPosition(*vit) << "):" << std::endl;
      
      Vec3d normal(0, 0, 0);
      double sum_areas = 0.0;
      bool label0 = false;
      bool label1 = false;
      bool label2 = false;
      for (VertexFaceIterator vfit = shellObj->vf_iter(*vit); vfit; ++vfit)
      {
        FaceVertexIterator fvit = shellObj->fv_iter(*vfit); assert(fvit);
        Vec3d x0 = shell->getVertexPosition(*fvit); ++fvit; assert(fvit);
        Vec3d x1 = shell->getVertexPosition(*fvit); ++fvit; assert(fvit);
        Vec3d x2 = shell->getVertexPosition(*fvit); ++fvit; assert(!fvit);
        
        Vec2i label = shell->getFaceLabel(*vfit);
        
        double area = (x1 - x0).cross(x2 - x0).norm() / 2;
        normal += (x1 - x0).cross(x2 - x0).normalized() * area * speeds[label.x()][label.y()];
        sum_areas += area;


        if(label.x() == 0 || label.y() == 0) label0 = true;
        if(label.x() == 1 || label.y() == 1) label1 = true;
        if(label.x() == 2 || label.y() == 2) label2 = true;
        
        //        std::cout << "face: label = " << label << " area = " << area << " normal = " << (x1 - x0).cross(x2 - x0).normalized() << " speed = " << speeds[label.x()][label.y()] << std::endl;
        //        std::cout << "  total normal = " << normal << std::endl;

      }
      if (normal.norm() > 0)
        normal.normalize();
      
      //if we're on the 3-way border, zero it out.
      //if(label0 && label1 && label2) normal = Vec3d(0,0,0);

      double speed = 0.1;
      velocities[*vit] = speed * normal;
    }
    
    //    double capped_dt = MeshSmoother::compute_max_timestep_quadratic_solve( surf.m_mesh.get_triangles(), surf.get_positions(), displacements, false );
    
    //    adaptive_dt = min( adaptive_dt, capped_dt );
    
    for (VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit)
    {
      //std::cout << "vertex " << shell->getVertexPosition(*vit) << " has velocity " << velocities[*vit] << std::endl;
      shell->setVertexVelocity(*vit, velocities[*vit]);
      shell->setVertexPosition(*vit, shell->getVertexPosition(*vit) + velocities[*vit] * dt);
    }
  }
  
  if (g_recording.isRecording())
  {
    ElTopo::SurfTrack * st = mesh2surftrack();
    g_recording.log() << "End time step" << std::endl;
    g_recording.recordSurfTrack(*st);
    delete st;
  }

}

void DoubleBubbleTest::s7_enright_velocity(double t, const Vec3d & pos, Vec3d & out)
{
  //
  //  code adapted from El Topo's Enright driver:
  //    https://github.com/tysonbrochu/eltopo/blob/master/talpa/drivers/enrightdriver.h
  //
  double x = pos[0]; 
  double y = pos[1]; 
  double z = pos[2];
  
  out = Vec3d( 2.0 * std::sin(M_PI*x) * std::sin(M_PI*x) * std::sin(2.0*M_PI*y) * std::sin(2.0*M_PI*z),
              -std::sin(2.0*M_PI*x) * std::sin(M_PI*y)*std::sin(M_PI*y) * std::sin(2.0*M_PI*z),
              -std::sin(2.0*M_PI*x) * std::sin(2.0*M_PI*y) * std::sin(M_PI*z) * std::sin(M_PI*z) );
  
  out *= sin(M_PI * t * 2 / 3);    // modulate with a period of 3
}

void DoubleBubbleTest::s12_zalesak_velocity(double t, const Vec3d & pos, Vec3d & out)
{
  double x = pos[0]; 
  double y = pos[1]; 
  double z = pos[2];
  
  out = Vec3d((z - 0.5), 0, -(x - 0.5));
}


void DoubleBubbleTest::AfterStep()
{
//  triangulation_added_vertices.clear();
//  triangulation_added_edges.clear();
//  triangulation_added_faces.clear();
//  svf->triangulateBBWalls(triangulation_added_vertices, triangulation_added_edges, triangulation_added_faces);
  
  if (m_active_scene == 7)
  {
    int nregion = 0;
    for (FaceIterator fit = shellObj->faces_begin(); fit != shellObj->faces_end(); ++fit)
    {
      Vec2i label = shell->getFaceLabel(*fit);
      if (label.x() + 1 > nregion) nregion = label.x() + 1;
      if (label.y() + 1 > nregion) nregion = label.y() + 1;
    }
    
    std::vector<Scalar> vol(nregion, 0);
    static std::vector<Scalar> init_vol(nregion, -1);
    
    Vec3d c = Vec3d(0, 0, 0);    
    for (FaceIterator fit = shellObj->faces_begin(); fit != shellObj->faces_end(); ++fit)
    {
      FaceVertexIterator fvit = shellObj->fv_iter(*fit); assert(fvit);
      Vec3d x0 = shell->getVertexPosition(*fvit); ++fvit; assert(fvit);
      Vec3d x1 = shell->getVertexPosition(*fvit); ++fvit; assert(fvit);
      Vec3d x2 = shell->getVertexPosition(*fvit); ++fvit; assert(!fvit);
      
      Vec2i label = shell->getFaceLabel(*fit);
      vol[label.x()] += (x0 - c).cross(x1 - c).dot(x2 - c);
      vol[label.y()] -= (x0 - c).cross(x1 - c).dot(x2 - c);
    }
    
    for (int i = 0; i < nregion; i++)
    {
      vol[i] /= 6;
      if (i == 0)
        vol[i] = -vol[i]; // this is used to compute total volume
      
      if (init_vol[i] < 0)
        init_vol[i] = vol[i];
      
      if (i == 0)
        std::cout << "Total volume = " << vol[i] << " error = " << fabs(vol[i] - init_vol[i]) * 100 / init_vol[i] << "%" << std::endl;
      else
        std::cout << "Region " << i << ": volume = " << vol[i] << " error = " << fabs(vol[i] - init_vol[i]) * 100 / init_vol[i] << "%" << std::endl;
    }
    
  }
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
  shell = new ElasticShell(shellObj, shellFaces, m_timestep, this);
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
  shell = new ElasticShell(shellObj, shellFaces, m_timestep, this);
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

void DoubleBubbleTest::setupScene3()
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
    shell = new ElasticShell(shellObj, shellFaces, m_timestep, this);
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
  shell = new ElasticShell(shellObj, shellFaces, m_timestep, this);
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
  int nsite = GetIntOpt("shell-x-resolution");
  srand(GetIntOpt("shell-y-resolution"));
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
  shell = new ElasticShell(shellObj, shellFaces, m_timestep, this);
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

void DoubleBubbleTest::setupScene7()
{
  //vertices
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);
  
  //create a sphere
  std::vector<VertexHandle> vertList;
  
  int N = GetIntOpt("shell-x-resolution");
  Scalar r = 0.1;
  Vec3d c = Vec3d(0.4, 0.6, 0.6);
  vertList.push_back(shellObj->addVertex());
  positions[vertList.back()] = Vec3d(c - Vec3d(0, 0, r));
  for (int j = 0; j < N - 1; j++)
  {
    for (int i = 0; i < N * 2; i++)
    {
      vertList.push_back(shellObj->addVertex());
      
      Scalar theta = (Scalar)i * 2 * M_PI / (N * 2);
      Scalar alpha = (Scalar)(j + 1) * M_PI / N - M_PI / 2;
      positions[vertList.back()] = c + r * Vec3d(cos(alpha) * cos(theta), cos(alpha) * sin(theta), sin(alpha));
    }
  }
  vertList.push_back(shellObj->addVertex());
  positions[vertList.back()] = Vec3d(c + Vec3d(0, 0, r));
  vertList.push_back(shellObj->addVertex());
  positions[vertList.back()] = Vec3d(c + Vec3d(0, 0, 0));
  
  for (int i = 0; i < shellObj->nv(); ++i)
  {
    velocities[vertList[i]] = Vec3d(0, 0, 0);
    undeformed[vertList[i]] = positions[vertList[i]];
  }
  
  std::vector<FaceHandle> faceList;
  FaceProperty<Vec2i> faceLabels(shellObj); //label face regions to do volume constrained bubbles  

  int Nsplit = 8;
  if (Nsplit == 2)
  {
    for (int j = 0; j < N; j++)
    {
      for (int i = 0; i < N * 2; i++)
      {
        int v0, v1, v2;
        v0 = (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
        v1 = (j == 0 ? 0 : 2 * N * (j - 1) + (i + 1) % (N * 2) + 1);
        v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
        if (!(v0 == v1 || v0 == v2 || v1 == v2))
        {
          faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
          faceLabels[faceList.back()] = Vec2i((i < N ? 1 : 2), 0);
        }
        
        v0 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
        v1 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + i + 1);
        v2 = (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
        if (!(v0 == v1 || v0 == v2 || v1 == v2))
        {
          faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
          faceLabels[faceList.back()] = Vec2i((i < N ? 1 : 2), 0);
        }
      }
    }
    
    for (int j = 0; j < N; j++)
    {
      int v0, v1, v2;
      v0 = (j == 0 ? 0 : 2 * N * (j - 1) + 0 + 1);
      v1 = 2 * (N - 1) * N + 2;
      v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + 0 + 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(2, 1);
      }
      
      v0 = (j == 0 ? 0 : 2 * N * (j - 1) + N + 1);
      v1 = 2 * (N - 1) * N + 2;
      v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + N + 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(1, 2);
      }
    }
    
  } else if (Nsplit == 4)
  {
    for (int j = 0; j < N; j++)
    {
      for (int i = 0; i < N * 2; i++)
      {
        int l;
        if (i < N / 2) l = 1;
        else if (i < N) l = 3;
        else if (i < N * 3 / 2) l = 4;
        else l = 2;
        
        int v0, v1, v2;
        v0 = (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
        v1 = (j == 0 ? 0 : 2 * N * (j - 1) + (i + 1) % (N * 2) + 1);
        v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
        if (!(v0 == v1 || v0 == v2 || v1 == v2))
        {
          faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
          faceLabels[faceList.back()] = Vec2i(l, 0);
        }
        
        v0 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
        v1 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + i + 1);
        v2 = (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
        if (!(v0 == v1 || v0 == v2 || v1 == v2))
        {
          faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
          faceLabels[faceList.back()] = Vec2i(l, 0);
        }
      }
    }
    
    for (int j = 0; j < N; j++)
    {
      int v0, v1, v2;
      v0 = (j == 0 ? 0 : 2 * N * (j - 1) + 0 + 1);
      v1 = 2 * (N - 1) * N + 2;
      v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + 0 + 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(2, 1);
      }
      
      v0 = (j == 0 ? 0 : 2 * N * (j - 1) + N + 1);
      v1 = 2 * (N - 1) * N + 2;
      v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + N + 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(3, 4);
      }

      v0 = (j == 0 ? 0 : 2 * N * (j - 1) + N / 2 + 1);
      v1 = 2 * (N - 1) * N + 2;
      v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + N / 2 + 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(1, 3);
      }

      v0 = (j == 0 ? 0 : 2 * N * (j - 1) + N * 3 / 2+ 1);
      v1 = 2 * (N - 1) * N + 2;
      v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + N * 3 / 2+ 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(4, 2);
      }
    }
    
  } else if (Nsplit == 8)
  {
    for (int j = 0; j < N; j++)
    {
      for (int i = 0; i < N * 2; i++)
      {
        int l;
        if (i < N / 2) l = 1;
        else if (i < N) l = 3;
        else if (i < N * 3 / 2) l = 4;
        else l = 2;
        
        if (j < N / 2) l += 4;
        
        int v0, v1, v2;
        v0 = (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
        v1 = (j == 0 ? 0 : 2 * N * (j - 1) + (i + 1) % (N * 2) + 1);
        v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
        if (!(v0 == v1 || v0 == v2 || v1 == v2))
        {
          faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
          faceLabels[faceList.back()] = Vec2i(l, 0);
        }
        
        v0 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
        v1 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + i + 1);
        v2 = (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
        if (!(v0 == v1 || v0 == v2 || v1 == v2))
        {
          faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
          faceLabels[faceList.back()] = Vec2i(l, 0);
        }
      }
    }
    
    for (int j = 0; j < N; j++)
    {
      int la = (j < N / 2 ? 4 : 0);
      
      int v0, v1, v2;
      v0 = (j == 0 ? 0 : 2 * N * (j - 1) + 0 + 1);
      v1 = 2 * (N - 1) * N + 2;
      v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + 0 + 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(2 + la, 1 + la);
      }
      
      v0 = (j == 0 ? 0 : 2 * N * (j - 1) + N + 1);
      v1 = 2 * (N - 1) * N + 2;
      v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + N + 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(3 + la, 4 + la);
      }
      
      v0 = (j == 0 ? 0 : 2 * N * (j - 1) + N / 2 + 1);
      v1 = 2 * (N - 1) * N + 2;
      v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + N / 2 + 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(1 + la, 3 + la);
      }
      
      v0 = (j == 0 ? 0 : 2 * N * (j - 1) + N * 3 / 2+ 1);
      v1 = 2 * (N - 1) * N + 2;
      v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + N * 3 / 2+ 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(4 + la, 2 + la);
      }
    }
    
    for (int i = 0; i < N * 2; i++)
    {
      int l;
      if (i < N / 2) l = 1;
      else if (i < N) l = 3;
      else if (i < N * 3 / 2) l = 4;
      else l = 2;
      
      int v0, v1, v2;
      v0 = 2 * N * (N / 2 - 1) + i + 1;
      v1 = 2 * (N - 1) * N + 2;
      v2 = 2 * N * (N / 2 - 1) + (i + 1) % (N * 2) + 1;
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(l, l + 4);
      }
    }
    
  } else  // default: Nsplit = 1
  {
    for (int j = 0; j < N; j++)
    {
      for (int i = 0; i < N * 2; i++)
      {
        int v0, v1, v2;
        v0 = (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
        v1 = (j == 0 ? 0 : 2 * N * (j - 1) + (i + 1) % (N * 2) + 1);
        v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
        if (!(v0 == v1 || v0 == v2 || v1 == v2))
        {
          faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
          faceLabels[faceList.back()] = Vec2i(1, 0);
        }
        
        v0 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
        v1 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + i + 1);
        v2 = (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
        if (!(v0 == v1 || v0 == v2 || v1 == v2))
        {
          faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
          faceLabels[faceList.back()] = Vec2i(1, 0);
        }
      }
    }
  }
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(shellObj); 
  DeformableObject::face_iter fIt;
  for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(shellObj, shellFaces, m_timestep, this);
  shellObj->addModel(shell);
  
  //positions
  //  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  shell->setFaceLabels(faceLabels);
}

void DoubleBubbleTest::setupScene8()
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
    
    int nv = 16;
    for(int i = 0; i < nv; ++i) 
        vertList.push_back(shellObj->addVertex());
    
    //create positions
    Scalar size = GetScalarOpt("shell-width");
    Vec3d center(0.5, 0.5, 0.5);
    Scalar d = 0.1464466094;
    positions[vertList[ 0]] = Vec3d(0, 0, 0.5);
    positions[vertList[ 1]] = Vec3d(d, 0, 0.5);
    positions[vertList[ 2]] = Vec3d(d, 0, 0.5) * size + center * (1 - size);
    positions[vertList[ 3]] = Vec3d(d, 1, 0.5) * size + center * (1 - size);
    positions[vertList[ 4]] = Vec3d(d, 1, 0.5);
    positions[vertList[ 5]] = Vec3d(0, 1, 0.5);
    
    positions[vertList[ 6]] = Vec3d(1, 0.5, 1);
    positions[vertList[ 7]] = Vec3d(1 - d, 0.5, 1);
    positions[vertList[ 8]] = Vec3d(1 - d, 0.5, 1) * size + center * (1 - size);
    positions[vertList[ 9]] = Vec3d(1 - d, 0.5, 0) * size + center * (1 - size);
    positions[vertList[10]] = Vec3d(1 - d, 0.5, 0);
    positions[vertList[11]] = Vec3d(1, 0.5, 0);
    
    positions[vertList[12]] = Vec3d(0.5, 0, 1);
    positions[vertList[13]] = Vec3d(0.5, 1, 1);
    positions[vertList[14]] = Vec3d(0.5, 1, 0);
    positions[vertList[15]] = Vec3d(0.5, 0, 0);
    
    int tet_corners[4];
    tet_corners[0] = 2;
    tet_corners[1] = 3;
    tet_corners[2] = 8;
    tet_corners[3] = 9;
    
    int res = GetIntOpt("shell-x-resolution");
    for (int i = 0; i < 4; i++)
    {
        for (int j = i + 1; j < 4; j++)
        {
            VertexHandle a = vertList[tet_corners[i]];
            VertexHandle b = vertList[tet_corners[j]];
            Vec3d xa = positions[a];
            Vec3d xb = positions[b];
            
            // res - 1 vertices created
            for (int k = 1; k < res; k++)
            {
                vertList.push_back(shellObj->addVertex());
                positions[vertList.back()] = xa + (xb - xa) * k / res;
            }
        }
    }

    for (int f = 0; f < 4; f++)
    {
        VertexHandle a, b, c, d;
        a = vertList[tet_corners[(f + 0) % 4]];
        b = vertList[tet_corners[(f + 1) % 4]];
        c = vertList[tet_corners[(f + 2) % 4]];
        d = vertList[tet_corners[(f + 3) % 4]];
        Vec3d xa, xb, xc, xd;
        xa = positions[a];
        xb = positions[b];
        xc = positions[c];
        xd = positions[d];
        
        // (res - 1) * (res - 2) vertices created
        for (int i = 2; i < res; i++)
        {
            for (int j = 1; j < i; j++)
            {
                vertList.push_back(shellObj->addVertex());
                positions[vertList.back()] = xc * (i - j) / res + xa * (res - i) / res + xb * j / res; 
            }
        }
    }
    
    for(int i = 0; i < shellObj->nv(); ++i)
    {
        undeformed[vertList[i]] = positions[vertList[i]];
        velocities[vertList[i]] = Vec3d(0, 0, 0);
    }

    std::vector<FaceHandle> faceList;
    FaceProperty<Vec2i> faceLabels(shellObj); //label face regions to do volume constrained bubbles  
    
    int regions[4] = { 0, 3, 2, 1 };    
    for (int f = 0; f < 4; f++)
    {
        VertexHandle a, b, c, d;
        a = vertList[tet_corners[(f + 0) % 4]];
        b = vertList[tet_corners[(f + 1) % 4]];
        c = vertList[tet_corners[(f + 2) % 4]];
        d = vertList[tet_corners[(f + 3) % 4]];
        Vec3d xa, xb, xc, xd;
        xa = positions[a];
        xb = positions[b];
        xc = positions[c];
        xd = positions[d];
        
        bool oriented = ((xa - center).cross(xb - center).dot(xc - center) > 0);
        Vec2i label = (oriented ? Vec2i(regions[f], 4) : Vec2i(4, regions[f]));

        for (int i = 0; i < res; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                Vec3d x0, x1, x2;
                x0 = xc * (i - j) / res + xa * (res - i) / res + xb * j / res; 
                x1 = xc * (i + 1 - j) / res + xa * (res - i - 1) / res + xb * j / res; 
                x2 = xc * (i - j) / res + xa * (res - i - 1) / res + xb * (j + 1) / res; 
                
                size_t v0, v1, v2;
                for (v0 = 0; v0 < vertList.size(); v0++)
                    if ((x0 - positions[vertList[v0]]).squaredNorm() < 1e-6)
                        break;
                assert(v0 < vertList.size());
                for (v1 = 0; v1 < vertList.size(); v1++)
                    if ((x1 - positions[vertList[v1]]).squaredNorm() < 1e-6)
                        break;
                assert(v1 < vertList.size());
                for (v2 = 0; v2 < vertList.size(); v2++)
                    if ((x2 - positions[vertList[v2]]).squaredNorm() < 1e-6)
                        break;
                assert(v2 < vertList.size());
                
                faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
                faceLabels[faceList.back()] = label;
            }
        }
        
        for (int i = 1; i < res; i++)
        {
            for (int j = 0; j < i; j++)
            {
                Vec3d x0, x1, x2;
                x0 = xc * (i - j) / res + xa * (res - i) / res + xb * j / res; 
                x1 = xc * (i - j) / res + xa * (res - i - 1) / res + xb * (j + 1) / res; 
                x2 = xc * (i - j - 1) / res + xa * (res - i) / res + xb * (j + 1) / res; 
                
                size_t v0, v1, v2;
                for (v0 = 0; v0 < vertList.size(); v0++)
                    if ((x0 - positions[vertList[v0]]).squaredNorm() < 1e-6)
                        break;
                assert(v0 < vertList.size());
                for (v1 = 0; v1 < vertList.size(); v1++)
                    if ((x1 - positions[vertList[v1]]).squaredNorm() < 1e-6)
                        break;
                assert(v1 < vertList.size());
                for (v2 = 0; v2 < vertList.size(); v2++)
                    if ((x2 - positions[vertList[v2]]).squaredNorm() < 1e-6)
                        break;
                assert(v2 < vertList.size());
                
                faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
                faceLabels[faceList.back()] = label;
            }
        }
    }
    
    faceList.push_back(shellObj->addFace(vertList[0], vertList[1], vertList[2])); faceLabels[faceList.back()] = Vec2i(1, 0);
    faceList.push_back(shellObj->addFace(vertList[3], vertList[4], vertList[5])); faceLabels[faceList.back()] = Vec2i(1, 0);
    faceList.push_back(shellObj->addFace(vertList[6], vertList[7], vertList[8])); faceLabels[faceList.back()] = Vec2i(3, 2);
    faceList.push_back(shellObj->addFace(vertList[9], vertList[10], vertList[11])); faceLabels[faceList.back()] = Vec2i(3, 2);

    faceList.push_back(shellObj->addFace(vertList[5], vertList[0], vertList[nv + (res - 1) * 0 + res / 2 - 1])); faceLabels[faceList.back()] = Vec2i(1, 0);
    faceList.push_back(shellObj->addFace(vertList[11], vertList[6], vertList[nv + (res - 1) * 5 + res / 2 - 1])); faceLabels[faceList.back()] = Vec2i(3, 2);
    
    faceList.push_back(shellObj->addFace(vertList[1], vertList[2], vertList[12])); faceLabels[faceList.back()] = Vec2i(0, 2);
    faceList.push_back(shellObj->addFace(vertList[12], vertList[8], vertList[7])); faceLabels[faceList.back()] = Vec2i(0, 2);
    faceList.push_back(shellObj->addFace(vertList[4], vertList[3], vertList[13])); faceLabels[faceList.back()] = Vec2i(3, 0);
    faceList.push_back(shellObj->addFace(vertList[13], vertList[8], vertList[7])); faceLabels[faceList.back()] = Vec2i(3, 0);
    faceList.push_back(shellObj->addFace(vertList[1], vertList[2], vertList[15])); faceLabels[faceList.back()] = Vec2i(2, 1);
    faceList.push_back(shellObj->addFace(vertList[15], vertList[9], vertList[10])); faceLabels[faceList.back()] = Vec2i(2, 1);
    faceList.push_back(shellObj->addFace(vertList[4], vertList[3], vertList[14])); faceLabels[faceList.back()] = Vec2i(1, 3);
    faceList.push_back(shellObj->addFace(vertList[14], vertList[9], vertList[10])); faceLabels[faceList.back()] = Vec2i(1, 3);
    
    for (int i = 0; i < res; i++)
    {
        int v0, v1, v2;
        v0 = (i < res / 2 ? 0 : 5);
        v1 = (i == 0       ? tet_corners[0] : (nv + (res - 1) * 0 + i - 1));
        v2 = (i == res - 1 ? tet_corners[1] : (nv + (res - 1) * 0 + i));        
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(1, 0);
        
        v0 = (i < res / 2 ? 6 : 11);
        v1 = (i == 0       ? tet_corners[2] : (nv + (res - 1) * 5 + i - 1));
        v2 = (i == res - 1 ? tet_corners[3] : (nv + (res - 1) * 5 + i));        
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(3, 2);
        
        v0 = 12;
        v1 = (i == 0       ? tet_corners[0] : (nv + (res - 1) * 1 + i - 1));
        v2 = (i == res - 1 ? tet_corners[2] : (nv + (res - 1) * 1 + i));        
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(0, 2);
        
        v0 = 13;
        v1 = (i == 0       ? tet_corners[1] : (nv + (res - 1) * 3 + i - 1));
        v2 = (i == res - 1 ? tet_corners[2] : (nv + (res - 1) * 3 + i));        
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(3, 0);
        
        v0 = 14;
        v1 = (i == 0       ? tet_corners[1] : (nv + (res - 1) * 4 + i - 1));
        v2 = (i == res - 1 ? tet_corners[3] : (nv + (res - 1) * 4 + i));        
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(1, 3);
        
        v0 = 15;
        v1 = (i == 0       ? tet_corners[0] : (nv + (res - 1) * 2 + i - 1));
        v2 = (i == res - 1 ? tet_corners[3] : (nv + (res - 1) * 2 + i));        
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(2, 1);
    }
    
    // push the vertices outside, inflating the tet
    for (int i = 0; i < 4; i++)
    {
        for (int j = i + 1; j < 4; j++)
        {
            VertexHandle a = vertList[tet_corners[i]];
            VertexHandle b = vertList[tet_corners[j]];
            Vec3d xa = positions[a];
            Vec3d xb = positions[b];
            
            VertexHandle c;
            VertexHandle d;
            for (int k = 0; k < 4; k++)
            {
                if (k != i && k != j)
                {
                    if (!c.isValid())
                        c = vertList[tet_corners[k]];
                    else
                        d = vertList[tet_corners[k]];
                }
            }
            Vec3d xc = positions[c];
            Vec3d xd = positions[d];
            
            // res - 1 vertices created
            for (int k = 1; k < res; k++)
            {
                Vec3d x = xa + (xb - xa) * k / res;
                
                size_t v0;
                for (v0 = 0; v0 < vertList.size(); v0++)
                    if ((x - positions[vertList[v0]]).squaredNorm() < 1e-6)
                        break;
                assert(v0 < vertList.size());
                
                positions[vertList[v0]] = (xc + xd) / 2 + (x - (xc + xd) / 2).normalized() * (xa - (xc + xd) / 2).norm();
            }
        }
    }
    
    for (int f = 0; f < 4; f++)
    {
        VertexHandle a, b, c, d;
        a = vertList[tet_corners[(f + 0) % 4]];
        b = vertList[tet_corners[(f + 1) % 4]];
        c = vertList[tet_corners[(f + 2) % 4]];
        d = vertList[tet_corners[(f + 3) % 4]];
        Vec3d xa, xb, xc, xd;
        xa = positions[a];
        xb = positions[b];
        xc = positions[c];
        xd = positions[d];
        
        // (res - 1) * (res - 2) vertices created
        for (int i = 2; i < res; i++)
        {
            for (int j = 1; j < i; j++)
            {
                Vec3d x = xc * (i - j) / res + xa * (res - i) / res + xb * j / res; 
                
                size_t v0;
                for (v0 = 0; v0 < vertList.size(); v0++)
                    if ((x - positions[vertList[v0]]).squaredNorm() < 1e-6)
                        break;
                assert(v0 < vertList.size());
                
                positions[vertList[v0]] = xd + (x - xd).normalized() * (xa - xd).norm();
            }
        }
    }
    

    
    //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
    FaceProperty<char> shellFaces(shellObj); 
    DeformableObject::face_iter fIt;
    for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
        shellFaces[*fIt] = true;
    
    //now create the physical model to hang on the mesh
    shell = new ElasticShell(shellObj, shellFaces, m_timestep, this);
    shellObj->addModel(shell);
    
    //positions
//  shell->setVertexUndeformed(undeformed);
    shell->setVertexPositions(positions);
    shell->setVertexVelocities(velocities);
    
    shell->setFaceLabels(faceLabels);
    
  
}

void DoubleBubbleTest::createIcoSphere(DeformableObject & mesh, Vec3d & center, Scalar r, int subdivision, std::vector<VertexHandle> & vertList, std::vector<FaceHandle> & faceList, VertexProperty<Vec3d> & positions)
{
  vertList.clear();
  faceList.clear();
  
  DeformableObject obj;
  VertexProperty<Vec3d> pos(&obj);
  
  // create icosahedron
  for (int i = 0; i < 12; i++)
    vertList.push_back(obj.addVertex());
  
  Scalar phi = (1.0 + sqrt(5.0)) / 2.0;
  Scalar len = Vec3d(0, 1, phi).norm();
  pos[*(vertList.rbegin() +  0)] = center + r * Vec3d(0,  1,  phi) / len; 
  pos[*(vertList.rbegin() +  1)] = center + r * Vec3d(0, -1,  phi) / len; 
  pos[*(vertList.rbegin() +  2)] = center + r * Vec3d(0,  1, -phi) / len; 
  pos[*(vertList.rbegin() +  3)] = center + r * Vec3d(0, -1, -phi) / len; 
  pos[*(vertList.rbegin() +  4)] = center + r * Vec3d( 1,  phi, 0) / len; 
  pos[*(vertList.rbegin() +  5)] = center + r * Vec3d(-1,  phi, 0) / len; 
  pos[*(vertList.rbegin() +  6)] = center + r * Vec3d( 1, -phi, 0) / len; 
  pos[*(vertList.rbegin() +  7)] = center + r * Vec3d(-1, -phi, 0) / len; 
  pos[*(vertList.rbegin() +  8)] = center + r * Vec3d( phi, 0,  1) / len; 
  pos[*(vertList.rbegin() +  9)] = center + r * Vec3d( phi, 0, -1) / len; 
  pos[*(vertList.rbegin() + 10)] = center + r * Vec3d(-phi, 0,  1) / len; 
  pos[*(vertList.rbegin() + 11)] = center + r * Vec3d(-phi, 0, -1) / len;
  
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  0), *(vertList.rbegin() +  1), *(vertList.rbegin() +  8)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  1), *(vertList.rbegin() +  0), *(vertList.rbegin() + 10)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  2), *(vertList.rbegin() +  3), *(vertList.rbegin() + 11)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  3), *(vertList.rbegin() +  2), *(vertList.rbegin() +  9)));
  
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  4), *(vertList.rbegin() +  5), *(vertList.rbegin() +  0)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  5), *(vertList.rbegin() +  4), *(vertList.rbegin() +  2)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  6), *(vertList.rbegin() +  7), *(vertList.rbegin() +  3)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  7), *(vertList.rbegin() +  6), *(vertList.rbegin() +  1)));
  
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  8), *(vertList.rbegin() +  9), *(vertList.rbegin() +  4)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  9), *(vertList.rbegin() +  8), *(vertList.rbegin() +  6)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() + 10), *(vertList.rbegin() + 11), *(vertList.rbegin() +  7)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() + 11), *(vertList.rbegin() + 10), *(vertList.rbegin() +  5)));
  
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  0), *(vertList.rbegin() +  8), *(vertList.rbegin() +  4)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  1), *(vertList.rbegin() +  6), *(vertList.rbegin() +  8)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  0), *(vertList.rbegin() +  5), *(vertList.rbegin() + 10)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  1), *(vertList.rbegin() + 10), *(vertList.rbegin() +  7)));
  
  faceList.push_back(obj.addFace(*(vertList.rbegin() + 11), *(vertList.rbegin() +  3), *(vertList.rbegin() +  7)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  5), *(vertList.rbegin() +  2), *(vertList.rbegin() + 11)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  6), *(vertList.rbegin() +  3), *(vertList.rbegin() +  9)));
  faceList.push_back(obj.addFace(*(vertList.rbegin() +  9), *(vertList.rbegin() +  2), *(vertList.rbegin() +  4)));
  
  for (int i = 0; i < subdivision; i++)
  {
    std::vector<Eigen::Matrix<VertexHandle, 3, 1> > faces_to_create;
    EdgeProperty<VertexHandle> midpoints(&obj);
    
    for (EdgeIterator eit = obj.edges_begin(); eit != obj.edges_end(); ++eit)
    {
      VertexHandle vh = obj.addVertex();
      vertList.push_back(vh);
      midpoints[*eit] = vh;
      pos[vh] = center + r * ((pos[obj.fromVertex(*eit)] + pos[obj.toVertex(*eit)]) / 2 - center).normalized();
    }
    
    for (FaceIterator fit = obj.faces_begin(); fit != obj.faces_end(); ++fit)
    {
      FaceVertexIterator fvit = obj.fv_iter(*fit); assert(fvit);
      VertexHandle v0 = *fvit; ++fvit; assert(fvit);
      VertexHandle v1 = *fvit; ++fvit; assert(fvit);
      VertexHandle v2 = *fvit; ++fvit; assert(!fvit);
      
      EdgeHandle e01 = findEdge(obj, v0, v1); assert(e01.isValid());
      EdgeHandle e12 = findEdge(obj, v1, v2); assert(e12.isValid());
      EdgeHandle e20 = findEdge(obj, v2, v0); assert(e20.isValid());

      VertexHandle m01 = midpoints[e01]; assert(m01.isValid());
      VertexHandle m12 = midpoints[e12]; assert(m12.isValid());
      VertexHandle m20 = midpoints[e20]; assert(m20.isValid());
      
      faces_to_create.push_back(Eigen::Matrix<VertexHandle, 3, 1>(v0, m01, m20));
      faces_to_create.push_back(Eigen::Matrix<VertexHandle, 3, 1>(v1, m12, m01));
      faces_to_create.push_back(Eigen::Matrix<VertexHandle, 3, 1>(v2, m20, m12));
      faces_to_create.push_back(Eigen::Matrix<VertexHandle, 3, 1>(m01, m12, m20));
      
    }

    std::vector<FaceHandle> faces_created;
    for (size_t j = 0; j < faces_to_create.size(); j++)
      faces_created.push_back(obj.addFace(faces_to_create[j].x(), faces_to_create[j].y(), faces_to_create[j].z()));

    for (size_t j = 0; j < faceList.size(); j++)
      obj.deleteFace(faceList[j], true);

    faceList = faces_created;
    
  }
  
  vertList.clear();
  faceList.clear();
  
  VertexProperty<VertexHandle> vert_map(&obj);
  for (VertexIterator vit = obj.vertices_begin(); vit != obj.vertices_end(); ++vit)
  {
    vertList.push_back(mesh.addVertex());
    positions[vertList.back()] = pos[*vit];
    vert_map[*vit] = vertList.back();
  }
  
  for (FaceIterator fit = obj.faces_begin(); fit != obj.faces_end(); ++fit)
  {
    FaceVertexIterator fvit = obj.fv_iter(*fit); assert(fvit);
    VertexHandle v0 = *fvit; ++fvit; assert(fvit);
    VertexHandle v1 = *fvit; ++fvit; assert(fvit);
    VertexHandle v2 = *fvit; ++fvit; assert(!fvit);
    
    faceList.push_back(mesh.addFace(vert_map[v0], vert_map[v1], vert_map[v2]));
  }
  
}

void DoubleBubbleTest::setupScene9() 
{
  //vertices
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);
  
  //create a sphere
  std::vector<VertexHandle> vertList;
  
  int N = GetIntOpt("shell-x-resolution");
  Scalar separate = GetScalarOpt("shell-width");
  Scalar r = GetScalarOpt("shell-height");
  Vec3d c1 = Vec3d(0.49 - separate, 0.51, 0.509765);
  Vec3d c2 = Vec3d(0.55 + separate, 0.53, 0.46876);
  
  std::vector<FaceHandle> faceList;
  FaceProperty<Vec2i> faceLabels(shellObj);  
  std::vector<VertexHandle> vs;
  std::vector<FaceHandle> fs;
  
  vs.clear();
  fs.clear();
  createIcoSphere(*shellObj, c1, r, N, vs, fs, positions);
  for (size_t i = 0; i < fs.size(); i++)
    faceLabels[fs[i]] = Vec2i(1, 0);
  vertList.insert(vertList.end(), vs.begin(), vs.end());
  faceList.insert(faceList.end(), fs.begin(), fs.end());
  
  vs.clear();
  fs.clear();
  createIcoSphere(*shellObj, c2, r, N, vs, fs, positions);
  for (size_t i = 0; i < fs.size(); i++)
    faceLabels[fs[i]] = Vec2i(2, 0);
  vertList.insert(vertList.end(), vs.begin(), vs.end());
  faceList.insert(faceList.end(), fs.begin(), fs.end());
  
  velocities.assign(Vec3d(0, 0, 0));
  undeformed = positions;
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(shellObj); 
  DeformableObject::face_iter fIt;
  for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(shellObj, shellFaces, m_timestep, this);
  shellObj->addModel(shell);
  
  //positions
  //  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  shell->setFaceLabels(faceLabels);

}

void DoubleBubbleTest::setupScene10() 
{
  //vertices
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);
  
  //create a sphere
  std::vector<VertexHandle> vertList;
  
  int N = GetIntOpt("shell-x-resolution");
  Scalar r = GetScalarOpt("shell-height");
  Vec3d c1 = Vec3d(0.5, 0.5, 0.5 - r * 0.707);
  vertList.push_back(shellObj->addVertex());
  positions[vertList.back()] = Vec3d(c1 + Vec3d(0, 0, r));
  for (int j = 0; j < N - 1; j++)
  {
    for (int i = 0; i < N * 2; i++)
    {
      vertList.push_back(shellObj->addVertex());
      
      Scalar theta = (Scalar)i * 2 * M_PI / (N * 2);
      Scalar alpha = (Scalar)(j + 1) * M_PI / N - M_PI / 2;
      positions[vertList.back()] = c1 + r * Vec3d(cos(alpha) * cos(theta), cos(alpha) * sin(theta), -sin(alpha));
    }
  }
  vertList.push_back(shellObj->addVertex());
  positions[vertList.back()] = Vec3d(c1 - Vec3d(0, 0, r));
  vertList.push_back(shellObj->addVertex());
  positions[vertList.back()] = Vec3d(0.5, 0.5, 0.5);
  
  Vec3d c2 = Vec3d(0.5, 0.5, 0.5 + r * 0.707);
  vertList.push_back(shellObj->addVertex());
  positions[vertList.back()] = Vec3d(c2 - Vec3d(0, 0, r));
  for (int j = 0; j < N - 1; j++)
  {
    for (int i = 0; i < N * 2; i++)
    {
      vertList.push_back(shellObj->addVertex());
      
      Scalar theta = (Scalar)i * 2 * M_PI / (N * 2);
      Scalar alpha = (Scalar)(j + 1) * M_PI / N - M_PI / 2;
      positions[vertList.back()] = c2 + r * Vec3d(cos(alpha) * cos(theta), cos(alpha) * sin(theta), sin(alpha));
    }
  }
  vertList.push_back(shellObj->addVertex());
  positions[vertList.back()] = Vec3d(c2 + Vec3d(0, 0, r));
  
  for (int i = 0; i < shellObj->nv(); ++i)
  {
    velocities[vertList[i]] = Vec3d(0, 0, 0);
    undeformed[vertList[i]] = positions[vertList[i]];
  }
  
  std::vector<FaceHandle> faceList;
  FaceProperty<Vec2i> faceLabels(shellObj); //label face regions to do volume constrained bubbles  
  
  for (int j = N / 4; j < N; j++)
  {
    for (int i = 0; i < N * 2; i++)
    {
      int v0, v1, v2;
      v0 = (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
      v1 = (j == 0 ? 0 : 2 * N * (j - 1) + (i + 1) % (N * 2) + 1);
      v2 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(0, 1);
      }
      
      v0 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
      v1 = (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + i + 1);
      v2 = (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(0, 1);
      }
    }
  }
  
  int offset = 2 * N * (N - 1) + 3;
  for (int j = N / 4 + 1; j < N; j++)
  {
    for (int i = 0; i < N * 2; i++)
    {
      int v0, v1, v2;
      v0 = offset + (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
      v1 = offset + (j == 0 ? 0 : 2 * N * (j - 1) + (i + 1) % (N * 2) + 1);
      v2 = offset + (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(2, 0);
      }
      
      v0 = offset + (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + (i + 1) % (N * 2) + 1);
      v1 = offset + (j == N - 1 ? 2 * (N - 1) * N + 1 : 2 * N * j + i + 1);
      v2 = offset + (j == 0 ? 0 : 2 * N * (j - 1) + i + 1);
      if (!(v0 == v1 || v0 == v2 || v1 == v2))
      {
        faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
        faceLabels[faceList.back()] = Vec2i(2, 0);
      }
    }
  }
  
  for (int i = 0; i < N * 2; i++)
  {
    int v0, v1, v2;
    v0 = 2 * N * (N / 4 - 1) + i + 1;
    v1 = 2 * N * (N / 4 - 1) + (i + 1) % (N * 2) + 1;
    v2 = offset + 2 * N * (N / 4) + (i + 1) % (N * 2) + 1;
    if (!(v0 == v1 || v0 == v2 || v1 == v2))
    {
      faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
      faceLabels[faceList.back()] = Vec2i(2, 0);
    }
    
    v0 = offset + 2 * N * (N / 4) + (i + 1) % (N * 2) + 1;
    v1 = offset + 2 * N * (N / 4) + i + 1;
    v2 = 2 * N * (N / 4 - 1) + i + 1;
    if (!(v0 == v1 || v0 == v2 || v1 == v2))
    {
      faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
      faceLabels[faceList.back()] = Vec2i(2, 0);
    }
  }
  
  for (int i = 0; i < N * 2; i++)
  {
    int v0, v1, v2;
    v0 = 2 * N * (N / 4 - 1) + i + 1;
    v1 = 2 * (N - 1) * N + 2;
    v2 = 2 * N * (N / 4 - 1) + (i + 1) % (N * 2) + 1;
    if (!(v0 == v1 || v0 == v2 || v1 == v2))
    {
      faceList.push_back(shellObj->addFace(vertList[v0], vertList[v1], vertList[v2]));
      faceLabels[faceList.back()] = Vec2i(2, 1);
    }
  }
  
  for (VertexIterator vit = shellObj->vertices_begin(); vit != shellObj->vertices_end(); ++vit)
    if (shellObj->vertexIncidentEdges(*vit) == 0)
      shellObj->deleteVertex(*vit);
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(shellObj); 
  DeformableObject::face_iter fIt;
  for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(shellObj, shellFaces, m_timestep, this);
  shellObj->addModel(shell);
  
  //positions
  //  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  shell->setFaceLabels(faceLabels);
  
}

void DoubleBubbleTest::setupScene11() 
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
  
  for(int i = 0; i < 9; ++i) 
  {
    vertList.push_back(shellObj->addVertex());
    velocities[vertList[i]] = Vec3d(0,0,0);
  }
  
  //create positions
  positions[vertList[ 0]] = Vec3d(0.1, 0, 0.1);
  positions[vertList[ 1]] = Vec3d(0.9, 0, 0.1);
  positions[vertList[ 2]] = Vec3d(0.9, 0, 0.9);
  positions[vertList[ 3]] = Vec3d(0.1, 0, 0.9);
  
  positions[vertList[ 4]] = Vec3d(0.1, 1, 0.1);
  positions[vertList[ 5]] = Vec3d(0.9, 1, 0.1);
  positions[vertList[ 6]] = Vec3d(0.9, 1, 0.9);
  positions[vertList[ 7]] = Vec3d(0.1, 1, 0.9);

  positions[vertList[ 8]] = Vec3d(0.5, 0.5, 0.5);
  
  for(int i = 0; i < 9; ++i)
    undeformed[vertList[i]] = positions[vertList[i]];
  
  std::vector<FaceHandle> faceList;
  FaceProperty<Vec2i> faceLabels(shellObj); //label face regions to do volume constrained bubbles  
  
  faceList.push_back(shellObj->addFace(vertList[ 0], vertList[ 1], vertList[ 8]));  faceLabels[faceList.back()] = Vec2i(0, 1);
  faceList.push_back(shellObj->addFace(vertList[ 1], vertList[ 2], vertList[ 8]));  faceLabels[faceList.back()] = Vec2i(0, 1);
  faceList.push_back(shellObj->addFace(vertList[ 2], vertList[ 3], vertList[ 8]));  faceLabels[faceList.back()] = Vec2i(0, 1);
  faceList.push_back(shellObj->addFace(vertList[ 3], vertList[ 0], vertList[ 8]));  faceLabels[faceList.back()] = Vec2i(0, 1);
  
  faceList.push_back(shellObj->addFace(vertList[ 4], vertList[ 5], vertList[ 8]));  faceLabels[faceList.back()] = Vec2i(2, 0);
  faceList.push_back(shellObj->addFace(vertList[ 5], vertList[ 6], vertList[ 8]));  faceLabels[faceList.back()] = Vec2i(2, 0);
  faceList.push_back(shellObj->addFace(vertList[ 6], vertList[ 7], vertList[ 8]));  faceLabels[faceList.back()] = Vec2i(2, 0);
  faceList.push_back(shellObj->addFace(vertList[ 7], vertList[ 4], vertList[ 8]));  faceLabels[faceList.back()] = Vec2i(2, 0);
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(shellObj); 
  DeformableObject::face_iter fIt;
  for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(shellObj, shellFaces, m_timestep, this);
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

void DoubleBubbleTest::serialize( std::ofstream& of )
{
  ElTopo::SurfTrack * st = mesh2surftrack();
  Recording::writeSurfTrack(of, *st);
  delete (st);
}

void DoubleBubbleTest::resumeFromfile( std::ifstream& ifs )
{
  ElTopo::SurfTrack * st = mesh2surftrack();
  Recording::readSurfTrack(ifs, *st);
  surftrack2mesh(*st);
  delete (st);
}

class Vec2iLess
{
public:
  bool operator () (const Vec2i & x, const Vec2i & y) const { return (x.x() < y.x() || (x.x() == y.x() && x.y() < y.y())); }
};

void DoubleBubbleTest::keyboard(unsigned char k, int x, int y)
{
  if (g_recording.isPlaybackOn())
  {
    if (k == '[')
    {
      int f = g_recording.currentFrame();
      g_recording.setCurrentFrame(f - 1);
    
      ElTopo::SurfTrack * st = mesh2surftrack();
      g_recording.loadRecording(*st);
      surftrack2mesh(*st);
      delete st;
      
      setTime(getDt() * (f - 1));
      glutPostRedisplay();

    } else if (k == ']')
    {
      int f = g_recording.currentFrame();
      g_recording.setCurrentFrame(f + 1);
      
      ElTopo::SurfTrack * st = mesh2surftrack();
      g_recording.loadRecording(*st);
      surftrack2mesh(*st);
      delete st;
      
      setTime(getDt() * (f + 1));
      glutPostRedisplay();
      
    } else if (k == '{')
    {
      int f = g_recording.currentFrame();
      g_recording.setCurrentFrame(f - 1);
      
      setTime(getDt() * (f - 1));
      glutPostRedisplay();
      
    } else if (k == '}')
    {
      int f = g_recording.currentFrame();
      g_recording.setCurrentFrame(f + 10);
      
      setTime(getDt() * (f + 10));
      glutPostRedisplay();
      
    } else if (k == '.')
    {
      ElTopo::SurfTrack * st = mesh2surftrack();
      g_recording.loadRecording(*st, 1);
      surftrack2mesh(*st);
      delete st;
      
      glutPostRedisplay();
    } else if (k == ',')
    {
      ElTopo::SurfTrack * st = mesh2surftrack();
      g_recording.loadRecording(*st, -1);
      surftrack2mesh(*st);
      delete st;
      
      glutPostRedisplay();
    } else if (k == '>')
    {
      int f = g_recording.currentStep();
      g_recording.setCurrentStep(f + 10);
      std::cout << "current step: " << g_recording.currentStep() << std::endl;
      glutPostRedisplay();
    } else if (k == '<')
    {
      int f = g_recording.currentStep();
      g_recording.setCurrentStep(f - 10);
      std::cout << "current step: " << g_recording.currentStep() << std::endl;
      glutPostRedisplay();
    }
    
  }
    
  if (k == 'n')
  {
    VertexHandle v = ViewController::singleton()->nearestVertex();
    EdgeHandle e = ViewController::singleton()->nearestEdge();
    FaceHandle f = ViewController::singleton()->nearestFace();
    
    DeformableObject & mesh = shell->getDefoObj();
    
    if (v.isValid())
    {
      std::cout << "VOI: " << v.idx() << " (" << shell->getVertexPosition(v) << ")" << std::endl;
      std::cout << "Region graph: " << std::endl;

      std::set<int> regionset;
      std::map<Vec2i, bool, Vec2iLess> graph;
      for (VertexFaceIterator vfit = mesh.vf_iter(v); vfit; ++vfit)
      {
        Vec2i label = shell->getFaceLabel(*vfit);
        if (label.y() < label.x()) std::swap(label.x(), label.y());
        regionset.insert(label.x());
        regionset.insert(label.y());
        graph[label] = true;
      }
      
      std::vector<int> regions;
      regions.assign(regionset.begin(), regionset.end());
      std::sort(regions.begin(), regions.end());
      std::cout << "   "; for (size_t i = 0; i < regions.size(); i++) std::cout << std::setfill(' ') << std::setw(2) << regions[i] << " "; std::cout << std::endl;
      for (size_t i = 0; i < regions.size(); i++)
      {
        std::cout << std::setfill(' ') << std::setw(2) << regions[i] << " ";
        for (size_t j = 0; j < regions.size(); j++)
        {
          Vec2i regionpair = (regions[i] < regions[j] ? Vec2i(regions[i], regions[j]) : Vec2i(regions[j], regions[i]));
          std::cout << " " << (graph[regionpair] ? "*" : " ") << " ";
        }
        std::cout << std::endl;
      }
    }
    if (e.isValid())
    {
      ElTopo::SurfTrack * st = mesh2surftrack();
      Vec3d x0 = shell->getVertexPosition(mesh.fromVertex(e));
      Vec3d x1 = shell->getVertexPosition(mesh.toVertex(e));
      ElTopo::Vec3d st_x0(x0.x(), x0.y(), x0.z());
      ElTopo::Vec3d st_x1(x1.x(), x1.y(), x1.z());
      size_t st_v0 = st->m_mesh.nv();
      size_t st_v1 = st->m_mesh.nv();
      for (size_t i = 0; i < st->m_mesh.nv(); i++)
      {
        if (mag(st->get_position(i) - st_x0) < 1e-6)
        {
          assert(st_v0 == st->m_mesh.nv());
          st_v0 = i;
        }
        if (mag(st->get_position(i) - st_x1) < 1e-6)
        {
          assert(st_v1 == st->m_mesh.nv());
          st_v1 = i;
        }
      }
      assert(st_v0 != st->m_mesh.nv());
      assert(st_v1 != st->m_mesh.nv());
      
      size_t st_e = st->m_mesh.get_edge_index(st_v0, st_v1);
      
      std::cout << "EOI: " << e.idx() << ": " << mesh.fromVertex(e).idx() << " (" << x0 << ") - " << mesh.toVertex(e).idx() << " (" << x1 << ") length = " << (x0 - x1).norm() << " feature = " << st->edge_is_feature(st_e) << std::endl;
      delete st;
    }
    if (f.isValid())
    {
      FaceVertexIterator fvit = mesh.fv_iter(f); assert(fvit);
      VertexHandle v0 = *fvit; ++fvit; assert(fvit);
      VertexHandle v1 = *fvit; ++fvit; assert(fvit);
      VertexHandle v2 = *fvit; ++fvit; assert(!fvit);
      std::cout << "FOI: " << f.idx() << ": " << v0.idx() << " (" << shell->getVertexPosition(v0) << "), " << v1.idx() << " (" << shell->getVertexPosition(v1) << "), " << v2.idx() << " (" << shell->getVertexPosition(v2) << ")" << std::endl;
    }
  }
  
}

std::vector<std::string> explode(const std::string & s, char delim) 
{
  std::vector<std::string> elems;
  std::stringstream ss(s);
  std::string item;
  while(std::getline(ss, item, delim)) 
    elems.push_back(item);
  return elems;
}

void DoubleBubbleTest::setupScene12() 
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

  std::vector<Vec3d> obj_vs;
  std::vector<Vec3i> obj_fs;
//  std::ifstream objfile("assets/doublebubbletest/zalesak.obj");
  std::ifstream objfile("assets/doublebubbletest/zalesak_multiphase.obj");
  
  // OBJ loader
  while (!objfile.eof())
  {
    std::string line;
    std::getline(objfile, line);
    std::stringstream liness(line);
    std::string ins;
    liness >> ins;
    if (ins == "v")
    {
      Vec3d v;
      liness >> v.x() >> v.y() >> v.z();
      obj_vs.push_back(v);
    } else if (ins == "f")
    {
      Vec3i f;
      for (int i = 0; i < 3; i++)
      {
        std::string s;
        liness >> s;
        std::vector<std::string> elems = explode(s, '/');
        assert(elems.size() == 1 || elems.size() == 3);
        std::stringstream ss(elems[0]);
        ss >> f[i];
        f[i]--;
      }
      obj_fs.push_back(f);
    }
  }
  objfile.close();
  
  std::cout << "OBJ load report: nv = " << obj_vs.size() << " nf = " << obj_fs.size() << std::endl;
  
  //create the mesh
  std::vector<VertexHandle> vertList;
  for (size_t i = 0; i < obj_vs.size(); i++) 
  {
    vertList.push_back(shellObj->addVertex());
    velocities[vertList[i]] = Vec3d(0, 0, 0);
    positions[vertList[i]] = obj_vs[i] / 200 + Vec3d(0.5, 0.5, 0.5);  // scale the model back to be within [0, 1]. note that remeshing resolution should be at least 0.1 to be high enough to preserve the sharp features in this model.
    undeformed[vertList[i]] = positions[vertList[i]];
  }
  
  std::vector<FaceHandle> faceList;
  FaceProperty<Vec2i> faceLabels(shellObj); //label face regions to do volume constrained bubbles  
  
  for (size_t i = 0; i < obj_fs.size(); i++)
  {
    faceList.push_back(shellObj->addFace(vertList[obj_fs[i].x()], vertList[obj_fs[i].y()], vertList[obj_fs[i].z()]));
    faceLabels[faceList.back()] = Vec2i(0, 1);  // placeholder here; real labels will be set below
  }
  
  // put a few seeds for each region, and decide face labels using visibility. for all convex regions, one seed for each region is enough.
  // regions unlabeled in the end will be assigned label 0.
  std::vector<std::pair<Vec3d, int> > seeds;
  seeds.push_back(std::pair<Vec3d, int>(Vec3d(0.5, 0.4, 0.3), 1)); 
  seeds.push_back(std::pair<Vec3d, int>(Vec3d(0.2, 0.4, 0.5), 2)); 
  seeds.push_back(std::pair<Vec3d, int>(Vec3d(0.5, 0.4, 0.7), 3)); 
  seeds.push_back(std::pair<Vec3d, int>(Vec3d(0.5, 0.6, 0.3), 4)); 
  seeds.push_back(std::pair<Vec3d, int>(Vec3d(0.2, 0.6, 0.5), 5)); 
  seeds.push_back(std::pair<Vec3d, int>(Vec3d(0.5, 0.6, 0.7), 6));
  
  // the loops here are unoptimized.
  for (size_t i = 0; i < obj_fs.size(); i++)
  {
    Vec3d & x0 = positions[vertList[obj_fs[i].x()]];
    Vec3d & x1 = positions[vertList[obj_fs[i].y()]];
    Vec3d & x2 = positions[vertList[obj_fs[i].z()]];
    Vec3d c = (x0 + x1 + x2) / 3;
    
    faceLabels[faceList[i]] = Vec2i(0, 0);
    
    for (size_t j = 0; j < seeds.size(); j++)
    {
      // reference for line-triangle intersection test:
      //  http://geomalgorithms.com/a06-_intersect-2.html
      
      Vec3d & seed = seeds[j].first;
      bool collision = false;
      
      for (size_t k = 0; k < obj_fs.size(); k++)
      {
        if (k == i)
          continue;
        
        Vec3d & v0 = positions[vertList[obj_fs[k].x()]];
        Vec3d & v1 = positions[vertList[obj_fs[k].y()]];
        Vec3d & v2 = positions[vertList[obj_fs[k].z()]];
        
        Vec3d u = v1 - v0;
        Vec3d v = v2 - v0;
        Vec3d n = u.cross(v);
      
        Scalar r = n.dot(v0 - c) / n.dot(seed - c);
        if (r > 1 || r < 0)
          continue;
        
        Vec3d p = c + (seed - c) * r;
        Vec3d w = p - v0;
        Scalar s = (u.dot(v) * w.dot(v) - v.dot(v) * w.dot(u)) / (u.dot(v) * u.dot(v) - u.dot(u) * v.dot(v));
        Scalar t = (u.dot(v) * w.dot(u) - u.dot(u) * w.dot(v)) / (u.dot(v) * u.dot(v) - u.dot(u) * v.dot(v));
        
        if (s >= 0 && s <= 1 && t >= 0 && t <= 1 && 1 - s - t >= 0)
        {
          collision = true;
          break;
        }
      }
      
      if (!collision)
      {
        faceLabels[faceList[i]][(x1 - x0).cross(x2 - x0).dot(seed - c) > 0 ? 1 : 0] = seeds[j].second;
      }
    }
  }
    
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(shellObj); 
  DeformableObject::face_iter fIt;
  for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(shellObj, shellFaces, m_timestep, this);
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

void DoubleBubbleTest::setupScene13()
{
  //vertices
  VertexProperty<Vec3d> undeformed(shellObj);
  VertexProperty<Vec3d> positions(shellObj);
  VertexProperty<Vec3d> velocities(shellObj);
  
  //edge properties
  EdgeProperty<Scalar> undefAngle(shellObj);
  EdgeProperty<Scalar> edgeAngle(shellObj);
  EdgeProperty<Scalar> edgeVel(shellObj);
  
  //create a sphere
  std::vector<VertexHandle> vertList;
  
  int N = GetIntOpt("shell-x-resolution");
  Scalar r = 0.1;
  Vec3d c = Vec3d(0.4, 0.6, 0.6);
  vertList.push_back(shellObj->addVertex());
  positions[vertList.back()] = Vec3d(c + Vec3d(0, 0, 0));
  vertList.push_back(shellObj->addVertex());
  positions[vertList.back()] = Vec3d(c + Vec3d(0, 0, r));  
  vertList.push_back(shellObj->addVertex());
  positions[vertList.back()] = Vec3d(c + Vec3d(0, r * 2 * sqrt(2) / 3, -r / 3));
  vertList.push_back(shellObj->addVertex());
  positions[vertList.back()] = Vec3d(c + Vec3d(r * sqrt(6) / 3, -r * sqrt(2) / 3, -r / 3));
  vertList.push_back(shellObj->addVertex());
  positions[vertList.back()] = Vec3d(c + Vec3d(-r * sqrt(6) / 3, -r * sqrt(2) / 3, -r / 3));
  
  for (int i = 0; i < shellObj->nv(); ++i)
  {
    velocities[vertList[i]] = Vec3d(0, 0, 0);
    undeformed[vertList[i]] = positions[vertList[i]];
  }
  
  std::vector<FaceHandle> faceList;
  FaceProperty<Vec2i> faceLabels(shellObj); //label face regions to do volume constrained bubbles  
  
  faceList.push_back(shellObj->addFace(vertList[ 0], vertList[ 1], vertList[ 2]));
  faceLabels[faceList.back()] = Vec2i(1, 2);
  faceList.push_back(shellObj->addFace(vertList[ 0], vertList[ 2], vertList[ 3]));
  faceLabels[faceList.back()] = Vec2i(1, 3);
  faceList.push_back(shellObj->addFace(vertList[ 0], vertList[ 3], vertList[ 1]));
  faceLabels[faceList.back()] = Vec2i(1, 4);
  faceList.push_back(shellObj->addFace(vertList[ 0], vertList[ 1], vertList[ 4]));
  faceLabels[faceList.back()] = Vec2i(2, 4);
  faceList.push_back(shellObj->addFace(vertList[ 0], vertList[ 2], vertList[ 4]));
  faceLabels[faceList.back()] = Vec2i(3, 2);
  faceList.push_back(shellObj->addFace(vertList[ 0], vertList[ 3], vertList[ 4]));
  faceLabels[faceList.back()] = Vec2i(4, 3);
  faceList.push_back(shellObj->addFace(vertList[ 1], vertList[ 2], vertList[ 3]));
  faceLabels[faceList.back()] = Vec2i(0, 1);
  faceList.push_back(shellObj->addFace(vertList[ 1], vertList[ 2], vertList[ 4]));
  faceLabels[faceList.back()] = Vec2i(2, 0);
  faceList.push_back(shellObj->addFace(vertList[ 2], vertList[ 3], vertList[ 4]));
  faceLabels[faceList.back()] = Vec2i(3, 0);
  faceList.push_back(shellObj->addFace(vertList[ 1], vertList[ 3], vertList[ 4]));
  faceLabels[faceList.back()] = Vec2i(0, 4);
  
  // subdivision
  for (int i = 0; i < N; i++)
  {
    EdgeProperty<VertexHandle> mid(shellObj);
    for (EdgeIterator eit = shellObj->edges_begin(); eit != shellObj->edges_end(); ++eit)
    {
      mid[*eit] = shellObj->addVertex();
      positions[mid[*eit]] = (positions[shellObj->fromVertex(*eit)] + positions[shellObj->toVertex(*eit)]) / 2;
      positions[mid[*eit]] = c + (positions[mid[*eit]] - c).normalized() * r;
    }
    
    std::vector<FaceHandle> oldfaces;
    for (FaceIterator fit = shellObj->faces_begin(); fit != shellObj->faces_end(); ++fit)
      oldfaces.push_back(*fit);
    
    for (size_t j = 0; j < oldfaces.size(); j++)
    {
      FaceHandle fh = oldfaces[j];
      FaceVertexIterator fvit = shellObj->fv_iter(fh); assert(fvit);
      VertexHandle v0 = *fvit; ++fvit; assert(fvit);
      VertexHandle v1 = *fvit; ++fvit; assert(fvit);
      VertexHandle v2 = *fvit; ++fvit; assert(!fvit);
      EdgeHandle e01 = findEdge(*shellObj, v0, v1); assert(e01.isValid());
      EdgeHandle e12 = findEdge(*shellObj, v1, v2); assert(e01.isValid());
      EdgeHandle e20 = findEdge(*shellObj, v2, v0); assert(e01.isValid());
      VertexHandle m01 = mid[e01];
      VertexHandle m12 = mid[e12];
      VertexHandle m20 = mid[e20];
      Vec2i label = faceLabels[fh];
      
      FaceHandle nf;
      if (v0 == vertList[0])
      {
        nf = shellObj->addFace(v0, v1, m12); faceLabels[nf] = label;
        nf = shellObj->addFace(v0, m12, v2); faceLabels[nf] = label;
      } else if (v1 == vertList[0])
      {
        nf = shellObj->addFace(v1, v2, m20); faceLabels[nf] = label;
        nf = shellObj->addFace(v1, m20, v0); faceLabels[nf] = label;
      } else if (v2 == vertList[0])
      {
        nf = shellObj->addFace(v2, v0, m01); faceLabels[nf] = label;
        nf = shellObj->addFace(v2, m01, v1); faceLabels[nf] = label;
      } else
      {
        nf = shellObj->addFace(v0, m01, m20); faceLabels[nf] = label;
        nf = shellObj->addFace(m01, v1, m12); faceLabels[nf] = label;
        nf = shellObj->addFace(m20, m12, v2); faceLabels[nf] = label;
        nf = shellObj->addFace(m12, m20, m01); faceLabels[nf] = label;
      }
      
      shellObj->deleteFace(fh, true);
    }
  }
  
  //create a face property to flag which of the faces are part of the object. (All of them, in this case.)
  FaceProperty<char> shellFaces(shellObj);
  DeformableObject::face_iter fIt;
  for(fIt = shellObj->faces_begin(); fIt != shellObj->faces_end(); ++fIt)
    shellFaces[*fIt] = true;
  
  //now create the physical model to hang on the mesh
  shell = new ElasticShell(shellObj, shellFaces, m_timestep, this);
  shellObj->addModel(shell);
  
  //positions
  //  shell->setVertexUndeformed(undeformed);
  shell->setVertexPositions(positions);
  shell->setVertexVelocities(velocities);
  
  shell->setFaceLabels(faceLabels);
}

void DoubleBubbleTest::collapse(const ElTopo::SurfTrack & st, size_t e)
{
//    std::cout << "collapse---" << std::endl;
    g_recording.log() << "Collapse edge " << e << std::endl;
    g_recording.recordSurfTrack(st);
}

void DoubleBubbleTest::split(const ElTopo::SurfTrack & st, size_t e)
{
//    std::cout << "split---" << std::endl;
    g_recording.log() << "Split edge " << e << std::endl;
    g_recording.recordSurfTrack(st);
}

void DoubleBubbleTest::flip(const ElTopo::SurfTrack & st, size_t e)
{
//    std::cout << "flip---" << std::endl;
    g_recording.log() << "Flip edge " << e << std::endl;
    g_recording.recordSurfTrack(st);
}

void DoubleBubbleTest::t1(const ElTopo::SurfTrack & st, size_t e)
{
//    std::cout << "t1---" << std::endl;
    g_recording.log() << "T1 pop vertex " << e << std::endl;
    g_recording.recordSurfTrack(st);
}

void DoubleBubbleTest::facesplit(const ElTopo::SurfTrack & st, size_t f)
{
//    std::cout << "face split---" << std::endl;
    g_recording.log() << "Split face " << f << std::endl;
    g_recording.recordSurfTrack(st);
}

void DoubleBubbleTest::snap(const ElTopo::SurfTrack & st, size_t v0, size_t v1)
{
//    std::cout << "snap---" << std::endl;
    g_recording.log() << "Snap vertex " << v0 << " to vertex " << v1 << std::endl;
    g_recording.recordSurfTrack(st);
}

void DoubleBubbleTest::smoothing(const ElTopo::SurfTrack & st)
{
//    std::cout << "smooth---" << std::endl;
    g_recording.log() << "Null space smoothing " << std::endl;
    g_recording.recordSurfTrack(st);
}

std::ostream & DoubleBubbleTest::log()
{
    return g_recording.log();
}
