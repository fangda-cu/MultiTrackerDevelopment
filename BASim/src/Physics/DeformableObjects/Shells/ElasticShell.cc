#include "BASim/src/Collisions/ElTopo/broadphasegrid.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShell.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"
#include "BASim/src/Collisions/ElTopo/ccd_wrapper.hh"
//#include "BASim/src/Physics/DeformableObjects/Shells/CSTMembraneForce.hh"
#include "BASim/src/Math/Math.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellVertexPointSpringForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellStickyRepulsionForce.hh"
#include "BASim/src/Collisions/ElTopo/collisionqueries.hh"

#include "BASim/src/Collisions/ElTopo/array3.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"

#include "surftrack.h"
#include "subdivisionscheme.h"

#include <algorithm>
#include <numeric>

namespace BASim {
  
ElasticShell::ElasticShell(DeformableObject* object, const FaceProperty<char>& shellFaces, Scalar timestep, SteppingCallback * stepping_callback, int doublebubble_scene) :
  PhysicalModel(*object), m_obj(object), 
    m_active_faces(shellFaces), 
    m_face_regions(object),
    m_vertex_constraint_labels(object),
    m_stepping_callback(stepping_callback),
    m_doublebubble_scene(doublebubble_scene)
{
}

ElasticShell::~ElasticShell() {
  
  //this includes deletion of member spring/repulsion forces
  for(unsigned int i = 0; i < m_shell_forces.size(); ++i)
    delete m_shell_forces[i];
  
}

void ElasticShell::computeConservativeForcesEnergy( VecXd& force , Scalar& energy)
{
  const std::vector<ElasticShellForce*>& forces = getForces();
  std::vector<ElasticShellForce*>::const_iterator fIt;
  
  VecXd curr_force(force.size());
  for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
    curr_force.setZero();
    (*fIt)->globalForce(curr_force);
    energy += (*fIt)->globalEnergy();
    force += curr_force;
  }

}

void ElasticShell::computeForces( VecXd& force )
{
  const std::vector<ElasticShellForce*>& forces = getForces();
  std::vector<ElasticShellForce*>::const_iterator fIt;
  
  VecXd curr_force(force.size());
  for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
    curr_force.setZero();
    (*fIt)->globalForce(curr_force);
    
    force += curr_force;
  }

}

void ElasticShell::computeJacobian( Scalar scale, MatrixBase& J )
{
  const std::vector<ElasticShellForce*>& forces = getForces();
  std::vector<ElasticShellForce*>::const_iterator fIt;

  for (fIt = forces.begin(); fIt != forces.end(); ++fIt)
    (*fIt)->globalJacobian(scale, J);
}

const std::vector<ElasticShellForce*>& ElasticShell::getForces() const
{
  return m_shell_forces;
}

void ElasticShell::addForce( ElasticShellForce* force )
{
  assert(force != NULL);

  m_shell_forces.push_back(force);
}

Vec3d ElasticShell::getFaceNormal(const FaceHandle& f) {
    std::vector<Vec3d> v;
    const DeformableObject& mesh = *m_obj;
    for( FaceVertexIterator fvit = mesh.fv_iter(f); fvit; ++fvit )
    {
        v.push_back(getVertexPosition(*fvit));
    }
    Vec3d n = (v[1] - v[0]).cross(v[2]-v[0]);
    n.normalize();
    return n;
}

void ElasticShell::getFaceNormals(FaceProperty<Vec3d> & fNormals) const{
    const DeformableObject& mesh = *m_obj;
    for( FaceIterator fit = mesh.faces_begin(); fit != mesh.faces_end(); ++fit ){
        std::vector<Vec3d> v;
        for( FaceVertexIterator fvit = mesh.fv_iter(*fit); fvit; ++fvit )
        {
          v.push_back(getVertexPosition(*fvit));
        }
        Vec3d n = (v[1] - v[0]).cross(v[2]-v[0]);
        n.normalize();
        fNormals[*fit] = n;
    }
}

void ElasticShell::getVertexNormals(VertexProperty<Vec3d> & vNormals) const{
    FaceProperty<Vec3d> fNormals(& getDefoObj());
    getFaceNormals(fNormals);
    DeformableObject& mesh = *m_obj;
    Scalar w = 0.0;
    for ( VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit){
        vNormals[*vit] = Vec3d(0.0,0.0,0.0);
        for ( VertexFaceIterator vfit = mesh.vf_iter(*vit); vfit; ++vfit){
            w = getArea(*vfit);
            vNormals[*vit] += w*fNormals[*vfit];
        }
        vNormals[*vit].normalize();
    }
}

Scalar ElasticShell::getArea(const FaceHandle& f, bool current) const  {
  FaceVertexIterator fvit = m_obj->fv_iter(f);
  VertexHandle v0_hnd = *fvit; ++fvit; assert(fvit);
  VertexHandle v1_hnd = *fvit; ++fvit; assert(fvit);
  VertexHandle v2_hnd = *fvit; ++fvit; assert(!fvit);

  //compute triangle areas
//  if(current) 
  assert(current);
  {
    Vec3d pos0 = m_obj->getVertexPosition(v0_hnd);
    Vec3d pos1 = m_obj->getVertexPosition(v1_hnd);
    Vec3d pos2 = m_obj->getVertexPosition(v2_hnd);
    
    Vec3d v0 = pos1 - pos0;
    Vec3d v1 = pos2 - pos0;
    Vec3d triVec = v0.cross(v1);
    return 0.5*triVec.norm();
  }
}

bool ElasticShell::isVertexActive( const VertexHandle& v ) const
{
  //determine if the vertex is on any active face
  VertexFaceIterator vf = m_obj->vf_iter(v);
  for(;vf; ++vf) {
    if(isFaceActive(*vf)) {
      return true;
    }
  }
  
  return false;
}

bool ElasticShell::isEdgeActive( const EdgeHandle& e) const {
  //if any adjacent face is active, we say this edge is active.
  EdgeFaceIterator ef = m_obj->ef_iter(e);
  for(;ef;++ef) {
    if(isFaceActive(*ef)) {
      return true;
    }
  }
  
  return false;
}

const Scalar& ElasticShell::getDof( const DofHandle& hnd ) const
{
  assert("No dof");
  //Just to satisfy my compiler - CB
#ifdef _MSC_VER
#pragma warning( disable : 4172 )
#endif
  Scalar dummy;
  return const_cast<Scalar&>(dummy);
}

void ElasticShell::setDof( const DofHandle& hnd, const Scalar& dof )
{
  assert("No dof");
}

const Scalar& ElasticShell::getVel( const DofHandle& hnd ) const
{
  assert("No dof");
  //Just to satisfy my compiler - CB
#ifdef _MSC_VER
#pragma warning( disable : 4172 )
#endif
  Scalar dummy;
  return const_cast<Scalar&>(dummy);
}

void ElasticShell::setVel( const DofHandle& hnd, const Scalar& vel )
{
  assert("No dof");
}

void ElasticShell::getScriptedDofs( IntArray& dofIndices, std::vector<Scalar>& dofValues, Scalar time ) const
{
}

void ElasticShell::startStep(Scalar time, Scalar timestep)
{
  //tell the forces to update anything they need to update
  const std::vector<ElasticShellForce*>& forces = getForces();
  for(unsigned int i = 0; i < forces.size(); ++i) {
    forces[i]->update();
  }
  
  m_stepping_callback->afterStartStep();
}

void ElasticShell::resolveCollisions(Scalar timestep)
{
  std::cout << "Resolving collisions with El Topo\n";
  //Convert the data to the form required by El Topo!
  std::vector<ElTopo::Vec3d> vert_new, vert_old;
  std::vector<ElTopo::Vec3st> tri_data;
  std::vector<ElTopo::Vec3d> masses;

  DeformableObject& mesh = getDefoObj();

  //Index mappings between us and El Topo
  VertexProperty<int> vert_numbers(&mesh);
  FaceProperty<int> face_numbers(&mesh);
  std::vector<VertexHandle> reverse_vertmap;
  std::vector<FaceHandle> reverse_trimap;

  //walk through vertices, create linear list, store numbering
  int id = 0;
  for(VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit) {
    VertexHandle vh = *vit;
    Vec3d vert = getVertexPosition(vh);
    Vec3d old_vert = getVertexDampingUndeformed(vh);
//    Scalar mass = getMass(vh);
    ElTopo::Vec3d mass(1, 1, 1);

    vert_new.push_back(ElTopo::Vec3d(vert[0], vert[1], vert[2]));
    vert_old.push_back(ElTopo::Vec3d(old_vert[0], old_vert[1], old_vert[2]));
    for (int i = 0; i < 3; i++)
      if (getDefoObj().isConstrainedInDirection(vh, i))
        mass[i] = numeric_limits<Scalar>::infinity();
    masses.push_back(mass);
    vert_numbers[vh] = id;
    reverse_vertmap.push_back(vh);

    ++id;
  }

  //walk through triangles, creating linear list, using the vertex numbering assigned above
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
    face_numbers[fh] = id;
    reverse_trimap.push_back(fh);
    ++id;
  }


  // build a DynamicSurface
  Scalar friction_coeff = 0;
  Scalar mean_edge_len = (m_et_remesh_edge_min_len + m_et_remesh_edge_max_len) / 2;
  ElTopo::DynamicSurface dynamic_surface( vert_old, tri_data, std::vector<ElTopo::Vec2i>(tri_data.size(), ElTopo::Vec2i(0, 0)), masses, m_et_collision_epsilon_fraction * mean_edge_len, friction_coeff, true, false );

  dynamic_surface.set_all_newpositions( vert_new );
    
//    for (size_t i = 0; i < vert_old.size(); i++)
//        std::cout << "old vertex " << i << ": " << vert_old[i] << std::endl;
//    for (size_t i = 0; i < vert_new.size(); i++)
//        std::cout << "new vertex " << i << ": " << vert_new[i] << std::endl;
  
  // advance by dt
  double actual_dt;
  dynamic_surface.integrate( timestep, actual_dt );
  if(actual_dt != timestep)
    std::cout << "XXXXXXXXX Failed to step the full length of the recommended step!XXXXX\n";
  // the dt used may be different than specified (if we cut the time step)
  
  //figure out what the actual velocities were, and update the mesh data
  for(unsigned int i = 0; i < vert_new.size(); ++i) {
    VertexHandle vh = reverse_vertmap[i];
    ElTopo::Vec3d pos = dynamic_surface.get_position(i);
    ElTopo::Vec3d vel = (pos - vert_old[i]) / timestep;
    Vec3d new_pos(pos[0], pos[1], pos[2]);
    Vec3d new_vel(vel[0], vel[1], vel[2]);
    setVertexPosition(vh, new_pos);
    setVertexVelocity(vh, new_vel);

    if(isnan(pos[0]) || isnan(pos[1]) || isnan(pos[2]))
      std::cout << "ElTopo Failed: NaN vertex\n";
    if(isinf(pos[0]) || isinf(pos[1]) || isinf(pos[2]))
      std::cout << "ElTopo Failed: Inf vertex\n";
  }

}


void ElasticShell::endStep(Scalar time, Scalar timestep) {
    
//    for (VertexIterator v = m_obj->vertices_begin(); v != m_obj->vertices_end(); ++v)
//    {
//        assert(getVertexPosition(*v) == getVertexPosition(*v));
//        std::cout << "vertex " << (*v).idx() << ": " << getVertexPosition(*v) << std::endl;
//    }

  if (m_stepping_callback)
    m_stepping_callback->beforeEndStep();
  
  // remove faces completely inside BB walls
  for (FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit)
  {
    FaceVertexIterator fvit = m_obj->fv_iter(*fit); assert(fvit);
    Vec3d x0 = getVertexPosition(*fvit); ++fvit; assert(fvit);
    Vec3d x1 = getVertexPosition(*fvit); ++fvit; assert(fvit);
    Vec3d x2 = getVertexPosition(*fvit); ++fvit; assert(!fvit);
    
    //    if (x0.y() == 0 || x1.y() == 0 || x2.y() == 0)
    //      std::cout << "face: " << x0 << " " << x1 << " " << x2 << std::endl;
    
    int w0 = onBBWall(x0);
    int w1 = onBBWall(x1);
    int w2 = onBBWall(x2);
    if (((w0 & w1) & w2) != 0)
    {
      //      std::cout << "face: " << x0 << " " << x1 << " " << x2 << std::endl;
      m_obj->deleteFace(*fit, false);
    }
  }
  
  // prune orphan edges and vertices
  for (EdgeIterator eit = m_obj->edges_begin(); eit != m_obj->edges_end(); ++eit)
    if (m_obj->edgeIncidentFaces(*eit) == 0)
      m_obj->deleteEdge(*eit, true);
  
  for (VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit)
    if (m_obj->vertexIncidentEdges(*vit) == 0)
      m_obj->deleteVertex(*vit);
  
  
  std::cout << "Starting endStep.\n";
  bool do_relabel = false;

  std::cout << "Vertex count: " << m_obj->nv() << std::endl;

  if (m_obj->nf() == 0) // ElTopo crashes if given a mesh with zero triangles, which is possible for example when we're only testing rods.
    return;
//  return;
  
  //El Topo collision processing.
  
//  resolveCollisions(timestep);

  //Remeshing
  if(m_do_remeshing) {
    std::cout << "Remeshing\n";
    remesh(timestep);

    //Relabel DOFs if necessary
    do_relabel = true;
  }
  
  std::cout << "nv = " << getDefoObj().nv() << " ne = " << getDefoObj().ne() << " nf = " << getDefoObj().nf() << " nt = " << getDefoObj().nt() << std::endl;

  static int only_twice = 0;

  if(do_relabel) {
    m_obj->computeDofIndexing();
  }

  //Update masses based on new areas/thicknesses
//  computeMasses();

  std::cout << "Completed endStep\n";

}

void ElasticShell::remesh(Scalar timestep, bool initial)
{
  // remove faces completely inside BB walls
  for (FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit)
  {
    FaceVertexIterator fvit = m_obj->fv_iter(*fit); assert(fvit);
    Vec3d x0 = getVertexPosition(*fvit); ++fvit; assert(fvit);
    Vec3d x1 = getVertexPosition(*fvit); ++fvit; assert(fvit);
    Vec3d x2 = getVertexPosition(*fvit); ++fvit; assert(!fvit);
    
    //    if (x0.y() == 0 || x1.y() == 0 || x2.y() == 0)
    //      std::cout << "face: " << x0 << " " << x1 << " " << x2 << std::endl;
    
    int w0 = onBBWall(x0);
    int w1 = onBBWall(x1);
    int w2 = onBBWall(x2);
    if (((w0 & w1) & w2) != 0)
    {
      //      std::cout << "face: " << x0 << " " << x1 << " " << x2 << std::endl;
      m_obj->deleteFace(*fit, false);
    }
  }
  
  // prune orphan edges and vertices
  for (EdgeIterator eit = m_obj->edges_begin(); eit != m_obj->edges_end(); ++eit)
    if (m_obj->edgeIncidentFaces(*eit) == 0)
      m_obj->deleteEdge(*eit, true);
  
  for (VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit)
    if (m_obj->vertexIncidentEdges(*vit) == 0)
      m_obj->deleteVertex(*vit);
  
    
  //Set up a SurfTrack, run remeshing, render the new mesh
  ElTopo::SurfTrackInitializationParameters construction_parameters;
//  construction_parameters.m_proximity_epsilon = m_collision_epsilon;
//  construction_parameters.m_merge_proximity_epsilon = 0.02 * m_remesh_edge_min_len;
//  construction_parameters.m_allow_vertex_movement_during_collapse = true;
//  construction_parameters.m_perform_smoothing = false;
//  construction_parameters.m_min_edge_length = m_remesh_edge_min_len;
//  construction_parameters.m_max_edge_length = m_remesh_edge_max_len;
//  construction_parameters.m_max_volume_change = 1e-4 * m_remesh_edge_min_len * m_remesh_edge_min_len * m_remesh_edge_min_len;
//  construction_parameters.m_min_triangle_angle = initial ? 0 : 3;
//  construction_parameters.m_max_triangle_angle = initial ? 180 : 177;
//  construction_parameters.m_large_triangle_angle_to_split = 160;
//  construction_parameters.m_min_triangle_area = 0.02*m_remesh_edge_min_len*m_remesh_edge_min_len;
//  construction_parameters.m_verbose = false;
//  construction_parameters.m_allow_non_manifold = true;
//  construction_parameters.m_allow_topology_changes = true;
//  construction_parameters.m_collision_safety = true;
//  construction_parameters.m_remesh_boundaries = true;
//  construction_parameters.m_t1_transition_enabled = m_remesh_t1transition;
//  construction_parameters.m_velocity_field_callback = NULL;
//  construction_parameters.m_pull_apart_distance = m_remesh_edge_min_len / 2;// (initial ? 0.1 : 0.02) * m_remesh_edge_min_len;
    
    Scalar mean_edge_len = (m_et_remesh_edge_min_len + m_et_remesh_edge_max_len) / 2;
    
    construction_parameters.m_proximity_epsilon = m_et_collision_epsilon_fraction * mean_edge_len;
    construction_parameters.m_merge_proximity_epsilon = m_et_merge_proximity_epsilon_fraction * mean_edge_len;
    construction_parameters.m_allow_vertex_movement_during_collapse = true;
    construction_parameters.m_perform_smoothing = m_et_perform_smoothing;
    construction_parameters.m_min_edge_length = m_et_remesh_edge_min_len;
    construction_parameters.m_max_edge_length = m_et_remesh_edge_max_len;
    construction_parameters.m_max_volume_change = m_et_max_volume_change_fraction * pow(mean_edge_len, 3);
    construction_parameters.m_min_triangle_angle = initial ? 0 : m_et_min_triangle_angle;
    construction_parameters.m_max_triangle_angle = initial ? 180 : m_et_max_triangle_angle;
    construction_parameters.m_large_triangle_angle_to_split = m_et_large_triangle_angle_to_split;
    construction_parameters.m_min_triangle_area = m_et_min_triangle_area_fraction * pow(mean_edge_len, 2);
    construction_parameters.m_verbose = false;
    construction_parameters.m_allow_non_manifold = true;
    construction_parameters.m_allow_topology_changes = true;
    construction_parameters.m_collision_safety = false;
    construction_parameters.m_remesh_boundaries = true;
    construction_parameters.m_t1_transition_enabled = m_et_t1_transition_enabled;
    construction_parameters.m_velocity_field_callback = NULL;
    construction_parameters.m_pull_apart_distance = m_et_t1_pull_apart_distance_fraction * mean_edge_len;
  
    if (m_doublebubble_scene == 14)
        construction_parameters.m_velocity_field_callback = this; // this is only turned on for specific scenes
  
    if (m_et_smooth_subdivision)
        construction_parameters.m_subdivision_scheme = new ElTopo::ModifiedButterflyScheme();
    else
        construction_parameters.m_subdivision_scheme = new ElTopo::MidpointScheme();
  //construction_parameters.m_subdivision_scheme = new ElTopo::QuadraticErrorMinScheme();
  //construction_parameters.m_subdivision_scheme = new ElTopo::ButterflyScheme();
  //construction_parameters.m_subdivision_scheme = new ElTopo::ModifiedButterflyScheme();

  construction_parameters.m_use_curvature_when_collapsing = false;
  construction_parameters.m_use_curvature_when_splitting = false;
  //construction_parameters.m_max_curvature_multiplier = 1000;
  //construction_parameters.m_min_curvature_multiplier = 1.0;

  std::vector<ElTopo::Vec3d> vert_data;
  std::vector<ElTopo::Vec3d> vert_vel;
  std::vector<ElTopo::Vec3st> tri_data;
  std::vector<ElTopo::Vec2i> tri_labels;
  std::vector<ElTopo::Vec3d> masses;

  DeformableObject& mesh = getDefoObj();
  
  //Index mappings between us and El Topo, used in remesh()
  VertexProperty<int> vert_numbers(&mesh);
  FaceProperty<int> face_numbers(&mesh);
  std::vector<VertexHandle> reverse_vertmap;
  std::vector<FaceHandle> reverse_trimap;
  
  reverse_vertmap.reserve(m_obj->nv());
  reverse_trimap.reserve(m_obj->nt());  
  
  vert_data.reserve(m_obj->nv());
  tri_data.reserve(m_obj->nt());
  tri_labels.reserve(m_obj->nt());
  masses.reserve(m_obj->nv());

  //walk through vertices, create linear list, store numbering
  int id = 0;
  for(VertexIterator vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit) {
    VertexHandle vh = *vit;
    Vec3d vert = getVertexPosition(vh);

    vert_data.push_back(ElTopo::Vec3d(vert[0], vert[1], vert[2]));
    Vec3d vel = getVertexVelocity(vh);
    vert_vel.push_back(ElTopo::Vec3d(vel[0], vel[1], vel[2]) * timestep);
      
    ElTopo::Vec3d mass(1, 1, 1);
    assert(getDefoObj().isConstrained(vh) == (getVertexConstraintLabel(vh) != 0));
    for (int i = 0; i < 3; i++)
      if (getDefoObj().isConstrainedInDirection(vh, i))
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
    tri_labels.push_back(ElTopo::Vec2i(m_face_regions[fh].x(), m_face_regions[fh].y()));
    face_numbers[fh] = id;
    reverse_trimap.push_back(fh);
    ++id;
  }

  std::cout << "Calling surface improvement\n";
  
  ElTopo::SurfTrack surface_tracker( vert_data, tri_data, tri_labels, masses, construction_parameters ); 
  surface_tracker.m_solid_vertices_callback = this;
  surface_tracker.m_mesheventcallback = m_mesheventcallback;
  surface_tracker.set_all_remesh_velocities(vert_vel);
    
//    for (size_t i = 0; i < vert_data.size(); i++)
//        std::cout << "vertex " << i << ": " << vert_data[i] << std::endl;
    
    //remove faces that are completely within a BB wall (equivalent to a flap face if BB walls are triangulated). these faces result in collision handling difficulties when they collide within BB walls.
    for (size_t i = 0; i < surface_tracker.m_mesh.nt(); i++)
    {
        ElTopo::Vec3st tri = surface_tracker.m_mesh.get_triangle(i);
        if (tri[0] == tri[1] && tri[0] == tri[2])
            continue;
        
        ElTopo::Vec3d v0 = surface_tracker.get_position(tri[0]);
        ElTopo::Vec3d v1 = surface_tracker.get_position(tri[1]);
        ElTopo::Vec3d v2 = surface_tracker.get_position(tri[2]);
        
        Vec3d x0(v0[0], v0[1], v0[2]);
        Vec3d x1(v1[0], v1[1], v1[2]);
        Vec3d x2(v2[0], v2[1], v2[2]);
        
        int onwall0 = onBBWall(x0);
        int onwall1 = onBBWall(x1);
        int onwall2 = onBBWall(x2);
        if ((onwall0 & onwall1 & onwall2) != 0)
        {
            surface_tracker.remove_triangle(i);
        }
    }
  
  for(int i = 0; i < m_remeshing_iters; ++i) {
    surface_tracker.topology_changes();
    surface_tracker.improve_mesh();
  }
  
  // copy ElTopo mesh back, instead of repeating the operation history incrementally.
  // this is possible because ElasticShell doesn't keep any other information that ElTopo doesn't have
  
  for (FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit)
    m_obj->deleteFace(*fit, true);

  assert(m_obj->nv() == 0);
  assert(m_obj->ne() == 0);
  assert(m_obj->nf() == 0);
  
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
    VertexHandle v = m_obj->addVertex();
    ElTopo::Vec3d x = surface_tracker.get_position(i);
    setVertexPosition(v, Vec3d(x[0], x[1], x[2]));
    setVertexVelocity(v, Vec3d(0, 0, 0));
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
    FaceHandle f = m_obj->addFace(reverse_vertmap[tri[0]], reverse_vertmap[tri[1]], reverse_vertmap[tri[2]]);
    setFaceLabel(f, Vec2i(surface_tracker.m_mesh.get_triangle_label(i)[0], surface_tracker.m_mesh.get_triangle_label(i)[1]));
    m_active_faces[f] = 1;
    
    face_numbers[f] = i;
    reverse_trimap[i] = f;
  }

  
  std::cout << "El Topo performed " << surface_tracker.m_mesh_change_history.size() << " improvement operations:\n";
  for(unsigned int j = 0; j < surface_tracker.m_mesh_change_history.size(); ++j) 
  {
    ElTopo::MeshUpdateEvent event = surface_tracker.m_mesh_change_history[j];
    std::cout << "Event type = " << event.m_type << std::endl;
  }

    double minangle = M_PI;
    double maxangle = 0;
    
    double minedge = 10;
    double maxedge = 0;
    
    for (size_t i = 0; i < surface_tracker.m_mesh.nt(); i++)
    {
        if (surface_tracker.m_mesh.get_triangle(i)[0] == surface_tracker.m_mesh.get_triangle(i)[1] && surface_tracker.m_mesh.get_triangle(i)[0] == surface_tracker.m_mesh.get_triangle(i)[2])
            continue;
        
        ElTopo::Vec3d v0 = surface_tracker.get_position(surface_tracker.m_mesh.get_triangle(i)[0]);
        ElTopo::Vec3d v1 = surface_tracker.get_position(surface_tracker.m_mesh.get_triangle(i)[1]);
        ElTopo::Vec3d v2 = surface_tracker.get_position(surface_tracker.m_mesh.get_triangle(i)[2]);
        
        ElTopo::Vec3d x0, x1, x2;
        double angle;
        
        x0 = v0; x1 = v1; x2 = v2;
        angle = acos((dot(x1 - x0, x1 - x0) + dot(x2 - x0, x2 - x0) - dot(x2 - x1, x2 - x1)) / 2 / mag(x1 - x0) / mag(x2 - x0));
        if (angle > maxangle) maxangle = angle;
        if (angle < minangle) minangle = angle;
        
        x0 = v1; x1 = v2; x2 = v0;
        angle = acos((dot(x1 - x0, x1 - x0) + dot(x2 - x0, x2 - x0) - dot(x2 - x1, x2 - x1)) / 2 / mag(x1 - x0) / mag(x2 - x0));
        if (angle > maxangle) maxangle = angle;
        if (angle < minangle) minangle = angle;
        
        x0 = v2; x1 = v0; x2 = v1;
        angle = acos((dot(x1 - x0, x1 - x0) + dot(x2 - x0, x2 - x0) - dot(x2 - x1, x2 - x1)) / 2 / mag(x1 - x0) / mag(x2 - x0));
        if (angle > maxangle) maxangle = angle;
        if (angle < minangle) minangle = angle;
        
        double e0 = mag(v2 - v1);
        double e1 = mag(v2 - v0);
        double e2 = mag(v1 - v0);
        
        if (e0 > maxedge) maxedge = e0;
        if (e0 < minedge) minedge = e0;
        if (e1 > maxedge) maxedge = e1;
        if (e1 < minedge) minedge = e1;
        if (e2 > maxedge) maxedge = e2;
        if (e2 < minedge) minedge = e2;
    }
    
    std::cout << "minangle = " << minangle * 180 / M_PI << " maxangle = " << maxangle * 180 / M_PI << " minedge = " << minedge << " maxedge = " << maxedge << std::endl;
    
  // remove faces completely inside BB walls
  for (FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit)
  {
    FaceVertexIterator fvit = m_obj->fv_iter(*fit); assert(fvit);
    Vec3d x0 = getVertexPosition(*fvit); ++fvit; assert(fvit);
    Vec3d x1 = getVertexPosition(*fvit); ++fvit; assert(fvit);
    Vec3d x2 = getVertexPosition(*fvit); ++fvit; assert(!fvit);
    
//    if (x0.y() == 0 || x1.y() == 0 || x2.y() == 0)
//      std::cout << "face: " << x0 << " " << x1 << " " << x2 << std::endl;
    
    int w0 = onBBWall(x0);
    int w1 = onBBWall(x1);
    int w2 = onBBWall(x2);
    if (((w0 & w1) & w2) != 0)
    {
//      std::cout << "face: " << x0 << " " << x1 << " " << x2 << std::endl;
      m_obj->deleteFace(*fit, false);
    }
  }
  
  // prune orphan edges and vertices
  for (EdgeIterator eit = m_obj->edges_begin(); eit != m_obj->edges_end(); ++eit)
    if (m_obj->edgeIncidentFaces(*eit) == 0)
      m_obj->deleteEdge(*eit, true);
  
  for (VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit)
    if (m_obj->vertexIncidentEdges(*vit) == 0)
      m_obj->deleteVertex(*vit);
  

}

int ElasticShell::onBBWall(const Vec3d & pos) const
{
  static const double WALL_THRESHOLD = 1e-6;
  
  int walls = 0;
  if (pos.x() < 0 + WALL_THRESHOLD)
    walls |= (1 << 0);
  if (pos.y() < 0 + WALL_THRESHOLD)
    walls |= (1 << 1);
  if (pos.z() < 0 + WALL_THRESHOLD)
    walls |= (1 << 2);
  if (pos.x() > 1 - WALL_THRESHOLD)
    walls |= (1 << 3);
  if (pos.y() > 1 - WALL_THRESHOLD)
    walls |= (1 << 4);
  if (pos.z() > 1 - WALL_THRESHOLD)
    walls |= (1 << 5);
  
  return walls;
}
  
Vec3d ElasticShell::enforceBBWallConstraint(const Vec3d & input, int constraints) const
{
  Vec3d output = input;
  if (constraints & (1 << 0))
    output.x() = 0;
  if (constraints & (1 << 1))
    output.y() = 0;
  if (constraints & (1 << 2))
    output.z() = 0;
  if (constraints & (1 << 3))
    output.x() = 1;
  if (constraints & (1 << 4))
    output.y() = 1;
  if (constraints & (1 << 5))
    output.z() = 1;
  
  return output;
}

ElTopo::Vec3d ElasticShell::sampleVelocity(ElTopo::Vec3d & pos)
{
  if (m_doublebubble_scene == 14)
  {
    return ElTopo::Vec3d(pos[0], -pos[1], 0);
  } else
  {
    return ElTopo::Vec3d(0, 0, 0);
  }
}
  
bool ElasticShell::sampleDirectionalDivergence(const ElTopo::Vec3d & pos, const ElTopo::Vec3d & dir, double & output)
{
  if (m_doublebubble_scene == 14)
  {
    output = dir[0] * dir[0] - dir[1] * dir[1];
    return true;
  } else
  {
    return false;
  }
}


bool ElasticShell::generate_collapsed_position(ElTopo::SurfTrack & st, size_t v0, size_t v1, ElTopo::Vec3d & pos)
{
  ElTopo::Vec3d x0 = st.get_position(v0);
  ElTopo::Vec3d x1 = st.get_position(v1);
  
  int label0 = onBBWall(Vec3d(x0[0], x0[1], x0[2]));
  int label1 = onBBWall(Vec3d(x1[0], x1[1], x1[2]));
  
  if (label0 == label1)
  {
    // on the same wall(s), prefer the one with higher max edge valence
    size_t maxedgevalence0 = 0;
    size_t maxedgevalence1 = 0;
    for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v0].size(); i++)
      if (st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v0][i]].size() > maxedgevalence0)
        maxedgevalence0 = st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v0][i]].size();
    for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v1].size(); i++)
      if (st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v1][i]].size() > maxedgevalence1)
        maxedgevalence1 = st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v1][i]].size();
    
    if (maxedgevalence0 == maxedgevalence1) // same max edge valence, use their midpoint
      pos = (x0 + x1) / 2;
    else if (maxedgevalence0 < maxedgevalence1)
      pos = x1;
    else
      pos = x0;

    return true;
    
  } else if ((label0 & ~label1) == 0)
  {
    // label0 is a proper subset of label1 (since label0 != label1)
    pos = x1;
    
    return true;
    
  } else if ((label1 & ~label0) == 0)
  {
    // label1 is a proper subset of label0
    pos = x0;
    
    return true;
    
  } else
  {
    // label0 and label1 are not subset of each other
    int newlabel = label0 | label1;
    assert(label0 != newlabel); // not subset of each other
    assert(label1 != newlabel);
    
    assert(!((label0 & (1 << 0)) != 0 && (label0 & (1 << 3)) != 0)); // can't have conflicting constraints in label0 and label1 already
    assert(!((label0 & (1 << 1)) != 0 && (label0 & (1 << 4)) != 0));
    assert(!((label0 & (1 << 2)) != 0 && (label0 & (1 << 5)) != 0));
    assert(!((label1 & (1 << 0)) != 0 && (label1 & (1 << 3)) != 0));
    assert(!((label1 & (1 << 1)) != 0 && (label1 & (1 << 4)) != 0));
    assert(!((label1 & (1 << 2)) != 0 && (label1 & (1 << 5)) != 0));
    
    bool conflict = false;
    if ((newlabel & (1 << 0)) != 0 && (newlabel & (1 << 3)) != 0) conflict = true;
    if ((newlabel & (1 << 1)) != 0 && (newlabel & (1 << 4)) != 0) conflict = true;
    if ((newlabel & (1 << 2)) != 0 && (newlabel & (1 << 5)) != 0) conflict = true;
    
    if (conflict)
    {
      // the two vertices are on opposite walls (conflicting constraints). Can't collapse this edge (which shouldn't have become a collapse candidate in the first place)
      return false;
    }
    
    pos = (x0 + x1) / 2;
    if (newlabel & (1 << 0))  pos[0] = 0;    // project the midpoint onto the constraint manifold (BB walls)
    if (newlabel & (1 << 1))  pos[1] = 0;
    if (newlabel & (1 << 2))  pos[2] = 0;
    if (newlabel & (1 << 3))  pos[0] = 1;
    if (newlabel & (1 << 4))  pos[1] = 1;
    if (newlabel & (1 << 5))  pos[2] = 1;
    
    return true;
  }
  
  return false;
}

bool ElasticShell::generate_split_position(ElTopo::SurfTrack & st, size_t v0, size_t v1, ElTopo::Vec3d & pos)
{
  pos = (st.get_position(v0) + st.get_position(v1)) / 2;
  
  return true;
}

ElTopo::Vec3c ElasticShell::generate_collapsed_solid_label(ElTopo::SurfTrack & st, size_t v0, size_t v1, const ElTopo::Vec3c & label0, const ElTopo::Vec3c & label1)
{
    ElTopo::Vec3d x0 = st.get_position(v0);
    ElTopo::Vec3d x1 = st.get_position(v1);
    
    int constraint0 = onBBWall(Vec3d(x0[0], x0[1], x0[2]));
    int constraint1 = onBBWall(Vec3d(x1[0], x1[1], x1[2]));
    
    assert(((constraint0 & (1 << 0)) || (constraint0 & (1 << 3))) == (bool)label0[0]);
    assert(((constraint0 & (1 << 1)) || (constraint0 & (1 << 4))) == (bool)label0[1]);
    assert(((constraint0 & (1 << 2)) || (constraint0 & (1 << 5))) == (bool)label0[2]);
    assert(((constraint1 & (1 << 0)) || (constraint1 & (1 << 3))) == (bool)label1[0]);
    assert(((constraint1 & (1 << 1)) || (constraint1 & (1 << 4))) == (bool)label1[1]);
    assert(((constraint1 & (1 << 2)) || (constraint1 & (1 << 5))) == (bool)label1[2]);
    
    ElTopo::Vec3c result;  // if either endpoint is constrained, the collapsed point shold be constrained. more specifically it should be on all the walls any of the two endpoints is on (implemented in generate_collapsed_position())
    int result_constraint = (constraint0 | constraint1);
    result[0] = ((result_constraint & (1 << 0)) || (result_constraint & (1 << 3)));
    result[1] = ((result_constraint & (1 << 1)) || (result_constraint & (1 << 4)));
    result[2] = ((result_constraint & (1 << 2)) || (result_constraint & (1 << 5)));
    
    return result;
}

ElTopo::Vec3c ElasticShell::generate_split_solid_label(ElTopo::SurfTrack & st, size_t v0, size_t v1, const ElTopo::Vec3c & label0, const ElTopo::Vec3c & label1)
{
    ElTopo::Vec3d x0 = st.get_position(v0);
    ElTopo::Vec3d x1 = st.get_position(v1);
    
    int constraint0 = onBBWall(Vec3d(x0[0], x0[1], x0[2]));
    int constraint1 = onBBWall(Vec3d(x1[0], x1[1], x1[2]));
    
    assert(((constraint0 & (1 << 0)) || (constraint0 & (1 << 3))) == (bool)label0[0]);
    assert(((constraint0 & (1 << 1)) || (constraint0 & (1 << 4))) == (bool)label0[1]);
    assert(((constraint0 & (1 << 2)) || (constraint0 & (1 << 5))) == (bool)label0[2]);
    assert(((constraint1 & (1 << 0)) || (constraint1 & (1 << 3))) == (bool)label1[0]);
    assert(((constraint1 & (1 << 1)) || (constraint1 & (1 << 4))) == (bool)label1[1]);
    assert(((constraint1 & (1 << 2)) || (constraint1 & (1 << 5))) == (bool)label1[2]);
    
    ElTopo::Vec3c result;  // the splitting midpoint has a positive constraint label only if the two endpoints are on a same wall (sharing a bit in their constraint bitfield representation)
    int result_constraint = (constraint0 & constraint1);
    result[0] = ((result_constraint & (1 << 0)) || (result_constraint & (1 << 3)));
    result[1] = ((result_constraint & (1 << 1)) || (result_constraint & (1 << 4)));
    result[2] = ((result_constraint & (1 << 2)) || (result_constraint & (1 << 5)));
    
    return result;
  
}

bool ElasticShell::generate_edge_popped_positions(ElTopo::SurfTrack & st, size_t oldv, const ElTopo::Vec2i & cut, ElTopo::Vec3d & pos_upper, ElTopo::Vec3d & pos_lower)
{
  ElTopo::Vec3d original_pos = st.get_position(oldv);
  int original_constraint = onBBWall(Vec3d(original_pos[0], original_pos[1], original_pos[2]));
  
  Vec3d new_pos_upper = enforceBBWallConstraint(Vec3d(pos_upper[0], pos_upper[1], pos_upper[2]), original_constraint);
  Vec3d new_pos_lower = enforceBBWallConstraint(Vec3d(pos_lower[0], pos_lower[1], pos_lower[2]), original_constraint);
  
  pos_upper = ElTopo::Vec3d(new_pos_upper.x(), new_pos_upper.y(), new_pos_upper.z());
  pos_lower = ElTopo::Vec3d(new_pos_lower.x(), new_pos_lower.y(), new_pos_lower.z());
  
  return true;
}

bool ElasticShell::generate_vertex_popped_positions(ElTopo::SurfTrack & st, size_t oldv, int A, int B, ElTopo::Vec3d & pos_a, ElTopo::Vec3d & pos_b)
{
  ElTopo::Vec3d original_pos = st.get_position(oldv);
  int original_constraint = onBBWall(Vec3d(original_pos[0], original_pos[1], original_pos[2]));
  
  Vec3d new_pos_a = enforceBBWallConstraint(Vec3d(pos_a[0], pos_a[1], pos_a[2]), original_constraint);
  Vec3d new_pos_b = enforceBBWallConstraint(Vec3d(pos_b[0], pos_b[1], pos_b[2]), original_constraint);
  
  pos_a = ElTopo::Vec3d(new_pos_a.x(), new_pos_a.y(), new_pos_a.z());
  pos_b = ElTopo::Vec3d(new_pos_b.x(), new_pos_b.y(), new_pos_b.z());
  
  return true;
}

bool ElasticShell::solid_edge_is_feature(const ElTopo::SurfTrack & st, size_t e)
{
  ElTopo::Vec3d x0 = st.get_position(st.m_mesh.m_edges[e][0]);
  ElTopo::Vec3d x1 = st.get_position(st.m_mesh.m_edges[e][1]);
  
  int constraint0 = onBBWall(Vec3d(x0[0], x0[1], x0[2]));
  int constraint1 = onBBWall(Vec3d(x1[0], x1[1], x1[2]));
  
  if (constraint0 & constraint1)  // edge is completely inside a wall
    return true;
  else
    return false;
}


} //namespace BASim
