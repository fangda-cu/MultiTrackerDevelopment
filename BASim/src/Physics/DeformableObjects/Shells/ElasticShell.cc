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
#include "Timer.h"

#include <algorithm>
#include <numeric>

namespace BASim {
  
ElasticShell::ElasticShell(DeformableObject* object, const FaceProperty<char>& shellFaces, Scalar timestep, SteppingCallback * stepping_callback) : 
  PhysicalModel(*object), m_obj(object), 
    m_active_faces(shellFaces), 
//    m_undef_xi(object),
//    m_damping_undef_xi(object),
//    m_vertex_masses(object),
//    m_edge_masses(object),
//    m_thicknesses(object),
    m_face_regions(object),
    m_vertex_constraint_labels(object),
//    m_volumes(object),
//    m_xi(object), 
//    m_xi_vel(object),
//    m_density(1),
    m_collision_epsilon(1e-5),
    m_vert_point_springs(NULL),
    m_repulsion_springs(NULL),
    m_sphere_collisions(false),
    m_object_collisions(false),
    m_ground_collisions(false),
    m_do_eltopo_collisions(false),
//    m_do_thickness_updates(true),
//    m_momentum_conserving_remesh(false)
    m_stepping_callback(stepping_callback)
{
  m_vert_point_springs = new ShellVertexPointSpringForce(*this, "VertPointSprings", timestep);
  m_repulsion_springs = new ShellStickyRepulsionForce(*this, "RepulsionSprings", timestep);

  addForce(m_vert_point_springs);
  addForce(m_repulsion_springs);
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

//void ElasticShell::setThickness( Scalar thickness )
//{
//  for(FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit) {
//    m_thicknesses[*fit] = thickness;
//    Scalar area = getArea(*fit, false);
//    m_volumes[*fit] = thickness * area;
//  }
//
//}
//
//void ElasticShell::setDensity(Scalar density) {
//  m_density = density;
//}

//void ElasticShell::setVertexUndeformed( const VertexProperty<Vec3d>& undef )
//{
//  m_undeformed_positions = undef;
//}

//void ElasticShell::setEdgeUndeformed( const EdgeProperty<Scalar>& undef )
//{
//  m_undef_xi = undef;
//}

//void ElasticShell::setVertexPositions( const VertexProperty<Vec3d>& positions )
//{
//  m_positions = positions;
//}

//void ElasticShell::setEdgeXis( const EdgeProperty<Scalar>& positions )
//{
//  m_xi = positions;
//}

//void ElasticShell::setVertexVelocities( const VertexProperty<Vec3d>& velocities) 
//{
//  m_velocities = velocities;
//}

//void ElasticShell::setEdgeVelocities(const EdgeProperty<Scalar>& velocities)
//{
//  m_xi_vel = velocities;
//}

//Scalar ElasticShell::getThickness(const VertexHandle& vh) const {
//  Scalar totalA = 0.0;
//  Scalar w;
//  Scalar total = 0.0;
//  for (VertexFaceIterator vfit = m_obj->vf_iter(vh); vfit; ++vfit){
//      w = getArea(*vfit);
//      totalA += w;
//      total += w*m_thicknesses[*vfit];
//  }
//
//  assert ( totalA > 0.);
//  assert ( total > 0.);
//  return total / totalA;
//}
//Scalar ElasticShell::getMaxThickness() const {
//  Scalar maxVal = -1000000;
//  for( FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit ){
//      if (m_thicknesses[*fit] > maxVal ) maxVal = m_thicknesses[*fit];
//  }
//  return maxVal;
//}
//Scalar ElasticShell::getMinThickness() const {
//  Scalar minVal = 1000000;
//  for( FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit ){
//      if (m_thicknesses[*fit] < minVal ) minVal = m_thicknesses[*fit];
//  }
//  return minVal;
//}

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


//void ElasticShell::getThickness(VertexProperty<Scalar> & vThickness) const{
//    for ( VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit){
//        vThickness[*vit] = getThickness(*vit);
//    }
//}

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
//  else 
//  {
//    Vec3d pos0 = m_obj->getVertexUndeformedPosition(v0_hnd);
//    Vec3d pos1 = m_obj->getVertexUndeformedPosition(v1_hnd);
//    Vec3d pos2 = m_obj->getVertexUndeformedPosition(v2_hnd);
//    
//    Vec3d v0 = pos1 - pos0;
//    Vec3d v1 = pos2 - pos0;
//    Vec3d triVec = v0.cross(v1);
//    return 0.5*triVec.norm();
//  }
}

//void ElasticShell::recomputeVertexMass(const VertexHandle& v) {
//  
//  //Compute vertex mass in a lumped fashion from incident faces.
//  Scalar volume = 0;
//  for (VertexFaceIterator vfit = m_obj->vf_iter(v); vfit; ++vfit) {
//    const FaceHandle& f = *vfit;
//    volume += getArea(f) * m_thicknesses[f] / 3.0;
//  }
//  m_vertex_masses[v] = volume * m_density;
//
//}

//void ElasticShell::computeMasses()
//{
//  //Compute vertex masses in a lumped mass way.
//
//  m_vertex_masses.assign(0);
//  m_edge_masses.assign(0);
//
//  Scalar area = 0;
//
//  //Iterate over all triangles active in this shell and accumulate vertex masses
//  for(FaceIterator f_iter = m_obj->faces_begin(); f_iter != m_obj->faces_end(); ++f_iter) {
//    FaceHandle& f_hnd = *f_iter;
//    if(m_active_faces[f_hnd]) {
//
//      //get the three vertices
//      FaceVertexIterator fvit = m_obj->fv_iter(f_hnd);
//      VertexHandle v0_hnd = *fvit; ++fvit; assert(fvit);
//      VertexHandle v1_hnd = *fvit; ++fvit; assert(fvit);
//      VertexHandle v2_hnd = *fvit; ++fvit; assert(!fvit);
//
//      //compute triangle areas
//      Vec3d v0 = getVertexPosition(v1_hnd) - getVertexPosition(v0_hnd);
//      Vec3d v1 = getVertexPosition(v2_hnd) - getVertexPosition(v0_hnd);
//      Vec3d triVec = v0.cross(v1);
//      Scalar area = 0.5*sqrt(triVec.dot(triVec)) / 3.0;
//      Scalar contribution = m_thicknesses[f_hnd] * m_density * area;
//      
//      //accumulate mass to the vertices
//      m_vertex_masses[v0_hnd] += contribution;
//      m_vertex_masses[v1_hnd] += contribution;
//      m_vertex_masses[v2_hnd] += contribution;
//
//      //set edge masses to zero, since we want to solve them quasistatically (they're derivative DOFs)
//      FaceEdgeIterator feit = m_obj->fe_iter(f_hnd);
//      EdgeHandle e0_hnd = *feit; ++feit; assert(feit);
//      EdgeHandle e1_hnd = *feit; ++feit; assert(feit);
//      EdgeHandle e2_hnd = *feit; ++feit; //assert(feit);
//
//      m_edge_masses[e0_hnd] = 0;
//      m_edge_masses[e1_hnd] = 0;
//      m_edge_masses[e2_hnd] = 0;
//
//      //store the current volumes
//      m_volumes[f_hnd] = 3*area*m_thicknesses[f_hnd];
//    }
//  }
// 
//  m_obj->updateVertexMasses();
//}

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
//  //they're all edge Dofs for a shell
//  assert(hnd.getType() == DofHandle::EDGE_DOF);
//
//  //return reference to the appropriate position in the vector
//  const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
  
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
//  //they're all vertex Dofs for a shell
//  assert(hnd.getType() == DofHandle::EDGE_DOF);
//
//  const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
//  m_xi[eh] = dof;
}

const Scalar& ElasticShell::getVel( const DofHandle& hnd ) const
{
  assert("No dof");
//  assert(hnd.getType() == DofHandle::EDGE_DOF);
//
//  //return reference to the appropriate position in the vector
//  const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
  
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
//  assert(hnd.getType() == DofHandle::EDGE_DOF);
//
//  const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
//  m_xi_vel[eh] = vel;
}

//const Scalar& ElasticShell::getMass( const DofHandle& hnd ) const
//{
//  assert("No dof");
////  assert(hnd.getType() == DofHandle::EDGE_DOF);
////
////  const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
////  return m_edge_masses[eh];
//}

void ElasticShell::getScriptedDofs( IntArray& dofIndices, std::vector<Scalar>& dofValues, Scalar time ) const
{
    // position dof scripting is moved to PositionDofsModel.
    
//    for(unsigned int i = 0; i < constrainedEdges.size(); ++i) {
//        int dofID = getEdgeDofBase(constrainedEdges[i]);
//        dofIndices.push_back(dofID);
//        dofValues.push_back(constrainedXiValues[i]);
//    }
}

void ElasticShell::startStep(Scalar time, Scalar timestep)
{
//
//  //update the damping "reference configuration" for computing viscous forces.
//  m_damping_undef_xi = m_xi;

  //tell the forces to update anything they need to update
  const std::vector<ElasticShellForce*>& forces = getForces();
  for(unsigned int i = 0; i < forces.size(); ++i) {
    forces[i]->update();
  }
  
  m_stepping_callback->afterStartStep();
  
}

void ElasticShell::resolveCollisions(Scalar timestep) {
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
  ElTopo::DynamicSurface dynamic_surface( vert_old, tri_data, std::vector<ElTopo::Vec2i>(tri_data.size(), ElTopo::Vec2i(0, 0)), masses, m_collision_epsilon, friction_coeff, true, false );

  dynamic_surface.set_all_newpositions( vert_new );
    
    for (size_t i = 0; i < dynamic_surface.m_mesh.nv(); i++)
    {
        assert(dynamic_surface.get_position(i)[0] == dynamic_surface.get_position(i)[0]);
        assert(dynamic_surface.get_position(i)[1] == dynamic_surface.get_position(i)[1]);
        assert(dynamic_surface.get_position(i)[2] == dynamic_surface.get_position(i)[2]);
        assert(dynamic_surface.get_newposition(i)[0] == dynamic_surface.get_newposition(i)[0]);
        assert(dynamic_surface.get_newposition(i)[1] == dynamic_surface.get_newposition(i)[1]);
        assert(dynamic_surface.get_newposition(i)[2] == dynamic_surface.get_newposition(i)[2]);
    }
    
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


void ElasticShell::addSelfCollisionForces() {
  
  //update the broad phase structure with the current mesh data
  Scalar collision_distance = m_collision_proximity;

  m_broad_phase.update_broad_phase_static(*m_obj, getVertexPositions(), collision_distance);
  
  //determine proximity of vertex triangle pairs and
  //add damped springs between them to handle new collisions
  
  //consider all vertices
  
  for(VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit) {
    VertexHandle vh = *vit;
    Vec3d vert_pos = getVertexPosition(vh);
    ElTopoCode::Vec3d vert_vel = ElTopoCode::toElTopo(getVertexVelocity(vh));
    ElTopoCode::Vec3d vertex_position = ElTopoCode::toElTopo(vert_pos);

    //construct bound box for the vertex, and find all triangles near it
    std::vector<unsigned int> overlapping_triangles;
    Vec3d low = vert_pos - collision_distance*Vec3d(1,1,1), high = vert_pos + collision_distance*Vec3d(1,1,1);
    
    m_broad_phase.get_potential_triangle_collisions(ElTopoCode::toElTopo(low), ElTopoCode::toElTopo(high), overlapping_triangles);
    for(unsigned int i = 0; i < overlapping_triangles.size(); ++i) {
      int tri_idx = overlapping_triangles[i];
      FaceHandle f(tri_idx);
      
      if(m_repulsion_springs->springExists(f, vh)) {
        continue;
      }

      ElTopoCode::Vec3d face_verts[3];
      ElTopoCode::Vec3d face_vels[3];
      int fv = 0;
      bool goodSpring = true;
      for(FaceVertexIterator fvit = m_obj->fv_iter(f); fvit; ++fvit) {
        face_verts[fv] = ElTopoCode::toElTopo(getVertexPosition(*fvit));
        face_vels[fv] = ElTopoCode::toElTopo(getVertexVelocity(*fvit));
        if(*fvit == vh)
          goodSpring = false;
        ++fv;
      }
      if(!goodSpring) {
        continue;
      }
      
      //check if the geometry is actually close enough to warrant a spring
      Vec3d barycoords;
      Scalar distance;
      ElTopoCode::Vec3d normal;
      check_point_triangle_proximity(vertex_position, face_verts[0], face_verts[1], face_verts[2], distance, barycoords[0], barycoords[1], barycoords[2], normal );
      
      //if such a spring doesn't already exist, add it
      ElTopoCode::Vec3d faceNormal = cross(face_verts[1] - face_verts[0], face_verts[2] - face_verts[0]);
      normalize(faceNormal);
      ElTopoCode::Vec3d offset = vertex_position - (barycoords[0]*face_verts[0] + barycoords[1]*face_verts[1] + barycoords[2]*face_verts[2]);
      Scalar normalDist = dot(offset,faceNormal);
      
      

      if(distance < collision_distance) {
        
        bool isNbrVert = false;
        if(barycoords[0]< 1e-5 || barycoords[1] < 1e-5 || barycoords[2] < 1e-5) {
          //the proposed connection point lies along an edge or at a vertex. This is troublesome, if any of the relevant edges
          //neighbours the "colliding" vertex (it means that the spring we would add is in the plane of the mesh!)
          int vertNum = 0;
          for(FaceVertexIterator fvit = m_obj->fv_iter(f); fvit && !isNbrVert; ++fvit, ++vertNum) {
            //grab the current vertex
            VertexHandle testVert = *fvit;

            //if it's not one of the zero values, check its neighbours for incidence
            if(barycoords[vertNum] > 1e-5) {
              for(VertexEdgeIterator veit = m_obj->ve_iter(testVert); veit; ++veit) {
                EdgeHandle curEdge = *veit;
                VertexHandle nbrVert;
                VertexHandle fromVert = m_obj->fromVertex(curEdge);
                if(fromVert == testVert)
                  nbrVert = m_obj->toVertex(curEdge);
                else
                  nbrVert = fromVert;

                if(nbrVert == vh) {
                  isNbrVert = true;
                  break;
                }
              }
            }
             //look at next vertex
          }
        }

        if(isNbrVert) 
          continue;

        //look at velocities of closest points to see if the geometry is approaching
        ElTopoCode::Vec3d close_point_vel = barycoords[0]*face_vels[0] + barycoords[1]*face_vels[1] + barycoords[2]*face_vels[2];
        ElTopoCode::Vec3d rel_vel = vert_vel - close_point_vel;
        //if(dot(rel_vel, offset) < 0) //add spring only if they are approaching
          m_repulsion_springs->addSpring(f, vh, barycoords, m_collision_spring_stiffness, m_collision_spring_damping, fabs(normalDist));
      }
    }
  }
  
  
}

void ElasticShell::getSpringList(std::vector<Vec3d>& start, std::vector<Vec3d>& end)  const {
  std::vector<VertexHandle> verts;
  std::vector<Vec3d> bary;
  std::vector<FaceHandle> faces;
  
  m_repulsion_springs->getSpringLists(verts, faces, bary);
  for(unsigned int i = 0; i < verts.size(); ++i) {
    start.push_back(getVertexPosition(verts[i]));
    
    FaceVertexIterator fvit = m_obj->fv_iter(faces[i]);
    Vec3d v0 = getVertexPosition(*fvit);++fvit;
    Vec3d v1 = getVertexPosition(*fvit);++fvit;
    Vec3d v2 = getVertexPosition(*fvit);
    Vec3d result = bary[i][0] * v0 + bary[i][1] * v1 + bary[i][2] * v2;
    end.push_back(result);
  }

}

void ElasticShell::setCollisionParams(Scalar proximity, Scalar epsilon, Scalar stiffness, Scalar damping) {
  m_collision_epsilon = epsilon;
  m_collision_spring_stiffness = stiffness;
  m_collision_spring_damping = damping;
  m_collision_proximity = proximity;
}

void ElasticShell::setGroundPlane(bool enabled, Scalar height, Scalar velocity) {
  m_ground_collisions = enabled;
  m_ground_height = height;
  m_ground_velocity = velocity;
}

void ElasticShell::setCollisionSphere(bool enabled, Scalar radius, Vec3d position, Vec3d velocity) {
  m_sphere_collisions = enabled;
  m_sphere_radius = radius;
  m_sphere_position = position;
  m_sphere_velocity = velocity;
}

void ElasticShell::setCollisionObject(bool enabled, const Vec3d& position, const Vec3d& velocity, const ElTopoCode::Array3f& phi_grid, 
                                      const Vec3d& origin, Scalar dx) {
  m_object_collisions = enabled;
  m_object_position = position;
  m_object_velocity = velocity;
  m_object_SDF = phi_grid;
  m_object_origin = origin;
  m_object_dx = dx;

}

void ElasticShell::setSelfCollision(bool enabled) {
  m_self_collisions = enabled;
}


void ElasticShell::endStep(Scalar time, Scalar timestep) {
  
  CSim::TimerMan::timer("endStep").start();
    
//    for (VertexIterator v = m_obj->vertices_begin(); v != m_obj->vertices_end(); ++v)
//    {
//        assert(getVertexPosition(*v) == getVertexPosition(*v));
//        std::cout << "vertex " << (*v).idx() << ": " << getVertexPosition(*v) << std::endl;
//    }

  CSim::TimerMan::timer("endStep/beforeEndStep").start();
  if (m_stepping_callback)
    m_stepping_callback->beforeEndStep();
  CSim::TimerMan::timer("endStep/beforeEndStep").stop();
  
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
  
  CSim::TimerMan::timer("endStep/resolveCollisions").start();
  if(m_do_eltopo_collisions)
    resolveCollisions(timestep);
  CSim::TimerMan::timer("endStep/resolveCollisions").stop();
  
  //Ground plane penalty force.
  if(m_ground_collisions) {

    //Hard constraints
    for(VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit) {
      Vec3d curPos = getVertexPosition(*(vit));
      if(curPos[1] < m_ground_height) {
        if(!getDefoObj().isConstrained(*vit)) {
          //constrainVertex(*vit, curPos);
          
          //Sinking
          //constrainVertex(*vit, new FixedVelocityConstraint(curPos, Vec3d(0, m_ground_velocity, 0), time));

          curPos[1] = m_ground_height ; //project constraint back to the ground plane. This might be a tad unsafe.

          //Conveying
          getDefoObj().constrainVertex(*vit, new FixedVelocityConstraint(curPos, Vec3d(0, 0, m_ground_velocity), time));
        }
      }
    }
  }

  if(m_sphere_collisions) {
    //Hard constraints 
    for(VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit) {
      Vec3d curPos = getVertexPosition(*(vit));
      Vec3d offset = curPos - (m_sphere_position+m_sphere_velocity*time);
      if(offset.norm() < m_sphere_radius) {
        if(!getDefoObj().isConstrained(*vit)) {

          //Sinking
          //constrainVertex(*vit, new FixedVelocityConstraint(curPos, Vec3d(0, m_ground_velocity, 0), time));

          //Conveying
          getDefoObj().constrainVertex(*vit, new FixedVelocityConstraint(curPos, m_sphere_velocity, time));
        }
      }
    }
  }

  if(m_object_collisions) {
    for(VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit) {
      Vec3d curPos = getVertexPosition(*(vit));
      Vec3d offset = (curPos - (m_object_origin + m_object_position)) / m_object_dx;
      offset[0] = clamp(offset[0], 0.0, m_object_SDF.ni-1.6);
      offset[1] = clamp(offset[1], 0.0, m_object_SDF.nj-1.6);
      offset[2] = clamp(offset[2], 0.0, m_object_SDF.nk-1.6);
      Scalar dist_value = m_object_SDF((int)offset[0], (int)offset[1], (int)offset[2]);      
      if(dist_value < 0) {
        if(!getDefoObj().isConstrained(*vit)) {

          //Fixed position
          getDefoObj().constrainVertex(*vit, curPos);
          //constrainVertex(*vit, new FixedVelocityConstraint(curPos, m_sphere_velocity, time));
        }
      }
    }
  }

  //apply penalty springs for self-collision
  if(m_self_collisions) {
    addSelfCollisionForces();
  }

//  //Adjust thicknesses based on area changes
//  if(m_do_thickness_updates)
//    updateThickness();

  if(m_inflow) {
    extendMesh(time);
    do_relabel = true;
  }

  if(m_delete_region) {
    deleteRegion();
    do_relabel = true;
  }

    if (m_stepping_callback)
        m_stepping_callback->endStep(1);
    
  //Remeshing
  if(m_do_remeshing) {
    std::cout << "Remeshing\n";
    CSim::TimerMan::timer("endStep/remesh").start();
    remesh(timestep);
    CSim::TimerMan::timer("endStep/remesh").stop();

      if (m_stepping_callback)
          m_stepping_callback->endStep(2);
      
    //Relabel DOFs if necessary
    do_relabel = true;
  }
  
  std::cout << "nv = " << getDefoObj().nv() << " ne = " << getDefoObj().ne() << " nf = " << getDefoObj().nf() << " nt = " << getDefoObj().nt() << std::endl;

  static int only_twice = 0;

  //Tearing processing
  if(m_tearing){
    std::cout << "Processing tearing. \n";
    //fracture();
    do_relabel = true;
    only_twice++;
  }

    if (m_stepping_callback)
        m_stepping_callback->endStep(3);
    
  if(do_relabel) {
    m_obj->computeDofIndexing();
  }

  //Update masses based on new areas/thicknesses
//  computeMasses();

  std::cout << "Completed endStep\n";

  if (m_stepping_callback)
    m_stepping_callback->endStep(4);
    
    CSim::TimerMan::timer("endStep").stop();

}


void ElasticShell::fracture() {

  //Set up a SurfTrack, run remeshing, render the new mesh
  ElTopo::SurfTrackInitializationParameters construction_parameters;
  construction_parameters.m_proximity_epsilon = m_collision_epsilon;
  construction_parameters.m_allow_vertex_movement_during_collapse = true;
  construction_parameters.m_perform_smoothing = false;
  construction_parameters.m_min_edge_length = m_remesh_edge_min_len;
  construction_parameters.m_max_edge_length = m_remesh_edge_max_len;
  construction_parameters.m_max_volume_change = numeric_limits<double>::max();   
  construction_parameters.m_min_triangle_angle = 5;
  construction_parameters.m_max_triangle_angle = 175;
  construction_parameters.m_verbose = false;

  construction_parameters.m_use_curvature_when_collapsing = false;
  construction_parameters.m_use_curvature_when_splitting = false;
  
  construction_parameters.m_subdivision_scheme = new ElTopo::ButterflyScheme();//ElTopo::MidpointScheme(); 

  construction_parameters.m_allow_non_manifold = true;
  construction_parameters.m_remesh_boundaries = true;
  construction_parameters.m_collision_safety = true;

  std::vector<ElTopo::Vec3d> vert_data;
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
    ElTopo::Vec3d mass(1, 1, 1);
      
    vert_data.push_back(ElTopo::Vec3d(vert[0], vert[1], vert[2]));
    for (int i = 0; i < 3; i++)
      if (getDefoObj().isConstrainedInDirection(vh, i))
        mass[i] = numeric_limits<Scalar>::infinity();
    masses.push_back(mass);
    vert_numbers[vh] = id;
    reverse_vertmap.push_back(vh);
      
    Vec3d pos = getVertexPosition(vh);

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

  //constrain all the vertices in faces that are used by collision springs to prevent remeshing there
  std::vector<VertexHandle> verts;
  std::vector<FaceHandle> faces;
  std::vector<Vec3d> coords;
  m_repulsion_springs->getSpringLists(verts, faces, coords);
  for(unsigned int i = 0; i < faces.size(); ++i) {
    FaceVertexIterator fvit = getDefoObj().fv_iter(faces[i]);
    for(;fvit; ++fvit) {
      VertexHandle vh = *fvit;
      masses[vert_numbers[vh]][0] = numeric_limits<Scalar>::infinity();
      masses[vert_numbers[vh]][1] = numeric_limits<Scalar>::infinity();
      masses[vert_numbers[vh]][2] = numeric_limits<Scalar>::infinity();
    }
  }

    ElTopo::SurfTrack surface_tracker( vert_data, tri_data, std::vector<ElTopo::Vec2i>(tri_data.size(), ElTopo::Vec2i(0, 0)), masses, construction_parameters ); 

  std::vector< std::pair<size_t,size_t> > edges_to_cut;
  for(EdgeIterator it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
    EdgeHandle eh = *it;
    if(shouldFracture(eh)) {
      VertexHandle vh0 = mesh.fromVertex(eh);
      VertexHandle vh1 = mesh.toVertex(eh);
      edges_to_cut.push_back(make_pair(vert_numbers[vh0],vert_numbers[vh1]));
    }
  }

  std::cout << "Doing cutting with El Topo.\n";
  surface_tracker.cut_mesh(edges_to_cut);

  
  std::cout << "Performing " << surface_tracker.m_mesh_change_history.size() << " Cutting Operations:\n";
  for(unsigned int j = 0; j < surface_tracker.m_mesh_change_history.size(); ++j) {
    ElTopo::MeshUpdateEvent event = surface_tracker.m_mesh_change_history[j];

    VertexHandle v0 = reverse_vertmap[event.m_v0];
    VertexHandle v1 = reverse_vertmap[event.m_v1];
    EdgeHandle eh = findEdge(mesh, v0, v1);

    assert(event.m_type == ElTopo::MeshUpdateEvent::EDGE_CUT);
      
    //Identify the edge based on its endpoint vertices instead

    //based on the type of cut, perform the appropriate operation
    bool aBound = m_obj->isBoundary(v0);
    bool bBound = m_obj->isBoundary(v1);
    
    std::vector<VertexHandle> srcVerts;
    if(aBound) srcVerts.push_back(v0);
    if(bBound) srcVerts.push_back(v1);
    if(!aBound & !bBound) {
      VertexHandle vShared = reverse_vertmap[event.m_v2];
      srcVerts.push_back(vShared);
    }
    std::vector<FaceHandle> newFaces;
    std::vector<FaceHandle> deletedFaces;
    std::vector<EdgeHandle> deleteEdges;
    std::vector<VertexHandle> newVerts;

    //Do exactly the same cut as El Topo
    for(unsigned int i = 0; i < event.m_created_verts.size(); ++i) {
      VertexHandle nv = mesh.addVertex();
      newVerts.push_back(nv);
        
      //Update the mapping
      vert_numbers[nv] = event.m_created_verts[i];
      if(reverse_vertmap.size() <= event.m_created_verts[i]) 
        reverse_vertmap.resize(event.m_created_verts[i]+1, VertexHandle(-1));
      reverse_vertmap[event.m_created_verts[i]] = nv;

      //Set the various data for the new vertex
      Vec3d new_pos(event.m_created_vert_data[i][0], 
                    event.m_created_vert_data[i][1], 
                    event.m_created_vert_data[i][2]);
      setVertexPosition(nv, new_pos);

      VertexHandle source = srcVerts[i];
      setVertexVelocity(nv, getVertexVelocity(source));
//      setUndeformedVertexPosition(nv, getVertexUndeformed(source));
    }
      
    // Add/update faces
    for(unsigned int i = 0; i < event.m_created_tri_data.size(); ++i) {
      ElTopo::Vec3st new_face = event.m_created_tri_data[i];

      //determine handles in our indexing
      VertexHandle v0 = reverse_vertmap[new_face[0]];
      VertexHandle v1 = reverse_vertmap[new_face[1]];
      VertexHandle v2 = reverse_vertmap[new_face[2]];
      assert(v0.isValid() && v1.isValid() && v2.isValid());

      FaceHandle newFaceHandle = mesh.addFace(v0, v1, v2); //build the face. 
      setFaceActive(newFaceHandle);

      face_numbers[newFaceHandle] = event.m_created_tris[i];
      if(reverse_trimap.size() <= event.m_created_tris[i]) 
        reverse_trimap.resize(event.m_created_tris[i]+1);
      reverse_trimap[event.m_created_tris[i]] = newFaceHandle;
        
      //copy data from the old version of the face
      FaceHandle oldFaceHandle = reverse_trimap[event.m_deleted_tris[i]];
//      m_thicknesses[newFaceHandle] = m_thicknesses[oldFaceHandle];
//      m_volumes[newFaceHandle] = m_volumes[oldFaceHandle];
        
    }
    
    //delete the dead faces, and clean up the mapping
    for(size_t i = 0; i < event.m_deleted_tris.size(); ++i) {
      size_t tri_to_delete = event.m_deleted_tris[i];
      FaceHandle localFace = reverse_trimap[tri_to_delete];
      mesh.deleteFace(localFace, true);

      reverse_trimap[event.m_deleted_tris[i]] = FaceHandle(-1);
    }
  }
  

}


bool ElasticShell::isInflow(const EdgeHandle & eh) const{
    bool isInflow = false;
    for(unsigned int i = 0; i < m_inflow_boundaries.size(); ++i) {
      if(std::find(m_inflow_boundaries[i].begin(), m_inflow_boundaries[i].end(),eh) != m_inflow_boundaries[i].end()) {
        isInflow = true;
        break;
      }
    }
    return isInflow;
}


bool ElasticShell::shouldFracture (const EdgeHandle & eh) const{
    
    //Ignore inflow edges
    if(isInflow(eh)) return false;

    //Ignore constrain edges
    if ( getDefoObj().isConstrained(m_obj->fromVertex(eh)) || getDefoObj().isConstrained(m_obj->toVertex(eh))) return false;

    //Ignore boundary edges
    if ( m_obj->isBoundary(eh) ) return false;

    //Get the average thickness weighted by areas for each edge
    Scalar thickness = 0.0;
    Scalar totalA = 0.0;
    
    bool both_thin = true; //just check if all faces are thin

    for ( EdgeFaceIterator efit = m_obj->ef_iter(eh); efit; ++efit){
        int face_nbr_count = 0;
        Scalar w = getArea(*efit);
//        thickness += w * getThickness(*efit);
//        totalA += w;
//        if(getThickness(*efit) > m_tear_thres)
        both_thin = false;
       
    }
    thickness /= totalA;
    Scalar p = (Scalar) rand() / (Scalar) RAND_MAX;
    
    return both_thin && (p <  m_tear_rand);
}


void ElasticShell::remesh(Scalar timestep, bool initial)
{
  CSim::TimerMan::timer("endStep/remesh/prep").start();
  
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
  construction_parameters.m_proximity_epsilon = m_collision_epsilon;
  construction_parameters.m_merge_proximity_epsilon = 0.02 * m_remesh_edge_min_len;
  construction_parameters.m_allow_vertex_movement_during_collapse = true;
  construction_parameters.m_perform_smoothing = false;
  construction_parameters.m_min_edge_length = m_remesh_edge_min_len;
  construction_parameters.m_max_edge_length = m_remesh_edge_max_len;
  construction_parameters.m_max_volume_change = numeric_limits<double>::max();   
  construction_parameters.m_min_triangle_angle = initial ? 0 : 3;
  construction_parameters.m_max_triangle_angle = initial ? 180 : 177;
  construction_parameters.m_large_triangle_angle_to_split = 160;
  construction_parameters.m_min_triangle_area = 0.0002*m_remesh_edge_min_len*m_remesh_edge_min_len;
  construction_parameters.m_verbose = false;
  construction_parameters.m_allow_non_manifold = true;
  construction_parameters.m_allow_topology_changes = true;
  construction_parameters.m_collision_safety = true;
  construction_parameters.m_remesh_boundaries = true;
  construction_parameters.m_t1_transition_enabled = m_remesh_t1transition;
  construction_parameters.m_velocity_field_callback = NULL;
    construction_parameters.m_pull_apart_distance = m_remesh_edge_min_len / 2;// (initial ? 0.1 : 0.02) * m_remesh_edge_min_len;
  
  ElTopo::SubdivisionScheme * mb = new ElTopo::ModifiedButterflyScheme();
  ElTopo::SubdivisionScheme * mp = new ElTopo::MidpointScheme();
  construction_parameters.m_subdivision_scheme = (m_remesh_smooth_subdivision ? mb : mp);
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

  CSim::TimerMan::timer("endStep/remesh/prep").stop();

  CSim::TimerMan::timer("endStep/remesh/copy1").start();
  
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

  //constrain all the vertices in faces that are used by collision springs to prevent remeshing there,
  //as well as the associated opposite vertex.
  std::vector<VertexHandle> verts;
  std::vector<FaceHandle> faces;
  std::vector<Vec3d> coords;
  m_repulsion_springs->getSpringLists(verts, faces, coords);
  for(unsigned int i = 0; i < faces.size(); ++i) {
    FaceVertexIterator fvit = getDefoObj().fv_iter(faces[i]);
    //vertices of the face...
    for(;fvit; ++fvit) {
      VertexHandle vh = *fvit;
      masses[vert_numbers[vh]][0] = numeric_limits<Scalar>::infinity();
      masses[vert_numbers[vh]][1] = numeric_limits<Scalar>::infinity();
      masses[vert_numbers[vh]][2] = numeric_limits<Scalar>::infinity();
    }
    //and the other vertex
    masses[vert_numbers[verts[i]]][0] = numeric_limits<Scalar>::infinity();
    masses[vert_numbers[verts[i]]][1] = numeric_limits<Scalar>::infinity();
    masses[vert_numbers[verts[i]]][2] = numeric_limits<Scalar>::infinity();
  }

  std::cout << "Calling surface improvement\n";
  
  CSim::TimerMan::timer("endStep/remesh/copy1").stop();
  
  CSim::TimerMan::timer("endStep/remesh/st_constuction").start();
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
  CSim::TimerMan::timer("endStep/remesh/st_constuction").stop();

  for(int i = 0; i < m_remeshing_iters; ++i) {
    CSim::TimerMan::timer("endStep/remesh/topology_change").start();
    surface_tracker.topology_changes();
    CSim::TimerMan::timer("endStep/remesh/topology_change").stop();
    CSim::TimerMan::timer("endStep/remesh/improve_mesh").start();
    surface_tracker.improve_mesh();
    CSim::TimerMan::timer("endStep/remesh/improve_mesh").stop();
  }
  
    if (m_mesheventcallback)
        m_mesheventcallback->flip(surface_tracker, 4);
    
  // copy ElTopo mesh back, instead of repeating the operation history incrementally.
  // this is possible because ElasticShell doesn't keep any other information that ElTopo doesn't have
  
  CSim::TimerMan::timer("endStep/remesh/copy2").start();
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
  CSim::TimerMan::timer("endStep/remesh/copy2").stop();
  
    if (m_mesheventcallback)
        m_mesheventcallback->flip(surface_tracker, 5);
    
    if (m_stepping_callback)
        m_stepping_callback->endStep(11);

  
  CSim::TimerMan::timer("endStep/remesh/clean_up").start();
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
    
    if (m_stepping_callback)
        m_stepping_callback->endStep(12);
    
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
  
    if (m_stepping_callback)
        m_stepping_callback->endStep(13);
  CSim::TimerMan::timer("endStep/remesh/clean_up").stop();

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
    
    if (!(((constraint0 & (1 << 0)) || (constraint0 & (1 << 3))) == (bool)label0[0]) ||
        !(((constraint0 & (1 << 1)) || (constraint0 & (1 << 4))) == (bool)label0[1]) ||
        !(((constraint0 & (1 << 2)) || (constraint0 & (1 << 5))) == (bool)label0[2]) ||
        !(((constraint1 & (1 << 0)) || (constraint1 & (1 << 3))) == (bool)label1[0]) ||
        !(((constraint1 & (1 << 1)) || (constraint1 & (1 << 4))) == (bool)label1[1]) ||
        !(((constraint1 & (1 << 2)) || (constraint1 & (1 << 5))) == (bool)label1[2]))
    {
        std::cout << "labels: " << (bool)label0[0] << (bool)label0[1] << (bool)label0[2] << " " << (bool)label1[0] << (bool)label1[1] << (bool)label1[2] << std::endl;
        std::cout << "consts: " << constraint0 << " " << constraint1 << std::endl;
        std::cout << "x0 = " << std::setprecision(15) << x0[0] << " " << x0[1] << " " << x0[2] << " x1 = " << x1[0] << " " << x1[1] << " " << x1[2] << std::endl;
    }
    
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

void ElasticShell::performSplit(const EdgeHandle& eh, const Vec3d& midpoint, VertexHandle& new_vert) {

  VertexHandle v0 = m_obj->fromVertex(eh);
  VertexHandle v1 = m_obj->toVertex(eh);
  assert(v0!=v1);
  Vec3d p0 = getVertexPosition(v0);
  Vec3d p1 = getVertexPosition(v1);

  //store momenta
//  Vec3d mom0 = getVertexVelocity(v0) * getMass(v0);
//  Vec3d mom1 = getVertexVelocity(v1) * getMass(v1);
  
  std::vector<FaceHandle> oldFaces;
  std::vector<Scalar> oldThicknesses;
  std::vector<Vec2i> oldRegions;
  std::vector<Vec3d> momenta;
  std::vector<VertexHandle> nbrVerts;
  for(EdgeFaceIterator efit = m_obj->ef_iter(eh); efit; ++efit) {
    FaceHandle f = *efit;
    oldFaces.push_back(f);
//    oldThicknesses.push_back(getThickness(f));
    oldRegions.push_back(getFaceLabel(f));

    //store the momenta of the associated vertex
    VertexHandle thirdVert;
    getFaceThirdVertex(*m_obj, f, eh, thirdVert);
    nbrVerts.push_back(thirdVert);
//    momenta.push_back(getVertexVelocity(thirdVert) * getMass(thirdVert));
  }

//  VertexHandle v2, v3;
//  getEdgeOppositeVertices(*m_obj, eh, v2, v3);/////////////

  ElTopoCode::Vec3d midpoint_ET = ElTopoCode::toElTopo(midpoint);

  //perform the actual split
  std::vector<FaceHandle> newFaces;
  VertexHandle v_new = splitEdge(*m_obj, eh, newFaces);

  Vec3d simple_midpoint = 0.5*(p1 + p0);
  
  //determine split fraction, and using it to lerp the vertex data
  Vec3d dx(p1-p0);
  double m2 = dx.squaredNorm();
  Scalar s = clamp((p1-midpoint).dot(dx)/m2, 0., 1.);
  
  Vec3d velocity = s*getVertexVelocity(v0) + (1-s)*getVertexVelocity(v1);
//  Vec3d undef = s*getVertexUndeformed(v0) + (1-s)*getVertexUndeformed(v1);
  setVertexVelocity(v_new, velocity);
//  setUndeformedVertexPosition(v_new, undef);

  //initially set to the exact midpoint for simplicity
  setVertexPosition(v_new, simple_midpoint);

  //set consistent volumes and thickness for new faces
  assert(oldFaces.size() == newFaces.size()/2);
  
  //Old way, ignores vertex movement, assumes simple (true) midpoint
  for(unsigned int i = 0; i < oldFaces.size(); ++i) {
    //copy thicknesses exactly
//    m_thicknesses[newFaces[i*2]] = oldThicknesses[i];
//    m_thicknesses[newFaces[i*2+1]] = oldThicknesses[i];
    
    //compute corresponding volumes
//    m_volumes[newFaces[i*2]] = oldThicknesses[i] * getArea(newFaces[i*2], true);
//    m_volumes[newFaces[i*2+1]] = oldThicknesses[i] * getArea(newFaces[i*2+1], true);
    
    //copy region labels
    m_face_regions[newFaces[i*2]] = oldRegions[i];
    m_face_regions[newFaces[i*2+1]] = oldRegions[i];
  }

  //Now update to reflect motion of the midpoint to the desired position
  setVertexPosition(v_new, midpoint);

  //Volumes haven't changed, but thicknesses have due to area change.
  for(unsigned int i = 0; i < oldFaces.size(); ++i) {
//    m_thicknesses[newFaces[i*2]] = m_volumes[newFaces[i*2]] / getArea(newFaces[i*2], true);
//    m_thicknesses[newFaces[i*2+1]] = m_volumes[newFaces[i*2+1]] / getArea(newFaces[i*2+1], true);
  }

//  if(m_momentum_conserving_remesh) {
//    //recompute and update the masses
//    recomputeVertexMass(v0);
//    m_obj->updateVertexMass(v0);
//    recomputeVertexMass(v1);
//    m_obj->updateVertexMass(v1);
//    for(unsigned int i = 0; i < nbrVerts.size(); ++i) {
//      recomputeVertexMass(nbrVerts[i]);
//      m_obj->updateVertexMass(nbrVerts[i]);
//    }
//
//    //Now adjust the velocity of the new vertex to reflect conservation of momentum
//    Vec3d momSumBefore = mom0 + mom1;
//    Vec3d momSumAfter = getMass(v0) * getVertexVelocity(v0) + getMass(v1) * getVertexVelocity(v1);
//    for(unsigned int i = 0; i < momenta.size(); ++i) {
//      momSumBefore += momenta[i];
//      momSumAfter += getMass(nbrVerts[i]) * getVertexVelocity(nbrVerts[i]);
//    }
//  
//    Vec3d momDif = momSumBefore - momSumAfter;
//  
//    //assign missing momentum to the new vertex, based on its mass
//    recomputeVertexMass(v_new);
//    m_obj->updateVertexMass(v_new);
//    setVertexVelocity(v_new, momDif / getMass(v_new));
//  }

  VertexFaceIterator vf_iter = m_obj->vf_iter(v_new);
  for(;vf_iter; ++vf_iter)
    setFaceActive(*vf_iter);

  new_vert = v_new;
}

void ElasticShell::performCollapse(const EdgeHandle& eh, const VertexHandle& vert_to_remove, const VertexHandle& vert_to_keep, const Vec3d& new_position) {
  
  //determine area of collapsing faces
  EdgeFaceIterator efit = m_obj->ef_iter(eh);
//  Scalar totalVolumeLoss = 0;
//  for(;efit; ++efit)
//    totalVolumeLoss += getVolume(*efit);

  //use the new position to determine how to lerp the vertex data.
  
  Vec3d p0 = getVertexPosition(vert_to_remove);
  Vec3d p1 = getVertexPosition(vert_to_keep);

  //store momenta for later
//  Vec3d momKeep = getVertexVelocity(vert_to_keep) * getMass(vert_to_keep), 
//        momRemove = getVertexVelocity(vert_to_remove) * getMass(vert_to_remove);
  
  //get data for the surrounding nbr verts too
//  std::vector<VertexHandle> nbrs;
//  std::vector<Vec3d> nbrMoms;
//  if(m_momentum_conserving_remesh) {
//    for(EdgeFaceIterator efit = m_obj->ef_iter(eh); efit; ++efit) {
//      FaceHandle fh = *efit;
//      VertexHandle thirdVert;
//      if(!getFaceThirdVertex(*m_obj, fh, eh, thirdVert))
//        std::cout << "Error, no 3rd vertex found.\n";
//
//      nbrs.push_back(thirdVert);
//      Vec3d mom = getVertexVelocity(thirdVert) * getMass(thirdVert);
//      nbrMoms.push_back(mom);
//    }
//  }

  

  Vec3d dx(p1-p0);
  double m2 = dx.squaredNorm();
  Scalar s = clamp((p1-new_position).dot(dx)/m2, 0., 1.);

  Vec3d newVelocity = s*getVertexVelocity(vert_to_remove) + (1-s)*getVertexVelocity(vert_to_keep);
//  Vec3d newUndef = s*getVertexUndeformed(vert_to_remove) + (1-s)*getVertexUndeformed(vert_to_keep);

  //do the collapse itself
  std::vector<EdgeHandle> deletedEdges;
  VertexHandle result = m_obj->collapseEdge(eh, vert_to_remove, deletedEdges);
  if(!result.isValid()) {
    std::cout << "\nError: Refused to perform edge collapse!\n\n";
  }
  
  //determine new positions
  setVertexPosition(vert_to_keep, new_position);
//  setUndeformedVertexPosition(vert_to_keep, newUndef);
  setVertexVelocity(vert_to_keep, newVelocity);
  //increment the thickness of the nearby faces to account for the lost volume

  //sum up the total area increases in the faces
//  VertexFaceIterator vfit = m_obj->vf_iter(vert_to_keep);
//  Scalar totalNewArea = 0;
//  for(;vfit; ++vfit) {
//    FaceHandle fh = *vfit;
//    totalNewArea += std::max(0.0, getArea(fh, true) - m_volumes[fh] / m_thicknesses[fh]); //only consider increases
//  }

  //add the volume losses in the surrounding faces
//  vfit = m_obj->vf_iter(vert_to_keep);
//  for(;vfit; ++vfit) {
//    FaceHandle fh = *vfit;
//    Scalar areaLoss = m_volumes[fh] / m_thicknesses[fh] - getArea(fh, true);
//    if(areaLoss > 0) {
//      totalVolumeLoss += areaLoss * m_thicknesses[fh];
//    }
//  }

//  vfit = m_obj->vf_iter(vert_to_keep);
//  for(;vfit; ++vfit) {
//    FaceHandle fh = *vfit;
//    Scalar newArea = getArea(fh) - m_volumes[fh] / m_thicknesses[fh];
//    if(newArea <= 0) { //face lost area
//      //keep the old thickness, assign the new volume according to the adjusted smaller area
//      m_volumes[fh] = getArea(fh)*m_thicknesses[fh];
//    }
//    else { //face gained area
//      //distribute this expanding face some of the volume
//      m_volumes[fh] += (newArea / totalNewArea) * totalVolumeLoss; //allocate extra volume proportional to area increases
//      m_thicknesses[fh] = m_volumes[fh] / getArea(fh);
//    }
//  }

  //Use the momenta to determine new velocities for the involved vertices
//  if(m_momentum_conserving_remesh) {
//    //First determine new masses.
//    recomputeVertexMass(vert_to_keep);
//    for(unsigned int i = 0; i < nbrs.size(); ++i)
//      recomputeVertexMass(nbrs[i]);
//  
//    //Update the vertex masses with the position dofs model
//    m_obj->updateVertexMass(vert_to_keep);
//    for(unsigned int i = 0; i < nbrs.size(); ++i)
//      m_obj->updateVertexMass(nbrs[i]);
//  
//    //Set new velocities accordingly, based on momentum
//    Scalar massNew0 = getMass(vert_to_keep);
//    setVertexVelocity(vert_to_keep, (momKeep + momRemove) / massNew0);
//  
//    for(unsigned int i = 0; i < nbrs.size(); ++i) {
//      Scalar newmass = getMass(nbrs[i]);
//      setVertexVelocity(nbrs[i], nbrMoms[i] / newmass);  
//    }
//  }

}


bool ElasticShell::performFlip(const EdgeHandle& eh, const FaceHandle f0, const FaceHandle& f1, FaceHandle & new_f0, FaceHandle & new_f1, EdgeHandle& newEdge) {

  
  VertexHandle v0 = m_obj->fromVertex(eh),
               v1 = m_obj->toVertex(eh);
  
  VertexHandle v2, v3;
  getFaceThirdVertex(*m_obj, f0, eh, v2);
  getFaceThirdVertex(*m_obj, f1, eh, v3);

  Vec3d x0 = getVertexPosition(v0),
        x1 = getVertexPosition(v1),
        x2 = getVertexPosition(v2),
        x3 = getVertexPosition(v3);

//  Scalar mass0 = getMass(v0);
//  Scalar mass1 = getMass(v1);
//  Scalar mass2 = getMass(v2);
//  Scalar mass3 = getMass(v3);

  Scalar edgeLen = (x0-x1).norm();
  Scalar edgeLen2 = edgeLen*edgeLen;
  //Get the perpendicular heights of the two triangles.
  Scalar height0 = getArea(f0) / edgeLen, 
         height1 = getArea(f1) / edgeLen;
  
  //Find the perpendicular closest point of each opposite vertex to the edge being flipped.
  Scalar s0 = (x1 - x0).dot(x2 - x0) / edgeLen2;
  Scalar s1 = (x1 - x0).dot(x3 - x0) / edgeLen2;
  
  //areas of the resulting triangles.
  Scalar areaNew0 = 0.5*(x2-x3).cross(x0-x3).norm();
  Scalar areaNew1 = 0.5*(x2-x3).cross(x1-x3).norm();

  /*Scalar thickNew0, thickNew1;
  Scalar volNew0, volNew1;*/

  if(s0 > 0 && s0 < 1 && s1 > 0 && s1 < 1) {
    //Do split-based volume redistribution
    //The idea is to find the split-point assuming it is on the straight geodesic line
    //on the existing triangles between the two opposing vertices.
    //This split point is used to conceptually subdivide the tris, and redistribute
    //the volume to the new triangles. The goal is to more closely approximate the original
    //thicknesses in the new geometry.

    //Construct the closest points.
    Vec3d cp0 = (1 - s0) * x0 + s0 * x1;
    Vec3d cp1 = (1 - s1) * x0 + s1 * x1;

    //By considering the similar triangles constructed from the triangle heights, we know
    //where the hypothetical split point should be.
    Scalar lerpFrac = height0 / (height0 + height1);
    Vec3d splitPoint = (1 - lerpFrac) * cp0 +  lerpFrac * cp1;

    //Compute 4 sub-areas based on the split point
    Scalar area0 = 0.5*(x2 - x0).cross(splitPoint-x0).norm();
    Scalar area1 = 0.5*(x2 - x1).cross(splitPoint-x1).norm();

    Scalar area2 = 0.5*(x3 - x0).cross(splitPoint-x0).norm();
    Scalar area3 = 0.5*(x3 - x1).cross(splitPoint-x1).norm();

    //determine the volumes of each sub-area
    Scalar vFrac0 = area0 / (area0+area1);
//    Scalar vol0 = vFrac0 * m_volumes[f0];
//    Scalar vol1 = (1-vFrac0) * m_volumes[f0];

    Scalar vFrac1 = area2 / (area2+area3);
//    Scalar vol2 = vFrac1 * m_volumes[f1];
//    Scalar vol3 = (1-vFrac1) * m_volumes[f1];

    //redistribute those volumes to the resulting triangles
//    volNew0 = vol0 + vol2;
//    volNew1 = vol1 + vol3;
    
    //now recover the thicknesses of those two triangles from their areas
//    thickNew0 = volNew0 / areaNew0; //associated to the "from" vertex, v0
//    thickNew1 = volNew1 / areaNew1; //associated to the "to" vertex, v1
  }
  else {
    //do simpler, smeared averaging redistribution (splitting-based version doesn't make sense for certain geometries)
//    Scalar totalVolume = m_volumes[f0] + m_volumes[f1];
//    thickNew0 = thickNew1 = totalVolume / (areaNew0 + areaNew1);
//    volNew0 = thickNew0 * areaNew0;
//    volNew1 = thickNew1 * areaNew1;
  }

  Vec2i oldLabels = m_face_regions[f0]; //assume f0 and f1 have the same labels, because if not this flip is nonsense anyway
  assert(oldLabels == m_face_regions[f1]);

  FaceHandle f0new, f1new;
  newEdge = flipEdge(*m_obj, eh, f0, f1, f0new, f1new);
  if(!newEdge.isValid()) {
    std::cout << "\nERROR: Edge flip failed for some reason...\n\n";
    return false;
  }
  
  setFaceActive(f0new);
  setFaceActive(f1new);

  //make sure we are assigning thicknesses to the correct fact
  VertexHandle testVert;
  getFaceThirdVertex(*m_obj, f0new, newEdge, testVert);
  if(testVert != v0)
    swap(f0new, f1new);

//  m_thicknesses[f0new] = thickNew0;
//  m_thicknesses[f1new] = thickNew1;
//  m_volumes[f0new] = volNew0;
//  m_volumes[f1new] = volNew1;

//  if(m_momentum_conserving_remesh) {
//    //now adjust velocities to achieve momentum conservation
//    //get old momenta
//    Vec3d mom0 = mass0*getVertexVelocity(v0), 
//          mom1 = mass1*getVertexVelocity(v1), 
//          mom2 = mass2*getVertexVelocity(v2),
//          mom3 = mass3*getVertexVelocity(v3);
//  
//    //update the mass of the vertex in the shell
//    recomputeVertexMass(v0); 
//    recomputeVertexMass(v1); 
//    recomputeVertexMass(v2); 
//    recomputeVertexMass(v3);
//  
//    //call up to pos_dofs_model to recompute the full mass of the node (may include rods and such)
//    m_obj->updateVertexMass(v0);
//    m_obj->updateVertexMass(v1);
//    m_obj->updateVertexMass(v2);
//    m_obj->updateVertexMass(v3);
//
//    Scalar massNew0 = getMass(v0), 
//      massNew1 = getMass(v1), 
//      massNew2 = getMass(v2), 
//      massNew3 = getMass(v3);
//
//    //set new velocities accordingly
//    setVertexVelocity(v0, mom0 / massNew0);
//    setVertexVelocity(v1, mom1 / massNew1);
//    setVertexVelocity(v2, mom2 / massNew2);
//    setVertexVelocity(v3, mom3 / massNew3);
//  }

  //keep previous region labels
  m_face_regions[f0new] = oldLabels;
  m_face_regions[f1new] = oldLabels;

  new_f0 = f0new;
  new_f1 = f1new;
  
  return true;
}
  
void ElasticShell::performZippering(EdgeHandle e0, EdgeHandle e1, const std::vector<FaceHandle> & faces_deleted, const std::vector<std::vector<VertexHandle> > & faces_to_create, const std::vector<Vec2i> & face_labels_to_create, const std::vector<std::pair<FaceHandle, Vec2i> > & face_labels_to_change, std::vector<FaceHandle> & faces_created)
{
  assert(faces_deleted.size() == 2);
  
  VertexHandle v00 = m_obj->fromVertex(e0);
  VertexHandle v01 = m_obj->toVertex(e0);
  VertexHandle v10 = m_obj->fromVertex(e1);
  VertexHandle v11 = m_obj->toVertex(e1);

  EdgeFaceIterator efit0 = m_obj->ef_iter(e0);  assert(efit0);
  FaceHandle f00 = *efit0;  ++efit0;  assert(efit0);
  FaceHandle f01 = *efit0;  ++efit0;  assert(!efit0);
  EdgeFaceIterator efit1 = m_obj->ef_iter(e1);  assert(efit1);
  FaceHandle f10 = *efit1;  ++efit1;  assert(efit1);
  FaceHandle f11 = *efit1;  ++efit1;  assert(!efit1);

  // the deleted faces
  assert((f00 == faces_deleted[0] && f01 == faces_deleted[1]) || (f00 == faces_deleted[1] && f01 == faces_deleted[0]));
//  assert((f10 == faces_deleted[2] && f11 == faces_deleted[3]) || (f10 == faces_deleted[3] && f11 == faces_deleted[2]));
  assert((f10 != faces_deleted[2] && f10 != faces_deleted[3]) || (f11 != faces_deleted[3] && f11 != faces_deleted[2]));
  
  VertexHandle v02, v03;
  getFaceThirdVertex(*m_obj, f00, e0, v02);
  getFaceThirdVertex(*m_obj, f01, e0, v03);
  VertexHandle v12, v13;
  getFaceThirdVertex(*m_obj, f10, e1, v12);
  getFaceThirdVertex(*m_obj, f11, e1, v13);

  Vec3d x00 = getVertexPosition(v00);
  Vec3d x01 = getVertexPosition(v01);
  Vec3d x02 = getVertexPosition(v02);
  Vec3d x03 = getVertexPosition(v03);
  Vec3d x10 = getVertexPosition(v10);
  Vec3d x11 = getVertexPosition(v11);
  Vec3d x12 = getVertexPosition(v12);
  Vec3d x13 = getVertexPosition(v13);
  
  for (size_t i = 0; i < face_labels_to_change.size(); i++)
  {
    m_face_regions[face_labels_to_change[i].first] = face_labels_to_change[i].second;
  }
  
  // no change to vertex mass/velocity for now; this doesn't conserve volume or momentum
  assert(faces_to_create.size() == face_labels_to_create.size());
  for (size_t i = 0; i < faces_to_create.size(); i++)
  {
    FaceHandle f = m_obj->addFace(faces_to_create[i][0], faces_to_create[i][1], faces_to_create[i][2]);
    assert(f.isValid());
    faces_created.push_back(f);

    m_face_regions[f] = Vec2i(face_labels_to_create[i]);
  }

  for (size_t i = 0; i < faces_deleted.size(); i++)
  {
    bool success = m_obj->deleteFace(faces_deleted[i], false);
    assert(success);
  }
  bool success = m_obj->deleteEdge(e0, false);
  assert(success);
    
}



void ElasticShell::updateThickness() {
  FaceIterator fit = m_obj->faces_begin();
  for(;fit != m_obj->faces_end(); ++fit) {
    FaceHandle f = *fit;
    if(!isFaceActive(f)) continue;

    //compute face area
    std::vector<Vec3d> verts(3);
    int i = 0;
    for(FaceVertexIterator fvit = m_obj->fv_iter(f); fvit; ++fvit, ++i)
      verts[i] = getVertexPosition(*fvit);
    Vec3d v0 = verts[1] - verts[0]; Vec3d v1 = verts[2] - verts[0];
    Scalar area = 0.5f*(v0.cross(v1)).norm();

    //figure out new thickness from previous volume, and update it
//    Scalar newThickness = m_volumes[f] / area;
//    m_thicknesses[f] = newThickness;
  }

}

void ElasticShell::deleteRegion() {
  //Delete faces that have entered the deletion zone.

  std::vector<FaceHandle> faces_to_remove;
  for(FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit) {
    //compute barycentre of the face, and if it's in the deletion region, kill it.
    FaceHandle fh = *fit;
    Vec3d barycentre(0,0,0);
    for(FaceVertexIterator fvit = m_obj->fv_iter(fh); fvit; ++fvit) {
      VertexHandle vh = *fvit;
      Vec3d pos = getVertexPosition(vh);
      barycentre += pos;

    }
    barycentre /= 3.0;
    
    //check if barycentre is in the deletion region
    if(barycentre[0] > m_delete_lower[0] && barycentre[1] >  m_delete_lower[1] && barycentre[2] > m_delete_lower[2] &&
       barycentre[0] < m_delete_upper[0] && barycentre[1] < m_delete_upper[1] && barycentre[2] < m_delete_upper[2]) {
      faces_to_remove.push_back(fh);
    }
  }
  for(unsigned int i = 0; i < faces_to_remove.size(); ++i) {
    removeFace(faces_to_remove[i]);
  }

}

void ElasticShell::removeFace(FaceHandle& f) {
  //safely remove a face
  //and all the constraints/springs that link to it.

  //collect the vertices that comprise the face
  VertexHandle faceVerts[3];
  FaceVertexIterator fvit = m_obj->fv_iter(f);
  int j = 0;
  for(;fvit;++fvit) {
    faceVerts[j] = *fvit;
    ++j;
  }
  
  //remove springs on this face
  m_repulsion_springs->clearSprings(f);

  //delete the face, and recursively, any unused vertices
  m_obj->deleteFace(f, true);

  //now remove the springs and constraints to any vertices deleted as a side effect
  for(int j = 0; j < 3; ++j) {
    if(!m_obj->vertexExists(faceVerts[j])) {
      getDefoObj().releaseVertex(faceVerts[j]);
      m_repulsion_springs->clearSprings(faceVerts[j]);
    }
  }
}

void ElasticShell::extendMesh(Scalar current_time) {

  for(unsigned int boundary = 0; boundary < m_inflow_boundaries.size(); ++boundary) {
    //do a simple (not very general) check to see if the inflow has moved far enough away
    EdgeHandle edge0 = m_inflow_boundaries[boundary][0];
    VertexHandle vfrom = m_obj->fromVertex(edge0);
    VertexHandle vto = m_obj->toVertex(edge0);
    Vec3d startPos = m_inflow_positions[boundary][0];
    Vec3d curPos = getVertexPosition(vfrom);
    Vec3d curPos2 = getVertexPosition(vto);
    
    //look at aspect ratio of this triangle, and if it's too bad, skip it this time around
    Scalar baseLength = (curPos - curPos2).norm();
    Scalar len1 = (curPos - startPos).norm();
    Scalar len2 = (curPos2 - startPos).norm();
    if(len1/baseLength < 0.7 || len2 / baseLength < 0.7) {
      continue;
    }


    int count = m_inflow_boundaries[boundary].size();
    int last = m_inflow_boundaries[boundary].size()-1;
    bool direction = m_inflow_lastdir[boundary];
    m_inflow_lastdir[boundary] = !m_inflow_lastdir[boundary]; //flip direction for next time

    //check if open or closed loop
    VertexHandle loopVertex = getSharedVertex(*m_obj, m_inflow_boundaries[boundary][0], m_inflow_boundaries[boundary][last]);

    //VertexHandle prevVert = m_obj->addVertex(); //loopVertex.isValid() ? loopVertex : m_obj->addVertex();
    VertexHandle prevVert = m_obj->addVertex(); //loopVertex.isValid() ? loopVertex : m_obj->addVertex();

    std::vector<VertexHandle> vertices;
    vertices.push_back(prevVert);
    std::vector<FaceHandle> faces;
    
    VertexHandle sharedVert = getSharedVertex(*m_obj, m_inflow_boundaries[boundary][0], m_inflow_boundaries[boundary][1]);
    VertexHandle prevLowerVert = getEdgesOtherVertex(*m_obj, m_inflow_boundaries[boundary][0], sharedVert);
    
    EdgeHandle prevEdge = m_obj->addEdge(prevLowerVert, prevVert);
    VertexHandle loopTopVertex = prevVert;
    EdgeHandle startEdge = prevEdge; //save this for wrapping around.
    
    std::vector<EdgeHandle> newList;
    for(unsigned int edge = 0; edge < m_inflow_boundaries[boundary].size(); ++edge) {
      
      getDefoObj().releaseVertex(prevLowerVert);

      EdgeHandle eh1 = m_inflow_boundaries[boundary][edge];
      EdgeHandle eh2 = m_inflow_boundaries[boundary][(edge+1)%count];

      VertexHandle sharedVert, otherVert; 
      if((int)edge < count - 1 || loopVertex.isValid()) { //a middle edge, or the last edge if things loop around
        sharedVert = getSharedVertex(*m_obj, eh1, eh2);
        otherVert = getEdgesOtherVertex(*m_obj, eh1, sharedVert);
      }
      else {
        otherVert = prevLowerVert;
        sharedVert = getEdgesOtherVertex(*m_obj, eh1, otherVert);
      }

      VertexHandle newVert;
      EdgeHandle newEdge3;
      if(edge == count-1 && loopVertex.isValid()) {
        newVert = loopTopVertex;
        newEdge3 = startEdge;
      }
      else {
        newVert = m_obj->addVertex();
        vertices.push_back(newVert);
        newEdge3 = m_obj->addEdge(sharedVert, newVert);
      }
      
      EdgeHandle newEdge4 = m_obj->addEdge(prevVert, newVert);

      EdgeHandle newEdge2;
      FaceHandle newFace1, newFace2;
      if( (edge % 2 == 0 && direction) || (edge % 2 == 1 && !direction)) { //A option
        newEdge2 = m_obj->addEdge(prevVert, sharedVert);

        newFace1 = m_obj->addFace(prevEdge, eh1, newEdge2);
        newFace2 = m_obj->addFace(newEdge2, newEdge3, newEdge4);
      }
      else { //B option
        newEdge2 = m_obj->addEdge(otherVert, newVert);

        newFace1 = m_obj->addFace(prevEdge, newEdge2, newEdge4);
        newFace2 = m_obj->addFace(eh1, newEdge3, newEdge2);
      }

      setFaceActive(newFace1);
      setFaceActive(newFace2);
      faces.push_back(newFace1);
      faces.push_back(newFace2);
      //Todo: mess with ordering to automatically determine the correct front/back faces
      //based on the faces we are connecting to.
      
      newList.push_back(newEdge4);

      //advance to the next edge
      prevVert = newVert;
      prevLowerVert = sharedVert;
      prevEdge = newEdge3;
    }
    getDefoObj().releaseVertex(prevLowerVert);
    

    m_inflow_boundaries[boundary] = newList;

    for(unsigned int i = 0; i < vertices.size(); ++i) {
      setVertexPosition(vertices[i], m_inflow_positions[boundary][i]);
//      setUndeformedVertexPosition(vertices[i], m_inflow_positions[boundary][i]);
      setVertexVelocity(vertices[i], m_inflow_velocities[boundary][i]);
//      m_vertex_masses[vertices[i]] = 0;
      m_obj->setVertexDampingUndeformedPosition(vertices[i], m_inflow_positions[boundary][i]);

      getDefoObj().constrainVertex(vertices[i], new FixedVelocityConstraint(m_inflow_positions[boundary][i], m_inflow_velocities[boundary][i], current_time));
    }


    for(unsigned int i = 0; i < faces.size(); ++i) {
//      m_thicknesses[faces[i]] = m_inflow_thickness;
//      m_volumes[faces[i]] = 0;
//      m_volumes[faces[i]] = getArea(faces[i])*m_inflow_thickness;
      FaceVertexIterator fvit = m_obj->fv_iter(faces[i]);
      for(;fvit;++fvit) {
        VertexHandle vh = *fvit;
//        m_vertex_masses[vh] += m_volumes[faces[i]] * m_density / 3.0;
      }
    }

  }
  
//  computeMasses();

}

void ElasticShell::setInflowSection(std::vector<EdgeHandle> edgeList, const Vec3d& vel, Scalar thickness) {
  m_inflow = true;
  m_inflow_boundaries.push_back(edgeList);
  //m_inflow_velocity.push_back(vel);
  m_inflow_thickness = thickness;
  m_inflow_lastdir.push_back(false);

  VertexHandle prevVert;
  std::vector<Vec3d> posList,velList;
  for(unsigned int edge = 0; edge < edgeList.size(); ++edge) {
    EdgeHandle eh1 = edgeList[edge];
    EdgeHandle eh2 = edgeList[(edge+1)%edgeList.size()];
    
    VertexHandle sharedVert, otherVert;
    if(edge == edgeList.size() - 1) {
      otherVert = prevVert;
      sharedVert = getEdgesOtherVertex(*m_obj, eh1, otherVert); 
    }
    else {
      sharedVert = getSharedVertex(*m_obj, eh1, eh2);
      otherVert = getEdgesOtherVertex(*m_obj, eh1, sharedVert);
      
    }

    Vec3d pos = getVertexPosition(otherVert);
    getDefoObj().constrainVertex(otherVert, new FixedVelocityConstraint(pos, vel, 0));
    posList.push_back(pos);
    velList.push_back(vel);
    
    prevVert = sharedVert;
    
  }
  
  VertexHandle wrapVert = getSharedVertex(*m_obj, edgeList[0], edgeList[edgeList.size()-1]);
  if(!wrapVert.isValid()) {
    Vec3d pos = getVertexPosition(prevVert);
    getDefoObj().constrainVertex(prevVert, new FixedVelocityConstraint(pos, vel, 0));
    posList.push_back(pos);
    velList.push_back(vel);
  } 

  m_inflow_positions.push_back(posList);
  m_inflow_velocities.push_back(velList);

}

void ElasticShell::setDeletionBox(const Vec3d& lowerBound, const Vec3d& upperBound) {
    m_delete_lower = lowerBound;
    m_delete_upper = upperBound;
    m_delete_region = true;
}


} //namespace BASim
