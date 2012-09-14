#include "BASim/src/Collisions/ElTopo/broadphasegrid.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShell.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"
#include "BASim/src/Collisions/ElTopo/ccd_wrapper.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/CSTMembraneForce.hh"
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

ElasticShell::ElasticShell(DeformableObject* object, const FaceProperty<char>& shellFaces, Scalar timestep) : 
  PhysicalModel(*object), m_obj(object), 
    m_active_faces(shellFaces), 
    m_undef_xi(object),
    m_damping_undef_xi(object),
    m_vertex_masses(object),
    m_edge_masses(object),
    m_thicknesses(object),
    m_face_regions(object),
    m_volumes(object),
    m_xi(object), 
    m_xi_vel(object),
    m_density(1),
    m_collision_epsilon(1e-5),
    m_vert_point_springs(NULL),
    m_repulsion_springs(NULL),
    m_sphere_collisions(false),
    m_object_collisions(false),
    m_ground_collisions(false),
    m_do_eltopo_collisions(false),
    m_do_thickness_updates(true)
{
  m_vert_point_springs = new ShellVertexPointSpringForce(*this, "VertPointSprings", timestep);
  m_repulsion_springs = new ShellStickyRepulsionForce(*this, "RepulsionSprings", timestep);

  addForce(m_vert_point_springs);
  addForce(m_repulsion_springs);
}

ElasticShell::~ElasticShell() {
  delete m_vert_point_springs;
  delete m_repulsion_springs;
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

void ElasticShell::setThickness( Scalar thickness )
{
  for(FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit) {
    m_thicknesses[*fit] = thickness;
    Scalar area = getArea(*fit, false);
    m_volumes[*fit] = thickness * area;
  }

}

void ElasticShell::setDensity(Scalar density) {
  m_density = density;
}

//void ElasticShell::setVertexUndeformed( const VertexProperty<Vec3d>& undef )
//{
//  m_undeformed_positions = undef;
//}

void ElasticShell::setEdgeUndeformed( const EdgeProperty<Scalar>& undef )
{
  m_undef_xi = undef;
}

//void ElasticShell::setVertexPositions( const VertexProperty<Vec3d>& positions )
//{
//  m_positions = positions;
//}

void ElasticShell::setEdgeXis( const EdgeProperty<Scalar>& positions )
{
  m_xi = positions;
}

//void ElasticShell::setVertexVelocities( const VertexProperty<Vec3d>& velocities) 
//{
//  m_velocities = velocities;
//}

void ElasticShell::setEdgeVelocities(const EdgeProperty<Scalar>& velocities)
{
  m_xi_vel = velocities;
}

void ElasticShell::setFaceLabels(const FaceProperty<Vec2i>& labels) {
  m_face_regions = labels;
}


Scalar ElasticShell::getThickness(const VertexHandle& vh) const {
  Scalar totalA = 0.0;
  Scalar w;
  Scalar total = 0.0;
  for (VertexFaceIterator vfit = m_obj->vf_iter(vh); vfit; ++vfit){
      w = getArea(*vfit);
      totalA += w;
      total += w*m_thicknesses[*vfit];
  }

  assert ( totalA > 0.);
  assert ( total > 0.);
  return total / totalA;
}
Scalar ElasticShell::getMaxThickness() const {
  Scalar maxVal = -1000000;
  for( FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit ){
      if (m_thicknesses[*fit] > maxVal ) maxVal = m_thicknesses[*fit];
  }
  return maxVal;
}
Scalar ElasticShell::getMinThickness() const {
  Scalar minVal = 1000000;
  for( FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit ){
      if (m_thicknesses[*fit] < minVal ) minVal = m_thicknesses[*fit];
  }
  return minVal;
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


void ElasticShell::getThickness(VertexProperty<Scalar> & vThickness) const{
    for ( VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit){
        vThickness[*vit] = getThickness(*vit);
    }
}

Scalar ElasticShell::getArea(const FaceHandle& f, bool current) const  {
  FaceVertexIterator fvit = m_obj->fv_iter(f);
  VertexHandle v0_hnd = *fvit; ++fvit; assert(fvit);
  VertexHandle v1_hnd = *fvit; ++fvit; assert(fvit);
  VertexHandle v2_hnd = *fvit; ++fvit; assert(!fvit);

  //compute triangle areas
  if(current) 
  {
    Vec3d pos0 = m_obj->getVertexPosition(v0_hnd);
    Vec3d pos1 = m_obj->getVertexPosition(v1_hnd);
    Vec3d pos2 = m_obj->getVertexPosition(v2_hnd);
    
    Vec3d v0 = pos1 - pos0;
    Vec3d v1 = pos2 - pos0;
    Vec3d triVec = v0.cross(v1);
    return 0.5*triVec.norm();
  }
  else 
  {
    Vec3d pos0 = m_obj->getVertexUndeformedPosition(v0_hnd);
    Vec3d pos1 = m_obj->getVertexUndeformedPosition(v1_hnd);
    Vec3d pos2 = m_obj->getVertexUndeformedPosition(v2_hnd);
    
    Vec3d v0 = pos1 - pos0;
    Vec3d v1 = pos2 - pos0;
    Vec3d triVec = v0.cross(v1);
    return 0.5*triVec.norm();
  }
}

void ElasticShell::computeMasses()
{
  //Compute vertex masses in a lumped mass way.

  m_vertex_masses.assign(0);
  m_edge_masses.assign(0);

  Scalar area = 0;

  //Iterate over all triangles active in this shell and accumulate vertex masses
  for(FaceIterator f_iter = m_obj->faces_begin(); f_iter != m_obj->faces_end(); ++f_iter) {
    FaceHandle& f_hnd = *f_iter;
    if(m_active_faces[f_hnd]) {

      //get the three vertices
      FaceVertexIterator fvit = m_obj->fv_iter(f_hnd);
      VertexHandle v0_hnd = *fvit; ++fvit; assert(fvit);
      VertexHandle v1_hnd = *fvit; ++fvit; assert(fvit);
      VertexHandle v2_hnd = *fvit; ++fvit; assert(!fvit);

      //compute triangle areas
      Vec3d v0 = getVertexPosition(v1_hnd) - getVertexPosition(v0_hnd);
      Vec3d v1 = getVertexPosition(v2_hnd) - getVertexPosition(v0_hnd);
      Vec3d triVec = v0.cross(v1);
      Scalar area = 0.5*sqrt(triVec.dot(triVec)) / 3.0;
      Scalar contribution = m_thicknesses[f_hnd] * m_density * area;
      
      //accumulate mass to the vertices
      m_vertex_masses[v0_hnd] += contribution;
      m_vertex_masses[v1_hnd] += contribution;
      m_vertex_masses[v2_hnd] += contribution;

      //set edge masses to zero, since we want to solve them quasistatically (they're derivative DOFs)
      FaceEdgeIterator feit = m_obj->fe_iter(f_hnd);
      EdgeHandle e0_hnd = *feit; ++feit; assert(feit);
      EdgeHandle e1_hnd = *feit; ++feit; assert(feit);
      EdgeHandle e2_hnd = *feit; ++feit; //assert(feit);

      m_edge_masses[e0_hnd] = 0;
      m_edge_masses[e1_hnd] = 0;
      m_edge_masses[e2_hnd] = 0;

      //store the current volumes
      m_volumes[f_hnd] = 3*area*m_thicknesses[f_hnd];
    }
  }
 
  m_obj->updateVertexMasses();
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
  //they're all edge Dofs for a shell
  assert(hnd.getType() == DofHandle::EDGE_DOF);

  //return reference to the appropriate position in the vector
  const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
  return const_cast<Scalar&>(m_xi[eh]);
}

void ElasticShell::setDof( const DofHandle& hnd, const Scalar& dof )
{
  //they're all vertex Dofs for a shell
  assert(hnd.getType() == DofHandle::EDGE_DOF);

  const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
  m_xi[eh] = dof;
}

const Scalar& ElasticShell::getVel( const DofHandle& hnd ) const
{
  assert(hnd.getType() == DofHandle::EDGE_DOF);

  //return reference to the appropriate position in the vector
  const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
  return const_cast<Scalar&>(m_xi_vel[eh]);
}

void ElasticShell::setVel( const DofHandle& hnd, const Scalar& vel )
{
  assert(hnd.getType() == DofHandle::EDGE_DOF);

  const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
  m_xi_vel[eh] = vel;
}

const Scalar& ElasticShell::getMass( const DofHandle& hnd ) const
{
  assert(hnd.getType() == DofHandle::EDGE_DOF);

  const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
  return m_edge_masses[eh];
}

void ElasticShell::getScriptedDofs( IntArray& dofIndices, std::vector<Scalar>& dofValues, Scalar time ) const
{
    // position dof scripting is moved to PositionDofsModel.
    
    for(unsigned int i = 0; i < constrainedEdges.size(); ++i) {
        int dofID = getEdgeDofBase(constrainedEdges[i]);
        dofIndices.push_back(dofID);
        dofValues.push_back(constrainedXiValues[i]);
    }
}

void ElasticShell::startStep(Scalar time, Scalar timestep)
{
  std::cout << "Starting startStep\n";

  //update the damping "reference configuration" for computing viscous forces.
  m_damping_undef_xi = m_xi;

  //tell the forces to update anything they need to update
  const std::vector<ElasticShellForce*>& forces = getForces();
  for(unsigned int i = 0; i < forces.size(); ++i) {
    forces[i]->update();
  }
  std::cout << "Done startStep\n";
}

void ElasticShell::resolveCollisions(Scalar timestep) {
  std::cout << "Resolving collisions with El Topo\n";
  //Convert the data to the form required by El Topo!
  std::vector<ElTopo::Vec3d> vert_new, vert_old;
  std::vector<ElTopo::Vec3st> tri_data;
  std::vector<Scalar> masses;

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
    Scalar mass = getMass(vh);

    vert_new.push_back(ElTopo::Vec3d(vert[0], vert[1], vert[2]));
    vert_old.push_back(ElTopo::Vec3d(old_vert[0], old_vert[1], old_vert[2]));
    if(getDefoObj().isConstrained(vh)) {
      masses.push_back(numeric_limits<Scalar>::infinity());
    }
    else {
      masses.push_back(mass);
    }
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
  ElTopo::DynamicSurface dynamic_surface( vert_old, tri_data, masses, m_collision_epsilon, friction_coeff, true, false );

  dynamic_surface.set_all_newpositions( vert_new );

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

  std::cout << "Starting endStep.\n";
  bool do_relabel = false;

  std::cout << "Vertex count: " << m_obj->nv() << std::endl;

  if (m_obj->nf() == 0) // ElTopo crashes if given a mesh with zero triangles, which is possible for example when we're only testing rods.
    return;
//  return;
  
  //El Topo collision processing.
  
  if(m_do_eltopo_collisions)
    resolveCollisions(timestep);
  
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

  //Adjust thicknesses based on area changes
  if(m_do_thickness_updates)
    updateThickness();

  if(m_inflow) {
    extendMesh(time);
    do_relabel = true;
  }

  if(m_delete_region) {
    deleteRegion();
    do_relabel = true;
  }

  //Remeshing
  if(m_do_remeshing) {
    std::cout << "Remeshing\n";
    remesh();
    std::cout << "Completed remeshing\n";

    //Relabel DOFs if necessary
    do_relabel = true;
  }

  static int only_twice = 0;

  //Tearing processing
  if(m_tearing){
    std::cout << "Processing tearing. \n";
    //fracture();
    fracture();
    do_relabel = true;
    only_twice++;
  }

  if(do_relabel) {
    std::cout << "Re-indexing\n";
    m_obj->computeDofIndexing();
  }

  //Update masses based on new areas/thicknesses
  computeMasses();

  std::cout << "Completed endStep\n";

}


void ElasticShell::fracture() {

  //Set up a SurfTrack, run remeshing, render the new mesh
  ElTopo::SurfTrackInitializationParameters construction_parameters;
  construction_parameters.m_proximity_epsilon = m_collision_epsilon;
  construction_parameters.m_allow_vertex_movement = false;
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
  std::vector<Scalar> masses;

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
    Scalar mass = getMass(vh);
    vert_data.push_back(ElTopo::Vec3d(vert[0], vert[1], vert[2]));
    if(getDefoObj().isConstrained(vh))
      masses.push_back(numeric_limits<Scalar>::infinity());
    else
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
      masses[vert_numbers[vh]] = numeric_limits<Scalar>::infinity();
    }
  }

  
  ElTopo::SurfTrack surface_tracker( vert_data, tri_data, masses, construction_parameters ); 

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
      setUndeformedVertexPosition(nv, getVertexUndeformed(source));
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
      m_thicknesses[newFaceHandle] = m_thicknesses[oldFaceHandle];
      m_volumes[newFaceHandle] = m_volumes[oldFaceHandle];
        
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
        thickness += w * getThickness(*efit);
        totalA += w;
        if(getThickness(*efit) > m_tear_thres)
          both_thin = false;
       
    }
    thickness /= totalA;
    Scalar p = (Scalar) rand() / (Scalar) RAND_MAX;
    
    return both_thin && (p <  m_tear_rand);
}


void ElasticShell::remesh()
{

  //Set up a SurfTrack, run remeshing, render the new mesh
  ElTopo::SurfTrackInitializationParameters construction_parameters;
  construction_parameters.m_proximity_epsilon = m_collision_epsilon;
  construction_parameters.m_merge_proximity_epsilon = 0.01;
  construction_parameters.m_allow_vertex_movement = false;
  construction_parameters.m_min_edge_length = m_remesh_edge_min_len;
  construction_parameters.m_max_edge_length = m_remesh_edge_max_len;
  construction_parameters.m_max_volume_change = numeric_limits<double>::max();   
  construction_parameters.m_min_triangle_angle = 15;
  construction_parameters.m_max_triangle_angle = 165;
  construction_parameters.m_verbose = false;
  construction_parameters.m_allow_non_manifold = true;
  construction_parameters.m_allow_topology_changes = true;
  construction_parameters.m_collision_safety = true;
  construction_parameters.m_remesh_boundaries = true;
  
  //construction_parameters.m_subdivision_scheme = new ElTopo::MidpointScheme();
  //construction_parameters.m_subdivision_scheme = new ElTopo::QuadraticErrorMinScheme();
  //construction_parameters.m_subdivision_scheme = new ElTopo::ButterflyScheme();
  construction_parameters.m_subdivision_scheme = new ElTopo::ModifiedButterflyScheme();

  construction_parameters.m_use_curvature_when_collapsing = false;
  construction_parameters.m_use_curvature_when_splitting = false;
  //construction_parameters.m_max_curvature_multiplier = 1000;
  //construction_parameters.m_min_curvature_multiplier = 1.0;

  std::vector<ElTopo::Vec3d> vert_data;
  std::vector<ElTopo::Vec3st> tri_data;
  std::vector<Scalar> masses;

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
    Scalar mass = getMass(vh);
    vert_data.push_back(ElTopo::Vec3d(vert[0], vert[1], vert[2]));
    if(getDefoObj().isConstrained(vh))
      masses.push_back(numeric_limits<Scalar>::infinity());
    else
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
      masses[vert_numbers[vh]] = numeric_limits<Scalar>::infinity();
    }
    //and the other vertex
    masses[vert_numbers[verts[i]]] = numeric_limits<Scalar>::infinity();
  }

  std::cout << "Calling surface improvement\n";
  ElTopo::SurfTrack surface_tracker( vert_data, tri_data, masses, construction_parameters ); 

  surface_tracker.improve_mesh();
  surface_tracker.topology_changes();

  std::cout << "Performing " << surface_tracker.m_mesh_change_history.size() << " Improvement Operations:\n";
  for(unsigned int j = 0; j < surface_tracker.m_mesh_change_history.size(); ++j) {
    ElTopo::MeshUpdateEvent event = surface_tracker.m_mesh_change_history[j];
    
 
    if(event.m_type == ElTopo::MeshUpdateEvent::FLAP_DELETE) {
      assert(event.m_deleted_tris.size() == 2);
      for(unsigned int t = 0; t < event.m_deleted_tris.size(); ++t) {
        int triNo = event.m_deleted_tris[t];
        FaceHandle faceToDelete = reverse_trimap[triNo];
        
        assert(m_obj->faceExists(faceToDelete));

        m_obj->deleteFace(faceToDelete, true);
        std::cout << "Deleted flap face\n";
      }
    }
    else if(event.m_type == ElTopo::MeshUpdateEvent::EDGE_COLLAPSE) {
      VertexHandle v0 = reverse_vertmap[event.m_v0];
      VertexHandle v1 = reverse_vertmap[event.m_v1];
      EdgeHandle eh = findEdge(mesh, v0, v1);

      assert(eh.isValid()); //ensure the desired edge exists

      //std::cout << "Collapse\n";
      Vec3d new_pos(event.m_vert_position[0], event.m_vert_position[1], event.m_vert_position[2]);
      VertexHandle dead_vert = reverse_vertmap[event.m_deleted_verts[0]];
      VertexHandle keep_vert = getEdgesOtherVertex(mesh, eh, dead_vert);
      int kept_vert_ET_id = vert_numbers[keep_vert];

      performCollapse(eh, dead_vert, keep_vert, new_pos);

      // Delete the vertex
      reverse_vertmap[event.m_deleted_verts[0]] = VertexHandle(-1); //mark vert as invalid
      
      // Update faces
      for(unsigned int i = 0; i < event.m_created_tri_data.size(); ++i) {
        ElTopo::Vec3st new_face = event.m_created_tri_data[i];

        //determine handles in our indexing
        VertexHandle v0 = reverse_vertmap[new_face[0]];
        VertexHandle v1 = reverse_vertmap[new_face[1]];
        VertexHandle v2 = reverse_vertmap[new_face[2]];
        assert(v0.isValid() && v1.isValid() && v2.isValid());
        
        bool face_matched = false;
        for(VertexFaceIterator vfit = mesh.vf_iter(keep_vert); vfit; ++vfit) {
          FaceHandle face_candidate = *vfit;
          if(isFaceMatch(mesh, face_candidate, v0, v1, v2)) {
            if(reverse_trimap.size() <= event.m_created_tris[i]) 
              reverse_trimap.resize(event.m_created_tris[i]+1);
            reverse_trimap[event.m_created_tris[i]] = face_candidate;
            face_numbers[face_candidate] = event.m_created_tris[i];
            face_matched = true;
            break;
          }
        }

        if(!face_matched) {
          std::cout << "Vertex indices: " << v0.idx() << " " << v1.idx() << " " << v2.idx() << std::endl;
          std::cout << "ERROR: Couldn't match the face - COLLAPSE.\n\n\n";
        }
      }
    }
    else if(event.m_type == ElTopo::MeshUpdateEvent::EDGE_SPLIT) {
      //Identify the edge based on its endpoint vertices instead
      //std::cout << "Split\n";

       VertexHandle v0 = reverse_vertmap[event.m_v0];
       VertexHandle v1 = reverse_vertmap[event.m_v1];
       EdgeHandle eh = findEdge(mesh, v0, v1);

       assert(eh.isValid()); //ensure the desired edge exists

      VertexHandle new_vert;
      Vec3d new_pos(event.m_vert_position[0], event.m_vert_position[1], event.m_vert_position[2]);
      
      //Do the same split as El Topo
      performSplit(eh, new_pos, new_vert);
      
      //now update the mapping between structures vertices
      vert_numbers[new_vert] = event.m_created_verts[0];
      if(reverse_vertmap.size() <= event.m_created_verts[0]) 
        reverse_vertmap.resize(event.m_created_verts[0]+1, VertexHandle(-1));
      reverse_vertmap[event.m_created_verts[0]] = new_vert; //the vertex will always be added at the end by El Topo
      
      // Update faces
      for(unsigned int i = 0; i < event.m_created_tri_data.size(); ++i) {
        ElTopo::Vec3st new_face = event.m_created_tri_data[i];
        
        //determine handles in our indexing
        VertexHandle v0 = reverse_vertmap[new_face[0]];
        VertexHandle v1 = reverse_vertmap[new_face[1]];
        VertexHandle v2 = reverse_vertmap[new_face[2]];
        assert(v0.isValid() && v1.isValid() && v2.isValid());

        bool face_matched = false;
        for(VertexFaceIterator vfit = mesh.vf_iter(new_vert); vfit; ++vfit) {
          FaceHandle face_candidate = *vfit;
          if(isFaceMatch(mesh, face_candidate, v0, v1, v2)) {
            if(reverse_trimap.size() <= event.m_created_tris[i]) 
              reverse_trimap.resize(event.m_created_tris[i]+1);
            reverse_trimap[event.m_created_tris[i]] = face_candidate;
            face_numbers[face_candidate] = event.m_created_tris[i];
            face_matched = true;
            break;
          }
        }
        
        if(!face_matched)
          std::cout << "ERROR: Couldn't match the face - SPLIT.\n\n\n";
      }
    }
    else if(event.m_type == ElTopo::MeshUpdateEvent::EDGE_FLIP) {
      
      //std::cout << "Flip\n";
      //Do the same flip as El Topo
      VertexHandle v0 = reverse_vertmap[event.m_v0];
      VertexHandle v1 = reverse_vertmap[event.m_v1];
      EdgeHandle eh = findEdge(mesh, v0, v1);

      assert(eh.isValid()); //ensure the desired edge exists

      assert(event.m_deleted_tris[i].size() == 2);
      FaceHandle f0 = reverse_trimap[event.m_deleted_tris[0]];
      FaceHandle f1 = reverse_trimap[event.m_deleted_tris[1]];
      EdgeHandle newEdge;
      performFlip(eh, f0, f1, newEdge);
      
      // Update faces
      for(unsigned int i = 0; i < event.m_created_tri_data.size(); ++i) {
        ElTopo::Vec3st new_face = event.m_created_tri_data[i];

        //determine handles in our indexing
        VertexHandle v0 = reverse_vertmap[new_face[0]];
        VertexHandle v1 = reverse_vertmap[new_face[1]];
        VertexHandle v2 = reverse_vertmap[new_face[2]];
        assert(v0.isValid() && v1.isValid() && v2.isValid());

        bool face_matched = false;
        for(EdgeFaceIterator efit = mesh.ef_iter(newEdge); efit; ++efit) {

          FaceHandle face_candidate = *efit;
          if(isFaceMatch(mesh, face_candidate, v0, v1, v2)) {

            if(reverse_trimap.size() <= event.m_created_tris[i]) 
              reverse_trimap.resize(event.m_created_tris[i]+1);

            reverse_trimap[event.m_created_tris[i]] = face_candidate;
            face_numbers[face_candidate] = event.m_created_tris[i];
            face_matched = true;
            break;
          }
        }

        if(!face_matched)
          std::cout << "ERROR: Couldn't match the face - FLIP.\n\n\n";
      }
      
    }
    else {
      std::cout << "ERROR: unknown remeshing operation";
    }
    
    //kill the reverse map elements for the deleted faces, since no longer valid.
    for(size_t i = 0; i < event.m_deleted_tris.size(); ++i) {
      reverse_trimap[event.m_deleted_tris[i]] = FaceHandle(-1);
    }
  }
  
}


void ElasticShell::performSplit(const EdgeHandle& eh, const Vec3d& midpoint, VertexHandle& new_vert) {

  VertexHandle v0 = m_obj->fromVertex(eh);
  VertexHandle v1 = m_obj->toVertex(eh);
  assert(v0!=v1);
  Vec3d p0 = getVertexPosition(v0);
  Vec3d p1 = getVertexPosition(v1);

  std::vector<FaceHandle> oldFaces;
  std::vector<Scalar> oldThicknesses;
  std::vector<Vec2i> oldRegions;
  for(EdgeFaceIterator efit = m_obj->ef_iter(eh); efit; ++efit) {
    FaceHandle f = *efit;
    oldFaces.push_back(f);
    oldThicknesses.push_back(getThickness(f));
    oldRegions.push_back(getFaceLabel(f));
  }

  VertexHandle v2, v3;
  getEdgeOppositeVertices(*m_obj, eh, v2, v3);

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
  Vec3d undef = s*getVertexUndeformed(v0) + (1-s)*getVertexUndeformed(v1);
  setVertexVelocity(v_new, velocity);
  setUndeformedVertexPosition(v_new, undef);

  //initially set to the exact midpoint for simplicity
  setVertexPosition(v_new, simple_midpoint);

  //set consistent volumes and thickness for new faces
  assert(oldFaces.size() == newFaces.size()/2);
  
  //Old way, ignores vertex movement, assumes simple (true) midpoint
  for(unsigned int i = 0; i < oldFaces.size(); ++i) {
    //copy thicknesses exactly
    m_thicknesses[newFaces[i*2]] = oldThicknesses[i];
    m_thicknesses[newFaces[i*2+1]] = oldThicknesses[i];
    
    //compute corresponding volumes
    m_volumes[newFaces[i*2]] = oldThicknesses[i] * getArea(newFaces[i*2], true);
    m_volumes[newFaces[i*2+1]] = oldThicknesses[i] * getArea(newFaces[i*2+1], true);
    
    //copy region labels
    m_face_regions[newFaces[i*2]] = oldRegions[i];
    m_face_regions[newFaces[i*2+1]] = oldRegions[i];
  }

  //Now update to reflect motion of the midpoint to the desired position
  setVertexPosition(v_new, midpoint);

  //Volumes haven't changed, but thicknesses have due to area change.
  for(unsigned int i = 0; i < oldFaces.size(); ++i) {
    m_thicknesses[newFaces[i*2]] = m_volumes[newFaces[i*2]] / getArea(newFaces[i*2], true);
    m_thicknesses[newFaces[i*2+1]] = m_volumes[newFaces[i*2+1]] / getArea(newFaces[i*2+1], true);
  }

  VertexFaceIterator vf_iter = m_obj->vf_iter(v_new);
  for(;vf_iter; ++vf_iter)
    setFaceActive(*vf_iter);

  new_vert = v_new;
}

void ElasticShell::performCollapse(const EdgeHandle& eh, const VertexHandle& vert_to_remove, const VertexHandle& vert_to_keep, const Vec3d& new_position) {
  
  //determine area of collapsing faces
  EdgeFaceIterator efit = m_obj->ef_iter(eh);
  Scalar totalVolumeLoss = 0;
  for(;efit; ++efit)
    totalVolumeLoss += getVolume(*efit);

  //use the new position to determine how to lerp the vertex data.
  
  Vec3d p0 = getVertexPosition(vert_to_remove);
  Vec3d p1 = getVertexPosition(vert_to_keep);

  Vec3d dx(p1-p0);
  double m2 = dx.squaredNorm();
  Scalar s = clamp((p1-new_position).dot(dx)/m2, 0., 1.);

  Vec3d newVelocity = s*getVertexVelocity(vert_to_remove) + (1-s)*getVertexVelocity(vert_to_keep);
  Vec3d newUndef = s*getVertexUndeformed(vert_to_remove) + (1-s)*getVertexUndeformed(vert_to_keep);

  //do the collapse itself
  std::vector<EdgeHandle> deletedEdges;
  VertexHandle result = m_obj->collapseEdge(eh, vert_to_remove, deletedEdges);
  if(!result.isValid())
    std::cout << "Refused to perform edge collapse!\n";

  //determine 
  setVertexPosition(vert_to_keep, new_position);

  setVertexVelocity(vert_to_keep, newVelocity);
  setUndeformedVertexPosition(vert_to_keep, newUndef);

  //increment the thickness of the nearby faces to account for the lost volume

  //sum up the total area increases in the faces
  VertexFaceIterator vfit = m_obj->vf_iter(vert_to_keep);
  Scalar totalNewArea = 0;
  for(;vfit; ++vfit) {
    FaceHandle fh = *vfit;
    totalNewArea += std::max(0.0, getArea(fh, true) - m_volumes[fh] / m_thicknesses[fh]); //only consider increases
  }

  //add the volume losses in the surrounding faces
  vfit = m_obj->vf_iter(vert_to_keep);
  for(;vfit; ++vfit) {
    FaceHandle fh = *vfit;
    Scalar areaLoss = m_volumes[fh] / m_thicknesses[fh] - getArea(fh, true);
    if(areaLoss > 0) {
      totalVolumeLoss += areaLoss * m_thicknesses[fh];
    }
  }

  vfit = m_obj->vf_iter(vert_to_keep);
  for(;vfit; ++vfit) {
    FaceHandle fh = *vfit;
    Scalar newArea = getArea(fh) - m_volumes[fh] / m_thicknesses[fh];
    if(newArea <= 0) { //face lost area
      //keep the old thickness, assign the new volume according to the adjusted smaller area
      m_volumes[fh] = getArea(fh)*m_thicknesses[fh];
    }
    else { //face gained area
      //distribute this expanding face some of the volume
      m_volumes[fh] += (newArea / totalNewArea) * totalVolumeLoss; //allocate extra volume proportional to area increases
      m_thicknesses[fh] = m_volumes[fh] / getArea(fh);
    }
  }
}


bool ElasticShell::performFlip(const EdgeHandle& eh, const FaceHandle f0, const FaceHandle& f1, EdgeHandle& newEdge) {
      
    //determine volume of the region being flipped
    Scalar totalVolume = m_volumes[f0] + m_volumes[f1];
    
    Vec2i oldLabels = m_face_regions[f0]; //assume f0 and f1 have the same labels, because if not this flip is nonsense

    newEdge = flipEdge2(*m_obj, eh, f0, f1);
    if(!newEdge.isValid()) {//couldn't flip the edge, such an edge already existed
      std::cout << "Edge flip failed for some reason...\n";
      return false;
    }

    FaceHandle f0new, f1new;
    getEdgeFacePair(*m_obj, newEdge, f0new, f1new);
    setFaceActive(f0new);
    setFaceActive(f1new);

    //assign new thicknesses
    Scalar f0newArea = getArea(f0new, true);
    Scalar f1newArea = getArea(f1new, true);
    Scalar newThickness = totalVolume / (f0newArea + f1newArea);
    m_thicknesses[f0new] = newThickness;
    m_thicknesses[f1new] = newThickness;
    m_volumes[f0new] = f0newArea*newThickness;
    m_volumes[f1new] = f1newArea*newThickness;
    
    //keep previous labels
    m_face_regions[f0new] = oldLabels;
    m_face_regions[f1new] = oldLabels;
    
    return true;
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
    Scalar newThickness = m_volumes[f] / area;
    m_thicknesses[f] = newThickness;
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
      setUndeformedVertexPosition(vertices[i], m_inflow_positions[boundary][i]);
      setVertexVelocity(vertices[i], m_inflow_velocities[boundary][i]);
      m_vertex_masses[vertices[i]] = 0;
      m_obj->setVertexDampingUndeformedPosition(vertices[i], m_inflow_positions[boundary][i]);

      getDefoObj().constrainVertex(vertices[i], new FixedVelocityConstraint(m_inflow_positions[boundary][i], m_inflow_velocities[boundary][i], current_time));
    }


    for(unsigned int i = 0; i < faces.size(); ++i) {
      m_thicknesses[faces[i]] = m_inflow_thickness;
      m_volumes[faces[i]] = 0;
      m_volumes[faces[i]] = getArea(faces[i])*m_inflow_thickness;
      FaceVertexIterator fvit = m_obj->fv_iter(faces[i]);
      for(;fvit;++fvit) {
        VertexHandle vh = *fvit;
        m_vertex_masses[vh] += m_volumes[faces[i]] * m_density / 3.0;
      }
    }

  }
  
  computeMasses();

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
