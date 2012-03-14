#include "BASim/src/Collisions/ElTopo/broadphasegrid.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShell.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"
#include "BASim/src/Collisions/ElTopo/ccd_wrapper.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/CSTMembraneForce.hh"
#include "BASim/src/Math/Math.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellVertexTriSpringForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellVertexPointSpringForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ShellStickyRepulsionForce.hh"
#include "BASim/src/Collisions/ElTopo/collisionqueries.hh"

#include "BASim/src/Collisions/ElTopo/array3.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"

#include "surftrack.h"


#include <algorithm>
#include <numeric>

namespace BASim {

ElasticShell::ElasticShell(DeformableObject* object, const FaceProperty<char>& shellFaces, Scalar timestep) : 
  PhysicalModel(*object), m_obj(object), 
    m_active_faces(shellFaces), 
    m_undeformed_positions(object), 
    m_undef_xi(object),
    m_damping_undeformed_positions(object), 
    m_damping_undef_xi(object),
    m_vertex_masses(object),
    m_edge_masses(object),
    m_thicknesses(object),
    m_volumes(object),
    m_positions(object), 
    m_xi(object), 
    m_velocities(object),
    m_xi_vel(object),
    m_density(1),
    m_proximity_epsilon(1e-5),
    m_vert_point_springs(NULL),
    m_repulsion_springs(NULL),
    m_sphere_collisions(false),
    m_object_collisions(false),
    m_ground_collisions(false)
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

void ElasticShell::setVertexUndeformed( const VertexProperty<Vec3d>& undef )
{
  m_undeformed_positions = undef;
}

void ElasticShell::setEdgeUndeformed( const EdgeProperty<Scalar>& undef )
{
  m_undef_xi = undef;
}

void ElasticShell::setVertexPositions( const VertexProperty<Vec3d>& positions )
{
  m_positions = positions;
}

void ElasticShell::setEdgeXis( const EdgeProperty<Scalar>& positions )
{
  m_xi = positions;
}

void ElasticShell::setVertexVelocities( const VertexProperty<Vec3d>& velocities) 
{
  m_velocities = velocities;
}

void ElasticShell::setEdgeVelocities(const EdgeProperty<Scalar>& velocities)
{
  m_xi_vel = velocities;
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
  if(current) {
    Vec3d v0 = m_positions[v1_hnd] - m_positions[v0_hnd];
    Vec3d v1 = m_positions[v2_hnd] - m_positions[v0_hnd];
    Vec3d triVec = v0.cross(v1);
    return 0.5*triVec.norm();
  }
  else {
    Vec3d v0 = m_undeformed_positions[v1_hnd] - m_undeformed_positions[v0_hnd];
    Vec3d v1 = m_undeformed_positions[v2_hnd] - m_undeformed_positions[v0_hnd];
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
      Vec3d v0 = m_positions[v1_hnd] - m_positions[v0_hnd];
      Vec3d v1 = m_positions[v2_hnd] - m_positions[v0_hnd];
      Vec3d triVec = v0.cross(v1);
      Scalar area = 0.5*sqrt(triVec.dot(triVec)) / 3.0;
      Scalar contribution = m_thicknesses[f_hnd] * m_density * area;
      
      //accumulate mass to the vertices
      m_vertex_masses[v0_hnd] += contribution;
      m_vertex_masses[v1_hnd] += contribution;
      m_vertex_masses[v2_hnd] += contribution;

      //also accumulate mass to the edges (this mass computation is probably not consistent with what we want)
      FaceEdgeIterator feit = m_obj->fe_iter(f_hnd);
      EdgeHandle e0_hnd = *feit; ++feit; assert(feit);
      EdgeHandle e1_hnd = *feit; ++feit; assert(feit);
      EdgeHandle e2_hnd = *feit; ++feit; //assert(feit);

      m_edge_masses[e0_hnd] += contribution;
      m_edge_masses[e1_hnd] += contribution;
      m_edge_masses[e2_hnd] += contribution;

      //store the current volumes
      m_volumes[f_hnd] = 3*area*m_thicknesses[f_hnd];
    }
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
  //they're all vertex Dofs for a shell
  assert(hnd.getType() == DofHandle::VERTEX_DOF || hnd.getType() == DofHandle::EDGE_DOF);

  //return reference to the appropriate position in the vector
  if(hnd.getType() == DofHandle::VERTEX_DOF) {
    const VertexHandle& vh = static_cast<const VertexHandle&>(hnd.getHandle());
    return const_cast<Vec3d&>(m_positions[vh])[hnd.getNum()];
  }
  else {
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    return const_cast<Scalar&>(m_xi[eh]);
  }
}

void ElasticShell::setDof( const DofHandle& hnd, const Scalar& dof )
{
  //they're all vertex Dofs for a shell
  assert(hnd.getType() == DofHandle::VERTEX_DOF || hnd.getType() == DofHandle::EDGE_DOF);

  if(hnd.getType() == DofHandle::VERTEX_DOF) {
    const VertexHandle& vh = static_cast<const VertexHandle&>(hnd.getHandle());
    m_positions[vh][hnd.getNum()] = dof;
  }
  else {
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    m_xi[eh] = dof;
  }
}

const Scalar& ElasticShell::getVel( const DofHandle& hnd ) const
{
  assert(hnd.getType() == DofHandle::VERTEX_DOF || hnd.getType() == DofHandle::EDGE_DOF);

  //return reference to the appropriate position in the vector
  if(hnd.getType() == DofHandle::VERTEX_DOF) {
    const VertexHandle& vh = static_cast<const VertexHandle&>(hnd.getHandle());
    return const_cast<Vec3d&>(m_velocities[vh])[hnd.getNum()];
  }
  else{
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    return const_cast<Scalar&>(m_xi_vel[eh]);
  }
}

void ElasticShell::setVel( const DofHandle& hnd, const Scalar& vel )
{
  assert(hnd.getType() == DofHandle::VERTEX_DOF || hnd.getType() == DofHandle::EDGE_DOF);

  if(hnd.getType() == DofHandle::VERTEX_DOF) {
    const VertexHandle& vh = static_cast<const VertexHandle&>(hnd.getHandle());
    m_velocities[vh][hnd.getNum()] = vel;
  }
  else{
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    m_xi_vel[eh] = vel;
  }
}

const Scalar& ElasticShell::getMass( const DofHandle& hnd ) const
{
  assert(hnd.getType() == DofHandle::VERTEX_DOF || hnd.getType() == DofHandle::EDGE_DOF);

  if(hnd.getType() == DofHandle::VERTEX_DOF) {
    const VertexHandle& vh = static_cast<const VertexHandle&>(hnd.getHandle());
    return m_vertex_masses[vh];
  }
  else {
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    return m_edge_masses[eh];
  }
}

void ElasticShell::getScriptedDofs( IntArray& dofIndices, std::vector<Scalar>& dofValues, Scalar time ) const
{
  for(unsigned int i = 0; i < m_constrained_vertices.size(); ++i) {
    
    int dofBase = getVertexDofBase(m_constrained_vertices[i]);
    Vec3d pos = m_constraint_positions[i]->operator()(time);
    dofIndices.push_back(dofBase); dofValues.push_back(pos[0]);
    dofIndices.push_back(dofBase+1); dofValues.push_back(pos[1]);
    dofIndices.push_back(dofBase+2); dofValues.push_back(pos[2]);
  }
}

void ElasticShell::constrainVertex( const VertexHandle& v, const Vec3d& pos )
{
  m_constrained_vertices.push_back(v);
  PositionConstraint* c = new FixedPositionConstraint(pos);
  m_constraint_positions.push_back(c);
}

void ElasticShell::constrainVertex( const VertexHandle& v, PositionConstraint* c )
{
  m_constrained_vertices.push_back(v);
  m_constraint_positions.push_back(c);
}

void ElasticShell::releaseVertex( const VertexHandle& v)
{
 
  bool deletedVertex = true;
  while(deletedVertex) {
    deletedVertex = false;
    //can only have one constraint or things get broken anyways(right?), so no need to search for multiple
    int index = -1;
    for(unsigned int i = 0; i < m_constraint_positions.size(); ++i) {
      if(m_constrained_vertices[i] == v) {
        index = i;
        deletedVertex = true;
        break;
      }
    }

    //remove the constraint
    if(index != -1) {
      delete m_constraint_positions[index];
      m_constraint_positions.erase(m_constraint_positions.begin()+index);
      m_constrained_vertices.erase(m_constrained_vertices.begin()+index);
    }
  }
  
  
}


void ElasticShell::startStep(Scalar time, Scalar timestep)
{
  std::cout << "Starting startStep\n";
  /* //Debugging forces
   const std::vector<ElasticShellForce*>& forces = getForces();
   std::vector<ElasticShellForce*>::const_iterator fIt;

   
   for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
      if(typeid(*(*fIt)) == typeid(CSTMembraneForce)) {
         Scalar potEnergy = (*fIt)->globalEnergy();
         //std::cout << "Energy " << potEnergy;
         VecXd curr_force(m_obj->nv()*3);
         curr_force.setZero();
         (*fIt)->setDebug(true);
         (*fIt)->globalForce(curr_force);
         (*fIt)->setDebug(false);

         ////Now sum forces on all the vertices along the circumference
         Vec3d forceSum;
         forceSum.setZero();
         for(VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit) {
            Vec3d vertPos = getVertexPosition(*vit);
            if(vertPos[2] < -1e-5 || vertPos[2] > 1e-5) continue; //skip non-circumferential vertices

            int startIndex = getVertexDofBase(*vit);
            Vec3d vertForce(curr_force[startIndex], curr_force[startIndex+1], curr_force[startIndex+2]);
            forceSum += vertForce;
         }
         std::cout << "Force sum is: " << forceSum << std::endl;
      }
   }
   */



  //update the damping "reference configuration" for computing viscous forces.
  m_damping_undeformed_positions = m_positions;
  m_damping_undef_xi = m_xi;

  //tell the forces to update anything they need to update
  const std::vector<ElasticShellForce*>& forces = getForces();
  for(unsigned int i = 0; i < forces.size(); ++i) {
    forces[i]->update();
  }
  std::cout << "Done startStep\n";
}

void ElasticShell::resolveCollisions(Scalar timestep) {

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
    Vec3d old_vert = m_damping_undeformed_positions[vh];
    Scalar mass = getMass(vh);

    vert_new.push_back(ElTopo::Vec3d(vert[0], vert[1], vert[2]));
    vert_old.push_back(ElTopo::Vec3d(old_vert[0], old_vert[1], old_vert[2]));
    if(isConstrained(vh)) {
      masses.push_back(numeric_limits<Scalar>::infinity());
    }
    else {
      masses.push_back(mass);
    }
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


  // build a DynamicSurface
  Scalar friction_coeff = 0;
  ElTopo::DynamicSurface dynamic_surface( vert_old, tri_data, masses, m_proximity_epsilon, friction_coeff, true, false );

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

  m_broad_phase.update_broad_phase_static(*m_obj, m_positions, collision_distance);
  
  //determine proximity of vertex triangle pairs and
  //add damped springs between them to handle new collisions
  
  //consider all vertices
  
  for(VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit) {
    VertexHandle vh = *vit;
    Vec3d vert_pos = m_positions[vh];
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
      int fv = 0;
      bool goodSpring = true;
      for(FaceVertexIterator fvit = m_obj->fv_iter(f); fvit; ++fvit) {
        face_verts[fv] = ElTopoCode::toElTopo(m_positions[*fvit]);
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
    start.push_back(m_positions[verts[i]]);
    
    FaceVertexIterator fvit = m_obj->fv_iter(faces[i]);
    Vec3d v0 = m_positions[*fvit];++fvit;
    Vec3d v1 = m_positions[*fvit];++fvit;
    Vec3d v2 = m_positions[*fvit];
    Vec3d result = bary[i][0] * v0 + bary[i][1] * v1 + bary[i][2] * v2;
    end.push_back(result);
  }

}

void ElasticShell::setCollisionParams(Scalar proximity, Scalar stiffness, Scalar damping) {
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


  //El Topo collision processing.
  std::cout << "Resolving collisions\n";
  resolveCollisions(timestep);
  std::cout << "Finished resolving collisions.\n";

  //Ground plane penalty force.
  if(m_ground_collisions) {

    //Hard constraints
    for(VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit) {
      Vec3d curPos = getVertexPosition(*(vit));
      if(curPos[1] < m_ground_height) {
        if(!isConstrained(*vit)) {
          //constrainVertex(*vit, curPos);

          //Sinking
          //constrainVertex(*vit, new FixedVelocityConstraint(curPos, Vec3d(0, m_ground_velocity, 0), time));

          //Conveying
          constrainVertex(*vit, new FixedVelocityConstraint(curPos, Vec3d(0, 0, m_ground_velocity), time));
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
        if(!isConstrained(*vit)) {

          //Sinking
          //constrainVertex(*vit, new FixedVelocityConstraint(curPos, Vec3d(0, m_ground_velocity, 0), time));

          //Conveying
          constrainVertex(*vit, new FixedVelocityConstraint(curPos, m_sphere_velocity, time));
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
        if(!isConstrained(*vit)) {

          //Fixed position
          constrainVertex(*vit, curPos);
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
    //for(int i = 0; i < m_remeshing_iters; ++i)
    //  remesh(m_remesh_edge_length);  
    remesh_new();
    std::cout << "Completed remeshing\n";

    //Relabel DOFs if necessary
    do_relabel = true;
  }

  static int only_twice = 0;

  //Tearing processing
  if(m_tearing){
    std::cout << "Processing tearing. \n";
    //fracture();
    fracture_new();
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


void ElasticShell::fracture_new() {

  //Set up a SurfTrack, run remeshing, render the new mesh
  ElTopo::SurfTrackInitializationParameters construction_parameters;
  construction_parameters.m_proximity_epsilon = m_proximity_epsilon;
  construction_parameters.m_allow_vertex_movement = false;
  construction_parameters.m_min_edge_length = 0.5*m_remesh_edge_length;
  construction_parameters.m_max_edge_length = 1.5*m_remesh_edge_length;
  construction_parameters.m_max_volume_change = numeric_limits<double>::max();   
  construction_parameters.m_min_triangle_angle = 5;
  construction_parameters.m_max_triangle_angle = 175;
  construction_parameters.m_verbose = false;
  construction_parameters.m_use_curvature_when_collapsing = false;
  construction_parameters.m_use_curvature_when_splitting = false;
  construction_parameters.m_allow_non_manifold = false;
  construction_parameters.m_collision_safety = true;
  //TODO If we try using other subdivision schemes for splitting edges
  //note that the manner in which the edge split operation is performed
  //will need to be adjusted to properly preserve total volume.

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
    if(isConstrained(vh))
      masses.push_back(numeric_limits<Scalar>::infinity());
    else
      masses.push_back(mass);
    vert_numbers[vh] = id;
    reverse_vertmap.push_back(vh);

    Vec3d pos = getVertexPosition(vh);

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

  std::cout << "Collecting desired cuts\n";
  std::vector< std::pair<size_t,size_t> > edges_to_cut;
  for(EdgeIterator it = mesh.edges_begin(); it != mesh.edges_end(); ++it) {
    EdgeHandle eh = *it;
    VertexHandle vh0 = mesh.fromVertex(eh);
    VertexHandle vh1 = mesh.toVertex(eh);
    Vec3d pos0 = getVertexPosition(vh0);
    Vec3d pos1 = getVertexPosition(vh1);
    if(mesh.isBoundary(eh)) continue;
    if(pos0[1] > 2.7 && pos0[1] < 3.3 && pos1[1] > 2.7 && pos1[1] < 3.3) {
      edges_to_cut.push_back(make_pair(vert_numbers[vh0],vert_numbers[vh1]));
    }
  }

  std::cout << "Doing cutting with El Topo\n";
  surface_tracker.cut_mesh(edges_to_cut);

  
  std::cout << "Performing " << surface_tracker.m_mesh_change_history.size() << " Cutting Operations:\n";
  for(unsigned int j = 0; j < surface_tracker.m_mesh_change_history.size(); ++j) {
    ElTopo::MeshUpdateEvent event = surface_tracker.m_mesh_change_history[j];

    VertexHandle v0 = reverse_vertmap[event.m_v0];
    VertexHandle v1 = reverse_vertmap[event.m_v1];
    EdgeHandle eh = findEdge(mesh, v0, v1);

    assert(event.m_type == ElTopo::MeshUpdateEvent::EDGE_CUT);
      
    //Identify the edge based on its endpoint vertices instead
    std::cout << "Cut\n";

    //based on the type of cut, perform the appropriate operation
    bool aBound = m_obj->isBoundary(v0);
    bool bBound = m_obj->isBoundary(v1);
    
    std::vector<VertexHandle> srcVerts;
    if(aBound) srcVerts.push_back(v0);
    if(bBound) srcVerts.push_back(v1);

    std::vector<FaceHandle> newFaces;
    std::vector<FaceHandle> deletedFaces;
    std::vector<EdgeHandle> deleteEdges;
    std::vector<VertexHandle> newVerts;

    //Do exactly the same cut as El Topo
    //if(aBound && ! bBound || !aBound && bBound) { //one vertex boundary
      
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

      FaceHandle newFaceHandle = mesh.addFace(v0, v1, v2); //build the face. do we need to correct the ordering?
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

void ElasticShell::fracture(){
    //Here comes all the fracture code


    //Figure out which edges will be taken down
    std::vector<EdgeHandle> edgesToFrac;
    getDesiredFractures(edgesToFrac);

    std::cout << "Fracturing " << edgesToFrac.size() << " edges" << std::endl;

//    if ( edgesToFrac.size() != 0 && m_obj->edgeExists(edgesToFrac[0]) && shouldFracture(edgesToFrac[0]))
//        performTearing(edgesToFrac[0]);
//    All of these are fractured interior, now figure out which will become boundaries
    for(unsigned int i = 0; i < edgesToFrac.size(); ++i){
//    for(unsigned int i = 0; i < 1 && i < edgesToFrac.size(); ++i){
        //skip constraints

//        if (fromBound[i] || toBound[i]){

            //Triple check that it will be ok to tear
            if ( m_obj->edgeExists(edgesToFrac[i]) && shouldFracture(edgesToFrac[i])){
                performTearing(edgesToFrac[i]);
            }

//        }
    }
//    for ( int i = 0 ; i < 5; ++i){
//        for(EdgeIterator eit = m_obj->edges_begin(); eit != m_obj->edges_end(); ++eit){
//
//        }
//    }
}





void ElasticShell::performTearing(const EdgeHandle & eh){

#ifndef NDEBUG
            int facesBef = m_obj->nf();
            int edgesBef = m_obj->ne();
            int vertsBef = m_obj->nv();
#endif

    //Check for self-intersections being induced...
    VertexHandle v0 = m_obj->fromVertex(eh);
    VertexHandle v1 = m_obj->toVertex(eh);
    bool aBound = m_obj->isBoundary(v0);
    bool bBound = m_obj->isBoundary(v1);
    assert(v0!=v1);

    std::vector<FaceHandle> oldFaces;
    std::vector<EdgeHandle> oldEdges;
    std::vector<FaceHandle> newFaces;
    VertexHandle newVerta, newVertb;

//    ElTopo::Vec3d midpoint_ET = ElTopo::toElTopo(midpoint);
    //if(edgeSplitCausesCollision(midpoint_ET, midpoint_ET, eh))
    //  return false;

//    //remove the edge and surrounding faces from the collision structure
//    m_broad_phase.remove_edge(eh.idx());
//    for(unsigned int i = 0; i < oldFaces.size(); ++i)
//      m_broad_phase.remove_triangle(oldFaces[i].idx());

    //perform the actual split
    if ( aBound && bBound){
        std::cout << "\tSeparating geom" << std::endl;
        tearEdge(*m_obj, eh, v0, v1, newVerta, newVertb, newFaces, oldFaces, oldEdges);
        setVertexVelocity(newVerta, getVertexVelocity(v0));
        setVertexPosition(newVerta, getVertexPosition(v0));
        setUndeformedVertexPosition(newVerta, getVertexUndeformed(v0));

        setVertexVelocity(newVertb, getVertexVelocity(v1));
        setVertexPosition(newVertb, getVertexPosition(v1));
        setUndeformedVertexPosition(newVertb, getVertexUndeformed(v1));



    } else if ( aBound && !bBound){
        std::cout << "\tSeparating from" << std::endl;
        tearVertexAlong(*m_obj, eh, v0, newVerta, newFaces, oldFaces, oldEdges);
        setVertexVelocity(newVerta, getVertexVelocity(v0));
        setVertexPosition(newVerta, getVertexPosition(v0));
        setUndeformedVertexPosition(newVerta, getVertexUndeformed(v0));

    } else if ( !aBound && bBound){
        std::cout << "\tSeparating to" << std::endl;
        tearVertexAlong(*m_obj, eh, v1, newVertb, newFaces, oldFaces, oldEdges);
        setVertexVelocity(newVertb, getVertexVelocity(v1));
        setVertexPosition(newVertb, getVertexPosition(v1));
        setUndeformedVertexPosition(newVertb, getVertexUndeformed(v1));
    } else{
        std::cout << "\t******************************Updated bounds make this edge non-fracturable" << std::endl;
        
        VHList newVerts;
        tearInteriorEdge(*m_obj, eh, v0, v1, newVerts, newFaces, oldFaces, oldEdges);

        assert(newVerts.size() == 4);
        assert(oldFaces.size() == 2);
        assert(newFaces.size() == 4);
        //Copy the values to the new vertices
        Vec3d newVel = 0.5f * (getVertexVelocity(v0) + getVertexVelocity(v1));
        Vec3d newPos = 0.5f * (getVertexPosition(v0) + getVertexPosition(v1));
        Vec3d newUnd = 0.5f * (getVertexUndeformed(v0) + getVertexUndeformed(v1));
        for ( unsigned int i = 0; i < newVerts.size(); ++i){
            setVertexVelocity(newVerts[i], newVel);
            setVertexPosition(newVerts[i], newPos);
            setUndeformedVertexPosition(newVerts[i], newUnd);
        }
        
    }

    //set consistent volumes and thickness for new faces
//    std::cout << "\tFaces added: " << newFaces.size() << std::endl;
//    std::cout << "\tEdges at this point(before deleting): " << m_obj->ne() << std::endl;
//    std::cout << "\tFaces that will be deleted: " << oldFaces.size() << std::endl;
//    std::cout << "\tEdges that will be deleted: " << oldEdges.size() << std::endl;

    if ( aBound || bBound ){//Update the thicknesses and volumes for noninterior
        assert(oldFaces.size() == newFaces.size());
        for(unsigned int i = 0; i < oldFaces.size(); ++i) {
          m_thicknesses[newFaces[i]] = m_thicknesses[oldFaces[i]];
          m_volumes[newFaces[i]] = m_volumes[oldFaces[i]];
        }
    }else {//and for interior
        assert(oldFaces.size() == 2);
        assert(newFaces.size() == 4);

        m_thicknesses[newFaces[0]] = m_thicknesses[oldFaces[0]];
        m_thicknesses[newFaces[1]] = m_thicknesses[oldFaces[0]];
        m_thicknesses[newFaces[2]] = m_thicknesses[oldFaces[1]];
        m_thicknesses[newFaces[3]] = m_thicknesses[oldFaces[1]];

        m_volumes[newFaces[0]] = m_volumes[newFaces[1]] = m_volumes[oldFaces[0]]/2.0;
        m_volumes[newFaces[2]] = m_volumes[newFaces[3]] = m_volumes[oldFaces[1]]/2.0;

    }

    //Time to delete the extra tris
    //Now delete all the extra things
    // -all faces that were added to the newverts
    // -all edges that were added to the newverts
    for (int i = 0; i < (int)oldFaces.size(); ++i){
        m_obj->deleteFace(oldFaces[i], false);
    }
    for (int i = 0; i < (int)oldEdges.size(); ++i){
        m_obj->deleteEdge(oldEdges[i], false);
    }

#ifndef NDEBUG
        //Check that the correct number of things was created
//            std::cout << "\tFaces after: " << m_obj->nf() << std::endl;
//            std::cout << "\tEdges after: " << m_obj->ne() << std::endl;
//            std::cout << "\tVerts after: " << m_obj->nv() << std::endl;
            if ( aBound && bBound){
                assert( (facesBef - m_obj->nf()) == 0);
                assert( (edgesBef - m_obj->ne()) == -1);
                assert( (vertsBef - m_obj->nv()) == -2);
            }
            else if ( aBound || bBound ){
                assert( (facesBef - m_obj->nf()) == 0);
                assert( (edgesBef - m_obj->ne()) == -1);
                assert( (vertsBef - m_obj->nv()) == -1);
            }
            else {
                //Check that the correct number of things was created
                assert( (facesBef - m_obj->nf()) == -2);
                assert( (edgesBef - m_obj->ne()) == -7);
                assert( (vertsBef - m_obj->nv()) == -4);
            }
#endif
//    std::cout << "\tEdges at this point(after deleting): " << m_obj->ne() << std::endl;

    for(int i = 0; i < (int) newFaces.size(); ++i){
        setFaceActive(newFaces[i]);
    }

    //update collision data structures

//    //add vertices
//    m_broad_phase.add_vertex(v_new.idx(), ElTopo::toElTopo(midpoint), m_proximity_epsilon);
//
//    //add edges
//    VertexEdgeIterator ve_iter = m_obj->ve_iter(v_new);
//    for(; ve_iter; ++ve_iter) {
//      std::vector<ElTopo::Vec3d> edge_verts;
//      EdgeVertexIterator ev_iter = m_obj->ev_iter(*ve_iter);
//      for(; ev_iter; ++ev_iter)
//        edge_verts.push_back(ElTopo::toElTopo(getVertexPosition(*ev_iter)));
//      m_broad_phase.add_edge((*ve_iter).idx(), edge_verts[0], edge_verts[1], m_proximity_epsilon);
//    }
//
//    //add tris
//    vf_iter = m_obj->vf_iter(v_new);
//    for(;vf_iter; ++vf_iter) {
//      std::vector<ElTopo::Vec3d> tri_verts;
//      FaceVertexIterator fv_iter = m_obj->fv_iter(*vf_iter);
//      for(; fv_iter; ++fv_iter)
//        tri_verts.push_back(ElTopo::toElTopo(getVertexPosition(*fv_iter)));
//      m_broad_phase.add_triangle((*vf_iter).idx(), tri_verts[0], tri_verts[1], tri_verts[2], m_proximity_epsilon);
//    }

//    new_vert = v_new;

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
void ElasticShell::getDesiredFractures(std::vector<EdgeHandle> & edges ){
    for ( EdgeIterator eit = m_obj->edges_begin(); eit != m_obj->edges_end(); ++eit){
        //Now look if it exceeds the threshold
        if ( shouldFracture(*eit) ){
           edges.push_back(*eit);
        }
    }
}
bool ElasticShell::shouldFracture (const EdgeHandle & eh) const{
    //Ignore inflow edges
    if(isInflow(eh)) return false;
    //Ignore constrain edges
    if ( isConstrained(m_obj->fromVertex(eh)) || isConstrained(m_obj->toVertex(eh))) return false;

    //Ignore boundary edges
    if ( m_obj->isBoundary(eh) ) return false;

    //Get the average thickness weighted by areas for each edge
    Scalar thickness = 0.0;
    Scalar totalA = 0.0;
    for ( EdgeFaceIterator efit = m_obj->ef_iter(eh); efit; ++efit){
        Scalar w = getArea(*efit);
        thickness += w * getThickness(*efit);
        totalA += w;
    }
    thickness /= totalA;
    Scalar p = (Scalar) rand() / (Scalar) RAND_MAX;
//    return (thickness < m_tear_thres) && (m_obj->isBoundary(m_obj->fromVertex(eh)) || m_obj->isBoundary(m_obj->toVertex(eh)))
//            && ( p <  m_tear_rand);
    return (thickness < m_tear_thres)  && ( p <  m_tear_rand);
}


void ElasticShell::remesh_new()
{

  //Set up a SurfTrack, run remeshing, render the new mesh
  ElTopo::SurfTrackInitializationParameters construction_parameters;
  construction_parameters.m_proximity_epsilon = m_proximity_epsilon;
  construction_parameters.m_allow_vertex_movement = false;
  construction_parameters.m_min_edge_length = 0.5*m_remesh_edge_length;
  construction_parameters.m_max_edge_length = 1.5*m_remesh_edge_length;
  construction_parameters.m_max_volume_change = numeric_limits<double>::max();   
  construction_parameters.m_min_triangle_angle = 5;
  construction_parameters.m_max_triangle_angle = 175;
  construction_parameters.m_verbose = false;
  construction_parameters.m_use_curvature_when_collapsing = false;
  construction_parameters.m_use_curvature_when_splitting = false;
  construction_parameters.m_allow_non_manifold = false;
  construction_parameters.m_collision_safety = true;
  //TODO If we try using other subdivision schemes for splitting edges
  //note that the manner in which the edge split operation is performed
  //will need to be adjusted to properly preserve total volume.

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
    if(isConstrained(vh))
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

  std::cout << "Calling surface improvement\n";
  ElTopo::SurfTrack surface_tracker( vert_data, tri_data, masses, construction_parameters ); 

  surface_tracker.improve_mesh();

  std::cout << "Performing " << surface_tracker.m_mesh_change_history.size() << " Improvement Operations:\n";
  for(unsigned int j = 0; j < surface_tracker.m_mesh_change_history.size(); ++j) {
    ElTopo::MeshUpdateEvent event = surface_tracker.m_mesh_change_history[j];
    
    VertexHandle v0 = reverse_vertmap[event.m_v0];
    VertexHandle v1 = reverse_vertmap[event.m_v1];
    EdgeHandle eh = findEdge(mesh, v0, v1);

    if(event.m_type == ElTopo::MeshUpdateEvent::EDGE_COLLAPSE) {
      std::cout << "Collapse\n";
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

        if(!face_matched)
          std::cout << "ERROR: Couldn't match the face.\n\n\n";
      }
    }
    else if(event.m_type == ElTopo::MeshUpdateEvent::EDGE_SPLIT) {
      //Identify the edge based on its endpoint vertices instead
      std::cout << "Split\n";

      VertexHandle new_vert;
      Vec3d new_pos(event.m_vert_position[0], event.m_vert_position[1], event.m_vert_position[2]);
      
      //Do the same split as El Topo
      performSplitET(eh, new_pos, new_vert);
      
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
          std::cout << "ERROR: Couldn't match the face.\n\n\n";
      }
    }
    else if(event.m_type == ElTopo::MeshUpdateEvent::EDGE_FLIP) {
      
      std::cout << "Flip\n";
      //Do the same flip as El Topo
      EdgeHandle newEdge;
      performFlip(eh, newEdge);
      
      
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
          std::cout << "ERROR: Couldn't match the face.\n\n\n";
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


void ElasticShell::remesh( Scalar desiredEdge )
{

  //m_broad_phase.update_broad_phase_static(*m_obj, m_positions, m_collision_proximity);
  
  //Parameters adapted from Jiao et al. "Anisotropic Mesh Adaptation for Evolving Triangulated Surfaces"
  Scalar ratio_R = 0.5;
  Scalar ratio_r = 0.1;
  
  Scalar minEdge = ratio_R*desiredEdge; //L in the jiao paper
  Scalar maxEdge = 1.5*desiredEdge; //S in the jiao paper

  Scalar minAngle = 15.0*M_PI/180.0;
  Scalar maxAngle = 160.0*M_PI/180.0;
  
  flipEdges();
  splitEdges(desiredEdge, maxEdge, maxAngle);
  collapseEdges(minAngle, desiredEdge, ratio_R, ratio_r, minEdge);
  
}

void ElasticShell::flipEdges() {

  EdgeIterator e_it = m_obj->edges_begin();
  for(;e_it != m_obj->edges_end(); ++e_it) {
    EdgeHandle eh = *e_it;

    bool isInflow = this->isInflow(eh);

    if(isInflow) continue;

    //if either of the faces has a spring stuck to it, don't flip!
    bool springAttached = false;
    for(EdgeFaceIterator efit = m_obj->ef_iter(eh); efit; ++efit) {
      FaceHandle fh = *efit;
      if(m_repulsion_springs->isFaceInUse(fh))
        springAttached = true;    
    }
    if(springAttached) continue;

    VertexHandle v0 = m_obj->fromVertex(eh);
    VertexHandle v1 = m_obj->toVertex(eh);
    assert(v0!=v1);
    Vec3d p0 = getVertexPosition(v0);
    Vec3d p1 = getVertexPosition(v1);
    VertexHandle v2, v3;
    Scalar edgeLength = (p1-p0).norm();
    bool success = getEdgeOppositeVertices(*m_obj, eh, v2, v3);
    if(!success) continue;


    Vec3d p2 = getVertexPosition(v2);
    Vec3d p3 = getVertexPosition(v3);
    Scalar oppEdgeLength = (p3-p2).norm();

    Vec3d dir0 = p0-p2, dir1 = p1 - p2;
    Scalar angle0 = acos(dir0.dot(dir1));

    dir0 = p0-p3; dir1 = p1 - p3;
    Scalar angle1 = acos(dir0.dot(dir1));
    
    if(angle0 + angle1 > M_PI) {

      //check area of proposed tris
      Scalar area0 = (p2-p0).cross(p3-p0).norm();
      if(area0 < 1e-12) continue;
      Scalar area1 = (p2-p1).cross(p3-p1).norm();
      if(area1 < 1e-12) continue;
      EdgeHandle newEdge;
      performFlip(eh, newEdge);
    }
  }

}



struct SortableEdge
{
  EdgeHandle edge_hnd;
  double edge_length;

  SortableEdge( EdgeHandle eh, double el ) : edge_hnd(eh), edge_length(el) {}

  bool operator<( const SortableEdge& other ) const
  {
    return (this->edge_length < other.edge_length);
  }
};
bool ElasticShell::splitEdges( double desiredEdge, double maxEdge, double maxAngle) {

  //sort edges in order of length, so we split long ones first
  std::vector<SortableEdge> sortable_edges_to_try;
  
  EdgeIterator e_it = m_obj->edges_begin();
  for(;e_it != m_obj->edges_end(); ++e_it) {   
    
    EdgeHandle eh = *e_it;

    VertexHandle vertex_a = m_obj->fromVertex(eh);
    VertexHandle vertex_b = m_obj->toVertex(eh);

    if ( !m_obj->edgeExists(eh) || !isEdgeActive(eh))   { continue; }     // skip inactive/non-existent edges
    
    //skip edges that are flowing into the domain
    bool isConstrained = isInflow(eh);

    //don't split constrained edges. (alternatively, we could split them, and add the new vert to the constrained list.)
    bool aConstrained = this->isConstrained(vertex_a);
    bool bConstrained = this->isConstrained(vertex_b);


//    //don't split constrained edges. (alternatively, we could split them, and add the new vert to the constrained list somehow.)
//    bool aConstrained = false, bConstrained = false;
//
//    for(unsigned int i = 0; i < m_constrained_vertices.size(); ++i) {
//      if(m_constrained_vertices[i] == vertex_a)
//        aConstrained = true;
//      if(m_constrained_vertices[i] == vertex_b)
//        bConstrained = true;
//    }

    if(aConstrained && bConstrained || isConstrained) continue;

    //don't split faces that have springs attached (for now).
    bool facesHaveSprings = false;
    for(EdgeFaceIterator efit = m_obj->ef_iter(eh); efit; ++efit) {
      FaceHandle fh = *efit;
      if(m_repulsion_springs->isFaceInUse(fh)) {
        facesHaveSprings = true;
        break;
      }
    }
    if(facesHaveSprings) continue;
    
    assert( m_obj->vertexExists(vertex_a) );
    assert( m_obj->vertexExists(vertex_b) );

    double current_length = (m_positions[ vertex_a ] - m_positions[ vertex_b ]).norm();

    if ( current_length > maxEdge ) {
      sortable_edges_to_try.push_back( SortableEdge( eh, current_length ) );
    }
  }

  
  //sort into ascending order
  std::sort( sortable_edges_to_try.begin(), sortable_edges_to_try.end() );
  bool anySplits = false;

  //iterate from biggest to smallest (end to start)
  for(int i = sortable_edges_to_try.size()-1; i >= 0; --i) {

    EdgeHandle eh = sortable_edges_to_try[i].edge_hnd;
    
    bool isConstrained = false;
    for(unsigned int i = 0; i < m_inflow_boundaries.size(); ++i) {
      if(std::find(m_inflow_boundaries[i].begin(), m_inflow_boundaries[i].end(),eh) != m_inflow_boundaries[i].end()) {
        isConstrained = true;
        break;
      }
    }
    if(isConstrained) continue;

    if ( !m_obj->edgeExists(eh) || !isEdgeActive(eh))   { continue; }
    if(isSplitDesired(eh, maxEdge, desiredEdge, maxAngle)) {
    
      VertexHandle newVert;
      bool splitSuccess = performSplit(eh, newVert);
      
      //if the edge was constrained, constrain this vertex
      anySplits |= splitSuccess;
    }
  }
  
  return anySplits;
}

bool ElasticShell::isConstrained(const VertexHandle& v) const {
  for(unsigned int i = 0; i < m_constrained_vertices.size(); ++i)
    if(m_constrained_vertices[i] == v)
      return true;
  return false;
}


bool ElasticShell::isSplitDesired(const EdgeHandle& eh, double maxEdge, double desiredEdge, double maxAngle) {

  Scalar sEdge = 0.25*desiredEdge;

  VertexHandle v0 = m_obj->fromVertex(eh);
  VertexHandle v1 = m_obj->toVertex(eh);
  //don't split constrained edges. (alternatively, we could split them, and add the new vert to the constrained list.)
  
  bool aConstrained = false, bConstrained = false;
  for(unsigned int i = 0; i < m_constrained_vertices.size(); ++i) {
    if(m_constrained_vertices[i] == v0)
      aConstrained = true;
    if(m_constrained_vertices[i] == v1)
      bConstrained = true;
  }

  if(aConstrained && bConstrained) return false;

  assert(v0!=v1);
  Vec3d p0 = getVertexPosition(v0);
  Vec3d p1 = getVertexPosition(v1);
  Scalar edgeLength = (p1-p0).norm();
  bool doSplit = false;

  if(edgeLength > maxEdge) {
    //"Absolute longness"

    //only do the split if the adjacent faces have no edges longer than this one.
    bool edgeLongest = true;
    for(EdgeFaceIterator ef_it = m_obj->ef_iter(eh); ef_it && edgeLongest; ++ef_it) {
      for(FaceEdgeIterator fe_it = m_obj->fe_iter(*ef_it); fe_it; ++fe_it) {
        EdgeHandle otherEdge = *fe_it;
        if(otherEdge == eh) continue;

        Scalar otherLen = (getVertexPosition(m_obj->fromVertex(otherEdge)) - getVertexPosition(m_obj->toVertex(otherEdge))).norm();
        if(otherLen > edgeLength) {
          edgeLongest = false;
          break;
        }
      }
    }
    if(!edgeLongest) return false;

    doSplit = true;
  }

  
  if(!doSplit && edgeLength > desiredEdge) {

    //Additional criterion from Jiao et al. - Brochu uses just the above

    //"Relative longness"
    VertexHandle opp0, opp1;
    getEdgeOppositeVertices(*m_obj, eh, opp0, opp1);
    bool bothAnglesSmall = true;
    if(opp0.isValid()) {
      Vec3d p2 = getVertexPosition(opp0);
      Vec3d dir0 = p0-p2, dir1 = p1 - p2;
      Scalar angle0 = acos(dir0.dot(dir1));
      if(angle0 > maxAngle)
        bothAnglesSmall = false;
    }
    if(opp1.isValid()) {
      Vec3d p3 = getVertexPosition(opp1);
      Vec3d dir0 = p0-p3, dir1 = p1 - p3;
      Scalar angle1 = acos(dir0.dot(dir1));
      if(angle1 > maxAngle)
        bothAnglesSmall = false;
    }

    if(!bothAnglesSmall) {
      //check if shortest edge is above a threshold
      Scalar shortestEdge = 100*maxEdge;
      for(EdgeFaceIterator ef_it = m_obj->ef_iter(eh); ef_it; ++ef_it) {
        for(FaceEdgeIterator fe_it = m_obj->fe_iter(*ef_it); fe_it; ++fe_it) {
          EdgeHandle otherEdge = *fe_it;

          if(!otherEdge.isValid())continue;

          Scalar otherLen = (getVertexPosition(m_obj->fromVertex(otherEdge)) - getVertexPosition(m_obj->toVertex(otherEdge))).norm();
          if(otherLen < shortestEdge) {
            shortestEdge = otherLen;
          }
        }
      }
      if(shortestEdge > sEdge) {
        doSplit = true;
      }
    }
  }
  

  return doSplit;
}

bool ElasticShell::performSplit(const EdgeHandle& eh, VertexHandle& new_vert) {
  
  VertexHandle v0 = m_obj->fromVertex(eh);
  VertexHandle v1 = m_obj->toVertex(eh);
  assert(v0!=v1);
  Vec3d p0 = getVertexPosition(v0);
  Vec3d p1 = getVertexPosition(v1);
  
  Vec3d midpoint = 0.5f*(p0+p1);

  std::vector<FaceHandle> oldFaces;
  std::vector<Scalar> oldThicknesses;
  for(EdgeFaceIterator efit = m_obj->ef_iter(eh); efit; ++efit) {
    FaceHandle f = *efit;
    oldFaces.push_back(f);
    oldThicknesses.push_back(getThickness(f));
  }
  
  VertexHandle v2, v3;
  getEdgeOppositeVertices(*m_obj, eh, v2, v3);

  ElTopoCode::Vec3d midpoint_ET = ElTopoCode::toElTopo(midpoint);

  //perform the actual split
  std::vector<FaceHandle> newFaces;
  VertexHandle v_new = splitEdge(*m_obj, eh, newFaces);

  Vec3d velocity = 0.5f*(getVertexVelocity(v0) + getVertexVelocity(v1));
  Vec3d undef = 0.5f*(getVertexUndeformed(v0) + getVertexUndeformed(v1));
  setVertexVelocity(v_new, velocity);
  setVertexPosition(v_new, midpoint);
  setUndeformedVertexPosition(v_new, undef);

  //set consistent volumes and thickness for new faces
  assert(oldFaces.size() == newFaces.size()/2);
  Scalar newVolume = 0;
  for(unsigned int i = 0; i < oldFaces.size(); ++i) {
    m_thicknesses[newFaces[i*2]] = oldThicknesses[i];
    m_thicknesses[newFaces[i*2+1]] = oldThicknesses[i];
    m_volumes[newFaces[i*2]] = oldThicknesses[i] * getArea(newFaces[i*2], true);
    m_volumes[newFaces[i*2+1]] = oldThicknesses[i] * getArea(newFaces[i*2+1], true);
  }

  VertexFaceIterator vf_iter = m_obj->vf_iter(v_new);
  for(;vf_iter; ++vf_iter)
    setFaceActive(*vf_iter);

  //add edges
  VertexEdgeIterator ve_iter = m_obj->ve_iter(v_new);
  for(; ve_iter; ++ve_iter) {
    std::vector<ElTopoCode::Vec3d> edge_verts; 
    EdgeVertexIterator ev_iter = m_obj->ev_iter(*ve_iter);
    for(; ev_iter; ++ev_iter)
      edge_verts.push_back(ElTopoCode::toElTopo(getVertexPosition(*ev_iter)));
  }

  //add tris
  vf_iter = m_obj->vf_iter(v_new);
  for(;vf_iter; ++vf_iter) {
    std::vector<ElTopoCode::Vec3d> tri_verts; 
    FaceVertexIterator fv_iter = m_obj->fv_iter(*vf_iter);
    for(; fv_iter; ++fv_iter)
      tri_verts.push_back(ElTopoCode::toElTopo(getVertexPosition(*fv_iter)));
  }
  
  new_vert = v_new;
  return true;
}

void ElasticShell::performSplitET(const EdgeHandle& eh, const Vec3d& midpoint, VertexHandle& new_vert) {

  VertexHandle v0 = m_obj->fromVertex(eh);
  VertexHandle v1 = m_obj->toVertex(eh);
  assert(v0!=v1);
  Vec3d p0 = getVertexPosition(v0);
  Vec3d p1 = getVertexPosition(v1);

  std::vector<FaceHandle> oldFaces;
  std::vector<Scalar> oldThicknesses;
  for(EdgeFaceIterator efit = m_obj->ef_iter(eh); efit; ++efit) {
    FaceHandle f = *efit;
    oldFaces.push_back(f);
    oldThicknesses.push_back(getThickness(f));
  }

  VertexHandle v2, v3;
  getEdgeOppositeVertices(*m_obj, eh, v2, v3);

  ElTopoCode::Vec3d midpoint_ET = ElTopoCode::toElTopo(midpoint);

  //perform the actual split
  std::vector<FaceHandle> newFaces;
  VertexHandle v_new = splitEdge(*m_obj, eh, newFaces);

  
  //determine split fraction, and using it to lerp the vertex data
  Vec3d dx(p1-p0);
  double m2 = dx.squaredNorm();
  Scalar s = clamp((p1-midpoint).dot(dx)/m2, 0., 1.);
  
  Vec3d velocity = s*getVertexVelocity(v0) + (1-s)*getVertexVelocity(v1);
  Vec3d undef = s*getVertexUndeformed(v0) + (1-s)*getVertexUndeformed(v1);
  setVertexVelocity(v_new, velocity);
  setVertexPosition(v_new, midpoint);
  setUndeformedVertexPosition(v_new, undef);

  //set consistent volumes and thickness for new faces
  //TODO If the new point (chosen by ET) is allowed to be off the original edge, areas also change and this
  //code will need to be adjusted accordingly.
  assert(oldFaces.size() == newFaces.size()/2);
  Scalar newVolume = 0;
  for(unsigned int i = 0; i < oldFaces.size(); ++i) {
    m_thicknesses[newFaces[i*2]] = oldThicknesses[i];
    m_thicknesses[newFaces[i*2+1]] = oldThicknesses[i];
    m_volumes[newFaces[i*2]] = oldThicknesses[i] * getArea(newFaces[i*2], true);
    m_volumes[newFaces[i*2+1]] = oldThicknesses[i] * getArea(newFaces[i*2+1], true);
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
  m_obj->collapseEdge(eh, vert_to_remove, deletedEdges);

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

bool ElasticShell::performFlip(const EdgeHandle& eh, EdgeHandle& newEdge) {
      
  //determine volume of the region being flipped
    FaceHandle f0, f1;
    getEdgeFacePair(*m_obj, eh, f0, f1);
    Scalar totalVolume = m_volumes[f0] + m_volumes[f1];

    //EdgeHandle newEdge(-1);
    newEdge = flipEdge(*m_obj, eh);
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
    
    return true;
}

void ElasticShell::updateBroadPhaseStatic(const VertexHandle& vertex_a) 
{

  //update the broad phase grid to reflect the fact that two of the vertices are moving
  //i.e. change their bounding boxes to include their whole motion, while all other data
  //has just static bound boxes.

  ElTopoCode::Vec3d offset(m_proximity_epsilon, m_proximity_epsilon, m_proximity_epsilon);
  ElTopoCode::Vec3d low, high;
  
  //update vertices bounding box to include pseudomotion
  ElTopoCode::Vec3d pos = ElTopoCode::toElTopo(m_positions[vertex_a]);
  m_broad_phase.update_vertex( vertex_a.idx(), pos + offset, pos - offset);

  for ( VertexEdgeIterator veit = m_obj->ve_iter(vertex_a); veit; ++veit)
  {
    //get current edge pos
    EdgeHandle eh = *veit;
    VertexHandle vha = m_obj->fromVertex(eh);
    VertexHandle vhb = m_obj->toVertex(eh);
    ElTopoCode::Vec3d pa = ElTopoCode::toElTopo(m_positions[vha]);
    ElTopoCode::Vec3d pb = ElTopoCode::toElTopo(m_positions[vhb]);
    ElTopoCode::minmax(pa, pb, low, high);

    m_broad_phase.update_edge( (*veit).idx(), low, high );
  }

  for ( VertexFaceIterator vfit = m_obj->vf_iter(vertex_a); vfit; ++vfit )
  {
    //get current tri pos
    FaceHandle fh = *vfit;
    int c = 0;
    for(FaceVertexIterator fvit = m_obj->fv_iter(fh); fvit; ++fvit, ++c) {
      VertexHandle vha = *fvit;
      ElTopoCode::Vec3d pa = ElTopoCode::toElTopo(m_positions[vha]);
      if(c == 0)
        low = high = pa;
      else
        ElTopoCode::update_minmax(pa, low, high);
    }

    m_broad_phase.update_triangle( (*vfit).idx(), low, high );
  }
  
}

void ElasticShell::updateBroadPhaseForCollapse(const VertexHandle& vertex_a, const ElTopoCode::Vec3d& new_pos_a, 
                                               const VertexHandle& vertex_b, const ElTopoCode::Vec3d& new_pos_b) 
{

  //update the broad phase grid to reflect the fact that two of the vertices are moving
  //i.e. change their bounding boxes to include their whole motion, while all other data
  //has just static bound boxes.

  ElTopoCode::Vec3d offset(m_proximity_epsilon, m_proximity_epsilon, m_proximity_epsilon);
  ElTopoCode::Vec3d low, high;

  VertexHandle verts[2] = {vertex_a, vertex_b};
  ElTopoCode::Vec3d verts_pos[2] = {new_pos_a, new_pos_b};

  for(int i = 0; i < 2; ++i) {
    //update vertices bounding box to include pseudomotion
    ElTopoCode::minmax(ElTopoCode::toElTopo(m_positions[verts[i]]), verts_pos[i], low, high);
    low -= offset; high += offset;
    m_broad_phase.update_vertex( verts[i].idx(), low, high );

    for ( VertexEdgeIterator veit = m_obj->ve_iter(verts[i]); veit; ++veit)
    {
      //get current edge pos
      EdgeHandle eh = *veit;
      VertexHandle vha = m_obj->fromVertex(eh);
      VertexHandle vhb = m_obj->toVertex(eh);
      ElTopoCode::Vec3d pa = ElTopoCode::toElTopo(m_positions[vha]);
      ElTopoCode::Vec3d pb = ElTopoCode::toElTopo(m_positions[vhb]);
      ElTopoCode::minmax(pa, pb, low, high);
      
      //check if either vertex moved, and if so, add the position to the bound box
      for(int j = 0; j < 2; ++j) {
        if(vha == verts[j]) {
          ElTopoCode::update_minmax(verts_pos[j], low, high);
        }
        if(vhb == verts[j]) {
          ElTopoCode::update_minmax(verts_pos[j], low, high);
        }
      }
      
      m_broad_phase.update_edge( (*veit).idx(), low, high );
    }

    for ( VertexFaceIterator vfit = m_obj->vf_iter(verts[i]); vfit; ++vfit )
    {
      //get current tri pos
      FaceHandle fh = *vfit;
      int c = 0;
      for(FaceVertexIterator fvit = m_obj->fv_iter(fh); fvit; ++fvit, ++c) {
        VertexHandle vha = *fvit;
        ElTopoCode::Vec3d pa = ElTopoCode::toElTopo(m_positions[vha]);
        if(c == 0)
          low = high = pa;
        else
          ElTopoCode::update_minmax(pa, low, high);
      }

      //check if any vertex moved, and if so, add that position to the bound box
      for(int j = 0; j < 2; ++j) {
        for(FaceVertexIterator fvit = m_obj->fv_iter(fh); fvit; ++fvit) {
          VertexHandle vha = *fvit;
          if(vha == verts[j]) {
            ElTopoCode::update_minmax(verts_pos[j], low, high);
          }
        }
      }

      m_broad_phase.update_triangle( (*vfit).idx(), low, high );
    }
  }
}

// --------------------------------------------------------
///
/// Continuous collision detection between two triangles.  Duplicates collision detection, so 
/// should be used rarely, and only when inconvenient to run edge-edge and point-tri tests individually.
///
// --------------------------------------------------------

bool ElasticShell::checkTriangleVsTriangleCollisionForCollapse( const FaceHandle& triangle_a, const FaceHandle& triangle_b,  
                                                               const VertexHandle& source_vert, const VertexHandle& dest_vert,
                                                               ElTopoCode::Vec3d new_position)
{
  // --------------------------------------------------------
  // Point-triangle
  // --------------------------------------------------------

  // one point from triangle_a vs triangle_b

  for ( FaceVertexIterator fvit = m_obj->fv_iter(triangle_a); fvit; ++fvit)
  {
    VertexHandle vertex_0 = *fvit;
    
    VertexHandle vert_ids[3];
    ElTopoCode::Vec3d vert_pos[3];
    ElTopoCode::Vec3d new_vert_pos[3];

    int c = 0;
    bool shared_vert = false;
    for(FaceVertexIterator fvit2 = m_obj->fv_iter(triangle_b); fvit2; ++fvit2, ++c) {
      VertexHandle v = *fvit2;
      if(v == vertex_0) {
        shared_vert = true;
        break;
      }
      vert_ids[c] = v;
      vert_pos[c] = ElTopoCode::toElTopo(m_positions[vert_ids[c]]);
      new_vert_pos[c] = (v == source_vert || v == dest_vert)? new_position : vert_pos[c];
    }

    if(shared_vert) continue;
  
    ElTopoCode::Vec3d vert_0_pos, vert_0_newpos;
    vert_0_pos = ElTopoCode::toElTopo(m_positions[ vertex_0 ]);
    vert_0_newpos = (vertex_0 == source_vert || vertex_0 == dest_vert)? new_position : vert_0_pos;

    if ( ElTopoCode::point_triangle_collision(  
      vert_0_pos, vert_0_newpos, vertex_0.idx(),
      vert_pos[0], new_vert_pos[0], vert_ids[0].idx(), 
      vert_pos[1], new_vert_pos[1], vert_ids[1].idx(),
      vert_pos[2], new_vert_pos[2], vert_ids[2].idx() ) )

    {
      return true;
    }
  }

  // one point from triangle_b vs triangle_a

  for ( FaceVertexIterator fvit = m_obj->fv_iter(triangle_b); fvit; ++fvit)
  {
    VertexHandle vertex_0 = *fvit;

    VertexHandle vert_ids[3];
    ElTopoCode::Vec3d vert_pos[3];
    ElTopoCode::Vec3d new_vert_pos[3];

    int c = 0;
    bool shared_vert = false;
    for(FaceVertexIterator fvit2 = m_obj->fv_iter(triangle_a); fvit2; ++fvit2, ++c) {
      VertexHandle v = *fvit2;
      if(v == vertex_0) {
        shared_vert = true;
        break;
      }
      vert_ids[c] = v;
      vert_pos[c] = ElTopoCode::toElTopo(m_positions[vert_ids[c]]);
      new_vert_pos[c] = (v == source_vert || v == dest_vert)? new_position : vert_pos[c];
    }

    if(shared_vert) continue;

    ElTopoCode::Vec3d vert_0_pos, vert_0_newpos;
    vert_0_pos = ElTopoCode::toElTopo(m_positions[ vertex_0 ]);
    vert_0_newpos = (vertex_0 == source_vert || vertex_0 == dest_vert)? new_position : vert_0_pos;

    if ( ElTopoCode::point_triangle_collision(  
      vert_0_pos, vert_0_newpos, vertex_0.idx(),
      vert_pos[0], new_vert_pos[0], vert_ids[0].idx(), 
      vert_pos[1], new_vert_pos[1], vert_ids[1].idx(),
      vert_pos[2], new_vert_pos[2], vert_ids[2].idx() ) )

    {
      return true;
    }
  }



  // --------------------------------------------------------
  // edge-edge
  // --------------------------------------------------------
  
 
  for ( FaceEdgeIterator feit = m_obj->fe_iter(triangle_a); feit; ++feit)
  {

    EdgeHandle edge_0 = *feit;
    
    // one edge
    VertexHandle vertex_0 = m_obj->fromVertex(edge_0);
    VertexHandle vertex_1 = m_obj->toVertex(edge_0);

    for ( FaceEdgeIterator feit2 = m_obj->fe_iter(triangle_b); feit2; ++feit2)
    {
      
      EdgeHandle edge_1 = *feit2;

      // another edge
      VertexHandle vertex_2 = m_obj->fromVertex(edge_1);
      VertexHandle vertex_3 = m_obj->toVertex(edge_1);

      if ( vertex_0 == vertex_2 || vertex_0 == vertex_3 || vertex_1 == vertex_2 || vertex_1 == vertex_3 )
      {
        continue;
      }

      VertexHandle vert_ids[4] = {vertex_0, vertex_1, vertex_2, vertex_3};
      ElTopoCode::Vec3d vert_pos[4], new_vert_pos[4];
      for(int c = 0; c < 4; ++c) {
        VertexHandle v = vert_ids[c];
        vert_pos[c] = ElTopoCode::toElTopo(m_positions[vert_ids[c]]);
        new_vert_pos[c] = (v == source_vert || v == dest_vert)? new_position : vert_pos[c];
      }

      if ( ElTopoCode::segment_segment_collision(  
        vert_pos[0], new_vert_pos[0], vertex_0.idx(), 
        vert_pos[1], new_vert_pos[1], vertex_1.idx(),
        vert_pos[2], new_vert_pos[2], vertex_2.idx(),
        vert_pos[3], new_vert_pos[3], vertex_3.idx() ) )               
      {
        return true;
      }
    }
  }

  return false;
}

void ElasticShell::collapseEdges(double minAngle, double desiredEdge, double ratio_R, double ratio_r, double minEdge) {
  int count = 0;

  EdgeIterator e_it = m_obj->edges_begin();
  for(;e_it != m_obj->edges_end(); ++e_it) {

    EdgeHandle eh = *e_it;

    VertexHandle v0 = m_obj->fromVertex(eh);
    VertexHandle v1 = m_obj->toVertex(eh);
    assert(v0!=v1);
    Vec3d p0 = getVertexPosition(v0);
    Vec3d p1 = getVertexPosition(v1);
    Scalar edgeLength = (p1-p0).norm();

    bool doCollapse = false;

    //"Absolute small angle", "Relative shortness", "Absolute small angle", "Relative small triangle"
    bool smallAngle = false;
    bool relativeShortness = false;
    bool isSmallest = true;
    Scalar longestEdgeAll = 0;
    for(EdgeFaceIterator ef_iter = m_obj->ef_iter(eh); ef_iter; ++ef_iter) {
      Scalar longestEdge = edgeLength;
      for(FaceEdgeIterator fe_iter = m_obj->fe_iter(*ef_iter); fe_iter; ++fe_iter) {
        EdgeHandle edge_other = *fe_iter;
        if(edge_other == eh) continue;

        VertexHandle fromV = m_obj->fromVertex(edge_other);
        VertexHandle toV = m_obj->toVertex(edge_other);
        Scalar thisLength = (getVertexPosition(fromV) - getVertexPosition(toV)).norm();
        longestEdge = std::max(longestEdge, thisLength);
        if(thisLength < edgeLength)
          isSmallest = false;
        if(smallAngle) continue;

        VertexHandle oppV = fromV == v0 || fromV == v1 ? toV : fromV;
        Vec3d p2 = getVertexPosition(oppV);
        Vec3d dir0 = p0-p2, dir1 = p1 - p2;
        Scalar angle0 = acos(dir0.dot(dir1));
        if(angle0 < minAngle) {
          smallAngle = true;
        }
      }

      longestEdgeAll = std::max(longestEdgeAll, longestEdge);

      if(smallAngle && longestEdge < desiredEdge) //absolute shortness
        doCollapse = true;

    }

    if(longestEdgeAll < desiredEdge && edgeLength < ratio_R*longestEdgeAll) //Relative small triangle
      doCollapse = true;

    if(edgeLength < ratio_r*longestEdgeAll) //relative shortness
      doCollapse = true;

    if(isSmallest && longestEdgeAll < minEdge) //absolute small angle
      doCollapse = true;


    if(doCollapse) {
      
      bool isInflow = false;
      for(unsigned int i = 0; i < m_inflow_boundaries.size(); ++i) {
        if(std::find(m_inflow_boundaries[i].begin(), m_inflow_boundaries[i].end(),eh) != m_inflow_boundaries[i].end()) {
          isInflow = true;
          break;
        }
      }
      if(isInflow) continue;
      
      
      //don't collapse faces that have springs attached (for now).
      bool facesHaveSprings = false;
      for(EdgeFaceIterator efit = m_obj->ef_iter(eh); efit; ++efit) {
        FaceHandle fh = *efit;
        if(m_repulsion_springs->isFaceInUse(fh)) {
          facesHaveSprings = true;
          break;
        }
      }
      if(facesHaveSprings) continue;

      
      //don't collapse a constrained vertex
      bool v0_pinned = false,v1_pinned = false;
      for(unsigned int i = 0; i < m_constrained_vertices.size(); ++i) {
        if(v0 == m_constrained_vertices[i])
          v0_pinned = true;
        if(v1 == m_constrained_vertices[i])
          v1_pinned = true;
      }

      
      //check if either point is on the boundary
      bool v0_bdry, v1_bdry;
      v0_bdry = isVertexOnBoundary(*m_obj, v0);
      v1_bdry = isVertexOnBoundary(*m_obj, v1);

      Vec3d newPoint, newVel, newUndef;
      if(v0_bdry && v1_bdry) continue; //both verts on the boundary, don't collapse!
      
      if(v0_pinned && v1_pinned) continue; //both verts pinned, don't collapse
      
      //TODO What about when one is boundary, and the other is pinned?

      VertexHandle vert_to_remove, vert_to_keep;
      
      std::vector<VertexHandle> options_remove_point; options_remove_point.reserve(3);
      std::vector<VertexHandle> options_keep_point; options_keep_point.reserve(3);
      std::vector<Vec3d> options_new_pos; options_new_pos.reserve(3);
      std::vector<Vec3d> options_new_vel; options_new_vel.reserve(3);
      std::vector<Vec3d> options_new_undef; options_new_undef.reserve(3);
      if(v0_bdry || v0_pinned) {
        options_remove_point.push_back(v1);
        options_keep_point.push_back(v0);
        options_new_pos.push_back(p0);
        options_new_vel.push_back(getVertexVelocity(v0));
        options_new_undef.push_back(getVertexUndeformed(v0));
      }
      else if(v1_bdry || v1_pinned) {
        vert_to_remove = v0;
        vert_to_keep = v1;
        options_remove_point.push_back(v0);
        options_keep_point.push_back(v1);
        options_new_pos.push_back(p1);
        options_new_vel.push_back(getVertexVelocity(v1));
        options_new_undef.push_back(getVertexUndeformed(v1));
      }
      else { //either one works, so consider several options in case one or more is bad!
        //use average point
        options_remove_point.push_back(v0);
        options_keep_point.push_back(v1);
        options_new_pos.push_back(0.5f*(p0+p1));
        options_new_vel.push_back(0.5*(getVertexVelocity(v0) + getVertexVelocity(v1)));
        options_new_undef.push_back(0.5*(getVertexUndeformed(v0) + getVertexUndeformed(v1)));

        //use point 0
        options_remove_point.push_back(v0);
        options_keep_point.push_back(v1);
        options_new_pos.push_back(p1);
        options_new_vel.push_back(getVertexVelocity(v1));
        options_new_undef.push_back(getVertexUndeformed(v1));

        //use point 1
        options_remove_point.push_back(v1);
        options_keep_point.push_back(v0);
        options_new_pos.push_back(p0);
        options_new_vel.push_back(getVertexVelocity(v0));
        options_new_undef.push_back(getVertexUndeformed(v0));
      }

      bool collapseOkay = false;
      int choice = 0;
      /*
      while(!collapseOkay && choice < (int)options_remove_point.size()) {
        
        //assume good until proven bad
        collapseOkay = true;
        
        //compute expected areas and normals of all faces involved in the collapse,
        //to ensure no areas go near zero and no normals get badly flipped
        VertexHandle firstVert = options_keep_point[choice];
        VertexHandle secondVert = options_remove_point[choice];

        for(VertexFaceIterator vfit = m_obj->vf_iter(firstVert); vfit; ++vfit) {
          FaceHandle faceToCheck = *vfit;
          Vec3d faceVerts[3];
          Vec3d faceVertsNew[3];
          int vNo = 0;
          bool collapsingTri = false;
          for(FaceVertexIterator fvit = m_obj->fv_iter(faceToCheck); fvit; ++fvit) {
            VertexHandle curVert = *fvit;
            if(curVert == secondVert) {
              collapsingTri = true;
              break;
            }
            faceVerts[vNo] = m_positions[curVert];
            faceVertsNew[vNo] = (curVert == firstVert) ? options_new_pos[choice] : m_positions[curVert];
            ++vNo;
          }
          if(collapsingTri) continue;

          Vec3d normalOld = (faceVerts[2] - faceVerts[0]).cross(faceVerts[1]-faceVerts[0]);
          Scalar areaOld = fabs(normalOld.norm())/2;

          Vec3d normalNew = (faceVertsNew[2] - faceVertsNew[0]).cross(faceVertsNew[1]-faceVertsNew[0]);
          Scalar areaNew = fabs(normalNew.norm())/2;
          if(areaNew < 0.0001 * square(m_remesh_edge_length)) {
            collapseOkay = false;
            std::cout << "Prevented small area collapse\n";
            break;
          }
          if(normalOld.dot(normalNew) <= 0) {//direction flip, don't do it!
            collapseOkay = false;
            std::cout << "Prevented direction flip collapse\n";
            break;
          }
        }

        if(!collapseOkay) {
          ++choice;
          continue;
        }

        //TODO Refactor redundant code.
        //compute expected areas and normals of all faces involved in the collapse,
        //to ensure no areas go near zero and no normals get badly flipped
        for(VertexFaceIterator vfit = m_obj->vf_iter(secondVert); vfit; ++vfit) {
          FaceHandle faceToCheck = *vfit;
          Vec3d faceVerts[3];
          Vec3d faceVertsNew[3];
          int vNo = 0;
          bool collapsingTri = false;
          for(FaceVertexIterator fvit = m_obj->fv_iter(faceToCheck); fvit; ++fvit) {
            VertexHandle curVert = *fvit;
            if(curVert == firstVert) {
              collapsingTri = true;
              break;
            }
            faceVerts[vNo] = m_positions[curVert];
            faceVertsNew[vNo] = (curVert == secondVert) ? options_new_pos[choice] : m_positions[curVert];
            ++vNo;
          }
          if(collapsingTri) continue;
          Vec3d normalOld = (faceVerts[2] - faceVerts[0]).cross(faceVerts[1]-faceVerts[0]);
          Scalar areaOld = fabs(normalOld.norm())/2;

          Vec3d normalNew = (faceVertsNew[2] - faceVertsNew[0]).cross(faceVertsNew[1]-faceVertsNew[0]);
          Scalar areaNew = fabs(normalNew.norm())/2;

          if(areaNew < 0.0001 * square(m_remesh_edge_length)) {
            collapseOkay = false;
            std::cout << "Prevented small area collapse\n";
            break;
          }
          if(normalOld.dot(normalNew) <= 0) {//direction flip, don't do it!
            collapseOkay = false;
            std::cout << "Prevented direction flip collapse\n";
            break;
          }
        }

        if(!collapseOkay) {
          ++choice;
          continue;
        }
      }
      
      if(choice == options_keep_point.size()) continue; //none of the options was tolerable, so skip it.
      */

      //pick out the one we liked.
      vert_to_keep = options_keep_point[choice];
      vert_to_remove = options_remove_point[choice];
      newPoint = options_new_pos[choice];
      newVel = options_new_vel[choice];
      newUndef = options_new_undef[choice];
     

      performCollapse(eh, vert_to_remove, vert_to_keep, newPoint);


    }
  }
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
      Vec3d pos = m_positions[vh];
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
      releaseVertex(faceVerts[j]);
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
    Vec3d curPos = m_positions[vfrom];
    Vec3d curPos2 = m_positions[vto];
    
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
      
      releaseVertex(prevLowerVert);

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
    releaseVertex(prevLowerVert);
    

    m_inflow_boundaries[boundary] = newList;

    for(unsigned int i = 0; i < vertices.size(); ++i) {
      setVertexPosition(vertices[i], m_inflow_positions[boundary][i]);
      setUndeformedVertexPosition(vertices[i], m_inflow_positions[boundary][i]);
      setVertexVelocity(vertices[i], m_inflow_velocities[boundary][i]);
      m_vertex_masses[vertices[i]] = 0;
      m_damping_undeformed_positions[vertices[i]] = m_inflow_positions[boundary][i];

      constrainVertex(vertices[i], new FixedVelocityConstraint(m_inflow_positions[boundary][i], m_inflow_velocities[boundary][i], current_time));
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
    constrainVertex(otherVert, new FixedVelocityConstraint(pos, vel, 0));
    posList.push_back(pos);
    velList.push_back(vel);
    
    prevVert = sharedVert;
    
  }
  
  VertexHandle wrapVert = getSharedVertex(*m_obj, edgeList[0], edgeList[edgeList.size()-1]);
  if(!wrapVert.isValid()) {
    Vec3d pos = getVertexPosition(prevVert);
    constrainVertex(prevVert, new FixedVelocityConstraint(pos, vel, 0));
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
/*
void ElasticShell::setMass( const DofHandle& hnd, const Scalar& mass )
{
  assert(hnd.getType() == DofHandle::VERTEX_DOF);

  const VertexHandle vh = static_cast<const VertexHandle&>(hnd.getHandle());
  m_vertex_masses[vh] = mass;
}
*/


} //namespace BASim
