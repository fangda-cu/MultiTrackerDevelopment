/**
 * \file ElasticShell.h
 *
 * \author batty@cs.columbia.edu
 * \date April 18, 2011
 */
#ifndef ELASTICSHELL_H
#define ELASTICSHELL_H

#include "BASim/src/Physics/DeformableObjects/PhysicalModel.hh"
#include "BASim/src/Core/TopologicalObject/TopObjProperty.hh"
#include "BASim/src/Collisions/ElTopo/broadphasegrid.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "surftrack.h"

namespace BASim {

const int ELASTIC_SHELL_DOFS_PER_VERTEX = 0; //nodal position vectors
const int ELASTIC_SHELL_DOFS_PER_EDGE = 0; //mid-edge normal bending DOFs (Grinspun et al. 2006)

class DeformableObject;
class ElasticShellForce;
class ShellVertexPointSpringForce;
class ShellStickyRepulsionForce;

class ElasticShell : public PhysicalModel, public ElTopo::SurfTrack::ConstrainedVerticesCollapsingCallback {

public:
  ElasticShell(DeformableObject* object, const FaceProperty<char>& shellFaces, Scalar timestep);
  ~ElasticShell();

  //*Inherited from PhysicalModel
  void computeForces(VecXd& force);
  void computeJacobian(Scalar scale, MatrixBase& J);
  void computeConservativeForcesEnergy(VecXd& f, Scalar& energy);

  const Scalar& getDof(const DofHandle& hnd) const;
  void setDof(const DofHandle& hnd, const Scalar& dof);
  const Scalar& getVel(const DofHandle& hnd) const;
  void setVel(const DofHandle& hnd, const Scalar& vel);
//  const Scalar& getMass(const DofHandle& hnd) const;
 
  int numVertexDofs() const { return ELASTIC_SHELL_DOFS_PER_VERTEX; }
  int numEdgeDofs() const { return ELASTIC_SHELL_DOFS_PER_EDGE; }
  int numFaceDofs() const { return 0; }
  int numTetDofs() const { return 0; }

  bool isVertexActive(const VertexHandle& v) const;
  bool isEdgeActive(const EdgeHandle& e) const; 
  bool isFaceActive(const FaceHandle& f) const { return m_active_faces[f] != 0; }
  bool isTetActive(const TetHandle& t) const { return false; }

  void getScriptedDofs(IntArray& dofIndices, std::vector<Scalar>& dofValues, Scalar time) const;

  void startStep(Scalar time, Scalar timestep);
  void endStep(Scalar time, Scalar timestep);

  //*Elastic Shell-specific
  void setFaceActive(const FaceHandle& f) {m_active_faces[f] = true; }

  const std::vector<ElasticShellForce*>& getForces() const;
  void addForce(ElasticShellForce* force);

//  void constrainEdgeXi(const EdgeHandle& eh, Scalar xiValue) {
//      constrainedEdges.push_back(eh);
//      constrainedXiValues.push_back(xiValue);
//  }

  void setRemeshing(bool enable, Scalar min_rez, Scalar max_rez, int iterations) {
    m_do_remeshing = enable;
    m_remesh_edge_max_len = max_rez;
    m_remesh_edge_min_len = min_rez;
    m_remeshing_iters = iterations;
  }
  
//  void setThicknessUpdating(bool enable) {
//    m_do_thickness_updates = enable;
//  }

  void setElTopoCollisions(bool enable) {
    m_do_eltopo_collisions = enable;
  }

  //All DOFs at once
  // these methods should have be removed because position access is now provided by DeformableObject; but 
  // too much code in other parts of the codebase need to change because they depend on this, so these
  // methods are kept and implemented to redirect the calls
  void setVertexPositions(const VertexProperty<Vec3d>& positions) { m_obj->setVertexPositions(positions); }
  void setVertexVelocities(const VertexProperty<Vec3d>& velocities) { m_obj->setVertexVelocities(velocities); }
//  void setVertexUndeformed(const VertexProperty<Vec3d>& undef) { m_obj->setVertexUndeformedPositions(undef); }
  
//  void setEdgeXis(const EdgeProperty<Scalar>& xi);
//  void setEdgeUndeformed(const EdgeProperty<Scalar>& undef);
//  void setEdgeVelocities(const EdgeProperty<Scalar>& vels);
  
  void setFaceLabels(const FaceProperty<Vec2i>& labels) { m_face_regions = labels; }
  Vec2i getFaceLabel(const FaceHandle& face) const { return m_face_regions[face]; }
  void setFaceLabel(const FaceHandle & face, const Vec2i & label) { m_face_regions[face] = label; }
  
  void setVertexConstraintLabel(const VertexProperty<int> & labels) { m_vertex_constraint_labels = labels; }
  VertexProperty<int> & getVertexConstraintLabels() { return m_vertex_constraint_labels; }
  int & getVertexConstraintLabel(const VertexHandle & vertex) { return m_vertex_constraint_labels[vertex]; }
  
  //Individual DOFs
//  Vec3d getVertexUndeformed(const VertexHandle& v) const { return m_obj->getVertexUndeformedPosition(v); }
  Vec3d getVertexPosition(const VertexHandle& v) const { return m_obj->getVertexPosition(v); }
  Vec3d getVertexVelocity(const VertexHandle& v) const { return m_obj->getVertexVelocity(v); }
  Vec3d getVertexDampingUndeformed(const VertexHandle& v) const { return m_obj->getVertexDampingUndeformedPosition(v); }

  const VertexProperty<Vec3d>& getVertexPositions() const{ return m_obj->getVertexPositions(); }

//  void setUndeformedVertexPosition(const VertexHandle& v, const Vec3d& pos) { m_obj->setVertexUndeformedPosition(v, pos); }
  void setVertexPosition(const VertexHandle& v, const Vec3d& pos) { m_obj->setVertexPosition(v, pos); }
  void setVertexVelocity(const VertexHandle& v, const Vec3d& vel) { m_obj->setVertexVelocity(v, vel); }

//  Scalar getEdgeUndeformedXi(const EdgeHandle& eh) const { return m_undef_xi[eh]; }
//  Scalar getEdgeXi(const EdgeHandle& eh) const { return m_xi[eh]; }
//  Scalar getEdgeVelocity(const EdgeHandle& eh) const { return m_xi_vel[eh]; }
//  Scalar getDampingUndeformedXi(const EdgeHandle& eh) const { return m_damping_undef_xi[eh]; }
  
//  const VertexProperty<Scalar> & getVertexMasses() const { return m_vertex_masses; }
//  const Scalar getModelVertexMass(const VertexHandle& vh) const {return m_vertex_masses[vh]; }

//  void recomputeVertexMass(const VertexHandle& v);
//  void computeMasses();

//  void setDensity(Scalar density);
//  void setThickness(Scalar thickness);

//  Scalar getMass(const VertexHandle& v) const { return m_obj->getVertexMass(v); }
//  Scalar getMass(const EdgeHandle& e) const { return m_edge_masses[e]; }
//  Scalar getThickness(const FaceHandle& f) const { return m_thicknesses[f]; }
//  void setThickness(const FaceHandle& f, Scalar thick) { m_thicknesses[f] = thick; m_volumes[f] = getArea(f) * thick; }
//  Scalar getThickness(const VertexHandle& vh) const;
//  Scalar getMaxThickness () const;
//  Scalar getMinThickness () const;
//  Scalar getVolume(const FaceHandle& f) const {return m_volumes[f]; }
  Scalar getArea(const FaceHandle& f, bool current = true) const;

  Vec3d getFaceNormal(const FaceHandle& f);
  void getFaceNormals(FaceProperty<Vec3d> & fNormals) const;
  void getVertexNormals(VertexProperty<Vec3d> & vNormals) const;
//  void getThickness(VertexProperty<Scalar> & vThickness) const;

  void addVertexPointSpring(const VertexHandle& v, const Vec3d& pos, Scalar stiffness, Scalar damping, Scalar length);
  void addVertexTriSpring(const FaceHandle& f, const VertexHandle& v, const Vec3d& pos, Scalar stiffness, Scalar damping, Scalar length);

  void setCollisionParams(Scalar proximity, Scalar epsilon, Scalar stiffness, Scalar damping);
  void setGroundPlane(bool enabled, Scalar height, Scalar velocity);
  void setSelfCollision(bool enabled);

  void setCollisionSphere(bool enabled, Scalar radius, Vec3d position, Vec3d velocity);
  void setCollisionObject(bool enabled, const Vec3d& position, const Vec3d& velocity, const ElTopoCode::Array3f& grid_data, const Vec3d& origin, Scalar dx);

  void setInflowSection(std::vector<EdgeHandle> edgeList, const Vec3d& vel, Scalar thickness);
  void setDeletionBox(const Vec3d& lowerBound, const Vec3d& upperBound);

  void remesh();
  bool generate_collapsed_position(ElTopo::SurfTrack & st, size_t v0, size_t v1, ElTopo::Vec3d & pos);  // SurfTrack::ConstrainedVerticesCollapsingCallback method

  void extendMesh(Scalar current_time);
  void deleteRegion();

  void removeFace(FaceHandle& f);

  void getSpringList(std::vector<Vec3d>& start, std::vector<Vec3d>& end) const;

  //Fracture functions
  typedef std::vector<VertexHandle> VHList;
  typedef std::vector<bool> BoolList;

  bool shouldFracture (const EdgeHandle & eh) const;
  bool isInflow(const EdgeHandle & eh) const;
  void setTearing(bool tearing, Scalar thres, Scalar rand){
      m_tearing = tearing;
      m_tear_thres = thres;
      m_tear_rand = rand;
  }

  void fracture();

  void getCollisionSphere(Vec3d& position, Scalar& radius) const {
    position = m_sphere_position;
    radius = m_sphere_radius;
  }


protected:

  void performTearing(const EdgeHandle & eh);

  void resolveCollisions(Scalar timestep);
  void updateThickness();

  void performSplit(const EdgeHandle& eh, const Vec3d& newpos, VertexHandle& newVert);
  void performCollapse(const EdgeHandle& eh, const VertexHandle& vert_to_remove, const VertexHandle& vert_to_keep, const Vec3d& new_position);
  bool performFlip(const EdgeHandle& eh, const FaceHandle f0, const FaceHandle& f1, EdgeHandle& newEdge);
  void performZippering(EdgeHandle e0, EdgeHandle e1, const std::vector<FaceHandle> & faces_deleted, const std::vector<std::vector<VertexHandle> > & faces_to_create, const std::vector<Vec2i> & face_labels_to_create, const std::vector<std::pair<FaceHandle, Vec2i> > & face_labels_to_change, std::vector<FaceHandle> & face_created);
  
  void addSelfCollisionForces();
  

  //Various shell data
//  EdgeProperty<Scalar> m_undef_xi;
//
//  EdgeProperty<Scalar> m_xi;
//
//  EdgeProperty<Scalar> m_xi_vel;
//  
//  VertexProperty<Scalar> m_vertex_masses;
//  EdgeProperty<Scalar> m_edge_masses;
  
  //"undeformed" configuration that is updated at each step to support Rayleigh damping/viscosity
  //This is also used as the "start of step" configuration for eltopo collision resolution
//  EdgeProperty<Scalar> m_damping_undef_xi;

//  FaceProperty<Scalar> m_thicknesses;
//  FaceProperty<Scalar> m_volumes;

  bool m_do_remeshing;
  Scalar m_remesh_edge_max_len;
  Scalar m_remesh_edge_min_len;
  int m_remeshing_iters;

//  bool m_momentum_conserving_remesh;

//  bool m_do_thickness_updates;
  bool m_do_eltopo_collisions;
  
//  std::vector<EdgeHandle> constrainedEdges;
//  std::vector<Scalar> constrainedXiValues;

//  Scalar m_density;

  FaceProperty<char> m_active_faces; //list of faces to which this model is applied
  //Note: this should ideally use booleans, but std::vector<bool> doesn't support references, which we need. (vector<bool> isn't technically a container)

  //For each face, a number used to identify distinct regions of the mesh
  //This is used to implement per-region volume constraints.
  FaceProperty<Vec2i> m_face_regions; 
  VertexProperty<int> m_vertex_constraint_labels;

  //The base object, and the list of forces
  DeformableObject* m_obj;
  std::vector<ElasticShellForce*> m_shell_forces;

  ShellVertexPointSpringForce* m_vert_point_springs;
  ShellStickyRepulsionForce* m_repulsion_springs;

  //To handle continually inflowing regions
  Scalar m_inflow_thickness;
  bool m_inflow;
  std::vector<std::vector<EdgeHandle> > m_inflow_boundaries;
  std::vector<std::vector<Vec3d> > m_inflow_positions;
  std::vector<std::vector<Vec3d> > m_inflow_velocities;
  std::vector<bool> m_inflow_lastdir;
  
  
  bool m_delete_region;
  Vec3d m_delete_lower, m_delete_upper;

  Scalar m_collision_epsilon; //epsilon tolerance for El Topo collisions
  Scalar m_collision_proximity; //distance at which to trigger spring penalty forces

  Scalar m_ground_height;
  bool m_self_collisions;
  bool m_ground_collisions;
  Scalar m_ground_velocity;
  Scalar m_collision_spring_stiffness, m_collision_spring_damping;
  Scalar m_sphere_radius;
  Vec3d m_sphere_position;
  bool m_sphere_collisions;
  Vec3d m_sphere_velocity;

  bool m_object_collisions;
  Vec3d m_object_position;
  Vec3d m_object_velocity;
  ElTopoCode::Array3f m_object_SDF;
  Vec3d m_object_origin;
  Scalar m_object_dx;

  ElTopoCode::BroadPhaseGrid m_broad_phase;

  //Fracture properties
  bool m_tearing;
  Scalar m_tear_thres;
  Scalar m_tear_rand;
};

}


#endif //ELASTICSHELL_H
