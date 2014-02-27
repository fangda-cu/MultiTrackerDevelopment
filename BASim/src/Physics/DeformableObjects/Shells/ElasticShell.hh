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

typedef Eigen::Matrix<int, 2, 2> Mat2i;
  
const int ELASTIC_SHELL_DOFS_PER_VERTEX = 0; //nodal position vectors
const int ELASTIC_SHELL_DOFS_PER_EDGE = 0; //mid-edge normal bending DOFs (Grinspun et al. 2006)

class DeformableObject;
class ElasticShellForce;
class ShellVertexPointSpringForce;
class ShellStickyRepulsionForce;

class ElasticShell : public PhysicalModel, public ElTopo::SurfTrack::SolidVerticesCallback, public ElTopo::T1Transition::VelocityFieldCallback {

public:
  class SteppingCallback
  {
  public:
    virtual void beforeEndStep() = 0;
    virtual void afterStartStep() = 0;
  };
  
public:
  ElasticShell(DeformableObject* object, const FaceProperty<char>& shellFaces, Scalar timestep, SteppingCallback * stepping_callback = NULL, int doublebubble_scene = -1);
  ~ElasticShell();

  //*Inherited from PhysicalModel
  void computeForces(VecXd& force);
  void computeJacobian(Scalar scale, MatrixBase& J);
  void computeConservativeForcesEnergy(VecXd& f, Scalar& energy);

  const Scalar& getDof(const DofHandle& hnd) const;
  void setDof(const DofHandle& hnd, const Scalar& dof);
  const Scalar& getVel(const DofHandle& hnd) const;
  void setVel(const DofHandle& hnd, const Scalar& vel);
 
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

  void setRemeshing(bool enable, int iterations)
  {
    m_do_remeshing = enable;
    m_remeshing_iters = iterations;
  }
    
    void setElTopoParams(
                         Scalar collision_epsilon_fraction,
                         Scalar merge_proximity_epsilon_fraction,
                         Scalar remesh_edge_min_len,
                         Scalar remesh_edge_max_len,
                         bool perform_smoothing,
                         Scalar max_volume_change_fraction,
                         Scalar min_triangle_angle,
                         Scalar max_triangle_angle,
                         Scalar large_triangle_angle_to_split,
                         Scalar min_triangle_area_fraction,
                         bool t1_transition_enabled,
                         Scalar t1_pull_apart_distance_fraction,
                         bool smooth_subdivision
                         )
    {
        m_et_collision_epsilon_fraction = collision_epsilon_fraction;
        m_et_merge_proximity_epsilon_fraction = merge_proximity_epsilon_fraction;
        m_et_remesh_edge_min_len = remesh_edge_min_len;
        m_et_remesh_edge_max_len = remesh_edge_max_len;
        m_et_perform_smoothing = perform_smoothing;
        m_et_max_volume_change_fraction = max_volume_change_fraction;
        m_et_min_triangle_angle = min_triangle_angle;
        m_et_max_triangle_angle = max_triangle_angle;
        m_et_large_triangle_angle_to_split = large_triangle_angle_to_split;
        m_et_min_triangle_area_fraction = min_triangle_area_fraction;
        m_et_t1_transition_enabled = t1_transition_enabled;
        m_et_t1_pull_apart_distance_fraction = t1_pull_apart_distance_fraction;
        m_et_smooth_subdivision = smooth_subdivision;
    }
  
  //All DOFs at once
  // these methods should have be removed because position access is now provided by DeformableObject; but
  // too much code in other parts of the codebase need to change because they depend on this, so these
  // methods are kept and implemented to redirect the calls
  void setVertexPositions(const VertexProperty<Vec3d>& positions) { m_obj->setVertexPositions(positions); }
  void setVertexVelocities(const VertexProperty<Vec3d>& velocities) { m_obj->setVertexVelocities(velocities); }
  
  void setFaceLabels(const FaceProperty<Vec2i>& labels) { m_face_regions = labels; }
  Vec2i getFaceLabel(const FaceHandle& face) const { return m_face_regions[face]; }
  void setFaceLabel(const FaceHandle & face, const Vec2i & label) { m_face_regions[face] = label; }
  
  void setVertexConstraintLabel(const VertexProperty<int> & labels) { m_vertex_constraint_labels = labels; }
  VertexProperty<int> & getVertexConstraintLabels() { return m_vertex_constraint_labels; }
  int & getVertexConstraintLabel(const VertexHandle & vertex) { return m_vertex_constraint_labels[vertex]; }
  
  //Individual DOFs
  Vec3d getVertexPosition(const VertexHandle& v) const { return m_obj->getVertexPosition(v); }
  Vec3d getVertexPosition(const VertexHandle& v, const VertexHandle& v0) const { return m_obj->getVertexPosition(v, v0); }   // refer to ElTopo::DynamicSurface::get_position()
  Vec3d getVertexVelocity(const VertexHandle& v) const { return m_obj->getVertexVelocity(v); }
  Vec3d getVertexDampingUndeformed(const VertexHandle& v) const { return m_obj->getVertexDampingUndeformedPosition(v); }

  const VertexProperty<Vec3d>& getVertexPositions() const{ return m_obj->getVertexPositions(); }

  void setVertexPosition(const VertexHandle& v, const Vec3d& pos) { m_obj->setVertexPosition(v, pos); }
  void setVertexVelocity(const VertexHandle& v, const Vec3d& vel) { m_obj->setVertexVelocity(v, vel); }

  Scalar getArea(const FaceHandle& f, bool current = true) const;

  Vec3d getFaceNormal(const FaceHandle& f);
  void getFaceNormals(FaceProperty<Vec3d> & fNormals) const;
  void getVertexNormals(VertexProperty<Vec3d> & vNormals) const;

  void remesh(Scalar timestep, bool initial = false);
  
  // SurfTrack::ConstrainedVerticesCollapsingCallback method  
  bool generate_collapsed_position(ElTopo::SurfTrack & st, size_t v0, size_t v1, ElTopo::Vec3d & pos);
  bool generate_split_position(ElTopo::SurfTrack & st, size_t v0, size_t v1, ElTopo::Vec3d & pos);
  ElTopo::Vec3c generate_collapsed_solid_label(ElTopo::SurfTrack & st, size_t v0, size_t v1, const ElTopo::Vec3c & label0, const ElTopo::Vec3c & label1);
  ElTopo::Vec3c generate_split_solid_label(ElTopo::SurfTrack & st, size_t v0, size_t v1, const ElTopo::Vec3c & label0, const ElTopo::Vec3c & label1);
  bool generate_edge_popped_positions(ElTopo::SurfTrack & st, size_t oldv, const ElTopo::Vec2i & cut, ElTopo::Vec3d & pos_upper, ElTopo::Vec3d & pos_lower);
  bool generate_vertex_popped_positions(ElTopo::SurfTrack & st, size_t oldv, int A, int B, ElTopo::Vec3d & pos_a, ElTopo::Vec3d & pos_b);
  bool solid_edge_is_feature(const ElTopo::SurfTrack & st, size_t e);
  
  // T1Transition::VelocityFieldCallback methods
  ElTopo::Vec3d sampleVelocity(ElTopo::Vec3d & pos);
  bool sampleDirectionalDivergence(const ElTopo::Vec3d & pos, const ElTopo::Vec3d & dir, double & output);
  
  //Fracture functions
  typedef std::vector<VertexHandle> VHList;
  typedef std::vector<bool> BoolList;

    bool m_remesh_t1transition;
    bool m_remesh_smooth_subdivision;

    void setMeshEventCallback(ElTopo::SurfTrack::MeshEventCallback * cb) { m_mesheventcallback = cb; }
    
protected:

  void resolveCollisions(Scalar timestep);

  int onBBWall(const Vec3d & pos) const;
  Vec3d enforceBBWallConstraint(const Vec3d & input, int constraints) const;
  

  FaceProperty<char> m_active_faces; //list of faces to which this model is applied
  //Note: this should ideally use booleans, but std::vector<bool> doesn't support references, which we need. (vector<bool> isn't technically a container)

  //For each face, a number used to identify distinct regions of the mesh
  //This is used to implement per-region volume constraints.
  FaceProperty<Vec2i> m_face_regions; 
  VertexProperty<int> m_vertex_constraint_labels;

  //The base object, and the list of forces
  DeformableObject* m_obj;
  std::vector<ElasticShellForce*> m_shell_forces;

  ElTopoCode::BroadPhaseGrid m_broad_phase;
  
  // stepping callback
  SteppingCallback * m_stepping_callback;

  // other callbacks
  ElTopo::SurfTrack::MeshEventCallback * m_mesheventcallback;
  
  // scene ID for double bubble tests, temporary
  int m_doublebubble_scene;

    // ElTopo parameters
    Scalar  m_et_collision_epsilon_fraction;        // epsilon for collisions (expressed as percentage of mean edge length)
    Scalar  m_et_merge_proximity_epsilon_fraction;  // proximity epsilon for triggering snapping (expressed as percentage of mean edge length)
    Scalar  m_et_remesh_edge_max_len;
    Scalar  m_et_remesh_edge_min_len;
    bool    m_et_perform_smoothing;     // old hard-coded value: no
    Scalar  m_et_max_volume_change_fraction;    // old hard-coded value: 1e-4
    Scalar  m_et_min_triangle_angle;    // old hard-coded value: 3
    Scalar  m_et_max_triangle_angle;    // old hard-coded value: 177
    Scalar  m_et_large_triangle_angle_to_split; // old hard-coded value: 160
    Scalar  m_et_min_triangle_area_fraction;    // old hard-coded value: 0.02
    bool    m_et_t1_transition_enabled;
    Scalar  m_et_t1_pull_apart_distance_fraction;   // old hard-coded value: 0.5
    bool    m_et_smooth_subdivision;
    
    // remesh parameters
    bool m_do_remeshing;
    int m_remeshing_iters;
    
};

}


#endif //ELASTICSHELL_H
