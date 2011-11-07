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

namespace BASim {

const int ELASTIC_SHELL_DOFS_PER_VERTEX = 3;
const int ELASTIC_SHELL_DOFS_PER_EDGE = 0; //for mid-edge normal bending discretization

class DeformableObject;
class ElasticShellForce;

class ElasticShell : public PhysicalModel {

public:
  ElasticShell(DeformableObject* object, const FaceProperty<char>& shellFaces);
  ~ElasticShell();

  //*Inherited from PhysicalModel
  void computeForces(VecXd& force);
  void computeJacobian(Scalar scale, MatrixBase& J);
  
  const Scalar& getDof(const DofHandle& hnd) const;
  void setDof(const DofHandle& hnd, const Scalar& dof);
  const Scalar& getVel(const DofHandle& hnd) const;
  void setVel(const DofHandle& hnd, const Scalar& vel);
  const Scalar& getMass(const DofHandle& hnd) const;
 
  int numVertexDofs() const { return ELASTIC_SHELL_DOFS_PER_VERTEX; }
  int numEdgeDofs() const { return ELASTIC_SHELL_DOFS_PER_EDGE; }
  int numFaceDofs() const { return 0; }
  int numTetDofs() const { return 0; }

  bool isVertexActive(const VertexHandle& v) const;
  bool isEdgeActive(const EdgeHandle& e) const; 
  bool isFaceActive(const FaceHandle& f) const { return m_active_faces[f] != 0; }
  bool isTetActive(const TetHandle& t) const { return false; }

  void setFaceActive(const FaceHandle& f) {m_active_faces[f] = true; }

  void getScriptedDofs(IntArray& dofIndices, std::vector<Scalar>& dofValues) const;

  //*Elastic Shell-specific
  const std::vector<ElasticShellForce*>& getForces() const;
  void addForce(ElasticShellForce* force);

  void setThickness(Scalar thickness);
  void setDensity(Scalar density);
  void setVertexPositions(const VertexProperty<Vec3d>& positions);
  void setVertexVelocities(const VertexProperty<Vec3d>& velocities);
  void setUndeformedConfig(const VertexProperty<Vec3d>& undef);
  void computeMasses();

  Vec3d getDampingUndeformedPosition(const VertexHandle& v) const { return m_damping_undeformed_positions[v]; }
  Vec3d getUndeformedPosition(const VertexHandle& v) const { return m_undeformed_positions[v]; }
  
  Vec3d getVertexPosition(const VertexHandle& v) const { return m_positions[v]; }
  Vec3d getVertexVelocity(const VertexHandle& v) const { return m_velocities[v]; }
  
  
  Scalar getXi(const EdgeHandle& eh);
  Scalar getUndeformedXi(const EdgeHandle& eh);
  Scalar getDampingUndeformedXi(const EdgeHandle& eh);

  Scalar getMass(const VertexHandle& v) const { return m_vertex_masses[v]; }
  Scalar getThickness(const FaceHandle& f) const { return m_thicknesses[f]; }
  Scalar getVolume(const FaceHandle& f) const {return m_volumes[f]; }
  Scalar getArea(const FaceHandle& f, bool current = true) const;
  void setVertexPosition(const VertexHandle& v, const Vec3d& pos) { m_positions[v] = pos; }
  void setUndeformedVertexPosition(const VertexHandle& v, const Vec3d& pos) { m_undeformed_positions[v] = pos; }
  void setVertexVelocity(const VertexHandle& v, const Vec3d& vel) { m_velocities[v] = vel; }

  void constrainVertex(const VertexHandle& v, const Vec3d& pos);

  void setInflowSection(std::vector<EdgeHandle> edgeList, const Vec3d& vel);

  void remesh(Scalar desiredEdge );
  void extendMesh();

  void startStep();
  void endStep();

protected:

  void updateThickness();

  bool splitEdges(double desiredEdge, double maxEdge, double maxAngle);
  void collapseEdges(double minAngle, double desiredEdge, double ratio_R, double ratio_r, double minEdge);
  void flipEdges();

  void updateBroadPhaseStatic(const VertexHandle& vertex_a);

  bool isSplitDesired(const EdgeHandle& eh, double maxEdge, double desiredEdge, double maxAngle);
  bool edgeSplitCausesCollision( const ElTopo::Vec3d& new_vertex_position, const ElTopo::Vec3d& new_vertex_smooth_position, EdgeHandle edge);
  bool performSplit(const EdgeHandle& eh);

  bool performCollapse(const EdgeHandle& eh);
  void updateBroadPhaseForCollapse(const VertexHandle& vertex_a, const ElTopo::Vec3d& new_pos_a, const VertexHandle& vertex_b, const ElTopo::Vec3d& new_pos_b);
  bool edgeCollapseCausesCollision(const VertexHandle& source_vertex, const VertexHandle& destination_vertex, const EdgeHandle& edge_index, const ElTopo::Vec3d& vertex_new_position );
  bool checkTriangleVsTriangleCollisionForCollapse( const FaceHandle& triangle_a, const FaceHandle& triangle_b, const VertexHandle& source_vert, const VertexHandle& dest_vert, ElTopo::Vec3d new_position);

  bool performFlip(const EdgeHandle& eh);
  bool edgeFlipCausesCollision( const EdgeHandle& edge_index, const VertexHandle& new_end_a, const VertexHandle& new_end_b);
  

  VertexProperty<Vec3d> m_positions;
  EdgeProperty<Scalar> m_xi;
  
  VertexProperty<Vec3d> m_velocities;
  EdgeProperty<Scalar> m_xi_vel;
  
  VertexProperty<Scalar> m_vertex_masses;
  EdgeProperty<Scalar> m_edge_masses;
  
  std::vector<VertexHandle> m_constrained_vertices;
  std::vector<Vec3d> m_constraint_positions;
  
  std::vector<std::vector<EdgeHandle> > m_inflow_boundaries;
  std::vector<std::vector<Vec3d> > m_inflow_positions;
  std::vector<Vec3d> m_inflow_velocity;
  std::vector<Scalar> m_inflow_thickness;

  VertexProperty<Vec3d> m_undeformed_positions;
  EdgeProperty<Scalar> m_undef_xi;
  
  VertexProperty<Vec3d> m_damping_undeformed_positions; //"undeformed" configuration that is updated at each step to support damping/viscosity
  EdgeProperty<Scalar> m_damping_undef_xi;

  FaceProperty<Scalar> m_thicknesses;
  FaceProperty<Scalar> m_volumes;
  Scalar m_thickness;
  Scalar m_density;
  
  FaceProperty<char> m_active_faces; //list of faces to which this model is applied
  //Note: this should ideally use booleans, but std::vector<bool> doesn't support references, which we need. (vector<bool> isn't technically a container)

  DeformableObject* m_obj;
  std::vector<ElasticShellForce*> m_shell_forces;
  
  //collision-safe remeshing stuff ->Move into subclass? Remesh-able shell?
  Scalar m_proximity_epsilon, m_improve_collision_epsilon;
  ElTopo::BroadPhaseGrid m_broad_phase;
};

}


#endif //ELASTICSHELL_H