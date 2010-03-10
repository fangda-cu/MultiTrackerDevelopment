/**
 * \file ElasticRod.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef ELASTICROD_HH
#define ELASTICROD_HH

#include "../PhysObject.hh"
#include "../../Collisions/CollisionObject.hh"

namespace BASim {

class RodForce;
class RodPenaltyForce;

/** Base class for rods. The degrees of freedom for rods are the
    vertex positions (3 dofs per vertex) and the angles between the
    reference frames and the material frames (1 dof per edge). The
    default indexing that is set up is to interleave the verrtex
    degrees of freedom with the edge degrees of freedom as follows:
    \f$\left<x_0,y_0,z_0,\theta_0,x_1,y_1,z_1,\theta_1,...\right>\f$.
*/
class ElasticRod : public PhysObject, public CollisionObject
{
public:

  typedef std::vector<RodForce*> RodForces;

  ElasticRod(int numVertices = 3, bool closed = false);
  virtual ~ElasticRod();

  /** \name Inherited from PhysObject */

  //@{

  virtual void setup();
  virtual void computeForces(VecXd& force);
  virtual void computeJacobian(MatrixBase& J);

  virtual const Scalar&
  getVertexDof(const vertex_handle& vh, int num) const;

  virtual void
  setVertexDof(const vertex_handle& vh, int num, const Scalar& dof);

  virtual const Scalar&
  getEdgeDof(const edge_handle& eh, int num) const;

  virtual void
  setEdgeDof(const edge_handle& eh, int num, const Scalar& dof);

  virtual const Scalar&
  getVertexVel(const vertex_handle& vh, int num) const;

  virtual void
  setVertexVel(const vertex_handle& vh, int num, const Scalar& vel);

  virtual const Scalar&
  getEdgeVel(const edge_handle& eh, int num) const;

  virtual void
  setEdgeVel(const edge_handle& eh, int num, const Scalar& vel);

  virtual const Scalar&
  getVertexMass(const vertex_handle& vh, int num) const;

  virtual void
  setVertexMass(const vertex_handle& vh, int num, const Scalar& mass);

  virtual const Scalar&
  getEdgeMass(const edge_handle& eh, int num) const;

  virtual void
  setEdgeMass(const edge_handle& eh, int num, const Scalar& mass);

  //@}

  edge_handle inEdge(const vertex_handle& vh) const;
  edge_handle outEdge(const vertex_handle& vh) const;

  const Scalar& getYoungsModulus() const;
  void setYoungsModulus(const Scalar& E);

  const Scalar& getShearModulus() const;
  void setShearModulus(const Scalar& G);

  const Scalar& getViscosity() const;
  void setViscosity(const Scalar& mu);

  virtual int vertIdx(const vertex_handle& vh, int coordinate) const;
  virtual int edgeIdx(const edge_handle& eh) const;

  virtual int vertIdx(int vertexNumber, int coordinate) const;
  virtual int edgeIdx(int edgeNumber) const;

  const Vec3d& getVertex(const vertex_handle& vh) const;
  void setVertex(const vertex_handle& vh, const Vec3d& v);

  const Vec3d& getFromVertex(const edge_handle& eh) const;
  const Vec3d& getToVertex(const edge_handle& eh) const;

  const Vec3d& getVertex(int i) const;
  void setVertex(int i, const Vec3d& v);

  const Scalar& getTheta(const edge_handle& eh) const;
  void setTheta(const edge_handle& eh, const Scalar& t);

  const Scalar& getTheta(int j) const;
  void setTheta(int j, const Scalar& t);

  const Vec3d& getVelocity(int i) const;
  void setVelocity(int i, const Vec3d& v);

  const Scalar& getThetaDot(int j) const;
  void setThetaDot(int j, const Scalar& td);

  const Vec3d& getEdge(const edge_handle& eh) const;
  void setEdge(const edge_handle& eh, const Vec3d& edge);

  const Vec3d& getEdge(int j) const;
  void setEdge(int j, const Vec3d& edge);

  const Vec3d& getTangent(const edge_handle& eh) const;
  void setTangent(const edge_handle& eh, const Vec3d& tangent);

  const Vec3d& getTangent(int j) const;
  void setTangent(int j, const Vec3d& tangent);

  const Scalar& getVertexMass(int i) const;
  void setVertexMass(int i, const Scalar& m);

  const Scalar& getEdgeInertia(int j) const;
  void setEdgeInertia(int j, const Scalar& I);

  /** Sets a circular cross-section for all of the edges in the rod.

      \param[in] r The radius of the cross-section.
  */
  void setRadius(const Scalar& r);

  /** Sets an elliptical cross-section for all of the edges in the rod.

      \param[in] a The major axis of the ellipse.
      \param[in] b The minor axis of the ellipse.
  */
  void setRadius(const Scalar& a, const Scalar& b);

  /** Sets a circular cross-section for a specific edge in the rod.

      \param[in] j The edge number.
      \param[in] r The radius of the cross-section.
  */
  void setRadius(int j, const Scalar& r);
  void setRadius(const edge_handle& eh, const Scalar& r);

  /** Sets an elliptical cross-section for a specific edge in the rod.

      \param[in] j The edge number.
      \param[in] a The major radius of the ellipse.
      \param[in] b The minor radius of the ellipse.
  */
  void setRadius(int j, const Scalar& a, const Scalar& b);

  /** Gets the radius of the rod. It only makes sense to use this call
      if all of the edges in the rod have uniform, circular
      cross-sections.

      \return The radius of the rod.
  */
  const Scalar& radius() const;

  /** Gets the radius of an edge in the rod. It only makes sense to
      use this call if the edge has a circular cross-section.

      \param[in] j The edge number.
      \return The radius of the edge.
  */
  const Scalar& radius(int j) const;

  /** Gets the major radius of an edge in the rod.

      \param[in] j The edge number.
      \return The major radius of the edge.
  */
  const Scalar& radiusA(int j) const;
  const Scalar& radiusA(const edge_handle& eh) const;

  /** Gets the minor radius of an edge in the rod.

      \param[in] j The edge number.
      \return The minor radius of the edge.
  */
  const Scalar& radiusB(int j) const;
  const Scalar& radiusB(const edge_handle& eh) const;

  const Vec3d& getReferenceDirector1(const edge_handle& eh) const;
  void setReferenceDirector1(const edge_handle& eh, const Vec3d& u);
  const Vec3d& getReferenceDirector2(const edge_handle& eh) const;
  void setReferenceDirector2(const edge_handle& eh, const Vec3d& v);

  const Vec3d& getReferenceDirector1(int j) const;
  void setReferenceDirector1(int j, const Vec3d& u);
  const Vec3d& getReferenceDirector2(int j) const;
  void setReferenceDirector2(int j, const Vec3d& v);

  const Vec3d& getMaterial1(const edge_handle& eh) const;
  void setMaterial1(const edge_handle& eh, const Vec3d& m1);
  const Vec3d& getMaterial2(const edge_handle& eh) const;
  void setMaterial2(const edge_handle& eh, const Vec3d& m2);

  const Vec3d& getMaterial1(int j) const;
  void setMaterial1(int j, const Vec3d& m1);
  const Vec3d& getMaterial2(int j) const;
  void setMaterial2(int j, const Vec3d& m2);

  const Scalar& getEdgeLength(const edge_handle& eh) const;
  void setEdgeLength(const edge_handle& eh, const Scalar& len);

  const Scalar& getEdgeLength(int j) const;
  void setEdgeLength(int j, const Scalar& len);

  const Scalar& getVoronoiLength(int i) const;
  void setVoronoiLength(int i, const Scalar& len);

  const Scalar& getPhi(const vertex_handle& vh) const;
  void setPhi(const vertex_handle& vh, const Scalar& phi);

  const Scalar& getPhi(int i) const;
  void setPhi(int i, const Scalar& phi);

  const Vec3d& getCurvatureBinormal(const vertex_handle& vh) const;
  void setCurvatureBinormal(const vertex_handle& vh, const Vec3d& kb);

  const Vec3d& getCurvatureBinormal(int i) const;
  void setCurvatureBinormal(int i, const Vec3d& kb);

  bool vertFixed(const vertex_handle& vh) const;
  bool vertFixed(int i) const;
  void fixVert(int i);
  void unfixVert(int i);

  bool edgeFixed(const edge_handle& eh) const;
  bool edgeFixed(int j) const;
  void fixEdge(int j);
  void unfixEdge(int j);

  const IntArray& fixed() const;

  bool quasistatic() const;
  void setQuasistatic(bool q);

  enum RefFrameType { SpaceParallel, TimeParallel };
  RefFrameType refFrameType() const;
  void setRefFrameType(RefFrameType type);

  const Scalar& density() const { return property(m_density); }
  void setDensity(const Scalar& d) { property(m_density) = d; }

  const RodForces& getForces() const;
  RodForces& getForces();
  void addForce(RodForce* force);

  virtual void updateProperties();
  virtual void updateReferenceProperties();
  virtual void verifyProperties();

  ////////////////////////////////////////////////////////////////////////////////
  // 
  // Needed for collisions. Should rods be based off of CollisionObject
  // so there is a uniform interface for colliding objects?
  
  void collisionsBegin(Real dt)
  {
      for (int i=0; i<nv(); ++i)
      {
          // Start positions for this timestep are the end positions from last timestep
          //
          getStartPositions()[i] = getEndPositions()[i];
  
          // Candidate end positions are the current positions
          //
          getEndPositions()[i] = getVertex(i);
  
          // Average velocity
          //
          getVelocities()[i] = (getEndPositions()[i] - getStartPositions()[i]) / dt;
      }      
  }

  void collisionsEnd(Real dt)
  {
      // Compute final, end-of-timestep positions
      //
      updateEndPositions(dt);
  
      for (int i=0; i<nv(); ++i)
      {
          Vec3d velocityChange = (getEndPositions()[i] - getVertex(i)) / dt;
          setVelocity(i, getVelocity(i) + velocityChange);
          setVertex(i, getEndPositions()[i]);

       //   cerr << "Vertex " << i << ", velocity change = " << velocityChange << endl;
      }      
  }

  void updateEndPositions(Real dt)
  {
      // Update the end-of-timestep positions using the current velocity
      //
      for (int i=0; i<nv(); ++i)
          getEndPositions()[i] = getStartPositions()[i] + getVelocities()[i] * dt;
  }

  void setCollisionStartPositions()
  {
    // Collision code uses the start position data, so update it using the
    // current vertex positions
    //
    //getStartPositions().resize(nv());
    for (int i=0; i<nv(); ++i)
        getStartPositions()[i] = getVertex(i);
  }

  std::vector<Vec3d>& getStartPositions() { return m_previousPositions; }
  std::vector<Vec3d>& getEndPositions() { return m_currentPositions; }
  
  // Interface from CollisionObject base clasee
  virtual Positions& getPositions() { return m_previousPositions; }
  virtual Velocities& getVelocities() { return m_currentVelocities; }
  virtual Indices& getEdgeIndices() { return m_edgeIndices; }
  virtual Indices& getTriangleIndices()
  {
     // What are you doing trying to get
     // triangle indices from a rod?!?
     //
     assert(false);

     // Compiler doesn't like it if you
     // don't return something, even if
     // it will never get here
     //
     return m_edgeIndices;
  }

  virtual Real getMassInverse(uint vertex)
  {
    if (vertFixed(vertex))
      return 0.0;
    else
      return (1.0 / getVertexMass(vertex));
  }

  virtual Real getFrictionCoefficient() { return m_friction; }
  void setFrictionCoefficient(Real k) { m_friction = k; }
  virtual double getCoefficientOfRestitution() { return m_cor; }
  void setCoefficientOfRestitution(double cor) { m_cor = cor; }
  virtual double getSeparationStrength() { return m_separationStrength; }
  void setSeparationStrength(double k) { m_separationStrength = k; }
  // Again we assume rods are cylindrical for collision
  virtual Real getThickness() { return radius(); }

  RodPenaltyForce *getPenaltyForce()
  { return m_rodPenaltyForce; }

  void setPenaltyForce(RodPenaltyForce* force)
  { m_rodPenaltyForce = force; }
  
  ////////////////////////////////////////////////////////////////////////////////
  
  
protected:

  void computeEdges();
  void computeTangents();
  void computeEdgeLengths();
  void computeVoronoiLengths();
  void computeReferenceDirectors();
  void computeSpaceParallel();
  void computeTimeParallel();
  void computeCurvatureBinormals();
  void computePhi();
  void computeMaterialDirectors();
  void computeVertexMasses();
  void computeEdgeInertias();
  void setupDofIndices();

  /** Computes the mass of an elliptical cylinder, which is the
      generic representation of each edge of the rod.

      \param[in] density The volumetric density of the cylinder.
      \param[in] a The major radius of the cylinder.
      \param[in] b The minor radius of the cylinder.
      \param[in] h The height of the cylinder.
      \return The mass of the cylinder.
  */
  Scalar computeMass(Scalar density, Scalar a, Scalar b, Scalar h);

  ObjPropHandle<RodForces> m_forces; ///< forces acting on the rod
  ObjPropHandle<bool> m_quasistatic;
  ObjPropHandle<RefFrameType> m_refFrameType;
  ObjPropHandle<Scalar> m_density;
  ObjPropHandle<IntArray> m_fixed;
  ObjPropHandle<IntArray> m_fixedVerts;
  ObjPropHandle<IntArray> m_fixedEdges;
  ObjPropHandle<Scalar> m_YoungsModulus;
  ObjPropHandle<Scalar> m_ShearModulus;
  ObjPropHandle<Scalar> m_viscosity;

  VPropHandle<Vec3d> m_vertexPositions;
  VPropHandle<Vec3d> m_vertexVelocities;
  VPropHandle<Scalar> m_voronoiLengths;
  VPropHandle<Scalar> m_vertexMasses;
  VPropHandle<bool> m_vertexFixed;
  VPropHandle<Scalar> m_phi; ///< twist of the reference frame
  VPropHandle<Vec3d> m_curvatureBinormal;
  VPropHandle<int> m_vertIdx;

  EPropHandle<Scalar> m_theta;
  EPropHandle<Scalar> m_thetaDot;
  EPropHandle< Util::pair<Scalar, Scalar> > m_edgeRadius;
  EPropHandle<Scalar> m_edgeInertias;
  EPropHandle< Util::pair<Vec3d,Vec3d> > m_referenceDirectors;
  EPropHandle< Util::pair<Vec3d,Vec3d> > m_materialDirectors;
  EPropHandle<Vec3d> m_edges;
  EPropHandle<Vec3d> m_tangents;
  EPropHandle<Scalar> m_edgeLengths; ///< lengths of edges
  EPropHandle<bool> m_edgeFixed;
  EPropHandle<int> m_edgeIdx;
  
  ////////////////////////////////////////////////////////////////////////////////
  //
  // Needed for collisions. They should really be integrated with the the rest of
  // the data but to get this working I'm leaving them seperate as it's easier
  // to debug
  
  std::vector<Vec3d> m_currentPositions;
  std::vector<Vec3d> m_previousPositions;
  std::vector<Vec3d> m_currentVelocities;
  std::vector<uint> m_edgeIndices;
  
  // This should be a property like everything else
  double m_friction;
  double m_cor;
  double m_separationStrength;

  // This is silly. The force now lives in the time stepper but the collision code
  // only has the rods currently so needs to ask the rod for the penalty force
  // so it can apply a force.
  RodPenaltyForce* m_rodPenaltyForce;
  //  
  ////////////////////////////////////////////////////////////////////////////////
};

typedef std::vector<ElasticRod *> ElasticRods;
typedef std::vector<ElasticRod *>::iterator ElasticRodsIterator;

#include "ElasticRod.inl"

} // namespace BASim

#endif // ELASTICROD_HH
