/**
 * \file RodBoundaryCondition.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 12/09/2009
 */

#ifndef RODBOUNDARYCONDITION_HH
#define RODBOUNDARYCONDITION_HH

#include "../../Core/TopologicalObject/TopObjHandles.hh"

namespace BASim {

class ElasticRod;

/** Class for managing fixed/scripted vertices and edges of a rod */
// TODO: For consistency, methods in RodBoundaryCondition should acccept
//       edge and/or vertex iterators to be consistent with the rest of
//       BASim.
class RodBoundaryCondition
{
public:

  typedef std::vector<int> BCList;

  explicit RodBoundaryCondition(ElasticRod& rod);

  const BCList& scriptedVertices() const;
  bool isVertexScripted(int vertIdx) const;
  void setDesiredVertexPosition(int vertIdx, double t, const Vec3d& x, const Vec3d& v);
  void setDesiredVertexPosition(int vertIdx, const Vec3d& x);
  Vec3d getDesiredVertexPosition(int vertIdx, double t);
  void releaseVertex(int vertIdx);

  const BCList& scriptedEdges() const;
  bool isEdgeScripted(int edgeIdx) const;
  void setDesiredEdgeAngle(int edgeIdx, double t, const Scalar& theta, const Scalar& thetaDot);
  void setDesiredEdgeAngle(int edgeIdx, const Scalar& theta);
  Scalar getDesiredEdgeAngle(int edgeIdx, double t);
  void releaseEdge(int edgeIdx);

protected:

  ElasticRod& m_rod;

  ObjPropHandle<BCList> m_scriptedVerts;
  VPropHandle<double> m_positionTimeAnchor;
  VPropHandle<Vec3d> m_desiredPositions;
  VPropHandle<Vec3d> m_desiredVelocities;
  VPropHandle<bool> m_isVertexScripted;

  ObjPropHandle<BCList> m_scriptedEdges;
  VPropHandle<double> m_thetaTimeAnchor;
  EPropHandle<Scalar> m_desiredTheta;
  EPropHandle<Scalar> m_desiredThetaDot;
  EPropHandle<bool> m_isMaterialScripted;
};

} // namespace BASim

#endif // RODBOUNDARYCONDITION_HH
