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

  RodBoundaryCondition(ElasticRod& rod);

  const BCList& scriptedVertices() const;
  bool isVertexScripted(int vertIdx) const;
  void setDesiredVertexPosition(int vertIdx, const Vec3d& position);
  const Vec3d& getDesiredVertexPosition(int vertIdx);
  void releaseVertex(int vertIdx);

  const BCList& scriptedEdges() const;
  bool isEdgeScripted(int edgeIdx) const;
  void setDesiredEdgeAngle(int edgeIdx, const Scalar& theta);
  const Scalar& getDesiredEdgeAngle(int edgeIdx);
  void releaseEdge(int edgeIdx);

protected:

  ElasticRod& m_rod;

  ObjPropHandle<BCList> m_scriptedVerts;
  VPropHandle<Vec3d> m_desiredPositions;
  VPropHandle<bool> m_isVertexScripted;

  ObjPropHandle<BCList> m_scriptedEdges;
  EPropHandle<Scalar> m_desiredTheta;
  EPropHandle<bool> m_isMaterialScripted;
};

} // namespace BASim

#endif // RODBOUNDARYCONDITION_HH
