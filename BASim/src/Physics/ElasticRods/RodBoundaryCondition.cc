#include "ElasticRod.hh"
#include "RodBoundaryCondition.hh"

namespace BASim {

RodBoundaryCondition::RodBoundaryCondition(ElasticRod& rod)
  : m_rod(rod)
{
  m_rod.add_property(m_scriptedVerts, "list of scripted vertices");
  m_rod.add_property(m_positionTimeAnchor, "time at which desired positions/velocities are anchored");
  m_rod.add_property(m_desiredPositions, "desired vertex positions", Vec3d::Zero().eval());
  m_rod.add_property(m_desiredVelocities, "desired vertex velocity", Vec3d::Zero().eval());
  m_rod.add_property(m_isVertexScripted, "is vertex scripted", false);

  m_rod.add_property(m_scriptedEdges, "list of scripted edges");
  m_rod.add_property(m_thetaTimeAnchor, "time at which desired theta/thetaDot are anchored");
  m_rod.add_property(m_desiredTheta, "desired theta values", 0.0);
  m_rod.add_property(m_desiredThetaDot, "desired thetaDot values", 0.0);
  m_rod.add_property(m_isMaterialScripted, "is material scripted", false);
}

const RodBoundaryCondition::BCList&
RodBoundaryCondition::scriptedVertices() const
{
  return m_rod.property(m_scriptedVerts);
}

bool RodBoundaryCondition::isVertexScripted(int vertIdx) const
{
  return m_rod.property(m_isVertexScripted)[vertIdx];
}

void RodBoundaryCondition::setDesiredVertexPosition(int vertIdx, const Vec3d& x)
{
  setDesiredVertexPosition( vertIdx, 0, x, Vec3d::Zero() );
}

void RodBoundaryCondition::setDesiredVertexPosition(int vertIdx, double t, const Vec3d& x, const Vec3d& v)
{
  assert(vertIdx >= 0);
  assert(vertIdx < m_rod.nv());

  m_rod.property(m_isVertexScripted)[vertIdx] = true;

  BCList& verts = m_rod.property(m_scriptedVerts);
  BCList::iterator result
    = std::find(verts.begin(), verts.end(), vertIdx);
  if (result == verts.end()) verts.push_back(vertIdx);

  m_rod.property(m_positionTimeAnchor)[vertIdx] = t;
  m_rod.property(m_desiredPositions)[vertIdx] = x;
  m_rod.property(m_desiredVelocities)[vertIdx] = v;

  std::cout << "RodBoundaryCondition::setDesiredVertexPosition: Setting new linear function for vertexIdx = " << vertIdx << " t0 = " << t << " x = " << x << " v = " << v << " x("<<(t+0.004167)<<") = " << getDesiredVertexPosition( vertIdx, t + 0.004167 ) << std::endl;
}

Vec3d RodBoundaryCondition::getDesiredVertexPosition(int vertIdx, double t)
{
  double t0 = m_rod.property(m_positionTimeAnchor)[vertIdx];
  Vec3d x0 = m_rod.property(m_desiredPositions)[vertIdx];
  Vec3d v0 = m_rod.property(m_desiredVelocities)[vertIdx];

  Vec3d result = x0 + (t-t0)*v0;

  std::cout << "RodBoundaryCondition::getDesiredVertexPosition: vertIdx = " << vertIdx << " t = " << t << " t0 = " << t0 << " x0 = " << x0 << " v0 = " << v0 << " result = " << result << std::endl;

  return result;
}

void RodBoundaryCondition::releaseVertex(int vertIdx)
{
  assert(vertIdx >= 0);
  assert(vertIdx < m_rod.nv());

  m_rod.property(m_isVertexScripted)[vertIdx] = false;

  BCList& verts = m_rod.property(m_scriptedVerts);

//  std::remove(verts.begin(), verts.end(), vertIdx);   // this doesn't change the size of vector.
	for(int i=0; i<(int)verts.size(); i++) {
		if (verts[i] == vertIdx) {
			verts.erase(verts.begin() + i);
			break;
		}
	}

}

const RodBoundaryCondition::BCList& RodBoundaryCondition::scriptedEdges() const
{
  return m_rod.property(m_scriptedEdges);
}

bool RodBoundaryCondition::isEdgeScripted(int edgeIdx) const
{
  return m_rod.property(m_isMaterialScripted)[edgeIdx];
}

void RodBoundaryCondition::setDesiredEdgeAngle(int edgeIdx, const Scalar& theta)
{
  setDesiredEdgeAngle( edgeIdx, 0, theta, 0 );
}

void RodBoundaryCondition::setDesiredEdgeAngle(int edgeIdx, double t, const Scalar& theta, const Scalar& thetaDot)
{
  assert(edgeIdx >= 0);
  assert(edgeIdx < m_rod.ne());

  m_rod.property(m_isMaterialScripted)[edgeIdx] = true;

  BCList& edges = m_rod.property(m_scriptedEdges);
  BCList::iterator result
    = std::find(edges.begin(), edges.end(), edgeIdx);
  if (result == edges.end()) edges.push_back(edgeIdx);

  m_rod.property(m_thetaTimeAnchor)[edgeIdx] = t;
  m_rod.property(m_desiredTheta)[edgeIdx]    = theta;
  m_rod.property(m_desiredThetaDot)[edgeIdx] = thetaDot;
}

Scalar RodBoundaryCondition::getDesiredEdgeAngle(int edgeIdx, double t)
{
  double t0        = m_rod.property(m_thetaTimeAnchor)[edgeIdx];
  Scalar theta0    = m_rod.property(m_desiredTheta)[edgeIdx];
  Scalar thetaDot0 = m_rod.property(m_desiredThetaDot)[edgeIdx];
 
  return theta0 + (t-t0)*thetaDot0;
}

void RodBoundaryCondition::releaseEdge(int edgeIdx)
{
  assert(edgeIdx >= 0);
  assert(edgeIdx < m_rod.ne());

  m_rod.property(m_isMaterialScripted)[edgeIdx] = false;

  BCList& edges = m_rod.property(m_scriptedEdges);
//  std::remove(edges.begin(), edges.end(), edgeIdx);	// this doesn't change the size of vector.
  
	for(int i=0; i<(int)edges.size(); i++) {
		if (edges[i] == edgeIdx) {
			edges.erase(edges.begin() + i);
			break;
		}
	}
}

} // namespace BASim
