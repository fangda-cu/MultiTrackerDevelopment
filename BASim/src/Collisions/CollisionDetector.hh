/*
 * CollisionDetector.hh
 *
 *  Created on: 17/03/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef COLLISIONDETECTOR_HH_
#define COLLISIONDETECTOR_HH_

#include "BoundingBox.hh"
#include "Geometry.hh"
#include "BVH.hh"

namespace BASim
{

class CollisionDetector
{
	const GeometricData m_geodata;
	std::vector<const TopologicalElement*> m_elements;
	const double m_time_step;
	BVH m_bvh;
	std::list<ContinuousTimeCollision>* m_collisions;

public:
	// During construction, the BVH tree is created around the initial geometry.
	CollisionDetector(const VecXd& x, const VecXd& v, const std::vector<std::pair<int, int> >& edges,
			const std::vector<TriangularFace>& faces, const std::vector<double>& vert_radii, const std::vector<double>& masses,
			int obj_start, const double& timestep);

	virtual ~CollisionDetector();

	void getContinuousTimeCollisions(std::list<ContinuousTimeCollision>& cllsns);

	void getImplicitPenaltyCollisions(std::vector<EdgeEdgeProximityCollision>& edge_edge_collisions,
			std::vector<VertexFaceProximityCollision>& vertex_face_collisions);
	void updateContinuousTimeCollisions();

private:
	// Compute the collisions that happen during the time step between elements in the subtrees node_a and node_b
	void computeContinuousTimeCollisions(const BVHNode& node_a, const BVHNode& node_b);

	// Update the BVH tree starting from node, taking into account the evolution during the time step,
	// i.e. insert m_geodata.m_points+m_time_step*m_geodata.m_velocities.
	void updateBoundingBox(BVHNode& node);

	// Compute the collisions that happen during the time step between elements in the leaves node_a and node_b
	void intersectContent(const BVHNode& node_a, const BVHNode& node_b);

	void intersectContentSelf(const BVHNode& node);

	// Determine if the collision happens; if so append the CTC to m_collisions.
	void appendContinuousTimeIntersection(const TopologicalElement* obj_a, const TopologicalElement* obj_b);
	void appendContinuousTimeIntersection(const YAEdge* edge_a, const YAEdge* edge_b);
	void appendContinuousTimeIntersection(const YAEdge* edge, const YATriangle* triangle);
	void appendContinuousTimeIntersection(const YATriangle* triangle, const YAEdge* edge);
	void appendContinuousTimeIntersection(const YATriangle* triangle_a, const YATriangle* triangle_b);
	void appendContinuousTimeIntersection(int v_index, const YATriangle* triangle);
	// Do the actual collision computation.
	bool analyseCollision(EdgeEdgeContinuousTimeCollision& edgeXedge);
	bool analyseCollision(VertexFaceContinuousTimeCollision& vertexXface);

	bool isVertexFixed(int vert_idx) const;
	bool isRodVertex(int vert) const;
	Vec3d computeRelativeVelocity(const int& idxa0, const int& idxa1, const int& idxb0, const int& idxb1, const double& s,
			const double& t);
	Vec3d computeRelativeVelocity(const int& vrtidx, const int& fcidx0, const int& fcidx1, const int& fcidx2, const double& u,
			const double& v, const double& w);
};
}

#endif /* COLLISIONDETECTOR_HH_ */
