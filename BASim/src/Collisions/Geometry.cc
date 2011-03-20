/*
 * Geometry.cc
 *
 *  Created on: 17/03/2011
 *      Author: jaubry
 */

#include "Geometry.hh"

namespace BASim
{

static const double MARGIN = 1.0e-6;

BoundingBox<Scalar> YAEdge::GetBBox(const GeometricData& geodata, const double time_step)
{
	const Point<Scalar>& v0 = (Point<Scalar> ) (Vec3d) geodata.GetVelocity(m_edge.first);
	const Point<Scalar>& v1 = (Point<Scalar> ) (Vec3d) geodata.GetVelocity(m_edge.second);

	const Point<Scalar>& p0 = (Point<Scalar> ) (Vec3d) geodata.GetPoint(m_edge.first) + v0 * time_step;
	const Point<Scalar>& p1 = (Point<Scalar> ) (Vec3d) geodata.GetPoint(m_edge.second) + v1 * time_step;

	const double r0 = geodata.GetRadius(m_edge.first) + MARGIN;
	const double r1 = geodata.GetRadius(m_edge.second) + MARGIN;

	BoundingBox<Scalar> bbox(p0.x() - r0, p0.y() - r0, p0.z() - r0, p0.x() + r0, p0.y() + r0, p0.z() + r0);
	bbox.Insert(p1, r1);

	return bbox;
}

BoundingBox<Scalar> YATriangle::GetBBox(const GeometricData& geodata, const double time_step)
{
	const Point<Scalar>& v0 = (Point<Scalar> ) (Vec3d) geodata.GetVelocity(m_triangle.idx[0]);
	const Point<Scalar>& v1 = (Point<Scalar> ) (Vec3d) geodata.GetVelocity(m_triangle.idx[1]);
	const Point<Scalar>& v2 = (Point<Scalar> ) (Vec3d) geodata.GetVelocity(m_triangle.idx[2]);

	const Point<Scalar>& p0 = (Point<Scalar> ) (Vec3d) geodata.GetPoint(m_triangle.idx[0]) + v0 * time_step;
	const Point<Scalar>& p1 = (Point<Scalar> ) (Vec3d) geodata.GetPoint(m_triangle.idx[1]) + v1 * time_step;
	const Point<Scalar>& p2 = (Point<Scalar> ) (Vec3d) geodata.GetPoint(m_triangle.idx[2]) + v2 * time_step;

	const double r0 = geodata.GetRadius(m_triangle.idx[0]) + MARGIN;
	const double r1 = geodata.GetRadius(m_triangle.idx[1]) + MARGIN;
	const double r2 = geodata.GetRadius(m_triangle.idx[2]) + MARGIN;

	BoundingBox<Scalar> bbox(p0.x() - r0, p0.y() - r0, p0.z() - r0, p0.x() + r0, p0.y() + r0, p0.z() + r0);
	bbox.Insert(p1, r1);
	bbox.Insert(p2, r2);

	return bbox;
}

}
