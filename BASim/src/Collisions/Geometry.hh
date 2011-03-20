/*
 * Geometry.hh
 *
 *  Created on: 17/03/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef GEOMETRY_HH_
#define GEOMETRY_HH_

#include "../Core/Definitions.hh"
#include "BoundingBox.hh"
#include "BVHAABB.hh"

namespace BASim
{

// Holds the actual geometry
class GeometricData
{
	const VecXd& m_points;
	const VecXd& m_velocities;
	const std::vector<double>& m_radii;
	const std::vector<double>& m_masses;
	int m_obj_start;

public:
	GeometricData(const VecXd& points, const VecXd& velocities, const std::vector<double>& radii,
			const std::vector<double>& masses, int obj_start) :
		m_points(points), m_velocities(velocities), m_radii(radii), m_masses(masses), m_obj_start(obj_start)
	{
	}

	Vec3d GetPoint(int i) const
	{
		return m_points.segment<3> (3 * i);
	}

	Vec3d GetVelocity(int i) const
	{
		return m_velocities.segment<3> (3 * i);
	}

	double GetRadius(int i) const
	{
		return m_radii[i];
	}

	double GetMass(int i) const
	{
		return m_masses[i];
	}
	int GetObjStart() const
	{
		return m_obj_start;
	}

};

// A virtual class to abstract handling of edges and faces
class TopologicalElement
{

public:
	// Return the bounding box of the object after it has moved for time_step
	virtual BoundingBox<Scalar> GetBBox(const GeometricData& geodata, const double time_step = 0) = 0;

};

class YAEdge: public TopologicalElement
{
public:
	std::pair<int, int> m_edge;

	explicit YAEdge(std::pair<int, int> edge) :
		m_edge(edge)
	{
	}

	YAEdge(int vtx0, int vtx1) :
		m_edge(vtx0, vtx1)
	{
	}

	BoundingBox<Scalar> GetBBox(const GeometricData& geodata, const double time_step = 0);

};

class YATriangle: public TopologicalElement
{
public:
	TriangularFace m_triangle;

	explicit YATriangle(TriangularFace triangle) :
		m_triangle(triangle)
	{
	}

	BoundingBox<Scalar> GetBBox(const GeometricData& geodata, const double time_step = 0);

};

class GeometryBBoxFunctor
{
	std::vector<TopologicalElement*>& m_objects;
	const GeometricData& m_geodata;

public:
	GeometryBBoxFunctor(std::vector<TopologicalElement*>& objects, const GeometricData& geodata) :
		m_objects(objects), m_geodata(geodata)
	{
	}

	uint32_t size() const
	{
		return (uint32_t) m_objects.size();
	}

	BoundingBox<Scalar> operator[](const uint32_t i)
	{
		return m_objects[i]->GetBBox(m_geodata);
	}

	void swap(uint32_t i, uint32_t j)
	{
		std::swap(m_objects[i], m_objects[j]);
	}

};
}

#endif /* GEOMETRY_HH_ */
