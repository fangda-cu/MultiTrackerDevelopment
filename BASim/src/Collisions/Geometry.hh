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

namespace BASim
{

/**
 * Struct to store face indices.
 */
struct TriangularFace
{
    int idx[3];

    TriangularFace()
    {
    }

    TriangularFace(int un, int deux, int trois)
    {
        idx[0] = un;
        idx[1] = deux;
        idx[2] = trois;
    }

};

// Holds a handle to the actual geometry and related info
class GeometricData
{
    const VecXd& m_points;
    const VecXd& m_velocities;
    const std::vector<double>& m_radii;
    const std::vector<double>& m_masses;
    const std::vector<bool>& m_collision_immune;
    int& m_obj_start;
    const double& m_implicit_thickness;
    const double& m_vertex_face_penalty;

public:
    GeometricData(const VecXd& points, const VecXd& velocities, const std::vector<double>& radii,
            const std::vector<double>& masses, std::vector<bool>& collision_immune, int& obj_start, const double& implicit_thickness,
            const double& vertex_face_penalty) :
        m_points(points), m_velocities(velocities), m_radii(radii), m_masses(masses), m_collision_immune(collision_immune),
                m_obj_start(obj_start), m_implicit_thickness(implicit_thickness), m_vertex_face_penalty(vertex_face_penalty)
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

    double GetImplicitThickness() const
    {
        return m_implicit_thickness;
    }

    double GetVertexFacePenalty() const
    {
        return m_vertex_face_penalty;
    }

    /*
     Vec3d computeRelativeVelocity(const int& idxa0, const int& idxa1, const int& idxb0, const int& idxb1, const double& s,
     const double& t) const
     {
     const Vec3d& v0 = GetVelocity(idxa0);
     const Vec3d& v1 = GetVelocity(idxa1);
     const Vec3d& v2 = GetVelocity(idxb0);
     const Vec3d& v3 = GetVelocity(idxb1);

     return ((1.0 - t) * v2 + t * v3) - ((1.0 - s) * v0 + s * v1);
     }
     */
    /*
     Vec3d computeRelativeVelocity(const int& vrtidx, const int& fcidx0, const int& fcidx1, const int& fcidx2, const double& u,
     const double& v, const double& w) const
     {
     const Vec3d& vp = GetVelocity(vrtidx);
     const Vec3d& vt0 = GetVelocity(fcidx0);
     const Vec3d& vt1 = GetVelocity(fcidx1);
     const Vec3d& vt2 = GetVelocity(fcidx2);

     return vp - (u * vt0 + v * vt1 + w * vt2);
     }
     */

    bool isVertexFixed(int vert_idx) const
    {
        return GetMass(vert_idx) == std::numeric_limits<double>::infinity();
    }

    bool isRodVertex(int vert) const
    {
        // Is a vertex if index is less than start of object vertices in global array
        return vert < GetObjStart();
    }

    bool IsCollisionImmune(int vert) const
    {
        return m_collision_immune[vert];
    }

};

// A virtual class to abstract handling of edges and faces
class TopologicalElement
{
public:
    TopologicalElement()
    {
    }

    virtual ~TopologicalElement()
    {
    }

    // Return the bounding box of the object after it has moved for time_step
    virtual BoundingBox<Scalar> GetBBox(const GeometricData& geodata, const double time_step = 0) const = 0;

    virtual bool IsFixed(const GeometricData& geodata) = 0;
};

class YAEdge: public TopologicalElement
{
    std::pair<int, int> m_edge;

public:
    explicit YAEdge(const std::pair<int, int>& edge) :
        m_edge(edge)
    {
    }

    YAEdge(int vtx0, int vtx1) :
        m_edge(vtx0, vtx1)
    {
    }

    virtual ~YAEdge()
    {
    }

    int first() const
    {
        return m_edge.first;
    }

    int second() const
    {
        return m_edge.second;
    }

    BoundingBox<Scalar> GetBBox(const GeometricData& geodata, const double time_step = 0) const
    {
        // std::cerr << "Computing bounding box for edge " << first() << " " << second() << std::endl;

        const Point<Scalar>& v0 = (Point<Scalar> ) (Vec3d) geodata.GetVelocity(first());
        const Point<Scalar>& v1 = (Point<Scalar> ) (Vec3d) geodata.GetVelocity(second());

        const Point<Scalar>& p0 = (Point<Scalar> ) (Vec3d) geodata.GetPoint(m_edge.first) + v0 * time_step;
        const Point<Scalar>& p1 = (Point<Scalar> ) (Vec3d) geodata.GetPoint(m_edge.second) + v1 * time_step;

        const double r0 = geodata.GetRadius(m_edge.first);
        const double r1 = geodata.GetRadius(m_edge.second);

        BoundingBox<Scalar> bbox(p0.x() - r0, p0.y() - r0, p0.z() - r0, p0.x() + r0, p0.y() + r0, p0.z() + r0);
        bbox.Insert(p1, r1);

        return bbox;
    }

    bool IsFixed(const GeometricData& geodata)
    {
        return geodata.isVertexFixed(m_edge.first) && geodata.isVertexFixed(m_edge.second);
    }

    bool IsFree(const GeometricData& geodata)
    {
        return !geodata.isVertexFixed(m_edge.first) && !geodata.isVertexFixed(m_edge.second);
    }

};

class YATriangle: public TopologicalElement
{
    TriangularFace m_triangle;

public:
    explicit YATriangle(const TriangularFace& triangle) :
        m_triangle(triangle)
    {
    }

    YATriangle(int un, int deux, int trois) :
        m_triangle(un, deux, trois)
    {
    }

    virtual ~YATriangle()
    {
    }

    int first() const
    {
        return m_triangle.idx[0];
    }

    int second() const
    {
        return m_triangle.idx[1];
    }

    int third() const
    {
        return m_triangle.idx[2];
    }

    BoundingBox<Scalar> GetBBox(const GeometricData& geodata, const double time_step = 0) const
    {
        const Point<Scalar>& v0 = (Point<Scalar> ) (Vec3d) geodata.GetVelocity(m_triangle.idx[0]);
        const Point<Scalar>& v1 = (Point<Scalar> ) (Vec3d) geodata.GetVelocity(m_triangle.idx[1]);
        const Point<Scalar>& v2 = (Point<Scalar> ) (Vec3d) geodata.GetVelocity(m_triangle.idx[2]);

        const Point<Scalar>& p0 = (Point<Scalar> ) (Vec3d) geodata.GetPoint(m_triangle.idx[0]) + v0 * time_step;
        const Point<Scalar>& p1 = (Point<Scalar> ) (Vec3d) geodata.GetPoint(m_triangle.idx[1]) + v1 * time_step;
        const Point<Scalar>& p2 = (Point<Scalar> ) (Vec3d) geodata.GetPoint(m_triangle.idx[2]) + v2 * time_step;

        const double r0 = geodata.GetRadius(m_triangle.idx[0]);
        const double r1 = geodata.GetRadius(m_triangle.idx[1]);
        const double r2 = geodata.GetRadius(m_triangle.idx[2]);

        BoundingBox<Scalar> bbox(p0.x() - r0, p0.y() - r0, p0.z() - r0, p0.x() + r0, p0.y() + r0, p0.z() + r0);
        bbox.Insert(p1, r1);
        bbox.Insert(p2, r2);

        return bbox;
    }

    bool IsFixed(const GeometricData& geodata)
    {
        return geodata.isVertexFixed(m_triangle.idx[0]) && geodata.isVertexFixed(m_triangle.idx[1]) && geodata.isVertexFixed(
                m_triangle.idx[2]);
    }

};

class GeometryBBoxFunctor
{
    std::vector<const TopologicalElement*>& m_objects;
    const GeometricData& m_geodata;

public:
    GeometryBBoxFunctor(std::vector<const TopologicalElement*>& objects, const GeometricData& geodata) :
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
