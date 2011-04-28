/*
 * Collision.hh
 *
 *  Created on: 22/03/2011
 *      Author: jaubry
 */

#ifndef COLLISION_HH_
#define COLLISION_HH_

#include "../Core/Definitions.hh"
#include "Geometry.hh"
#include <iostream>

namespace BASim
{

class Collision
{
public:
    Collision(const GeometricData& geodata) :
        m_geodata(geodata)
    {
    }

    virtual bool analyseCollision(double time_step = 0) = 0;

protected:
    // Link to the actual geometric data (point coordinates, velocities etc.)
    const GeometricData& m_geodata;
    // Flag indicating whether the collision has been analysed and the previous variables properly set.
    bool m_analysed;

};

class EdgeFaceIntersection: public Collision
{
public:
    EdgeFaceIntersection(const GeometricData& geodata, const YAEdge* edge, const YATriangle* triangle) :
        Collision(geodata)
    {
        v0 = edge->first();
        v1 = edge->second();
        f0 = triangle->first();
        f1 = triangle->second();
        f2 = triangle->third();
    }

    virtual bool analyseCollision(double time_step = 0);

    // Edge vertices
    int v0;
    int v1;
    // Face vertices
    int f0;
    int f1;
    int f2;

    // Barycentric coordinates of the intersection point on the vertex
    double s;
    // Barycentric coordinates of the intersection point on the triangle
    double u;
    double v;
    double w;

};

class ProximityCollision: public Collision
{
public:
    ProximityCollision(const GeometricData& geodata) :
        Collision(geodata)
    {
    }
    virtual bool analyseCollision(double time_step = 0) = 0;

    // Penetration depth (sum of radii - minimum distance)
    double pen;
    // Collision normal
    Vec3d n;
};

/**
 * Struct to store information needed to resolve a "proximity" collision between two edges.
 */
class EdgeEdgeProximityCollision: public ProximityCollision
{
public:
    EdgeEdgeProximityCollision(const GeometricData& geodata) :
        ProximityCollision(geodata)
    {
    }
    virtual bool analyseCollision(double time_step)
    {
        return false;
    }

    // Vertices involved in collision
    int e0_v0;
    int e0_v1;
    int e1_v0;
    int e1_v1;
    // Radii of edge 0 and edge 1
    double r0;
    double r1;
    // Barycentric coordinates of closest points between edges
    double s;
    double t;
};

/**
 * Struct to store information needed to resolve a "proximity" collision between a vertex and a face.
 */
class VertexFaceProximityCollision: public ProximityCollision
{
public:
    VertexFaceProximityCollision(const GeometricData& geodata) :
        ProximityCollision(geodata)
    {
    }

    virtual bool analyseCollision(double time_step)
    {
        return false;
    }

    // Index of vertex
    int v0;
    // Index of face vertices
    int f0;
    int f1;
    int f2;
    // Radii of vertex and face
    double r0;
    double r1;
    // Barycentric coordinates of closest point on triangle
    double u;
    double v;
    double w;

    // Collision normal points TOWARDS the vertex
};

/**
 * Struct to store information needed to resolve a "implicity penalty" collision between a vertex and a face.
 */
class VertexFaceImplicitPenaltyCollision
{
public:
    // Index of vertex
    int v0;
    // Index of face vertices
    int f0;
    int f1;
    int f2;
    // Radii of vertex and face
    double r0;
    double r1;
    // Extra thickness
    double h;
    // Penalty stifness
    double k;
    // Collision normal, points OUTWARDS the object's face
    Vec3d n;
    // Contact (closest to the vertex) point on the face
    Vec3d cp;
};

// Virtual base class for continuous time collisions.
class CTCollision: public Collision
{
protected:
    // Time of the collision (scaled to be between 0 and 1)
    double m_time;
    // Collision normal
    Vec3d m_normal;
    // Relative velocity along the oriented normal
    double m_relative_velocity;

public:
    CTCollision(const GeometricData& geodata) :
        Collision(geodata)
    {
    }

    virtual ~CTCollision()
    {
    }

    void setTime(double time)
    {
        m_time = time;
    }

    double getTime() const
    {
        assert(m_analysed);
        return m_time;
    }

    Vec3d GetNormal() const
    {
        assert(m_analysed);
        return m_normal;
    }

    double GetRelativeVelocity() const
    {
        assert(m_analysed);
        return m_relative_velocity;
    }

    void ApplyRelativeVelocityKick()
    {
        assert(m_analysed);
        m_relative_velocity -= 1.0e6;
    }

    // From the initial collision data (vertices, velocities and time step) determine whether the collision happened, where and when.
    virtual bool analyseCollision(double time_step) = 0;
    virtual Vec3d computeInelasticImpulse() = 0;
    virtual bool IsFixed() = 0;
    virtual void Print(std::ostream& os) = 0;

    friend bool CompareTimes(const Collision* cllsnA, const Collision* cllsnB);

private:
    virtual double computeRelativeVelocity() const = 0;

};

inline bool CompareTimes(const Collision* cllsnA, const Collision* cllsnB)
{
    const CTCollision* ct_cllsnA = dynamic_cast<const CTCollision*> (cllsnA);
    const CTCollision* ct_cllsnB = dynamic_cast<const CTCollision*> (cllsnB);

    if (ct_cllsnA && ct_cllsnB)
        return ct_cllsnA->m_time < ct_cllsnB->m_time;
    else
        // Should never be there -> throw?
        return false;
}

/**
 * Struct to store information needed to resolve a continuous collision between two edges.
 */
class EdgeEdgeCTCollision: public CTCollision
{
public:
    // Vertices involved in collision
    int e0_v0;
    int e0_v1;
    int e1_v0;
    int e1_v1;
    // Barycentric coordinates of closest points between edges
    double s;
    double t;

    EdgeEdgeCTCollision(const GeometricData& geodata, const YAEdge* edge_a, const YAEdge* edge_b) :
        CTCollision(geodata)
    {
        e0_v0 = edge_a->first();
        e0_v1 = edge_a->second();
        e1_v0 = edge_b->first();
        e1_v1 = edge_b->second();
    }

    virtual ~EdgeEdgeCTCollision()
    {
    }

    bool IsRodRod()
    {
        return m_geodata.isRodVertex(e0_v0) && m_geodata.isRodVertex(e1_v0);
    }

    virtual bool analyseCollision(double time_step);
    virtual Vec3d computeInelasticImpulse();

    virtual bool IsFixed()
    {
        return m_geodata.isVertexFixed(e0_v0) && m_geodata.isVertexFixed(e0_v1) && m_geodata.isVertexFixed(e1_v0)
                && m_geodata.isVertexFixed(e1_v1);
    }

    virtual bool IsCollisionImmune()
    {
        return m_geodata.IsCollisionImmune(e0_v0) || m_geodata.IsCollisionImmune(e0_v1) || m_geodata.IsCollisionImmune(e1_v0)
                || m_geodata.IsCollisionImmune(e1_v1);
    }

    virtual void Print(std::ostream& os)
    {
        os << *this << std::endl;
    }

    friend std::ostream& operator<<(std::ostream& os, const EdgeEdgeCTCollision& eecol);
    //   virtual void Print(std::ostream& os); { os << *this << std::endl; };

private:
    virtual double computeRelativeVelocity() const;

};

/**
 * Struct to store information needed to resolve a continuous time collision between a vertex
 * and a face.
 */
class VertexFaceCTCollision: public CTCollision
{
public:
    // Index of vertex
    int v0;
    // Index of face vertices
    int f0;
    int f1;
    int f2;
    // Barycentric coordinates of closest point on triangle
    double u;
    double v;
    double w;

    VertexFaceCTCollision(const GeometricData& geodata, int v_index, const YATriangle* triangle) :
        CTCollision(geodata)
    {
        v0 = v_index;
        f0 = triangle->first();
        f1 = triangle->second();
        f2 = triangle->third();
    }

    virtual ~VertexFaceCTCollision()
    {
    }

    virtual bool analyseCollision(double time_step);
    virtual Vec3d computeInelasticImpulse();
    virtual bool IsFixed()
    {
        return m_geodata.isVertexFixed(v0) && m_geodata.isVertexFixed(f0) && m_geodata.isVertexFixed(f1)
                && m_geodata.isVertexFixed(f2);
    }
    virtual void Print(std::ostream& os)
    {
        os << *this << std::endl;
    }

    friend std::ostream& operator<<(std::ostream& os, const VertexFaceCTCollision& vfcol);
    //   virtual void Print(std::ostream& os) { os << *this << std::endl; };

private:
    virtual double computeRelativeVelocity() const;

};

// TODO: Move this out of here!
class mycomparison
{
public:
    bool operator()(const CTCollision& cllsnA, const CTCollision& cllsnB) const
    {
        return cllsnA.getTime() < cllsnB.getTime();
    }
};

}

#endif /* COLLISION_HH_ */
