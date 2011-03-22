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

namespace BASim
{

/**
 * Struct to store information needed to resolve a "proximity" collision between two edges.
 */
struct EdgeEdgeProximityCollision
{
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
    // Penetration depth (sum of radii - minimum distance)
    double pen;
    // Collision normal
    Vec3d n;
};

/**
 * Struct to store information needed to resolve a "proximity" collision between a vertex and a face.
 */
struct VertexFaceProximityCollision
{
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
    // Penetration depth (sum of radii - minimum distance)
    double pen;
    // Collision normal, points TOWARDS the vertex
    Vec3d n;
};

/**
 * Struct to store information needed to resolve a "implicity penalty" collision between a vertex and a face.
 */
struct VertexFaceImplicitPenaltyCollision
{
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




class CTCollision
{
    // Time of the collision (scaled to be between 0 and 1)
protected:
    double m_time;

public:
    // Collision normal
    Vec3d n;

public:

    virtual ~CTCollision()
    {
    }

    void setTime(double time)
    {
        m_time = time;
    }

    double getTime() const
    {
        return m_time;
    }

    // From the initial collision data (vertices, velocities and time step) determine whether the collision happened, where and when.
    virtual bool analyseCollision(const GeometricData& geodata, double time_step) = 0;
    virtual bool IsFixed(const GeometricData& geodata) = 0;

    friend bool CompareTimes(const CTCollision* cllsnA, const CTCollision* cllsnB) { return cllsnA->m_time < cllsnB->m_time; }
};

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

    EdgeEdgeCTCollision(const YAEdge* edge_a, const YAEdge* edge_b)
    {
        e0_v0 = edge_a->first();
        e0_v1 = edge_a->second();
        e1_v0 = edge_b->first();
        e1_v1 = edge_b->second();
    }

    virtual bool analyseCollision(const GeometricData& geodata, double time_step);

    bool IsRodRod(const GeometricData& geodata)
    {
        return geodata.isRodVertex(e0_v0) && geodata.isRodVertex(e1_v0);
    }

    bool IsFixed(const GeometricData& geodata)
    {
        return geodata.isVertexFixed(e0_v0) &&  geodata.isVertexFixed(e0_v1) && geodata.isVertexFixed(e1_v0) &&  geodata.isVertexFixed(e1_v1);
    }

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

    VertexFaceCTCollision(int v_index, const YATriangle* triangle)
    {
        v0 = v_index;
        f0 = triangle->first();
        f1 = triangle->second();
        f2 = triangle->third();
    }

    virtual bool analyseCollision(const GeometricData& geodata, double time_step);

    bool IsFixed(const GeometricData& geodata)
    {
        return geodata.isVertexFixed(v0) && geodata.isVertexFixed(f0) && geodata.isVertexFixed(f1) && geodata.isVertexFixed(f2);
    }
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
