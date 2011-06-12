/*
 * CTCollision.hh
 *
 *  Created on: 13/06/2011
 *      Author: jaubry
 */

#ifndef CTCOLLISION_HH_
#define CTCOLLISION_HH_

#include "Collision.hh"

namespace BASim
{

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

    double GetCachedRelativeVelocity() const
    {
        assert(m_analysed);
        return m_relative_velocity;
    }

    void ApplyRelativeVelocityKick()
    {
        assert(0); // BAD BAD BAD BAD BAD
        //        assert(m_analysed);
        //  m_relative_velocity -= 1.0e6;
    }

    // From the initial collision data (vertices, velocities and time step) determine whether the collision happened, where and when.
    virtual bool analyseCollision(double time_step) = 0;
    virtual Vec3d computeInelasticImpulse() = 0;
    virtual bool IsFixed() = 0;
    virtual void Print(std::ostream& os) = 0;
    virtual int GetRodVertex() const = 0;

    friend bool CompareTimes(const Collision* cllsnA, const Collision* cllsnB);

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

    virtual int GetRodVertex() const;

    friend std::ostream& operator<<(std::ostream& os, const EdgeEdgeCTCollision& eecol);
    //   virtual void Print(std::ostream& os); { os << *this << std::endl; };

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
    virtual int GetRodVertex() const;

    friend std::ostream& operator<<(std::ostream& os, const VertexFaceCTCollision& vfcol);
    //   virtual void Print(std::ostream& os) { os << *this << std::endl; };

    virtual double computeRelativeVelocity() const;

};

}

#endif /* CTCOLLISION_HH_ */
