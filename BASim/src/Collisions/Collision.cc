/*
 * Collisions.cc
 *
 *  Created on: 22/03/2011
 *      Author: jaubry
 */

#include "Collision.hh"
#include "CollisionUtils.hh"
#include "../Util/TextLog.hh"

namespace BASim
{
static const double SQ_TOLERANCE = 1e-12;

/**
 * Class EdgeFaceIntersection
 */
bool EdgeFaceIntersection::analyseCollision(double)
{
    const Vec3d offset = m_geodata.GetPoint(v0);

    const Vec3d p0 = m_geodata.GetPoint(v0) - offset;
    const Vec3d p1 = m_geodata.GetPoint(v1) - offset;
    const Vec3d q0 = m_geodata.GetPoint(f0) - offset;
    const Vec3d q1 = m_geodata.GetPoint(f1) - offset;
    const Vec3d q2 = m_geodata.GetPoint(f2) - offset;

    std::vector<double> times;
    std::vector<double> errors;
    getIntersectionPoint(p0, p1, q0, q1, q2, times, errors);
    assert(times.size() == errors.size());

    for (size_t j = 0; j < times.size(); ++j)
    {
        double dtime = times[j] * 1.0;

        // Determine if the collision actually happens
        const Vec3d& pcol = p0 + dtime * (p1 - p0);
        const Vec3d& f0col = q0;
        const Vec3d& f1col = q1;
        const Vec3d& f2col = q2;

        Vec3d cp = ClosestPtPointTriangle(pcol, f0col, f1col, f2col);

        // If the intersection point and the close point on the triangle are close and that point is interior, register a collision
        if ((pcol - cp).squaredNorm() <= SQ_TOLERANCE)
        {
            s = times[j];

            Barycentric(f0col, f1col, f2col, pcol, u, v, w);
            // Barycentric coords could be outside of [0,1] right now because we've extended the triangles a little bit
            assert(approxEq(u + v + w, 1.0));

            if (u > 0 && v > 0 && w > 0)
                return m_analysed = true;
        }
    }
    return false;
}

/**
 * Class EdgeEdgeProximityCollision
 */
bool EdgeEdgeProximityCollision::analyseCollision(double)
{
    return false;
}

double EdgeEdgeProximityCollision::computeRelativeVelocity() const // Assumes m_normal, s and t have been computed
{
    const Vec3d& v0 = m_geodata.GetVelocity(e0_v0);
    const Vec3d& v1 = m_geodata.GetVelocity(e0_v1);
    const Vec3d& v2 = m_geodata.GetVelocity(e1_v0);
    const Vec3d& v3 = m_geodata.GetVelocity(e1_v1);

    return (((1.0 - t) * v2 + t * v3) - ((1.0 - s) * v0 + s * v1)).dot(m_normal);
}

Vec3d EdgeEdgeProximityCollision::computeInelasticImpulse()
{
    assert(m_analysed);

    double ma0 = m_geodata.GetMass(e0_v0);
    double ma1 = m_geodata.GetMass(e0_v1);
    double mb0 = m_geodata.GetMass(e1_v0);
    double mb1 = m_geodata.GetMass(e1_v1);

    // Assumes negative relative velocity
    Vec3d numerator = -m_relative_velocity * m_normal;
    double denominator = (1 - s) * (1 - s) / ma0 + s * s / ma1 + (1 - t) * (1 - t) / mb0 + t * t / mb1;

    return numerator / denominator;
}

std::ostream& operator<<(std::ostream& os, const EdgeEdgeProximityCollision& eecol)
{
    os << "Edge edge collision!\n";
    os << "Normal: " << eecol.m_normal << '\n';
    os << "Relative velocity: " << eecol.m_relative_velocity << '\n';
    os << "Edge 0: " << eecol.e0_v0 << eecol.m_geodata.GetPoint(eecol.e0_v0) << ' ' << eecol.e0_v1 << eecol.m_geodata.GetPoint(
            eecol.e0_v1) << '\n';
    os << "Edge 1: " << eecol.e1_v0 << eecol.m_geodata.GetPoint(eecol.e1_v0) << ' ' << eecol.e1_v1 << eecol.m_geodata.GetPoint(
            eecol.e1_v1) << '\n';
    os << "Barycentric coordinates: " << eecol.s << ' ' << eecol.t;

    return os;
}

/**
 * Class VertexFaceProximityCollision
 */
bool VertexFaceProximityCollision::analyseCollision(double)
{
    assert(v0 != f0 && v0 != f1 && v0 != f2);

    // TODO: Add check for both having fixed vertices

    const Vec3d p1 = m_geodata.GetPoint(v0);
    const Vec3d t0 = m_geodata.GetPoint(f0);
    const Vec3d t1 = m_geodata.GetPoint(f1);
    const Vec3d t2 = m_geodata.GetPoint(f2);

    cp = ClosestPtPointTriangle(p1, t0, t1, t2);
    double sqrdist = (p1 - cp).squaredNorm();

    if (sqrdist < (m_geodata.GetImplicitThickness() + r0 + r1) * (m_geodata.GetImplicitThickness() + r0 + r1))
    {
        k = m_geodata.GetVertexFacePenalty();
        h = m_geodata.GetImplicitThickness();

        m_normal = (t1 - t0).cross(t2 - t0);
        assert(m_normal.norm() > 0.0);

        double nnorm = m_normal.norm();
        m_normal /= nnorm;

        assert(fabs(m_normal.norm() - 1.0) < 1.0e-6);

        if (nnorm == 0.0 || fabs(m_normal.norm() - 1.0) > 1.0e-6)
        {
            std::cerr << "WARNING, IGNORING COLLISION DUE TO DEGENERATE NORMAL" << std::endl;
            return false;
        }

        return m_analysed = true;
    }

    return false;
}

std::ostream& operator<<(std::ostream& os, const VertexFaceProximityCollision& vfcol)
{
    os << "Vertex face proximity!\n";
    os << "Vertex: " << vfcol.m_geodata.GetPoint(vfcol.v0) << '\n';
    os << "Close point: " << vfcol.cp << '\n';
    os << "Normal: " << vfcol.m_normal << '\n';

    return os;
}

}
