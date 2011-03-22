/*
 * Collisions.cc
 *
 *  Created on: 22/03/2011
 *      Author: jaubry
 */

#include "Collision.hh"
#include "CollisionUtils.hh"

namespace BASim
{

/**
 * Class EdgeEdgeCTCollision
 */

bool EdgeEdgeCTCollision::analyseCollision(const GeometricData& geodata, double time_step)
{
    const Vec3d& p0 = geodata.GetPoint(e0_v0);
    const Vec3d& q0 = geodata.GetPoint(e0_v1);
    const Vec3d& p1 = geodata.GetPoint(e1_v0);
    const Vec3d& q1 = geodata.GetPoint(e1_v1);

    if (p0 == p1 || q0 == q1 || q0 == p1 || p0 == q1)
        return false;

    const Vec3d& vp0 = geodata.GetVelocity(e0_v0);
    const Vec3d& vq0 = geodata.GetVelocity(e0_v1);
    const Vec3d& vp1 = geodata.GetVelocity(e1_v0);
    const Vec3d& vq1 = geodata.GetVelocity(e1_v1);

    std::vector<double> times;
    std::vector<double> errors;
    getCoplanarityTimes(p0, q0, p1, q1, p0 + time_step * vp0, q0 + time_step * vq0, p1 + time_step * vp1, q1 + time_step * vq1,
            times, errors);

    // Loop over the coplanarity times until we find a bona fide collision
    for (size_t j = 0; j < times.size(); ++j)
    {
        //   if (!isProperCollisionTime(times[j]))

        double dtime = times[j] * time_step;

        // Determine if the collision actually happens
        const Vec3d& p1col = p0 + dtime * vp0;
        const Vec3d& q1col = q0 + dtime * vq0;
        const Vec3d& p2col = p1 + dtime * vp1;
        const Vec3d& q2col = q1 + dtime * vq1;

        Vec3d c0;
        Vec3d c1;
        double sqrdist = ClosestPtSegmentSegment(p1col, q1col, p2col, q2col, s, t, c0, c1);

        // If, when they are coplanar, the objects are sufficiently close, register a collision
        if (sqrdist < 1.0e-12)
        {
            m_time = times[j];

            // Compute a collision normal at the time of the collision. For a first attempt, take the cross product of the edges.
            m_normal = (q1col - p1col).cross(q2col - p2col);
            double nnorm = m_normal.norm();

            // If the edges happen to be parallel
            if (nnorm == 0.0)
            {
                std::cerr << "WARNING, FALLING BACK ON OTHER CODE" << std::endl;
                // Use the pre-timestep positions of the collision points to generate a collision normal
                m_normal = ((1.0 - t) * p1 + t * q1).cross((1.0 - s) * p0 + s * q0);
                nnorm = m_normal.norm();

                // Crazy corner case! Handle it if it comes up :)
                if (nnorm == 0.0)
                {
                    bool lazy_person_implemented_this_normal_computation = false;
                    if (!lazy_person_implemented_this_normal_computation)
                        std::cerr << "YOU HIT AN UNSOPPORTED CODE PATH. REALLY DEGENERATE." << std::endl;
                    assert(lazy_person_implemented_this_normal_computation);
                    exit(1);
                }
            }
            m_normal /= nnorm;

            // Make sure the normal produces a negative relative velocity
            if (computeRelativeVelocity(geodata) > 0.0)
                m_normal = -m_normal;

            return true;
        }
    }
    return false;
}

double EdgeEdgeCTCollision::computeRelativeVelocity(const GeometricData& geodata) const
{
    const Vec3d& v0 = geodata.GetVelocity(e0_v0);
    const Vec3d& v1 = geodata.GetVelocity(e0_v1);
    const Vec3d& v2 = geodata.GetVelocity(e1_v0);
    const Vec3d& v3 = geodata.GetVelocity(e1_v1);

    return (((1.0 - t) * v2 + t * v3) - ((1.0 - s) * v0 + s * v1)).dot(m_normal);
}

Vec3d EdgeEdgeCTCollision::computeInelasticImpulse(const GeometricData& geodata, const double& relvel)
{
    double ma0 = geodata.GetMass(e0_v0);
    double ma1 = geodata.GetMass(e0_v1);
    double mb0 = geodata.GetMass(e1_v0);
    double mb1 = geodata.GetMass(e1_v1);

    // Assumes negative relative velocity
    Vec3d numerator = -relvel * m_normal;
    double denominator = (1 - s) * (1 - s) / ma0 + s * s / ma1 + (1 - t) * (1 - t) / mb0 + t * t / mb1;

    return numerator / denominator;
}

std::ostream& operator<<(std::ostream& os, const EdgeEdgeCTCollision& eecol)
{
    os << "Edge 0: " << eecol.e0_v0 << ' ' << eecol.e0_v1 << '\n';
    os << "Edge 1: " << eecol.e1_v0 << ' ' << eecol.e1_v1 << '\n';
    os << "Time: " << eecol.m_time << '\n';
    os << "Normal: " << eecol.m_normal << '\n';
    os << "Barycentric coordinates: " << eecol.s << ' ' << eecol.t;

    return os;
}

/**
 * Class VertexFaceCTCollision
 */

bool VertexFaceCTCollision::analyseCollision(const GeometricData& geodata, double time_step)
{
    const Vec3d& p = geodata.GetPoint(v0);
    const Vec3d& pf0 = geodata.GetPoint(f0);
    const Vec3d& pf1 = geodata.GetPoint(f1);
    const Vec3d& pf2 = geodata.GetPoint(f2);

    const Vec3d& vp = geodata.GetVelocity(v0);
    const Vec3d& vf0 = geodata.GetVelocity(f0);
    const Vec3d& vf1 = geodata.GetVelocity(f1);
    const Vec3d& vf2 = geodata.GetVelocity(f2);

    std::vector<double> times;
    std::vector<double> errors;
    getCoplanarityTimes(p, pf0, pf1, pf2, p + time_step * vp, pf0 + time_step * vf0, pf1 + time_step * vf1,
            pf2 + time_step * vf2, times, errors);
    assert(times.size() == errors.size());

    for (size_t j = 0; j < times.size(); ++j)
    {
        //   if (!isProperCollisionTime(times[j]))

        double dtime = times[j] * time_step;

        // TODO: Use barycentric coordinates or point-triangle closest point < epsilon here? closest point < epsilon really just extends the triangle a little bit.
        // Determine if the collision actually happens
        const Vec3d& pcol = p + dtime * vp;
        const Vec3d& f0col = pf0 + dtime * vf0;
        const Vec3d& f1col = pf1 + dtime * vf1;
        const Vec3d& f2col = pf2 + dtime * vf2;

        Vec3d cp = ClosestPtPointTriangle(pcol, f0col, f1col, f2col);

        // If, when they are coplanar, the objects are sufficiently close, register a collision
        if ((pcol - cp).squaredNorm() < 1.0e-12)
        {
            m_time = times[j];

            // Compute a collision normal at the time of the collision. For a first attempt,
            // take the cross product of two face edges.
            m_normal = (f2col - f0col).cross(f1col - f0col);
            double nnorm = m_normal.norm();

            // If the normal is of magnitude 0 (triangle degenerates into a line)
            if (nnorm == 0.0)
            {
                bool lazy_person_implemented_this_normal_computation = false;
                if (!lazy_person_implemented_this_normal_computation)
                    std::cerr
                            << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m YOU HIT AN UNSOPPORTED CODE PATH. REALLY DEGENERATE."
                            << std::endl;
                assert(lazy_person_implemented_this_normal_computation);
                exit(1);
            }

            m_normal /= nnorm;

            Barycentric(f0col, f1col, f2col, pcol, u, v, w);
            // Barycentric coords could be outside of [0,1] right now because we've extended the triangles a little bit
            assert(approxEq(u + v + w, 1.0));

            if (computeRelativeVelocity(geodata) > 0.0)
                m_normal = -m_normal;

            return true;
        }
    }
    return false;
}

double VertexFaceCTCollision::computeRelativeVelocity(const GeometricData& geodata) const
{
    const Vec3d& vp = geodata.GetVelocity(v0);
    const Vec3d& vt0 = geodata.GetVelocity(f0);
    const Vec3d& vt1 = geodata.GetVelocity(f1);
    const Vec3d& vt2 = geodata.GetVelocity(f2);

    return (vp - (u * vt0 + v * vt1 + w * vt2)).dot(m_normal);
}

Vec3d VertexFaceCTCollision::computeInelasticImpulse(const GeometricData& geodata, const double& relvel)
{
    double mvrt = geodata.GetMass(v0);
    double mfc0 = geodata.GetMass(f0);
    double mfc1 = geodata.GetMass(f1);
    double mfc2 = geodata.GetMass(f2);

    Vec3d numerator = -relvel * m_normal;
    double denominator = 1 / mvrt + u * u / mfc0 + v * v / mfc1 + w * w / mfc2;

    return numerator / denominator;
}

std::ostream& operator<<(std::ostream& os, const VertexFaceCTCollision& vfcol)
{
    os << "Vertex: " << vfcol.f0 << '\n';
    os << "Face: " << vfcol.f0 << ' ' << vfcol.f1 << ' ' << vfcol.f2 << '\n';
    os << "Time: " << vfcol.m_time << '\n';
    os << "Normal: " << vfcol.m_normal << '\n';
    os << "Barycentric coordinates: " << vfcol.u << ' ' << vfcol.v << ' ' << vfcol.w;

    return os;
}

}
