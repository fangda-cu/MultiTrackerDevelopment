/**
 * \file CollisionUtils.hh
 *
 * \author (The Internets + Books)
 * \date 04/17/2010
 */

#ifndef COLLISIONUTILS_HH
#define COLLISIONUTILS_HH

#ifdef WETA
#include "../Core/Definitions.hh"
#include "../Math/Math.hh"
#else
#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Math/Math.hh"
#endif

namespace BASim
{

inline Vec3d ClosestPtPointSegment(const Vec3d& point, const Vec3d& first, const Vec3d& last)
{
    if (approxEq(point, first))
        return first;
    if (approxEq(point, last))
        return last;

    const Scalar vx = last[0] - first[0];
    const Scalar vy = last[1] - first[1];
    const Scalar vz = last[2] - first[2];
    const Scalar inv = 1.0f / (vy * vy + vx * vx + vz * vz);
    const Scalar dtp = vx * point[0] + vy * point[1] + vz * point[2];
    const Scalar flx = first[2] * last[1] - first[1] * last[2];
    const Scalar fly = first[0] * last[2] - first[2] * last[0];
    const Scalar flz = first[1] * last[0] - first[0] * last[1];

    // First project on the line
    Vec3d projonline((vz * fly - vy * flz + vx * dtp) * inv, (vx * flz - vz * flx + vy * dtp) * inv,
            (vy * flx - vx * fly + vz * dtp) * inv);

    // Check if we are outside the interval
    if ((projonline - last).dot(first - last) <= 0)
        return last;
    else if ((projonline - first).dot(last - first) <= 0)
        return first;
    else
        return projonline;

}

// Closest point on a triangle to a vertex
Vec3d ClosestPtPointTriangle(const Vec3d& p, const Vec3d& a, const Vec3d& b, const Vec3d& c);

// Computes the squared distance between and closest points of two edges. 
double ClosestPtSegmentSegment(const Vec3d& p1, const Vec3d& q1, const Vec3d& p2, const Vec3d& q2, double& s, double& t,
        Vec3d& c1, Vec3d& c2);

// Computes the barycentric coordiantes of a point wrt a triangle
void Barycentric(const Vec3d& a, const Vec3d& b, const Vec3d& c, const Vec3d& p, double& u, double& v, double& w);

/////////////////////////////////////////////////////////////////
// Collision detection code adapted from Robert Bridson's website

void addUnique(std::vector<double>& a, double e);

double triple(const Vec3d &a, const Vec3d &b, const Vec3d &c);

double signed_volume(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2, const Vec3d& x3);

void getCoplanarityTimes(const Vec3d& x0, const Vec3d& x1, const Vec3d& x2, const Vec3d& x3, const Vec3d& xnew0,
        const Vec3d& xnew1, const Vec3d& xnew2, const Vec3d& xnew3, std::vector<double>& times, std::vector<double>& errors);

void getIntersectionPoint(const Vec3d& x0, const Vec3d& xnew0, const Vec3d& xnew1, const Vec3d& xnew2, const Vec3d& xnew3,
        std::vector<double>& times, std::vector<double>& errors);

}

#endif

// TODO:
//   o It would be nice to handle degenerate cases better in these methods.  
//     They all handle degenerate cases, but getting PREDICTABLE behavior out would rock!!!
