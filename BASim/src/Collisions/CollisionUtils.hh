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

void getIntersectionPoint(const Vec3d& x0,  const Vec3d& xnew0,
        const Vec3d& xnew1, const Vec3d& xnew2, const Vec3d& xnew3, std::vector<double>& times, std::vector<double>& errors);

}

#endif

// TODO:
//   o It would be nice to handle degenerate cases better in these methods.  
//     They all handle degenerate cases, but getting PREDICTABLE behavior out would rock!!!
