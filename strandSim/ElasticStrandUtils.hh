/*
 * ElasticStrandUtils.hh
 *
 *  Created on: 11/07/2011
 *      Author: jaubry
 */

#ifndef ELASTICSTRANDUTILS_HH_
#define ELASTICSTRANDUTILS_HH_

#include <limits>
#include "Definitions.hh"

namespace strandsim
{

template<typename NormableT>
inline bool isClose(const NormableT& x1, const NormableT& x2)
{
    return (x1 - x2).norm() < std::numeric_limits<ScalarT>::epsilon();
}

template<typename ScalarT>
inline bool isSmall(ScalarT x)
{
    return fabs(x) < std::numeric_limits<ScalarT>::epsilon();
}

inline Vec3d parallelTransport(const Vec3d& u, const Vec3d& t1, const Vec3d& t2)
{
    assert(isSmall(t1.norm() - 1));
    assert(isSmall(t2.norm() - 1));

    Vec3d b = t1.cross(t2);
    if (isSmall(b.norm())) // vectors are colinear
        return u;

    b.normalize();
    const Vec3d n1 = t1.cross(b);
    const Vec3d n2 = t2.cross(b);
    return u.dot(t1) * t2 + u.dot(n1) * n2 + u.dot(b) * b;
}

inline Vec3d findNormal(const Vec3d& u)
{
    assert(u.norm() != 0);
    Vec3d v;

    int maxCoordinate = 0;
    for (int i = 0; i < 3; ++i)
    {
        if (isSmall(u[i]))
        {
            v[i] = 1;
            goto finished;
        }
        if (fabs(u[i]) > fabs(u[maxCoordinate]))
            maxCoordinate = i;
    }
    {
        const int idx = (maxCoordinate + 1) % 3;
        v[idx] = u[maxCoordinate];
        v[maxCoordinate] = -u[idx];
    }

    finished: assert(isSmall(u.dot(v)));
    return v;
}

}
#endif /* ELASTICSTRANDUTILS_HH_ */
