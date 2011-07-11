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

inline Vec3d parallelTransport(const Vec3d& u, const Vec3d& t1, const Vec3d& t2)
{
    assert(fabs(t1.norm() - 1) < std::numeric_limits<Scalar>::epsilon());
    assert(fabs(t2.norm() - 1) < std::numeric_limits<Scalar>::epsilon());

    Vec3d b = t1.cross(t2);
    if (b.norm() < std::numeric_limits<Scalar>::epsilon()) // vectors are colinear
        return u;

    b.normalize();
    const Vec3d n1 = t1.cross(b);
    const Vec3d n2 = t2.cross(b);
    return u.dot(t1) * t2 + u.dot(n1) * n2 + u.dot(b) * b;
}

}
#endif /* ELASTICSTRANDUTILS_HH_ */
