/*
 * BoundingBox.hh
 *
 *  Created on: 17/03/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef BOUNDINGBOX_HH_
#define BOUNDINGBOX_HH_

#include <limits>
#include <algorithm>
#include <iostream>
#include "../Core/Definitions.hh"

namespace BASim
{

// typedef unsigned int uint;

template<typename ScalarT>
Eigen::Matrix<ScalarT, 3, 1> min(const Eigen::Matrix<ScalarT, 3, 1>& a, const Eigen::Matrix<ScalarT, 3, 1>& b)
{
    return Eigen::Matrix<ScalarT, 3, 1>(std::min(a[0], b[0]), std::min(a[1], b[1]), std::min(a[2], b[2]));
}

template<typename ScalarT>
Eigen::Matrix<ScalarT, 3, 1> max(const Eigen::Matrix<ScalarT, 3, 1>& a, const Eigen::Matrix<ScalarT, 3, 1>& b)
{
    return Eigen::Matrix<ScalarT, 3, 1>(std::max(a[0], b[0]), std::max(a[1], b[1]), std::max(a[2], b[2]));
}

template<typename ScalarT>
struct BoundingBox
{
    typedef Eigen::Matrix<ScalarT, 3, 1> PointType;
    PointType min;
    PointType max;

    BoundingBox() :
                min(std::numeric_limits<ScalarT>::max(), std::numeric_limits<ScalarT>::max(),
                        std::numeric_limits<ScalarT>::max()),
                max(-std::numeric_limits<ScalarT>::max(), -std::numeric_limits<ScalarT>::max(),
                        -std::numeric_limits<ScalarT>::max())
    {
    }

    explicit BoundingBox(PointType point) :
        min(point), max(point)
    {
    }

    BoundingBox(ScalarT minx, ScalarT miny, ScalarT minz, ScalarT maxx, ScalarT maxy, ScalarT maxz) :
        min(minx, miny, minz), max(maxx, maxy, maxz)
    {
    }

    bool IsValid() const
    {
        return max[0] >= min[0] && max[1] >= min[1] && max[2] >= min[2];
    }

    void Reset()
    {
        *this = BoundingBox();
    }

    // Expand the box to contain the given point
    inline void Insert(const PointType& p)
    {
        min[0] = std::min(min[0], p[0]);
        min[1] = std::min(min[1], p[1]);
        min[2] = std::min(min[2], p[2]);
        max[0] = std::max(max[0], p[0]);
        max[1] = std::max(max[1], p[1]);
        max[2] = std::max(max[2], p[2]);
    }

    // Expand the box to contain the given sphere
    inline void Insert(const PointType& p, const ScalarT& radius)
    {
        min[0] = std::min(min[0], p[0] - radius);
        min[1] = std::min(min[1], p[1] - radius);
        min[2] = std::min(min[2], p[2] - radius);
        max[0] = std::max(max[0], p[0] + radius);
        max[1] = std::max(max[1], p[1] + radius);
        max[2] = std::max(max[2], p[2] + radius);
    }

    // Expand the box to contain the given box
    inline void Insert(const BoundingBox& box)
    {
        Insert(box.min);
        Insert(box.max);
    }

    ScalarT Volume() const
    {
        return (max[0] - min[0]) * (max[1] - min[1]) * (max[2] - min[2]);
    }
};

template<typename ScalarT>
bool Intersect(const BoundingBox<ScalarT>& bbox_a, const BoundingBox<ScalarT>& bbox_b)
{
    if ((bbox_a.max[0] < bbox_b.min[0]) || (bbox_b.max[0] < bbox_a.min[0]) || (bbox_a.max[1] < bbox_b.min[1]) || (bbox_b.max[1]
            < bbox_a.min[1]) || (bbox_a.max[2] < bbox_b.min[2]) || (bbox_b.max[2] < bbox_a.min[2]))
        return false;
    else
        return true;
}

}

#endif /* BOUNDINGBOX_HH_ */
