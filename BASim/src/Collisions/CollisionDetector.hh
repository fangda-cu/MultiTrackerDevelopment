/*
 * CollisionDetector.hh
 *
 *  Created on: 17/03/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef COLLISIONDETECTOR_HH_
#define COLLISIONDETECTOR_HH_

#include "CollisionDetectorBase.hh"

namespace BASim
{

class CollisionDetector: public CollisionDetectorBase
{
    std::vector<const TopologicalElement*> m_elements;
    BVH m_bvh;

public:
    // During construction, the BVH tree is created around the initial geometry.
    // Parameter num_threads = -1 will cause the number of threads to be set equal to the number of available processors.
    CollisionDetector(const GeometricData& geodata, const std::vector<std::pair<int, int> >& edges,
            const std::vector<TriangularFace>& faces, const double& timestep, bool skip_rod_rod = true, int num_threads = -1);

    virtual ~CollisionDetector();

    void getCollisions(std::list<Collision*>& cllsns, CollisionFilter collision_filter);

    void buildBVH();

private:
    void computeCollisions(const BVHNode& node_a, const BVHNode& node_b);
};

}

#endif /* COLLISIONDETECTOR_HH_ */
