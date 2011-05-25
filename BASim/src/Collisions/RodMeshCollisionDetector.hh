/*
 * RodMeshCollisionDetector.hh
 *
 *  Created on: 25/05/2011
 *      Author: jaubry
 */

#ifndef RODMESHCOLLISIONDETECTOR_HH_
#define RODMESHCOLLISIONDETECTOR_HH_

#include "CollisionDetector.hh"

namespace BASim
{

class RodMeshCollisionDetector: public CollisionDetectorBase
{
    std::vector<const TopologicalElement*> m_rod_elements;
    std::vector<const TopologicalElement*> m_mesh_elements;

    BVH m_rod_bvh;
    BVH m_mesh_bvh;

public:
    RodMeshCollisionDetector(const GeometricData& geodata, const std::vector<std::pair<int, int> >& edges,
            const std::vector<TriangularFace>& faces, const double& timestep, bool skip_rod_rod = true, int num_threads = -1);

    virtual ~RodMeshCollisionDetector();

    void getCollisions(std::list<Collision*>& cllsns, CollisionFilter collision_filter);

    void buildBVH()
    {
        build_rod_BVH();
        build_mesh_BVH();
    }

    void build_rod_BVH();

    void build_mesh_BVH();

    friend class BVHParallelizer;

private:
    void computeCollisions(const BVHNode& node_a, const BVHNode& node_b);

};

}
#endif /* RODMESHCOLLISIONDETECTOR_HH_ */
