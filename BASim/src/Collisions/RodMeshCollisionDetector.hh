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

    void getCollisions(std::list<Collision*>& cllsns, CollisionFilter collision_filter, bool update_mesh_bbox = true);

    void buildBVH()
    {
        build_rod_BVH();
        build_mesh_BVH();
    }

    void build_rod_BVH();

    void build_mesh_BVH();

    void rebuildRodElements(const std::vector<std::pair<int, int> >& edges);

private:
    void computeCollisions(const BVHNodeType& node_a, const BVHNodeType& node_b);
};

}
#endif /* RODMESHCOLLISIONDETECTOR_HH_ */
