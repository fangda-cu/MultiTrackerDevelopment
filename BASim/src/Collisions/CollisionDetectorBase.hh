/*
 * CollisionDetectorBase.hh
 *
 *  Created on: 25/05/2011
 *      Author: jaubry
 */

#ifndef COLLISIONDETECTORBASE_HH_
#define COLLISIONDETECTORBASE_HH_

#include "Geometry.hh"
#include "Collision.hh"
#include "CTCollision.hh"
#include "BoundingBox.hh"
#include "BVH.hh"
#include <list>
#include <set>
#include "../Threads/MultithreadedStepper.hh"

namespace BASim
{

enum CollisionFilter
{
    ContinuousTime = 0, Proximity, EdgeFace
};

class CollisionDetectorBase
{
protected:
    const GeometricData& m_geodata;
    const double& m_time_step;
    bool m_skip_rod_rod;
    std::list<Collision*>* m_collisions_list;
    CollisionFilter m_collision_filter;
    threads::Mutex m_collisions_mutex;
    int m_num_threads;

public:
    int m_potential_collisions;

    CollisionDetectorBase(const GeometricData& geodata, const std::vector<std::pair<int, int> >& edges,
            const std::vector<TriangularFace>& faces, const double& timestep, bool skip_rod_rod = true, int num_threads = -1);

    virtual ~CollisionDetectorBase();

    virtual void getCollisions(std::list<Collision*>& cllsns, CollisionFilter collision_filter, bool) = 0;

    virtual void buildBVH() = 0;

    void updateContinuousTimeCollisions();

    void setSkipRodRodCollisions(bool skipRodRodCollisions)
    {
        m_skip_rod_rod = skipRodRodCollisions;
    }

protected:
    void getReady(std::list<Collision*>& cllsns, CollisionFilter collision_filter);

    // Update the BVH tree starting from node, taking into account the evolution during the time step,
    // i.e. insert m_geodata.m_points+m_time_step*m_geodata.m_velocities.
    void updateBoundingBox(BVH& bvh, const std::vector<const TopologicalElement*>& elements, BVHNode& node);

    // Collision detection
    virtual void computeCollisions(const BVHNode& node_a, const BVHNode& node_b) = 0;

    // Depending on m_collision_filter, determine and appends the relevant collision type between topological elements to m_collisions_list
    bool appendCollision(const TopologicalElement* obj_a, const TopologicalElement* obj_b);

    // Determine if the collision happens during the current time step; if so append the CTC to m_collisions_list.
    bool appendContinuousTimeCollision(const YAEdge* edge_a, const YAEdge* edge_b);
    bool appendContinuousTimeCollision(const YAEdge* edge, const YATriangle* triangle);
    bool appendContinuousTimeCollision(int v_index, const YATriangle* triangle);

    // Determine whether the edge intersects the triangle; if so append the VFI to m_collisions_list
    bool appendEdgeFaceIntersection(const YAEdge* edge, const YATriangle* triangle);

    // Determine if a close encounter has happened
    bool appendProximityCollision(const YAEdge* edge_a, const YAEdge* edge_b);
    bool appendProximityCollision(const YAEdge* edge, const YATriangle* triangle);
    bool appendProximityCollision(int v_index, const YATriangle* triangle);

    bool isVertexFixed(int vert_idx) const;
    bool isRodVertex(int vert) const;

    class BVHParallelizer
    {
        const BVHNode& m_node_a;
        const BVHNode& m_node_b;
        CollisionDetectorBase* m_coldet;

    public:
        BVHParallelizer(CollisionDetectorBase* coldet, const BVHNode& node_a, const BVHNode& node_b) :
            m_coldet(coldet), m_node_a(node_a), m_node_b(node_b)
        {
            // std::cerr << "Constructing BVHParallelizer" << std::endl;
        }

        ~BVHParallelizer()
        {
        }

        bool execute() const
        {
            m_coldet->computeCollisions(m_node_a, m_node_b);

            return true;
        }
    };
};

}

#endif /* COLLISIONDETECTORBASE_HH_ */
