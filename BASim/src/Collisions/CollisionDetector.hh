/*
 * CollisionDetector.hh
 *
 *  Created on: 17/03/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef COLLISIONDETECTOR_HH_
#define COLLISIONDETECTOR_HH_

#include "BoundingBox.hh"
#include "Geometry.hh"
#include "Collision.hh"
#include "../Threads/MultithreadedStepper.hh"
#include "BVH.hh"
#include <list>
#include <set>

namespace BASim
{

enum CollisionFilter
{
    ContinuousTime = 0, Proximity, EdgeFace
};

class CollisionDetector
{
    const GeometricData& m_geodata;
    std::vector<const TopologicalElement*> m_elements;
    const double& m_time_step;
    bool m_skip_rod_rod;
    BVH m_bvh;
    std::list<Collision*>* m_collisions;
    CollisionFilter m_collision_filter;
    threads::Mutex m_collisions_mutex;
    int m_num_threads;

public:
    // During construction, the BVH tree is created around the initial geometry.
    // Parameter num_threads = -1 will cause the number of threads to be set equal to the number of available processors.
    CollisionDetector(const GeometricData& geodata, const std::vector<std::pair<int, int> >& edges,
            const std::vector<TriangularFace>& faces, const double& timestep, bool skip_rod_rod = true, int num_threads = -1);

    virtual ~CollisionDetector();

    void getCollisions(std::list<Collision*>& cllsns, CollisionFilter collision_filter);

    void updateContinuousTimeCollisions();

    void skipRodRodCollisions(bool skipRodRodCollisions)
    {
        m_skip_rod_rod = skipRodRodCollisions;
    }

    friend class BVHParallelizer;

private:
    // Update the BVH tree starting from node, taking into account the evolution during the time step,
    // i.e. insert m_geodata.m_points+m_time_step*m_geodata.m_velocities.
    void updateBoundingBox(BVHNode& node);

    // Proximity collision detection
    void computeCollisions(const BVHNode& node_a, const BVHNode& node_b) ;

    // Depending on m_collision_filter, determine and appends the relevant collision type between topological elements to m_collisions
    void appendCollision(const TopologicalElement* obj_a, const TopologicalElement* obj_b);

    // Determine if the collision happens during the current time step; if so append the CTC to m_collisions.
    void appendContinuousTimeCollision(const TopologicalElement* obj_a, const TopologicalElement* obj_b);
    void appendContinuousTimeCollision(const YAEdge* edge_a, const YAEdge* edge_b);
    void appendContinuousTimeCollision(const YAEdge* edge, const YATriangle* triangle);
    void appendContinuousTimeCollision(const YATriangle* triangle_a, const YATriangle* triangle_b);
    void appendContinuousTimeCollision(int v_index, const YATriangle* triangle);

    // Determine whether the edge intersects the triangle; if so append the VFI to m_collisions
    void appendEdgeFaceIntersection(const YAEdge* edge, const YATriangle* triangle);

    bool isVertexFixed(int vert_idx) const;
    bool isRodVertex(int vert) const;
  //Vec3d computeRelativeVelocity(const int& idxa0, const int& idxa1, const int& idxb0, const int& idxb1, const double& s,
  //          const double& t);
  //  Vec3d computeRelativeVelocity(const int& vrtidx, const int& fcidx0, const int& fcidx1, const int& fcidx2, const double& u,
  //          const double& v, const double& w);
};

class BVHParallelizer
{
    const BVHNode& m_node_a;
    const BVHNode& m_node_b;
    CollisionDetector& m_coldet;

public:
    BVHParallelizer(CollisionDetector& coldet, const BVHNode& node_a, const BVHNode& node_b) :
        m_coldet(coldet), m_node_a(node_a), m_node_b(node_b)
    {
    }

    ~BVHParallelizer()
    {
    }

    bool execute() const
    {
        m_coldet.computeCollisions(m_node_a, m_node_b);

        return true;
    }
};

}

#endif /* COLLISIONDETECTOR_HH_ */
