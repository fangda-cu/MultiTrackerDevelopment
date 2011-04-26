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

class IntPair
{
public:

    IntPair(const int& i, const int& j)
    {
        m_i = std::min(i, j);
        m_j = std::max(i, j);
    }

    bool operator==(const IntPair& rhs) const
    {
        assert(m_i <= m_j);
        assert(rhs.m_i <= rhs.m_j);

        return (m_i == rhs.m_i && m_j == rhs.m_j);
    }

    bool operator!=(const IntPair& rhs) const
    {
        assert(m_i <= m_j);
        assert(rhs.m_i <= rhs.m_j);

        return !(*this == rhs);
    }

    bool operator<(const IntPair& rhs) const
    {
        assert(m_i <= m_j);
        assert(rhs.m_i <= rhs.m_j);

        if (m_i != rhs.m_i)
            return m_i < rhs.m_i;
        else
            return m_j < rhs.m_j;
    }

    int first() const
    {
        return m_i;
    }
    int second() const
    {
        return m_j;
    }

private:
    int m_i;
    int m_j;
};

enum CollisionType
{
    ContinuousTime = 0, Proximity, VertexFace
};

class CollisionDetector
{
    const GeometricData& m_geodata;
    std::vector<const TopologicalElement*> m_elements;
    const double& m_time_step;
    bool m_skip_rod_rod;
    BVH m_bvh;
    std::list<Collision*>* m_collisions;
    CollisionType m_ctype;
    //    std::list<ProximityCollision*>* m_prox_collisions;
    threads::Mutex m_collisions_mutex;
    int m_num_threads;

public:
    // During construction, the BVH tree is created around the initial geometry.
    // Parameter num_threads = -1 will cause the number of threads to be set equal to the number of available processors.
    CollisionDetector(const GeometricData& geodata, const std::vector<std::pair<int, int> >& edges,
            const std::vector<TriangularFace>& faces, const double& timestep, bool skip_rod_rod = true, int num_threads = -1);

    virtual ~CollisionDetector();

    void getContinuousTimeCollisions(std::list<Collision*>& cllsns);

    void getImplicitPenaltyCollisions(std::list<Collision*>& cllsns);

    void getVertexFaceIntersections(std::list<Collision*>& cllsns);

    void getCollisions(std::list<Collision*>& cllsns, CollisionType ctype);

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

    void appendVertexFaceIntersection(const YAEdge* edge, const YATriangle* triangle);

    void appendIntersection(const TopologicalElement* obj_a, const TopologicalElement* obj_b) ;

    // Compute the collisions that happen during the time step between elements in the subtrees node_a and node_b
//    void computeContinuousTimeCollisions(const BVHNode& node_a, const BVHNode& node_b);

    // Determine if the collision happens; if so append the CTC to m_collisions.
    void appendContinuousTimeIntersection(const TopologicalElement* obj_a, const TopologicalElement* obj_b);
    void appendContinuousTimeIntersection(const YAEdge* edge_a, const YAEdge* edge_b);
    void appendContinuousTimeIntersection(const YAEdge* edge, const YATriangle* triangle);
    void appendContinuousTimeIntersection(const YATriangle* triangle, const YAEdge* edge);
    void appendContinuousTimeIntersection(const YATriangle* triangle_a, const YATriangle* triangle_b);
    void appendContinuousTimeIntersection(int v_index, const YATriangle* triangle);

    bool isVertexFixed(int vert_idx) const;
    bool isRodVertex(int vert) const;
    Vec3d computeRelativeVelocity(const int& idxa0, const int& idxa1, const int& idxb0, const int& idxb1, const double& s,
            const double& t);
    Vec3d computeRelativeVelocity(const int& vrtidx, const int& fcidx0, const int& fcidx1, const int& fcidx2, const double& u,
            const double& v, const double& w);
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
