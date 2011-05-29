/*
 * CollisionDetector.cc
 *
 *  Created on: 17/03/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#include "CollisionDetector.hh"

namespace BASim
{
CollisionDetector::CollisionDetector(const GeometricData& geodata, const std::vector<std::pair<int, int> >& edges,
        const std::vector<TriangularFace>& faces, const double& timestep, bool skip_rod_rod, int num_threads) :
        CollisionDetectorBase(geodata, edges, faces, timestep, skip_rod_rod, num_threads)
{
    // std::cerr << "Constructing CollisionDetector" << std::endl;
    m_elements.reserve(edges.size() + faces.size());
    for (std::vector<std::pair<int, int> >::const_iterator i = edges.begin(); i != edges.end(); i++)
        m_elements.push_back(new YAEdge(*i));
    for (std::vector<TriangularFace>::const_iterator i = faces.begin(); i != faces.end(); i++)
        m_elements.push_back(new YATriangle(*i));

    buildBVH();
}

CollisionDetector::~CollisionDetector()
{
    for (std::vector<const TopologicalElement*>::iterator i = m_elements.begin(); i != m_elements.end(); i++)
        delete *i;
}

void CollisionDetector::getCollisions(std::list<Collision*>& cllsns, CollisionFilter collision_filter, bool)
{
    m_potential_collisions = 0;
    m_collision_filter = collision_filter;
    m_collisions_list = &cllsns;
    m_collisions_list->clear();
    std::vector<BVHParallelizer*> steppers;
    BVHNode& root = m_bvh.GetNode(0);
    updateBoundingBox(m_bvh, m_elements, root);

    if (root.IsLeaf()) // Can't really call this a tree, can we?
    {
        std::cerr << "Root is a leaf!" << std::endl;
        computeCollisions(root, root);
        return;
    }

    BVHNode& h = m_bvh.GetNode(root.ChildIndex());
    BVHNode& g = m_bvh.GetNode(root.ChildIndex() + 1);

    // If tree has depth 1, detect collisions at this level.
    if (h.IsLeaf() || g.IsLeaf())
    {
        steppers.reserve(3);
        steppers.push_back(new BVHParallelizer(this, h, h));
        steppers.push_back(new BVHParallelizer(this, h, g));
        steppers.push_back(new BVHParallelizer(this, g, g));
        MultithreadedStepper<std::vector<BVHParallelizer*> > (steppers, m_num_threads).Execute();
        return;
    }

    BVHNode& hh = m_bvh.GetNode(h.ChildIndex());
    BVHNode& hg = m_bvh.GetNode(h.ChildIndex() + 1);
    BVHNode& gh = m_bvh.GetNode(g.ChildIndex());
    BVHNode& gg = m_bvh.GetNode(g.ChildIndex() + 1);

    // If tree has depth 2, detect collisions at this level.
    if (hh.IsLeaf() || hg.IsLeaf() || gh.IsLeaf() || gg.IsLeaf())
    {
        steppers.reserve(10);
        steppers.push_back(new BVHParallelizer(this, hh, hh));
        steppers.push_back(new BVHParallelizer(this, hh, hg));
        steppers.push_back(new BVHParallelizer(this, hh, gh));
        steppers.push_back(new BVHParallelizer(this, hh, gg));
        steppers.push_back(new BVHParallelizer(this, hg, hg));
        steppers.push_back(new BVHParallelizer(this, hg, gh));
        steppers.push_back(new BVHParallelizer(this, hg, gg));
        steppers.push_back(new BVHParallelizer(this, gh, gh));
        steppers.push_back(new BVHParallelizer(this, gh, gg));
        steppers.push_back(new BVHParallelizer(this, gg, gg));
        MultithreadedStepper<std::vector<BVHParallelizer*> > (steppers, m_num_threads).Execute();
        return;
    }

    // If the tree is deep enough, launch recursive parallel collision detection from this level.
    BVHNode& hhh = m_bvh.GetNode(hh.ChildIndex());
    BVHNode& hhg = m_bvh.GetNode(hh.ChildIndex() + 1);
    BVHNode& hgh = m_bvh.GetNode(hg.ChildIndex());
    BVHNode& hgg = m_bvh.GetNode(hg.ChildIndex() + 1);
    BVHNode& ghh = m_bvh.GetNode(gh.ChildIndex());
    BVHNode& ghg = m_bvh.GetNode(gh.ChildIndex() + 1);
    BVHNode& ggh = m_bvh.GetNode(gg.ChildIndex());
    BVHNode& ggg = m_bvh.GetNode(gg.ChildIndex() + 1);
    steppers.reserve(36);
    steppers.push_back(new BVHParallelizer(this, hhh, hhh));
    steppers.push_back(new BVHParallelizer(this, hhh, hhg));
    steppers.push_back(new BVHParallelizer(this, hhh, hgh));
    steppers.push_back(new BVHParallelizer(this, hhh, hgg));
    steppers.push_back(new BVHParallelizer(this, hhh, ghh));
    steppers.push_back(new BVHParallelizer(this, hhh, ghg));
    steppers.push_back(new BVHParallelizer(this, hhh, ggh));
    steppers.push_back(new BVHParallelizer(this, hhh, ggg));
    steppers.push_back(new BVHParallelizer(this, hhg, hhg));
    steppers.push_back(new BVHParallelizer(this, hhg, hgh));
    steppers.push_back(new BVHParallelizer(this, hhg, hgg));
    steppers.push_back(new BVHParallelizer(this, hhg, ghh));
    steppers.push_back(new BVHParallelizer(this, hhg, ghg));
    steppers.push_back(new BVHParallelizer(this, hhg, ggh));
    steppers.push_back(new BVHParallelizer(this, hhg, ggg));
    steppers.push_back(new BVHParallelizer(this, hgh, hgh));
    steppers.push_back(new BVHParallelizer(this, hgh, hgg));
    steppers.push_back(new BVHParallelizer(this, hgh, ghh));
    steppers.push_back(new BVHParallelizer(this, hgh, ghg));
    steppers.push_back(new BVHParallelizer(this, hgh, ggh));
    steppers.push_back(new BVHParallelizer(this, hgh, ggg));
    steppers.push_back(new BVHParallelizer(this, hgg, hgg));
    steppers.push_back(new BVHParallelizer(this, hgg, ghh));
    steppers.push_back(new BVHParallelizer(this, hgg, ghg));
    steppers.push_back(new BVHParallelizer(this, hgg, ggh));
    steppers.push_back(new BVHParallelizer(this, hgg, ggg));
    steppers.push_back(new BVHParallelizer(this, ghh, ghh));
    steppers.push_back(new BVHParallelizer(this, ghh, ghg));
    steppers.push_back(new BVHParallelizer(this, ghh, ggh));
    steppers.push_back(new BVHParallelizer(this, ghh, ggg));
    steppers.push_back(new BVHParallelizer(this, ghg, ghg));
    steppers.push_back(new BVHParallelizer(this, ghg, ggh));
    steppers.push_back(new BVHParallelizer(this, ghg, ggg));
    steppers.push_back(new BVHParallelizer(this, ggh, ggh));
    steppers.push_back(new BVHParallelizer(this, ggh, ggg));
    steppers.push_back(new BVHParallelizer(this, ggg, ggg));
    MultithreadedStepper<std::vector<BVHParallelizer*> > (steppers, m_num_threads).Execute();

    for (std::vector<BVHParallelizer*>::iterator i = steppers.begin(); i != steppers.end(); i++)
        delete *i;

    //    STOP_TIMER("CollisionDetector::getContinuousTimeCollisions");
}

void CollisionDetector::buildBVH()
{
    GeometryBBoxFunctor bboxfunctor(m_elements, m_geodata);
    BVHBuilder bvh_builder;
    bvh_builder.build(bboxfunctor, &m_bvh);
}

void CollisionDetector::computeCollisions(const BVHNode& node_a, const BVHNode& node_b)
{
    // If the bounding volumes do not overlap, there are no possible collisions between their objects
    if (!Intersect(node_a.BBox(), node_b.BBox()))
        return;

    // If both bounding volumes are leaves, add their contents to list potential collisions
    if (node_a.IsLeaf() && node_b.IsLeaf())
    {
        if (&node_a != &node_b)
        {
            const uint32_t leaf_a_begin = node_a.LeafBegin();
            const uint32_t leaf_a_end = node_a.LeafEnd();
            const uint32_t leaf_b_begin = node_b.LeafBegin();
            const uint32_t leaf_b_end = node_b.LeafEnd();

            for (uint32_t i = leaf_a_begin; i < leaf_a_end; ++i)
                for (uint32_t j = leaf_b_begin; j < leaf_b_end; ++j)
                    appendCollision(m_elements[i], m_elements[j]);
        }

    }

    // If one bounding volume is a leaf, we must recurse on the other volume
    else if (node_a.IsLeaf())
    {
        computeCollisions(node_a, m_bvh.GetNode(node_b.ChildIndex()));
        computeCollisions(node_a, m_bvh.GetNode(node_b.ChildIndex() + 1));
    }
    else if (node_b.IsLeaf())
    {
        computeCollisions(m_bvh.GetNode(node_a.ChildIndex()), node_b);
        computeCollisions(m_bvh.GetNode(node_a.ChildIndex() + 1), node_b);
    }
    else
    {
        computeCollisions(m_bvh.GetNode(node_a.ChildIndex()), m_bvh.GetNode(node_b.ChildIndex()));
        computeCollisions(m_bvh.GetNode(node_a.ChildIndex() + 1), m_bvh.GetNode(node_b.ChildIndex()));
        if (&node_a != &node_b) // We need only to explore one side of the diagonal
            computeCollisions(m_bvh.GetNode(node_a.ChildIndex()), m_bvh.GetNode(node_b.ChildIndex() + 1));
        computeCollisions(m_bvh.GetNode(node_a.ChildIndex() + 1), m_bvh.GetNode(node_b.ChildIndex() + 1));
    }
}

} // namespace BASim
