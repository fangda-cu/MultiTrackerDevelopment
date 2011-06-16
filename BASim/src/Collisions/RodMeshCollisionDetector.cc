/*
 * RodMeshCollisionDetector.cc
 *
 *  Created on: 25/05/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#include "RodMeshCollisionDetector.hh"
#include "../Util/TextLog.hh"

namespace BASim
{

RodMeshCollisionDetector::RodMeshCollisionDetector(const GeometricData& geodata,
        const std::vector<std::pair<int, int> >& edges, const std::vector<TriangularFace>& faces, const double& timestep,
        bool skip_rod_rod, int num_threads) :
    CollisionDetectorBase(geodata, edges, faces, timestep, skip_rod_rod, num_threads)
{
    if (num_threads > 0)
        m_num_threads = num_threads;
    else
        m_num_threads = sysconf(_SC_NPROCESSORS_ONLN);

    m_rod_elements.reserve(edges.size());
    m_mesh_elements.reserve(faces.size());
    for (std::vector<std::pair<int, int> >::const_iterator i = edges.begin(); i != edges.end(); i++)
        m_rod_elements.push_back(new YAEdge(*i));
    for (std::vector<TriangularFace>::const_iterator i = faces.begin(); i != faces.end(); i++)
        m_mesh_elements.push_back(new YATriangle(*i));

    //  build_mesh_BVH();
    // build_rod_BVH();
}

RodMeshCollisionDetector::~RodMeshCollisionDetector()
{
    m_collisions_list = NULL;
    for (std::vector<const TopologicalElement*>::iterator i = m_rod_elements.begin(); i != m_rod_elements.end(); i++)
        delete *i;
    for (std::vector<const TopologicalElement*>::iterator i = m_mesh_elements.begin(); i != m_mesh_elements.end(); i++)
        delete *i;
}

void RodMeshCollisionDetector::build_rod_BVH()
{
    GeometryBBoxFunctor bboxfunctor(m_rod_elements, m_geodata);
    BVHBuilder bvh_builder;
    bvh_builder.build(bboxfunctor, &m_rod_bvh);
}

void RodMeshCollisionDetector::build_mesh_BVH()
{
    GeometryBBoxFunctor bboxfunctor(m_mesh_elements, m_geodata);
    BVHBuilder bvh_builder;
    bvh_builder.build(bboxfunctor, &m_mesh_bvh);
}

void RodMeshCollisionDetector::getCollisions(std::list<Collision*>& cllsns, CollisionFilter collision_filter,
        bool update_mesh_bbox)
{
    getReady(cllsns, collision_filter);

    std::vector<BVHParallelizer*> steppers;

    BVHNode& rod_root = m_rod_bvh.GetNode(0);
    DebugStream(g_log, "") << "Updating rods bounding box\n";
    updateBoundingBox(m_rod_bvh, m_rod_elements, rod_root);
    BVHNode& mesh_root = m_mesh_bvh.GetNode(0);
    if (update_mesh_bbox)
    {
        DebugStream(g_log, "") << "Updating rods bounding box\n";
        updateBoundingBox(m_mesh_bvh, m_mesh_elements, mesh_root);
    }

    if (mesh_root.IsLeaf() || rod_root.IsLeaf()) // Lazy!
    {
        computeCollisions(mesh_root, rod_root);
        return;
    }

    BVHNode& mesh_h = m_mesh_bvh.GetNode(mesh_root.ChildIndex());
    BVHNode& mesh_g = m_mesh_bvh.GetNode(mesh_root.ChildIndex() + 1);
    BVHNode& rod_h = m_rod_bvh.GetNode(rod_root.ChildIndex());
    BVHNode& rod_g = m_rod_bvh.GetNode(rod_root.ChildIndex() + 1);

    if (mesh_h.IsLeaf() || mesh_g.IsLeaf() || rod_h.IsLeaf() || rod_g.IsLeaf()) // Lazy!
    {
        steppers.push_back(new BVHParallelizer(this, mesh_h, rod_h));
        steppers.push_back(new BVHParallelizer(this, mesh_h, rod_g));
        steppers.push_back(new BVHParallelizer(this, mesh_g, rod_h));
        steppers.push_back(new BVHParallelizer(this, mesh_g, rod_g));
        MultithreadedStepper<std::vector<BVHParallelizer*> > (steppers, m_num_threads).Execute();
        return;
    }

    // else
    BVHNode& mesh_hh = m_mesh_bvh.GetNode(mesh_h.ChildIndex());
    BVHNode& mesh_hg = m_mesh_bvh.GetNode(mesh_h.ChildIndex() + 1);
    BVHNode& mesh_gh = m_mesh_bvh.GetNode(mesh_g.ChildIndex());
    BVHNode& mesh_gg = m_mesh_bvh.GetNode(mesh_g.ChildIndex() + 1);
    BVHNode& rod_hh = m_rod_bvh.GetNode(rod_h.ChildIndex());
    BVHNode& rod_hg = m_rod_bvh.GetNode(rod_h.ChildIndex() + 1);
    BVHNode& rod_gh = m_rod_bvh.GetNode(rod_g.ChildIndex());
    BVHNode& rod_gg = m_rod_bvh.GetNode(rod_g.ChildIndex() + 1);

    steppers.push_back(new BVHParallelizer(this, mesh_hh, rod_hh));
    steppers.push_back(new BVHParallelizer(this, mesh_hg, rod_hh));
    steppers.push_back(new BVHParallelizer(this, mesh_gh, rod_hh));
    steppers.push_back(new BVHParallelizer(this, mesh_gg, rod_hh));
    steppers.push_back(new BVHParallelizer(this, mesh_hh, rod_hg));
    steppers.push_back(new BVHParallelizer(this, mesh_hg, rod_hg));
    steppers.push_back(new BVHParallelizer(this, mesh_gh, rod_hg));
    steppers.push_back(new BVHParallelizer(this, mesh_gg, rod_hg));
    steppers.push_back(new BVHParallelizer(this, mesh_hh, rod_gh));
    steppers.push_back(new BVHParallelizer(this, mesh_hg, rod_gh));
    steppers.push_back(new BVHParallelizer(this, mesh_gh, rod_gh));
    steppers.push_back(new BVHParallelizer(this, mesh_gg, rod_gh));
    steppers.push_back(new BVHParallelizer(this, mesh_hh, rod_gg));
    steppers.push_back(new BVHParallelizer(this, mesh_hg, rod_gg));
    steppers.push_back(new BVHParallelizer(this, mesh_gh, rod_gg));
    steppers.push_back(new BVHParallelizer(this, mesh_gg, rod_gg));
    MultithreadedStepper<std::vector<BVHParallelizer*> > (steppers, m_num_threads).Execute();

    for (std::vector<BVHParallelizer*>::iterator i = steppers.begin(); i != steppers.end(); i++)
        delete *i;
}

void RodMeshCollisionDetector::computeCollisions(const BVHNode& mesh_node, const BVHNode& rod_node)
{
    // If the bounding volumes do not overlap, there are no possible collisions between their objects
    if (!Intersect(mesh_node.BBox(), rod_node.BBox()))
        return;

    // If both bounding volumes are leaves, add their contents to list potential collisions
    if (mesh_node.IsLeaf() && rod_node.IsLeaf())
    {
        const uint32_t mesh_leaf_begin = mesh_node.LeafBegin();
        const uint32_t mesh_leaf_end = mesh_node.LeafEnd();
        const uint32_t rod_leaf_begin = rod_node.LeafBegin();
        const uint32_t rod_leaf_end = rod_node.LeafEnd();
        for (uint32_t i = mesh_leaf_begin; i < mesh_leaf_end; ++i)
            for (uint32_t j = rod_leaf_begin; j < rod_leaf_end; ++j)
                appendCollision(m_mesh_elements[i], m_rod_elements[j]);
    }
    // If one bounding volume is a leaf, we must recurse on the other volume
    else if (mesh_node.IsLeaf())
    {
        computeCollisions(mesh_node, m_rod_bvh.GetNode(rod_node.ChildIndex()));
        computeCollisions(mesh_node, m_rod_bvh.GetNode(rod_node.ChildIndex() + 1));
    }
    else if (rod_node.IsLeaf())
    {
        computeCollisions(m_mesh_bvh.GetNode(mesh_node.ChildIndex()), rod_node);
        computeCollisions(m_mesh_bvh.GetNode(mesh_node.ChildIndex() + 1), rod_node);
    }
    else
    {
        computeCollisions(m_mesh_bvh.GetNode(mesh_node.ChildIndex()), m_rod_bvh.GetNode(rod_node.ChildIndex()));
        computeCollisions(m_mesh_bvh.GetNode(mesh_node.ChildIndex() + 1), m_rod_bvh.GetNode(rod_node.ChildIndex()));
        computeCollisions(m_mesh_bvh.GetNode(mesh_node.ChildIndex()), m_rod_bvh.GetNode(rod_node.ChildIndex() + 1));
        computeCollisions(m_mesh_bvh.GetNode(mesh_node.ChildIndex() + 1), m_rod_bvh.GetNode(rod_node.ChildIndex() + 1));
    }
}

void RodMeshCollisionDetector::rebuildRodElements(const std::vector<std::pair<int, int> >& edges)
{ // TODO: something smarter
    m_rod_elements.clear();

    for (std::vector<std::pair<int, int> >::const_iterator i = edges.begin(); i != edges.end(); i++)
        m_rod_elements.push_back(new YAEdge(*i));

}

}
