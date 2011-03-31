/*
 * CollisionDetector.cc
 *
 *  Created on: 17/03/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#include "CollisionDetector.hh"
#include "CollisionUtils.hh"
#include "../Core/Timer.hh"

#define SKIP_ROD_ROD true

namespace BASim
{
CollisionDetector::CollisionDetector(const GeometricData& geodata, const std::vector<std::pair<int, int> >& edges,
        const std::vector<TriangularFace>& faces, const double& timestep) :
    m_geodata(geodata), m_time_step(timestep), m_collisions(NULL), m_collisions_mutex()
{
    m_elements.reserve(edges.size() + faces.size());
    for (uint32_t i = 0; i < edges.size(); i++)
        m_elements.push_back(new YAEdge(edges[i]));
    for (uint32_t i = 0; i < faces.size(); i++)
        m_elements.push_back(new YATriangle(faces[i]));

    GeometryBBoxFunctor bboxfunctor(m_elements, m_geodata);
    MiddleBVHBuilder bvh_builder;
    bvh_builder.build(bboxfunctor, &m_bvh);
}

CollisionDetector::~CollisionDetector()
{
    m_collisions = NULL;
}

void CollisionDetector::getContinuousTimeCollisions(std::list<CTCollision*>& cllsns)
{
    //    START_TIMER("CollisionDetector::getContinuousTimeCollisions");
    m_collisions = &cllsns;
    m_collisions->clear();
    std::vector<BVHParallelizer*> steppers;
    const int num_threads = 8;
    BVHNode& root = m_bvh.GetNode(0);
    updateBoundingBox(root);

    if (root.IsLeaf()) // Can't really call this a tree, can we?
    {
        computeContinuousTimeCollisions(root, root);
        return;
    }

    BVHNode& h = m_bvh.GetNode(root.ChildIndex());
    BVHNode& g = m_bvh.GetNode(root.ChildIndex() + 1);

    // If tree has depth 1, detect collisions at this level.
    if (h.IsLeaf() || g.IsLeaf())
    {
        steppers.reserve(3);
        steppers.push_back(new BVHParallelizer(*this, h, h));
        steppers.push_back(new BVHParallelizer(*this, h, g));
        steppers.push_back(new BVHParallelizer(*this, g, g));
        MultithreadedStepper<std::vector<BVHParallelizer*> > (steppers, num_threads).Execute();
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
        steppers.push_back(new BVHParallelizer(*this, hh, hh));
        steppers.push_back(new BVHParallelizer(*this, hh, hg));
        steppers.push_back(new BVHParallelizer(*this, hh, gh));
        steppers.push_back(new BVHParallelizer(*this, hh, gg));
        steppers.push_back(new BVHParallelizer(*this, hg, hg));
        steppers.push_back(new BVHParallelizer(*this, hg, gh));
        steppers.push_back(new BVHParallelizer(*this, hg, gg));
        steppers.push_back(new BVHParallelizer(*this, gh, gh));
        steppers.push_back(new BVHParallelizer(*this, gh, gg));
        steppers.push_back(new BVHParallelizer(*this, gg, gg));
        MultithreadedStepper<std::vector<BVHParallelizer*> > (steppers, num_threads).Execute();
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
    steppers.push_back(new BVHParallelizer(*this, hhh, hhh));
    steppers.push_back(new BVHParallelizer(*this, hhh, hhg));
    steppers.push_back(new BVHParallelizer(*this, hhh, hgh));
    steppers.push_back(new BVHParallelizer(*this, hhh, hgg));
    steppers.push_back(new BVHParallelizer(*this, hhh, ghh));
    steppers.push_back(new BVHParallelizer(*this, hhh, ghg));
    steppers.push_back(new BVHParallelizer(*this, hhh, ggh));
    steppers.push_back(new BVHParallelizer(*this, hhh, ggg));
    steppers.push_back(new BVHParallelizer(*this, hhg, hhg));
    steppers.push_back(new BVHParallelizer(*this, hhg, hgh));
    steppers.push_back(new BVHParallelizer(*this, hhg, hgg));
    steppers.push_back(new BVHParallelizer(*this, hhg, ghh));
    steppers.push_back(new BVHParallelizer(*this, hhg, ghg));
    steppers.push_back(new BVHParallelizer(*this, hhg, ggh));
    steppers.push_back(new BVHParallelizer(*this, hhg, ggg));
    steppers.push_back(new BVHParallelizer(*this, hgh, hgh));
    steppers.push_back(new BVHParallelizer(*this, hgh, hgg));
    steppers.push_back(new BVHParallelizer(*this, hgh, ghh));
    steppers.push_back(new BVHParallelizer(*this, hgh, ghg));
    steppers.push_back(new BVHParallelizer(*this, hgh, ggh));
    steppers.push_back(new BVHParallelizer(*this, hgh, ggg));
    steppers.push_back(new BVHParallelizer(*this, hgg, hgg));
    steppers.push_back(new BVHParallelizer(*this, hgg, ghh));
    steppers.push_back(new BVHParallelizer(*this, hgg, ghg));
    steppers.push_back(new BVHParallelizer(*this, hgg, ggh));
    steppers.push_back(new BVHParallelizer(*this, hgg, ggg));
    steppers.push_back(new BVHParallelizer(*this, ghh, ghh));
    steppers.push_back(new BVHParallelizer(*this, ghh, ghg));
    steppers.push_back(new BVHParallelizer(*this, ghh, ggh));
    steppers.push_back(new BVHParallelizer(*this, ghh, ggg));
    steppers.push_back(new BVHParallelizer(*this, ghg, ghg));
    steppers.push_back(new BVHParallelizer(*this, ghg, ggh));
    steppers.push_back(new BVHParallelizer(*this, ghg, ggg));
    steppers.push_back(new BVHParallelizer(*this, ggh, ggh));
    steppers.push_back(new BVHParallelizer(*this, ggh, ggg));
    steppers.push_back(new BVHParallelizer(*this, ggg, ggg));
    MultithreadedStepper<std::vector<BVHParallelizer*> > (steppers, num_threads).Execute();

    //    STOP_TIMER("CollisionDetector::getContinuousTimeCollisions");
}

void CollisionDetector::updateContinuousTimeCollisions()
{
    for (std::list<CTCollision*>::iterator i = m_collisions->begin(); i != m_collisions->end(); i++)
        if (!(*i)->analyseCollision(m_geodata, m_time_step))
            m_collisions->erase(i--);
}

void CollisionDetector::getImplicitPenaltyCollisions(std::vector<EdgeEdgeProximityCollision>& edge_edge_collisions,
        std::vector<VertexFaceProximityCollision>& vertex_face_collisions)
{
    std::cerr << "\033[31mIMPLICIT PENALTY COLLISION DETECTION NOT IMPLEMENTED YET\033[0m"; // FIXME
}

void CollisionDetector::updateBoundingBox(BVHNode& node)
{
    BVHNode::BBoxType& bbox = node.BBox();
    bbox.Reset();
    if (node.IsLeaf()) // The leaf's bounding box contains the whole trajectory of its object during this time step.
    {
        const uint32_t leaf_begin = node.LeafBegin();
        const uint32_t leaf_end = node.LeafEnd();
        for (uint32_t i = leaf_begin; i < leaf_end; ++i)
        {
            bbox.Insert(m_elements[i]->GetBBox(m_geodata, 0.0));
            bbox.Insert(m_elements[i]->GetBBox(m_geodata, m_time_step));
        }
    }
    else // Update the children, then this node's bounding box
    {
        BVHNode& hansel = m_bvh.GetNode(node.ChildIndex());
        updateBoundingBox(hansel);
        bbox.Insert(hansel.BBox());
        BVHNode& gretel = m_bvh.GetNode(node.ChildIndex() + 1);
        updateBoundingBox(gretel);
        bbox.Insert(gretel.BBox());
    }
}

void CollisionDetector::computeContinuousTimeCollisions(const BVHNode& node_a, const BVHNode& node_b)
{
    if (// (&node_a > &node_b) ||
    (!Intersect(node_a.BBox(), node_b.BBox())))
        return;

    // If both nodes are leaves, intersect their contents
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
                    appendContinuousTimeIntersection(m_elements[i], m_elements[j]);
        }
    }
    // Else recurse on either the node that is not a leave, or both.
    else if (node_a.IsLeaf())
    {
        computeContinuousTimeCollisions(node_a, m_bvh.GetNode(node_b.ChildIndex()));
        computeContinuousTimeCollisions(node_a, m_bvh.GetNode(node_b.ChildIndex() + 1));
    }
    else if (node_b.IsLeaf())
    {
        computeContinuousTimeCollisions(m_bvh.GetNode(node_a.ChildIndex()), node_b);
        computeContinuousTimeCollisions(m_bvh.GetNode(node_a.ChildIndex() + 1), node_b);
    }
    else
    {
        computeContinuousTimeCollisions(m_bvh.GetNode(node_a.ChildIndex()), m_bvh.GetNode(node_b.ChildIndex()));
        computeContinuousTimeCollisions(m_bvh.GetNode(node_a.ChildIndex() + 1), m_bvh.GetNode(node_b.ChildIndex()));
        if (&node_a != &node_b) // We need only to explore one side of the diagonal
            computeContinuousTimeCollisions(m_bvh.GetNode(node_a.ChildIndex()), m_bvh.GetNode(node_b.ChildIndex() + 1));
        computeContinuousTimeCollisions(m_bvh.GetNode(node_a.ChildIndex() + 1), m_bvh.GetNode(node_b.ChildIndex() + 1));
    }
}

void CollisionDetector::appendContinuousTimeIntersection(const TopologicalElement* elem_a, const TopologicalElement* elem_b)
{
    const YAEdge* edge_a = dynamic_cast<const YAEdge*> (elem_a);
    const YATriangle* triangle_a = dynamic_cast<const YATriangle*> (elem_a);
    const YAEdge* edge_b = dynamic_cast<const YAEdge*> (elem_b);
    const YATriangle* triangle_b = dynamic_cast<const YATriangle*> (elem_b);

    if (edge_a && edge_b)
        appendContinuousTimeIntersection(edge_a, edge_b);
    else if (edge_a && triangle_b)
        appendContinuousTimeIntersection(edge_a, triangle_b);
    else if (triangle_a && edge_b)
        appendContinuousTimeIntersection(triangle_a, edge_b);
    else if (triangle_a && triangle_b)
        appendContinuousTimeIntersection(triangle_a, triangle_b);
}

void CollisionDetector::appendContinuousTimeIntersection(const YAEdge* edge_a, const YAEdge* edge_b)
{
    //    Timer::getTimer("CollisionDetector::appendContinuousTimeIntersection edge edge").start();

    EdgeEdgeCTCollision* edgeXedge = new EdgeEdgeCTCollision(edge_a, edge_b);

    if (SKIP_ROD_ROD && edgeXedge->IsRodRod(m_geodata)) // Detect rod-rod collisions and skip them.
        return;

    if (edgeXedge->analyseCollision(m_geodata, m_time_step))
    {
        m_collisions_mutex.Lock();
        m_collisions->push_back(edgeXedge);
        m_collisions_mutex.Unlock();
    }

    //    Timer::getTimer("CollisionDetector::appendContinuousTimeIntersection edge edge").stop();
}

void CollisionDetector::appendContinuousTimeIntersection(int v_index, const YATriangle* triangle)
{
    //    Timer::getTimer("CollisionDetector::appendContinuousTimeIntersection vertex face").start();

    VertexFaceCTCollision* vertexXface = new VertexFaceCTCollision(v_index, triangle);

    // If vertex is fixed, if face is fixed, nothing to do
    if (vertexXface->IsFixed(m_geodata))
        return;

    if (vertexXface->analyseCollision(m_geodata, m_time_step))
    {
        m_collisions_mutex.Lock();
        m_collisions->push_back(vertexXface);
        m_collisions_mutex.Unlock();
    }

    //    Timer::getTimer("CollisionDetector::appendContinuousTimeIntersection vertex face").stop();
}

void CollisionDetector::appendContinuousTimeIntersection(const YAEdge* edge, const YATriangle* triangle)
{

    YAEdge edge_2(triangle->first(), triangle->second());
    YAEdge edge_1(triangle->third(), triangle->first());
    YAEdge edge_0(triangle->second(), triangle->third());
    appendContinuousTimeIntersection(edge, &edge_0);
    appendContinuousTimeIntersection(edge, &edge_1);
    appendContinuousTimeIntersection(edge, &edge_2);

    appendContinuousTimeIntersection(edge->first(), triangle);
    appendContinuousTimeIntersection(edge->second(), triangle);
}

void CollisionDetector::appendContinuousTimeIntersection(const YATriangle* triangle, const YAEdge* edge)
{
    appendContinuousTimeIntersection(edge, triangle);
}

void CollisionDetector::appendContinuousTimeIntersection(const YATriangle* triangle_a, const YATriangle* triangle_b)
{
    // Do nothing.
    return;
}

}
