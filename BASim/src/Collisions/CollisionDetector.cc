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
    m_geodata(geodata), m_time_step(timestep), m_collisions(NULL)
{
    m_elements.reserve(edges.size() + faces.size());
    for (uint32_t i = 0; i < edges.size(); i++)
        m_elements.push_back(new YAEdge(edges[i]));
    for (uint32_t i = 0; i < faces.size(); i++)
        m_elements.push_back(new YATriangle(faces[i]));

    GeometryBBoxFunctor bboxfunctor(m_elements, m_geodata);
    MiddleBVHBuilder bvh_builder;
    bvh_builder.SetMaxLeafSize(1u);
    bvh_builder.build(bboxfunctor, &m_bvh);
}

CollisionDetector::~CollisionDetector()
{
    m_collisions = NULL;
}

void CollisionDetector::getContinuousTimeCollisions(std::list<CTCollision*>& cllsns)
{
    m_collisions = &cllsns;
    m_collisions->clear();
    BVHNode& root = m_bvh.GetNode(0);
    updateBoundingBox(root);
    computeContinuousTimeCollisions(root, root);
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
    std::cerr << "\033[31mIMPLICIT PENALTY COLLISION DETECTION NOT IMPLEMENTED YET\033[0m";
}

void CollisionDetector::updateBoundingBox(BVHNode& node)
{
    BVHNode::BBoxType& bbox = node.BBox();
    bbox.Reset();
    if (node.IsLeaf()) // Expand the leaf's bounding box to contain all its objects
    {
        const uint32_t leaf_begin = node.LeafBegin();
        const uint32_t leaf_end = node.LeafEnd();
        for (uint32_t i = leaf_begin; i < leaf_end; ++i)
        {
            bbox.Insert(m_elements[i]->GetBBox(m_geodata, 0.0));
            bbox.Insert(m_elements[i]->GetBBox(m_geodata, m_time_step));
        }
    }
    else
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
    if (!Intersect(node_a.BBox(), node_b.BBox()))
        return;

    // If both nodes are leaves, intersect their contents
    if (node_a.IsLeaf() && node_b.IsLeaf())
    {
        //	Timer::getTimer("Leaf leaf intersection").start();
        intersectContent(node_a, node_b);
        //	Timer::getTimer("Leaf leaf intersection").stop();
    }
    // Else recurse on either the node that is not a leave, or the largest volume one.
    else if (node_a.IsLeaf())
    {
        computeContinuousTimeCollisions(node_a, m_bvh.GetNode(node_b.ChildIndex()));
        computeContinuousTimeCollisions(node_a, m_bvh.GetNode(node_b.ChildIndex() + 1));
    }
    else if (node_b.IsLeaf() || node_a.BBox().Volume() >= node_b.BBox().Volume())
    {
        computeContinuousTimeCollisions(m_bvh.GetNode(node_a.ChildIndex()), node_b);
        computeContinuousTimeCollisions(m_bvh.GetNode(node_a.ChildIndex() + 1), node_b);
    }
    else
    {
        computeContinuousTimeCollisions(node_a, m_bvh.GetNode(node_b.ChildIndex()));
        computeContinuousTimeCollisions(node_a, m_bvh.GetNode(node_b.ChildIndex() + 1));
    }
}

void CollisionDetector::intersectContent(const BVHNode& node_a, const BVHNode& node_b)
{
    // First take care of what happens within each box
    intersectContentSelf(node_a);
    intersectContentSelf(node_b);

    // Now intersect between objects in node_a and objects in node_b
    const uint32_t leaf_a_begin = node_a.LeafBegin();
    const uint32_t leaf_a_end = node_a.LeafEnd();
    const uint32_t leaf_b_begin = node_b.LeafBegin();
    const uint32_t leaf_b_end = node_b.LeafEnd();

    for (uint32_t i = leaf_a_begin; i < leaf_a_end; ++i)
        for (uint32_t j = leaf_b_begin; j < leaf_b_end; ++j)
            appendContinuousTimeIntersection(m_elements[i], m_elements[j]);
}

void CollisionDetector::intersectContentSelf(const BVHNode& node)
{
    const uint32_t leaf_begin = node.LeafBegin();
    const uint32_t leaf_end = node.LeafEnd();

    for (uint32_t i = leaf_begin; i < leaf_end; ++i)
        for (uint32_t j = leaf_begin; j < i; ++j)
            appendContinuousTimeIntersection(m_elements[i], m_elements[j]);
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
    // else do nothing. Maybe it's time to introduce exception handling to this code.

}

void CollisionDetector::appendContinuousTimeIntersection(const YAEdge* edge_a, const YAEdge* edge_b)
{
    if (edge_a == edge_b)
        return;

    EdgeEdgeCTCollision* edgeXedge = new EdgeEdgeCTCollision(edge_a, edge_b);

    if (SKIP_ROD_ROD && edgeXedge->IsRodRod(m_geodata)) // Detect rod-rod collisions and skip them.
        return;

    if (edgeXedge->analyseCollision(m_geodata, m_time_step))
        m_collisions->push_back(edgeXedge);
}

void CollisionDetector::appendContinuousTimeIntersection(int v_index, const YATriangle* triangle)
{
    VertexFaceCTCollision* vertexXface = new VertexFaceCTCollision(v_index, triangle);

    // If vertex is fixed, if face is fixed, nothing to do
    if (vertexXface->IsFixed(m_geodata))
        return;

    if (vertexXface->analyseCollision(m_geodata, m_time_step))
        m_collisions->push_back(vertexXface);
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
