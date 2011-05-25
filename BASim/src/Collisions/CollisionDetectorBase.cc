/*
 * CollisionDetectorBase.cc
 *
 *  Created on: 25/05/2011
 *      Author: jaubry
 */

#include "CollisionDetectorBase.hh"
#include "CollisionUtils.hh"
#include "../Core/Timer.hh"

namespace BASim
{

CollisionDetectorBase::CollisionDetectorBase(const GeometricData& geodata, const std::vector<std::pair<int, int> >& edges,
        const std::vector<TriangularFace>& faces, const double& timestep, bool skip_rod_rod, int num_threads) :
    m_geodata(geodata), m_time_step(timestep), m_skip_rod_rod(skip_rod_rod), m_collisions_list(NULL), m_collisions_mutex()
{
    // std::cerr << "Constructing CollisionDetectorBase" << std::endl;
    if (num_threads > 0)
        m_num_threads = num_threads;
    else
        m_num_threads = sysconf(_SC_NPROCESSORS_ONLN);
}

CollisionDetectorBase::~CollisionDetectorBase()
{
    m_collisions_list = NULL;
}


void CollisionDetectorBase::updateContinuousTimeCollisions()
{
    for (std::list<Collision*>::iterator collision = m_collisions_list->begin(); collision != m_collisions_list->end(); collision++)
        if (!(*collision)->analyseCollision(m_time_step))
        {
            delete *collision;
            m_collisions_list->erase(collision--);
        }
}

void CollisionDetectorBase::appendCollision(const TopologicalElement* elem_a, const TopologicalElement* elem_b)
{
    const YAEdge* edge_a = dynamic_cast<const YAEdge*> (elem_a);
    const YATriangle* triangle_a = dynamic_cast<const YATriangle*> (elem_a);
    const YAEdge* edge_b = dynamic_cast<const YAEdge*> (elem_b);
    const YATriangle* triangle_b = dynamic_cast<const YATriangle*> (elem_b);

    switch (m_collision_filter)
    {
    case ContinuousTime:
        if (edge_a && edge_b)
            appendContinuousTimeCollision(edge_a, edge_b);
        else if (edge_a && triangle_b)
            appendContinuousTimeCollision(edge_a, triangle_b);
        else if (triangle_a && edge_b)
            appendContinuousTimeCollision(edge_b, triangle_a);
        break;

    case Proximity:
        if (edge_a && edge_b)
            appendProximityCollision(edge_a, edge_b);
        else if (edge_a && triangle_b)
            appendProximityCollision(edge_a, triangle_b);
        else if (triangle_a && edge_b)
            appendProximityCollision(edge_b, triangle_a);
        break;

    case EdgeFace:
        if (edge_a && triangle_b)
            appendEdgeFaceIntersection(edge_a, triangle_b);
        if (edge_b && triangle_a)
            appendEdgeFaceIntersection(edge_b, triangle_a);
        break;
    }
}

void CollisionDetectorBase::updateBoundingBox(BVH& bvh, const std::vector<const TopologicalElement*>& elements, BVHNode& node)
{
    BVHNode::BBoxType& bbox = node.BBox();
    bbox.Reset();
    if (node.IsLeaf()) // The leaf's bounding box contains the whole trajectory of its object during this time step.
    {
        const uint32_t leaf_begin = node.LeafBegin();
        const uint32_t leaf_end = node.LeafEnd();
        for (uint32_t i = leaf_begin; i < leaf_end; ++i)
        {
            bbox.Insert(elements[i]->GetBBox(m_geodata, 0.0));
            bbox.Insert(elements[i]->GetBBox(m_geodata, m_time_step));
        }
    }
    else // Update the children, then this node's bounding box
    {
        BVHNode& hansel = bvh.GetNode(node.ChildIndex());
        updateBoundingBox(bvh, elements, hansel);
        bbox.Insert(hansel.BBox());
        BVHNode& gretel = bvh.GetNode(node.ChildIndex() + 1);
        updateBoundingBox(bvh, elements, gretel);
        bbox.Insert(gretel.BBox());
    }
}

void CollisionDetectorBase::appendContinuousTimeCollision(const YAEdge* edge_a, const YAEdge* edge_b)
{
    //    Timer::getTimer("CollisionDetectorBase::appendContinuousTimeCollision edge edge").start();

    EdgeEdgeCTCollision* edgeXedge = new EdgeEdgeCTCollision(m_geodata, edge_a, edge_b);

    if ((m_skip_rod_rod && edgeXedge->IsRodRod()) || edgeXedge->IsCollisionImmune())
    { // Detect rod-rod collisions and skip them.
        //std::cout << "CollisionDetectorBase: Skipping rod-rod collision" << std::endl;
        delete edgeXedge;
        return;
    }

    if (edgeXedge->analyseCollision(m_time_step))
    {
        m_collisions_mutex.Lock();
        m_collisions_list->push_back(edgeXedge); // Will be deleted in BARodStepper::executeIterativeInelasticImpulseResponse()
        m_collisions_mutex.Unlock();
        // std::cout << "CollisionDetectorBase: Found edge-edge collision" << std::endl;
    }
    else
        delete edgeXedge;

    //    Timer::getTimer("CollisionDetectorBase::appendContinuousTimeCollision edge edge").stop();
}

void CollisionDetectorBase::appendContinuousTimeCollision(int v_index, const YATriangle* triangle)
{
    //    Timer::getTimer("CollisionDetectorBase::appendContinuousTimeCollision vertex face").start();

    VertexFaceCTCollision* vertexXface = new VertexFaceCTCollision(m_geodata, v_index, triangle);

    // If vertex is fixed, if face is fixed, nothing to do
    if (vertexXface->IsFixed() || m_geodata.IsCollisionImmune(v_index))
    {
        delete vertexXface;
        // std::cout << "CollisionDetectorBase: Skipping vertex " << v_index << " - face " << triangle->first() << "/" << triangle->second() << "/" << triangle->third() << " collision with fixed vertex and/or face" << std::endl;
        return;
    }

    if (vertexXface->analyseCollision(m_time_step))
    {
        m_collisions_mutex.Lock();
        m_collisions_list->push_back(vertexXface); // Will be deleted in BARodStepper::executeIterativeInelasticImpulseResponse()
        m_collisions_mutex.Unlock();
        // std::cout << "CollisionDetectorBase: Found vertex-face collision" << std::endl;
    }
    else
        delete vertexXface;

    //    Timer::getTimer("CollisionDetectorBase::appendContinuousTimeCollision vertex face").stop();
}

void CollisionDetectorBase::appendContinuousTimeCollision(const YAEdge* edge, const YATriangle* triangle)
{
    YAEdge edge_2(triangle->first(), triangle->second());
    YAEdge edge_1(triangle->third(), triangle->first());
    YAEdge edge_0(triangle->second(), triangle->third());
    appendContinuousTimeCollision(edge, &edge_0);
    appendContinuousTimeCollision(edge, &edge_1);
    appendContinuousTimeCollision(edge, &edge_2);

    appendContinuousTimeCollision(edge->first(), triangle);
    appendContinuousTimeCollision(edge->second(), triangle);
}

void CollisionDetectorBase::appendEdgeFaceIntersection(const YAEdge* edge_a, const YATriangle* triangle)
{
    EdgeFaceIntersection* edgeXface = new EdgeFaceIntersection(m_geodata, edge_a, triangle);

    if (edgeXface->analyseCollision())
    {
        m_collisions_mutex.Lock();
        m_collisions_list->push_back(edgeXface);
        m_collisions_mutex.Unlock();
    }
    else
        delete edgeXface;
}

void CollisionDetectorBase::appendProximityCollision(const YAEdge* edge_a, const YAEdge* edge_b)
{

    EdgeEdgeProximityCollision* edgeXedge = new EdgeEdgeProximityCollision(m_geodata, edge_a, edge_b);

    if ((m_skip_rod_rod && edgeXedge->IsRodRod()) || edgeXedge->IsCollisionImmune())
    { // Detect rod-rod collisions and skip them.
        //std::cout << "CollisionDetectorBase: Skipping rod-rod collision" << std::endl;
        delete edgeXedge;
        return;
    }

    if (edgeXedge->analyseCollision())
    {
        m_collisions_mutex.Lock();
        m_collisions_list->push_back(edgeXedge); // Will be deleted in BARodStepper::executeImplicitPenaltyResponse()
        m_collisions_mutex.Unlock();
        // std::cout << "CollisionDetectorBase: Found edge-edge collision" << std::endl;
    }
    else
        delete edgeXedge;

}

void CollisionDetectorBase::appendProximityCollision(const YAEdge* edge, const YATriangle* triangle)
{
    YAEdge edge_2(triangle->first(), triangle->second());
    YAEdge edge_1(triangle->third(), triangle->first());
    YAEdge edge_0(triangle->second(), triangle->third());
    appendProximityCollision(edge, &edge_0);
    appendProximityCollision(edge, &edge_1);
    appendProximityCollision(edge, &edge_2);

    appendProximityCollision(edge->first(), triangle);
    appendProximityCollision(edge->second(), triangle);
}

void CollisionDetectorBase::appendProximityCollision(int v_index, const YATriangle* triangle)
{

    VertexFaceProximityCollision* vertexXface = new VertexFaceProximityCollision(m_geodata, v_index, triangle);

    // If vertex is fixed, if face is fixed, nothing to do
    if (vertexXface->IsFixed() || m_geodata.IsCollisionImmune(v_index))
    {
        delete vertexXface;
        // std::cout << "CollisionDetectorBase: Skipping vertex " << v_index << " - face " << triangle->first() << "/" << triangle->second() << "/" << triangle->third() << " collision with fixed vertex and/or face" << std::endl;
        return;
    }

    if (vertexXface->analyseCollision())
    {
        m_collisions_mutex.Lock();
        m_collisions_list->push_back(vertexXface); // Will be deleted in BARodStepper::executeImplicitPenaltyResponse()
        m_collisions_mutex.Unlock();
        // std::cout << "CollisionDetectorBase: Found vertex-face collision" << std::endl;
    }
    else
        delete vertexXface;

    //    Timer::getTimer("CollisionDetectorBase::appendContinuousTimeCollision vertex face").stop();
}



}
