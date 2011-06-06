/*
 * CollisionDetectorBase.cc
 *
 *  Created on: 25/05/2011
 *      Author: jaubry
 */

#include "CollisionDetectorBase.hh"
#include "CollisionUtils.hh"
#include "../Core/Timer.hh"
#include "../Util/TextLog.hh"

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

void CollisionDetectorBase::getReady(std::list<Collision*>& cllsns, CollisionFilter collision_filter)
{
    m_potential_collisions = 0;
    m_collision_filter = collision_filter;
    assert(cllsns.empty());
    m_collisions_list = &cllsns;
}

void CollisionDetectorBase::updateContinuousTimeCollisions()
{
    for (std::list<Collision*>::iterator collision = m_collisions_list->begin(); collision != m_collisions_list->end(); collision++)
    {
        bool collisionDetected = (*collision)->analyseCollision(m_time_step);
        if (!collisionDetected)
        {
            delete *collision;
            m_collisions_list->erase(collision--);
        }
    }
}

bool CollisionDetectorBase::appendCollision(const TopologicalElement* elem_a, const TopologicalElement* elem_b)
{
    m_potential_collisions++;

    const YAEdge* edge_a = dynamic_cast<const YAEdge*> (elem_a);
    const YATriangle* triangle_a = dynamic_cast<const YATriangle*> (elem_a);
    const YAEdge* edge_b = dynamic_cast<const YAEdge*> (elem_b);
    const YATriangle* triangle_b = dynamic_cast<const YATriangle*> (elem_b);

    switch (m_collision_filter)
    {
    case ContinuousTime:
        if (edge_a && edge_b)
            return appendContinuousTimeCollision(edge_a, edge_b);
        else if (edge_a && triangle_b)
            return appendContinuousTimeCollision(edge_a, triangle_b);
        else if (triangle_a && edge_b)
            return appendContinuousTimeCollision(edge_b, triangle_a);
        break;

    case Proximity:
        if (edge_a && edge_b)
            return appendProximityCollision(edge_a, edge_b);
        else if (edge_a && triangle_b)
            return appendProximityCollision(edge_a, triangle_b);
        else if (triangle_a && edge_b)
            return appendProximityCollision(edge_b, triangle_a);
        break;

    case EdgeFace:
        if (edge_a && triangle_b)
            return appendEdgeFaceIntersection(edge_a, triangle_b);
        if (edge_b && triangle_a)
            return appendEdgeFaceIntersection(edge_b, triangle_a);
        break;

    default:
        return false;
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
            const YAEdge* edge = dynamic_cast<const YAEdge*> (elements[i]);
            if (!edge || !edge->IsCollisionImmune(m_geodata)) // Immune edges shouldn't be taken into account here.
            {
                bbox.Insert(elements[i]->GetBBox(m_geodata, 0.0));
                bbox.Insert(elements[i]->GetBBox(m_geodata, m_time_step));
            }
        }
        if (bbox.Volume() > 100)
        {
            WarningStream(g_log, "") << "LARGE BOUNDING BOX RESET TO ZERO\n";
            bbox = BBoxType();
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

bool CollisionDetectorBase::appendContinuousTimeCollision(const YAEdge* edge_a, const YAEdge* edge_b)
{
    //    Timer::getTimer("CollisionDetectorBase::appendContinuousTimeCollision edge edge").start();

    EdgeEdgeCTCollision* edgeXedge = new EdgeEdgeCTCollision(m_geodata, edge_a, edge_b);

    bool collisionDetected = edgeXedge->analyseCollision(m_time_step);

    if ((m_skip_rod_rod && edgeXedge->IsRodRod()) || edgeXedge->IsCollisionImmune() || !collisionDetected)
    {
        delete edgeXedge;
        return false;
    }
    else
    {
        m_collisions_mutex.Lock();
        m_collisions_list->push_back(edgeXedge); // Will be deleted in BARodStepper::executeIterativeInelasticImpulseResponse()
        m_collisions_mutex.Unlock();
        return true;
    }
    //    Timer::getTimer("CollisionDetectorBase::appendContinuousTimeCollision edge edge").stop();
}

bool CollisionDetectorBase::appendContinuousTimeCollision(int v_index, const YATriangle* triangle)
{
    //    Timer::getTimer("CollisionDetectorBase::appendContinuousTimeCollision vertex face").start();

    VertexFaceCTCollision* vertexXface = new VertexFaceCTCollision(m_geodata, v_index, triangle);

    bool collisionDetected = vertexXface->analyseCollision(m_time_step);

    // If vertex is fixed, if face is fixed, nothing to do
    if (vertexXface->IsFixed() || m_geodata.IsCollisionImmune(v_index) || !collisionDetected)
    {
        delete vertexXface;
        return false;
    }
    else
    {
        m_collisions_mutex.Lock();
        m_collisions_list->push_back(vertexXface); // Will be deleted in BARodStepper::executeIterativeInelasticImpulseResponse()
        m_collisions_mutex.Unlock();
        return true;
    }
    //    Timer::getTimer("CollisionDetectorBase::appendContinuousTimeCollision vertex face").stop();
}

bool CollisionDetectorBase::appendContinuousTimeCollision(const YAEdge* edge, const YATriangle* triangle)
{
    YAEdge edge_2(triangle->first(), triangle->second());
    YAEdge edge_1(triangle->third(), triangle->first());
    YAEdge edge_0(triangle->second(), triangle->third());

    bool ee0 = appendContinuousTimeCollision(edge, &edge_0);
    bool ee1 = appendContinuousTimeCollision(edge, &edge_1);
    bool ee2 = appendContinuousTimeCollision(edge, &edge_2);

    bool vf0 = appendContinuousTimeCollision(edge->first(), triangle);
    bool vf1 = appendContinuousTimeCollision(edge->second(), triangle);

    return ee0 || ee1 || ee2 || vf0 || vf1;
}

bool CollisionDetectorBase::appendEdgeFaceIntersection(const YAEdge* edge_a, const YATriangle* triangle)
{
    EdgeFaceIntersection* edgeXface = new EdgeFaceIntersection(m_geodata, edge_a, triangle);

    bool collisionDetected = edgeXface->analyseCollision();

    if (collisionDetected)
    {
        m_collisions_mutex.Lock();
        m_collisions_list->push_back(edgeXface);
        m_collisions_mutex.Unlock();
        return true;
    }
    else
    {
        delete edgeXface;
        return false;
    }
}

bool CollisionDetectorBase::appendProximityCollision(const YAEdge* edge_a, const YAEdge* edge_b)
{

    EdgeEdgeProximityCollision* edgeXedge = new EdgeEdgeProximityCollision(m_geodata, edge_a, edge_b);

    bool collisionDetected = edgeXedge->analyseCollision();

    if ((m_skip_rod_rod && edgeXedge->IsRodRod()) || edgeXedge->IsCollisionImmune() || !collisionDetected)
    {
        delete edgeXedge;
        return false;
    }
    else
    {
        m_collisions_mutex.Lock();
        m_collisions_list->push_back(edgeXedge); // will be deleted at the end of BARodStepper::step()
        m_collisions_mutex.Unlock();
        return true;
    }
}

bool CollisionDetectorBase::appendProximityCollision(const YAEdge* edge, const YATriangle* triangle)
{
    YAEdge edge_2(triangle->first(), triangle->second());
    YAEdge edge_1(triangle->third(), triangle->first());
    YAEdge edge_0(triangle->second(), triangle->third());
    bool ee0 = appendProximityCollision(edge, &edge_0);
    bool ee1 = appendProximityCollision(edge, &edge_1);
    bool ee2 = appendProximityCollision(edge, &edge_2);

    bool vf0 = appendProximityCollision(edge->first(), triangle);
    bool vf1 = appendProximityCollision(edge->second(), triangle);

    return ee0 || ee1 || ee2 || vf0 || vf1;
}

bool CollisionDetectorBase::appendProximityCollision(int v_index, const YATriangle* triangle)
{

    VertexFaceProximityCollision* vertexXface = new VertexFaceProximityCollision(m_geodata, v_index, triangle);

    bool collisionDetected = vertexXface->analyseCollision();

    // If vertex is fixed, if face is fixed, nothing to do
    if (vertexXface->IsFixed() || m_geodata.IsCollisionImmune(v_index) || !collisionDetected)
    {
        delete vertexXface;
        return false;
    }
    else
    {
        m_collisions_mutex.Lock();
        m_collisions_list->push_back(vertexXface); // will be deleted at the end of BARodStepper::step()
        m_collisions_mutex.Unlock();
        return true;
    }
    //    Timer::getTimer("CollisionDetectorBase::appendContinuousTimeCollision vertex face").stop();
}

}
