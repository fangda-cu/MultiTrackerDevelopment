/*
 * CollisionDetectorBase.cc
 *
 *  Created on: 25/05/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
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

template void CollisionDetectorBase::updateCollisions<std::list<Collision*> >(std::list<Collision*>& collisions);

template<typename CollisionContainerT>
void CollisionDetectorBase::updateCollisions(CollisionContainerT& collisions)
{
    for (typename CollisionContainerT::iterator collision = collisions.begin(); collision != collisions.end(); collision++)
    {
        bool collisionDetected = (*collision)->analyseCollision(m_time_step);
        if (!collisionDetected)
        {
            delete *collision;
            collisions.erase(collision--);
        }
    }
}

void CollisionDetectorBase::getReady(std::list<Collision*>& cllsns, CollisionFilter collision_filter)
{
    m_potential_collisions = 0;
    m_collision_filter = collision_filter;
    assert(cllsns.empty());
    m_collisions_list = &cllsns;
}

template<>
bool CollisionDetectorBase::appendCollision<EdgeFace>(const YAEdge* edge_a, const YATriangle* triangle)
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

bool CollisionDetectorBase::appendCollision(const TopologicalElement* elem_a, const TopologicalElement* elem_b)
{
    m_potential_collisions++;

    const YAEdge* edge_a = dynamic_cast<const YAEdge*> (elem_a);
    const YATriangle* triangle_a = dynamic_cast<const YATriangle*> (elem_a);
    const YAEdge* edge_b = dynamic_cast<const YAEdge*> (elem_b);
    const YATriangle* triangle_b = dynamic_cast<const YATriangle*> (elem_b);

    switch (m_collision_filter)
    // Pending full templatization of the class
    {
    case ContinuousTime:
        if (edge_a && edge_b)
            return appendCollision<ContinuousTime> (edge_a, edge_b);
        else if (edge_a && triangle_b)
            return appendCollision<ContinuousTime> (edge_a, triangle_b);
        else if (triangle_a && edge_b)
            return appendCollision<ContinuousTime> (edge_b, triangle_a);
        break;

    case Proximity:
        if (edge_a && edge_b)
            return appendCollision<Proximity> (edge_a, edge_b);
        else if (edge_a && triangle_b)
            return appendCollision<Proximity> (edge_a, triangle_b);
        else if (triangle_a && edge_b)
            return appendCollision<Proximity> (edge_b, triangle_a);
        break;

    case EdgeFace:
        if (edge_a && triangle_b)
            return appendCollision<EdgeFace> (edge_a, triangle_b);
        if (edge_b && triangle_a)
            return appendCollision<EdgeFace> (edge_b, triangle_a);
        break;

    default:
        return false;
    }
}

void CollisionDetectorBase::updateBoundingBox(BVH& bvh, const std::vector<const TopologicalElement*>& elements, BVHNodeType& node)
{
    BVHNodeType::BBoxType& bbox = node.BBox();
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
        if (bbox.Volume() > 1e5)
        {
            WarningStream(g_log, "") << "LARGE BOUNDING BOX RESET TO ZERO (volume = " << bbox.Volume() << ")\n";
            DebugStream(g_log, "") << "Bounding box coordinates are: " << bbox.min << " " << bbox.max << '\n';
            DebugStream(g_log, "") << "Number of elements = " << leaf_end - leaf_begin << '\n';

            bbox = BBoxType();
        }
    }
    else // Update the children, then this node's bounding box
    {
        BVHNodeType& hansel = bvh.GetNode(node.ChildIndex());
        updateBoundingBox(bvh, elements, hansel);
        bbox.Insert(hansel.BBox());
        BVHNodeType& gretel = bvh.GetNode(node.ChildIndex() + 1);
        updateBoundingBox(bvh, elements, gretel);
        bbox.Insert(gretel.BBox());
    }
}

template<CollisionFilter CF>
bool CollisionDetectorBase::appendCollision(const YAEdge* edge, const YATriangle* triangle)
{
    YAEdge edge_2(triangle->first(), triangle->second());
    YAEdge edge_1(triangle->third(), triangle->first());
    YAEdge edge_0(triangle->second(), triangle->third());

    bool ee0 = appendCollision<CF> (edge, &edge_0);
    bool ee1 = appendCollision<CF> (edge, &edge_1);
    bool ee2 = appendCollision<CF> (edge, &edge_2);

    bool vf0 = appendCollision<CF> (edge->first(), triangle);
    bool vf1 = appendCollision<CF> (edge->second(), triangle);

    return ee0 || ee1 || ee2 || vf0 || vf1;
}

template<CollisionFilter CF>
bool CollisionDetectorBase::appendCollision(const YAEdge* edge_a, const YAEdge* edge_b)
{
    typedef typename CollisionTraits<CF>::EdgeEdgeCollisionType EdgeEdgeCollisionType;
    EdgeEdgeCollisionType* edgeXedge = new EdgeEdgeCollisionType(m_geodata, edge_a, edge_b);

    bool collisionDetected = edgeXedge->analyseCollision(m_time_step);

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

template<CollisionFilter CF>
bool CollisionDetectorBase::appendCollision(int v_index, const YATriangle* triangle)
{
    typedef typename CollisionTraits<CF>::VertexFaceCollisionType VertexFaceCollisionType;

    VertexFaceCollisionType* vertexXface = new VertexFaceCollisionType(m_geodata, v_index, triangle);

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
        m_collisions_list->push_back(vertexXface); // will be deleted at the end of BARodStepper::step()
        m_collisions_mutex.Unlock();
        return true;
    }
}

}
