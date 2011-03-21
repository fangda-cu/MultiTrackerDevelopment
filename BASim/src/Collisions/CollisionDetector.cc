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
CollisionDetector::CollisionDetector(const VecXd& x, const VecXd& v, const std::vector<std::pair<int, int> >& edges,
        const std::vector<TriangularFace>& faces, const std::vector<double>& vertex_radii, const std::vector<double>& masses,
        int obj_start, const double& timestep) :
    m_geodata(x, v, vertex_radii, masses, obj_start), m_time_step(timestep), m_collisions(NULL)
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

void CollisionDetector::getContinuousTimeCollisions(std::list<ContinuousTimeCollision>& cllsns)
{
    Timer::getTimer("Collision detector").start();

    m_collisions = &cllsns;
    m_collisions->clear();
    BVHNode& root = m_bvh.GetNode(0);

    Timer::getTimer("Expanding boxes").start();
    updateBoundingBox(root);
    Timer::getTimer("Expanding boxes").stop();

    Timer::getTimer("Computing collisions").start();
    computeContinuousTimeCollisions(root, root);
    Timer::getTimer("Computing collisions").stop();

    Timer::getTimer("Collision detector").stop();
}

void CollisionDetector::updateContinuousTimeCollisions()
{
    for (std::list<ContinuousTimeCollision>::iterator i = m_collisions->begin(); i != m_collisions->end(); i++)
    {
        ContinuousTimeCollision& col = *i;
        if (col.isEdgeEdge())
        {
            if (!analyseCollision(col.getEdgeEdge()))
                m_collisions->erase(i--);
        }
        else if (col.isVertexFace())
        {
            if (!analyseCollision(col.getVertexFace()))
                m_collisions->erase(i--);
        }
    }
}

void CollisionDetector::getImplicitPenaltyCollisions(std::vector<EdgeEdgeProximityCollision>& edge_edge_collisions,
        std::vector<VertexFaceProximityCollision>& vertex_face_collisions)
{
    std::cerr << "IMPLICIT PENALTY COLLISION DETECTION NOT IMPLEMENTED YET";
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
        {
            const TopologicalElement* object_a = m_elements[i];
            const TopologicalElement* object_b = m_elements[j];
            appendContinuousTimeIntersection(object_a, object_b);
        }
}

void CollisionDetector::CollisionDetector::intersectContentSelf(const BVHNode& node)
{
    const uint32_t leaf_begin = node.LeafBegin();
    const uint32_t leaf_end = node.LeafEnd();
    for (uint32_t i = leaf_begin; i < leaf_end; ++i)
        for (uint32_t j = leaf_begin; j < i; ++j)
        {
            const TopologicalElement* object_a = m_elements[i];
            const TopologicalElement* object_b = m_elements[j];
            appendContinuousTimeIntersection(object_a, object_b);
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
    // else do nothing. Maybe it's time to introduce exception handling to this code.

}

void CollisionDetector::appendContinuousTimeIntersection(const YAEdge* edge_a, const YAEdge* edge_b)
{
    if (edge_a == edge_b)
        return;

    EdgeEdgeContinuousTimeCollision edgeXedge;
    edgeXedge.e0_v0 = edge_a->first();
    edgeXedge.e0_v1 = edge_a->second();
    edgeXedge.e1_v0 = edge_b->first();
    edgeXedge.e1_v1 = edge_b->second();

    if (SKIP_ROD_ROD && isRodVertex(edgeXedge.e0_v0) && isRodVertex(edgeXedge.e1_v0)) // FIXME: Detect rod-rod collisions and skip them.
        return;

    if (analyseCollision(edgeXedge))
        m_collisions->push_back(ContinuousTimeCollision(edgeXedge));
}

void CollisionDetector::appendContinuousTimeIntersection(int v_index, const YATriangle* triangle)
{
    VertexFaceContinuousTimeCollision vertexXface;
    vertexXface.v0 = v_index;
    vertexXface.f0 = triangle->first();
    vertexXface.f1 = triangle->second();
    vertexXface.f2 = triangle->third();

    // If vertex is fixed, if face is fixed, nothing to do
    if (isVertexFixed(vertexXface.v0) && isVertexFixed(vertexXface.f0) && isVertexFixed(vertexXface.f1) && isVertexFixed(
            vertexXface.f2))
        return;

    if (analyseCollision(vertexXface))
        m_collisions->push_back(ContinuousTimeCollision(vertexXface));
}

void CollisionDetector::appendContinuousTimeIntersection(const YAEdge* edge, const YATriangle* triangle)
{
    /*
     YAEdge edge_2(triangle->first(), triangle->second());
     YAEdge edge_1(triangle->third(), triangle->first());
     YAEdge edge_0(triangle->second(), triangle->third());
     appendContinuousTimeIntersection(edge, &edge_0);
     appendContinuousTimeIntersection(edge, &edge_1);
     appendContinuousTimeIntersection(edge, &edge_2);
     */
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

bool CollisionDetector::analyseCollision(EdgeEdgeContinuousTimeCollision& edgeXedge)
{
    const Vec3d& p0 = m_geodata.GetPoint(edgeXedge.e0_v0);
    const Vec3d& q0 = m_geodata.GetPoint(edgeXedge.e0_v1);
    const Vec3d& p1 = m_geodata.GetPoint(edgeXedge.e1_v0);
    const Vec3d& q1 = m_geodata.GetPoint(edgeXedge.e1_v1);

    if (p0 == p1 || q0 == q1 || q0 == p1 || p0 == q1)
        return false;

    const Vec3d& vp0 = m_geodata.GetVelocity(edgeXedge.e0_v0);
    const Vec3d& vq0 = m_geodata.GetVelocity(edgeXedge.e0_v1);
    const Vec3d& vp1 = m_geodata.GetVelocity(edgeXedge.e1_v0);
    const Vec3d& vq1 = m_geodata.GetVelocity(edgeXedge.e1_v1);

    std::vector<double> times;
    std::vector<double> errors;
    getCoplanarityTimes(p0, q0, p1, q1, p0 + m_time_step * vp0, q0 + m_time_step * vq0, p1 + m_time_step * vp1,
            q1 + m_time_step * vq1, times, errors);

    // Loop over the coplanarity times until we find a bona fide collision
    for (size_t j = 0; j < times.size(); ++j)
    {
        /*
         if (!isProperCollisionTime(times[j]))
         continue;
         */

        double dtime = times[j] * m_time_step;

        // Determine if the collision actually happens
        const Vec3d& p1col = p0 + dtime * vp0;
        const Vec3d& q1col = q0 + dtime * vq0;
        const Vec3d& p2col = p1 + dtime * vp1;
        const Vec3d& q2col = q1 + dtime * vq1;

        Vec3d c0;
        Vec3d c1;
        double sqrdist = ClosestPtSegmentSegment(p1col, q1col, p2col, q2col, edgeXedge.s, edgeXedge.t, c0, c1);

        // If, when they are coplanar, the objects are sufficiently close, register a collision
        if (sqrdist < 1.0e-12)
        {
            edgeXedge.time = times[j];

            // Compute a collision normal at the time of the collision. For a first attempt, take the cross product of the edges.
            edgeXedge.n = (q1col - p1col).cross(q2col - p2col);
            double nnorm = edgeXedge.n.norm();

            // If the edges happen to be parallel
            if (nnorm == 0.0)
            {
                std::cerr << "WARNING, FALLING BACK ON OTHER CODE" << std::endl;
                // Use the pre-timestep positions of the collision points to generate a collision normal
                edgeXedge.n = ((1.0 - edgeXedge.t) * p1 + edgeXedge.t * q1).cross((1.0 - edgeXedge.s) * p0 + edgeXedge.s * q0);
                nnorm = edgeXedge.n.norm();

                // Crazy corner case! Handle it if it comes up :)
                if (edgeXedge.n.norm() == 0.0)
                {
                    bool lazy_person_implemented_this_normal_computation = false;
                    if (!lazy_person_implemented_this_normal_computation)
                        std::cerr << "YOU HIT AN UNSOPPORTED CODE PATH. REALLY DEGENERATE." << std::endl;
                    assert(lazy_person_implemented_this_normal_computation);
                    exit(1);
                }
            }

            edgeXedge.n /= nnorm;

            // Make sure the normal produces a negative relative velocity
            Vec3d relvel = computeRelativeVelocity(edgeXedge.e0_v0, edgeXedge.e0_v1, edgeXedge.e1_v0, edgeXedge.e1_v1,
                    edgeXedge.s, edgeXedge.t);
            if (edgeXedge.n.dot(relvel) > 0.0)
                edgeXedge.n *= -1.0;

            return true;
        }
    }
    return false;
}

bool CollisionDetector::analyseCollision(VertexFaceContinuousTimeCollision& vertexXface)
{
    const Vec3d& p = m_geodata.GetPoint(vertexXface.v0);
    const Vec3d& f0 = m_geodata.GetPoint(vertexXface.f0);
    const Vec3d& f1 = m_geodata.GetPoint(vertexXface.f1);
    const Vec3d& f2 = m_geodata.GetPoint(vertexXface.f2);

    const Vec3d& vp = m_geodata.GetVelocity(vertexXface.v0);
    const Vec3d& vf0 = m_geodata.GetVelocity(vertexXface.f0);
    const Vec3d& vf1 = m_geodata.GetVelocity(vertexXface.f1);
    const Vec3d& vf2 = m_geodata.GetVelocity(vertexXface.f2);

    std::vector<double> times;
    std::vector<double> errors;
    getCoplanarityTimes(p, f0, f1, f2, p + m_time_step * vp, f0 + m_time_step * vf0, f1 + m_time_step * vf1,
            f2 + m_time_step * vf2, times, errors);
    assert(times.size() == errors.size());

    for (size_t j = 0; j < times.size(); ++j)
    {
        /*
         if (!isProperCollisionTime(times[j]))
         continue;
         */

        double dtime = times[j] * m_time_step;

        // TODO: Use barycentric coordinates or point-triangle closest point < epsilon here? closest point < epsilon really just extends the triangle a little bit.
        // Determine if the collision actually happens
        const Vec3d& pcol = p + dtime * vp;
        const Vec3d& f0col = f0 + dtime * vf0;
        const Vec3d& f1col = f1 + dtime * vf1;
        const Vec3d& f2col = f2 + dtime * vf2;

        Vec3d cp = ClosestPtPointTriangle(pcol, f0col, f1col, f2col);

        // If, when they are coplanar, the objects are sufficiently close, register a collision
        if ((pcol - cp).squaredNorm() < 1.0e-12)
        {
            vertexXface.time = times[j];

            // Compute a collision normal at the time of the collision. For a first attempt,
            // take the cross product of two face edges.
            vertexXface.n = (f2col - f0col).cross(f1col - f0col);
            double nnorm = vertexXface.n.norm();

            // If the normal is of magnitude 0 (triangle degenerates into a line)
            if (nnorm == 0.0)
            {
                bool lazy_person_implemented_this_normal_computation = false;
                if (!lazy_person_implemented_this_normal_computation)
                    std::cerr
                            << "\033[31;1mWARNING IN BRIDSON STEPPER:\033[m YOU HIT AN UNSOPPORTED CODE PATH. REALLY DEGENERATE."
                            << std::endl;
                assert(lazy_person_implemented_this_normal_computation);
                exit(1);
            }

            vertexXface.n /= nnorm;
            assert(fabs(vertexXface.n.norm() - 1.0) < 1.0e-6);

            Barycentric(f0col, f1col, f2col, pcol, vertexXface.u, vertexXface.v, vertexXface.w);
            // Barycentric coords could be outside of [0,1] right now because we've extended the triangles a little bit
            assert(approxEq(vertexXface.u + vertexXface.v + vertexXface.w, 1.0));

            Vec3d relvel = computeRelativeVelocity(vertexXface.v0, vertexXface.f0, vertexXface.f1, vertexXface.f2,
                    vertexXface.u, vertexXface.v, vertexXface.w);
            if (vertexXface.n.dot(relvel) > 0.0)
                vertexXface.n *= -1.0;

            return true;
        }
    }
    return false;
}

bool CollisionDetector::isVertexFixed(int vert_idx) const
{
    return m_geodata.GetMass(vert_idx) == std::numeric_limits<double>::infinity();
}

bool CollisionDetector::isRodVertex(int vert) const
{
    // Is a vertex if index is less than start of object vertices in global array
    return vert < m_geodata.GetObjStart();
}

Vec3d CollisionDetector::computeRelativeVelocity(const int& idxa0, const int& idxa1, const int& idxb0, const int& idxb1,
        const double& s, const double& t)
{
    const Vec3d& v0 = m_geodata.GetVelocity(idxa0);
    const Vec3d& v1 = m_geodata.GetVelocity(idxa1);
    const Vec3d& v2 = m_geodata.GetVelocity(idxb0);
    const Vec3d& v3 = m_geodata.GetVelocity(idxb1);

    return ((1.0 - t) * v2 + t * v3) - ((1.0 - s) * v0 + s * v1);
}

Vec3d CollisionDetector::computeRelativeVelocity(const int& vrtidx, const int& fcidx0, const int& fcidx1, const int& fcidx2,
        const double& u, const double& v, const double& w)
{
    const Vec3d& vp = m_geodata.GetVelocity(vrtidx);
    const Vec3d& vt0 = m_geodata.GetVelocity(fcidx0);
    const Vec3d& vt1 = m_geodata.GetVelocity(fcidx1);
    const Vec3d& vt2 = m_geodata.GetVelocity(fcidx2);

    return vp - (u * vt0 + v * vt1 + w * vt2);
}

}
