/*
 * BVH.cc
 *
 *  Created on: 23/05/2011
 *      Author: jaubry
 */

#include "BVH.hh"

namespace BASim
{

void BVHBuilder::build(GeometryBBoxFunctor& bboxes, BVH* bvh)
{
    const uint32_t n = bboxes.size();

    {
        BBoxType node_bbox = compute_bbox(bboxes, 0u, n);

        m_stack.push(StackNode(0u, // node index
                0u, // begin index
                n, // end index
                0u, // node depth
                node_bbox)); // node bbox
    }

    std::vector<BVHNode>& nodes = bvh->GetNodeVector();
    nodes.resize(1u);

    while (m_stack.empty() == false)
    {
        StackNode node = m_stack.top();
        m_stack.pop();

        const BBoxType node_bbox = compute_bbox(bboxes, node.m_begin, node.m_end);

        if (node.m_end - node.m_begin <= m_max_leaf_size)
        {
            nodes[node.m_node_index] = BVHNode(node_bbox, node.m_begin, node.m_end);
            continue;
        }
        else
        {
            BBoxType kd_bbox = presplit(node_bbox, node.m_kd_bbox);

            Vector_Type edge = kd_bbox.max - kd_bbox.min;

            uint32_t split_dim = uint32_t(std::max_element(&edge[0], &edge[0] + 3) - &edge[0]);

            Scalar split_plane = (kd_bbox.max[split_dim] + kd_bbox.min[split_dim]) * 0.5;

            uint32_t split_index = partition(bboxes, node.m_begin, node.m_end, split_dim, split_plane);

            if (split_index == node.m_begin || split_index == node.m_end)
            {
                edge = node_bbox.max - node_bbox.min;
                split_dim = uint32_t(std::max_element(&edge[0], &edge[0] + 3) - &edge[0]);
                split_plane = (node_bbox.max[split_dim] + node_bbox.min[split_dim]) * 0.5;
                split_index = partition(bboxes, node.m_begin, node.m_end, split_dim, split_plane);
            }

            if (split_index == node.m_begin || split_index == node.m_end)
            {
                edge = node_bbox.max - node_bbox.min;
                split_dim = uint32_t(std::max_element(&edge[0], &edge[0] + 3) - &edge[0]);

                Scalar mean = 0.0;
                for (uint32_t i = node.m_begin; i < node.m_end; i++)
                    mean += (bboxes[i].min[split_dim] + bboxes[i].max[split_dim]) * 0.5;

                split_plane = mean / Scalar(node.m_begin - node.m_end);
                split_index = partition(bboxes, node.m_begin, node.m_end, split_dim, split_plane);
            }

            if (split_index == node.m_begin || split_index == node.m_end)
            {
                nodes[node.m_node_index] = BVHNode(node_bbox, node.m_begin, node.m_end);
                continue;
            }

            const uint32_t child_index = uint32_t(nodes.size());
            nodes.resize(nodes.size() + 2u);

            nodes[node.m_node_index] = BVHNode(node_bbox, child_index);

            BBoxType right_bbox = node.m_kd_bbox;
            right_bbox.min[split_dim] = split_plane;
            StackNode right_node(child_index + 1u, split_index, node.m_end, node.m_depth + 1u, right_bbox);
            m_stack.push(right_node);

            BBoxType left_bbox = node.m_kd_bbox;
            left_bbox.max[split_dim] = split_plane;
            StackNode left_node(child_index, node.m_begin, split_index, node.m_depth + 1u, left_bbox);
            m_stack.push(left_node);
        }
    }
}

BBoxType BVHBuilder::presplit(const BBoxType& node_bbox, const BBoxType& kd_bbox)
{
    int tests[3];
    for (uint32_t i = 0; i < 3; i++)
        tests[i] = 8;

    BBoxType out_bbox = kd_bbox;
    while (tests[0] && tests[1] && tests[2])
    {
        const Vector_Type edge = out_bbox.max - out_bbox.min;
        const uint32_t split_dim = uint32_t(std::max_element(&edge[0], &edge[0] + 3) - &edge[0]);
        const Scalar split_plane = (out_bbox.min[split_dim] + out_bbox.max[split_dim]) * 0.5;

        tests[split_dim]--;
        if (split_plane < node_bbox.min[split_dim])
        {
            out_bbox.min[split_dim] = split_plane;
            continue;
        }
        if (split_plane > node_bbox.max[split_dim])
        {
            out_bbox.max[split_dim] = split_plane;
            continue;
        }
        return out_bbox;
    }
    return out_bbox;
}

void swap(BVH& a, BVH& b)
{
    std::swap(a.GetNodeVector(), b.GetNodeVector());
}

bool do_overlap(const BBoxType& bbox1, const BBoxType& bbox2)
{
    if (bbox1.max.x() < bbox2.min.x() || bbox2.max.x() < bbox1.min.x() || bbox1.max.y() < bbox2.min.y() || bbox2.max.y()
            < bbox1.min.y() || bbox1.max.z() < bbox2.min.z() || bbox2.max.z() < bbox1.min.z())
        return false;
    else
        return true;
}

bool is_contained(const BBoxType& bbox1, const BBoxType& bbox2, const Scalar tol)
{
    if (bbox1.max.x() > bbox2.max.x() + tol || bbox1.min.x() < bbox2.min.x() - tol || bbox1.max.y() > bbox2.max.y() + tol
            || bbox1.min.y() < bbox2.min.y() - tol || bbox1.max.z() > bbox2.max.z() + tol || bbox1.min.z() < bbox2.min.z()
            - tol)
        return false;
    else
        return true;
}

BBoxType intersection(const BBoxType& bbox1, const BBoxType& bbox2)
{
    BBoxType bb;

    bb.min = max(bbox1.min, bbox2.min);
    bb.max = min(bbox1.max, bbox2.max);

    return bb;
}

bool is_left(const BBoxType& bbox, const uint32_t axis, const Scalar pivot)
{
    switch (axis)
    {
    case 0:
        return bbox.min.x() + bbox.max.x() < pivot * 2.0;
    case 1:
        return bbox.min.y() + bbox.max.y() < pivot * 2.0;
    case 2:
        return bbox.min.z() + bbox.max.z() < pivot * 2.0;
    }
    return false; // Never reached, this is just to suppress a warning.
}

void insert(BBoxType& bbox, const BBoxType& bbox2)
{
    bbox.min = min(bbox.min, bbox2.min);
    bbox.max = max(bbox.max, bbox2.max);

}

BBoxType merge(const BBoxType& bbox1, const BBoxType& bbox2)
{
    BBoxType bb;
    bb.min = min(bbox1.min, bbox2.min);
    bb.max = max(bbox1.max, bbox2.max);
    return bb;
}

uint32_t partition(GeometryBBoxFunctor& bboxes, const uint32_t begin, const uint32_t end, const uint32_t axis,
        const Scalar pivot)
{
    uint32_t i = begin;
    uint32_t j = end - 1;
    while (i != j)
    {
        if (is_left(bboxes[i], axis, pivot) == false)
        {
            bboxes.swap(i, j);
            --j;
        }
        else
            ++i;
    }
    if (is_left(bboxes[i], axis, pivot) == true)
        return ++i;
    else
        return i;
}

BBoxType compute_bbox(GeometryBBoxFunctor& bboxes, const uint32_t begin, const uint32_t end)
{
    BBoxType node_bbox;
    for (uint32_t i = begin; i < end; i++)
        insert(node_bbox, bboxes[i]);

    return node_bbox;
}

}
