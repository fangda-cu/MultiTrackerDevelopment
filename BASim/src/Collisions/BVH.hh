/*
 * BVH.cc
 *
 *  Created on: 23/05/2011
 *      Author: jaubry
 */

#ifndef BVH_HH
#define BVH_HH

#include <vector>
#include <stack>
#include <limits>
#include <algorithm>

#include "BVHNode.hh"
#include "Geometry.hh"

namespace BASim
{

typedef BoundingBox<Scalar> BBoxType;
typedef BVHNode<BBoxType> BVHNodeType;

class BVH
{
public:
    typedef BVHNodeType Node_Type; // for pantaray::kNN<BVH>

    /// empty constructor
    BVH()
    {
    }

    /// returns the size of this object in bytes
    size_t ByteSize() const
    {
        return sizeof(BVHNodeType) * m_nodes.size() + sizeof(BVH);
    }

    /// get node vector
    const std::vector<BVHNodeType>& GetNodeVector() const
    {
        return m_nodes;
    }

    /// get node vector
    std::vector<BVHNodeType>& GetNodeVector()
    {
        return m_nodes;
    }

    /// get nodes pointer
    const BVHNodeType* GetNodes() const
    {
        return &m_nodes[0];
    }

    /// get nodes pointer
    BVHNodeType* GetNodes()
    {
        return &m_nodes[0];
    }

    /// get the i-th node
    const BVHNodeType& GetNode(const unsigned int i) const
    {
        return m_nodes[i];
    }

    /// get the i-th node
    BVHNodeType& GetNode(const unsigned int i)
    {
        return m_nodes[i];
    }

private:
    std::vector<BVHNodeType> m_nodes; ///< bvh nodes
};

void swap(BVH& a, BVH& b);

template<typename BBoxFunctorT>
class BVHBuilder
{
public:
    typedef BBoxType::PointType PointType;

    /// empty constructor
    BVHBuilder() :
        m_max_leaf_size(1u)
    {
    }

    void build(BBoxFunctorT& bboxes, BVH* bvh);

private:
    BBoxType presplit(const BBoxType& node_bbox, const BBoxType& kd_bbox);

    struct StackNode
    {
        StackNode()
        {
        }
        StackNode(const unsigned int node, const unsigned int begin, const unsigned int end, const unsigned int depth, const BBoxType& kd_bbox) :
            m_node_index(node), m_begin(begin), m_end(end), m_depth(depth), m_kd_bbox(kd_bbox)
        {
        }

        unsigned int m_node_index;
        unsigned int m_begin;
        unsigned int m_end;
        unsigned int m_depth;
        BBoxType m_kd_bbox;
    };

    BVH* m_bvh; ///< output bvh
    std::stack<StackNode> m_stack; ///< internal stack
    const unsigned int m_max_leaf_size;///< maximum leaf size
};

bool do_overlap(const BBoxType& bbox1, const BBoxType& bbox2);

bool is_contained(const BBoxType& bbox1, const BBoxType& bbox2, const Scalar tol = 0.0f);

BBoxType intersection(const BBoxType& bbox1, const BBoxType& bbox2);

typedef BBoxType BBoxType;

bool is_left(const BBoxType& bbox, const unsigned int axis, const Scalar pivot);

void insert(BBoxType& bbox, const BBoxType& bbox2);

BBoxType merge(const BBoxType& bbox1, const BBoxType& bbox2);

template<typename BBoxFunctorT>
unsigned int partition(BBoxFunctorT& bboxes, const unsigned int begin, const unsigned int end, const unsigned int axis, const Scalar pivot);

template<typename BBoxFunctorT>
BBoxType compute_bbox(BBoxFunctorT& bboxes, const unsigned int begin, const unsigned int end);

}

#endif
