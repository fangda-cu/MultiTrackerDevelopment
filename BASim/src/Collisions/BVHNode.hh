#ifndef TREES_BVHNODE_HH
#define TREES_BVHNODE_HH 1

#include "BoundingBox.hh"

namespace BASim
{

typedef double Scalar;

struct BVHNode
{
    typedef BoundingBox<Scalar> BBoxType;
    typedef BBoxType::PointType PointType;

    BVHNode() :
        m_end(uint32_t(-1))
    {
    }

    BVHNode(const BBoxType& bbox, const uint32_t child_index) :
        m_bbox(bbox), m_index(child_index), m_end(uint32_t(-1))
    {
    }

    BVHNode(const BBoxType& bbox, const uint32_t leaf_begin, const uint32_t leaf_end) :
        m_bbox(bbox), m_index(leaf_begin), m_end(leaf_end)
    {
    }

    bool IsLeaf() const
    {
        return m_end != uint32_t(-1);
    }

    uint32_t LeafBegin() const
    {
        return m_index;
    }

    uint32_t LeafEnd() const
    {
        return m_end;
    }

    uint32_t ChildIndex() const
    {
        return m_index;
    }

    const BBoxType& BBox() const
    {
        return m_bbox;
    }

    BBoxType& BBox()
    {
        return m_bbox;
    }

    void SetBBox(const BBoxType& bbox)
    {
        m_bbox = bbox;
    }

    BBoxType m_bbox; ///< the node's bbox
    uint32_t m_index; ///< if an INNER node: the node's first child index, if LEAF: the leaf begin index
    uint32_t m_end; ///< if an INNER node: -1, if a LEAF: the leaf end index

};

}

#endif
