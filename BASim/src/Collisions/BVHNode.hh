// Based on PantaRay's BVHNode.hh

#ifndef TREES_BVHNODE_HH
#define TREES_BVHNODE_HH 1

#include "BoundingBox.hh"

namespace BASim
{

typedef double Scalar;

struct BVHNode
{
	typedef Point<Scalar> Vector_Type;
	typedef BoundingBox<Scalar> BBoxType;

	/// empty constructor
	BVHNode() :
		m_end(uint32_t(-1))
	{
	}

	/// inner node constructor
	BVHNode(const BBoxType& bbox, const uint32_t child_index) :
		m_bbox(bbox), m_index(child_index), m_end(uint32_t(-1))
	{
	}

	/// leaf constructor
	BVHNode(const BBoxType& bbox, const uint32_t leaf_begin, const uint32_t leaf_end) :
		m_bbox(bbox), m_index(leaf_begin), m_end(leaf_end)
	{
	}

	/// is a leaf?
	bool IsLeaf() const
	{
		return m_end != uint32_t(-1);
	}

	/// leaf begin
	uint32_t LeafBegin() const
	{
		return m_index;
	}

	/// leaf end
	uint32_t LeafEnd() const
	{
		return m_end;
	}

	/// child index
	uint32_t ChildIndex() const
	{
		return m_index;
	}

	/// node bbox
	const BBoxType& BBox() const
	{
		return m_bbox;
	}

	/// node bbox
	BBoxType& BBox()
	{
		return m_bbox;
	}

	/// set bbox
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
