// Based on PantaRay's BVH.hh

#ifndef TREES_BVH_HH
#define TREES_BVH_HH 1

#include <vector>
#include <stack>
#include <limits>
#include <algorithm>

#include "BVHNode.hh"

namespace BASim
{

typedef BoundingBox<Scalar> BBoxType;

class BVH
{
public:
	typedef BVHNode Node_Type;
	typedef std::vector<Node_Type> NodeVector_Type;

	/// empty constructor
	BVH()
	{
	}

	/// returns the size of this object in bytes
	size_t ByteSize() const
	{
		return sizeof(Node_Type) * m_nodes.size() + sizeof(BVH);
	}

	/// get node vector
	const std::vector<Node_Type>& GetNodeVector() const
	{
		return m_nodes;
	}

	/// get node vector
	std::vector<Node_Type>& GetNodeVector()
	{
		return m_nodes;
	}

	/// get nodes pointer
	const Node_Type* GetNodes() const
	{
		return &m_nodes[0];
	}

	/// get nodes pointer
	Node_Type* GetNodes()
	{
		return &m_nodes[0];
	}

	/// get the i-th node
	const Node_Type& GetNode(const uint32_t i) const
	{
		return m_nodes[i];
	}

	/// get the i-th node
	Node_Type& GetNode(const uint32_t i)
	{
		return m_nodes[i];
	}

private:
	std::vector<Node_Type> m_nodes; ///< bvh nodes
};

//template <uint32_t DIM>
inline void swap(BVH& a, BVH& b)
{
	std::swap(a.GetNodeVector(), b.GetNodeVector());
}

/// A fast middle BVH builder
class MiddleBVHBuilder
{
public:
	typedef Point<Scalar> Vector_Type;
	typedef BVH BVH_Type;
	typedef BVHNode BVHNode_Type;

	/// empty constructor
	MiddleBVHBuilder() :
		m_max_leaf_size(1u)
	{
	}

	/// set maximum leaf size
	void SetMaxLeafSize(const uint32_t max_leaf_size)
	{
		m_max_leaf_size = max_leaf_size;
	}

	/// build a bvh
	///
	/// The BboxVectorT template argument must obey the following interface:
	///
	///     uint32_t BboxVectorT::size() const;
	///     Bbox BboxVectorT::operator[] (const uint32_t i);
	///     void BboxVectorT::swap(uint32_t i, uint32_t j);
	///
	/// \param bboxes       a BboxVectorT, representing the bboxes to be sorted
	/// \param bvh          the output bvh
	///
	template<typename BBoxFunctorT>
	void build(BBoxFunctorT& bboxes, BVH_Type* bvh);

private:
	BBoxType presplit(const BBoxType& node_bbox, const BBoxType& kd_bbox);

	struct StackNode
	{
		StackNode()
		{
		}
		StackNode(const uint32_t node, const uint32_t begin, const uint32_t end, const uint32_t depth, const BBoxType& kd_bbox) :
			m_node_index(node), m_begin(begin), m_end(end), m_depth(depth), m_kd_bbox(kd_bbox)
		{
		}

		uint32_t m_node_index;
		uint32_t m_begin;
		uint32_t m_end;
		uint32_t m_depth;
		BBoxType m_kd_bbox;
	};

	BVH_Type* m_bvh; ///< output bvh
	std::stack<StackNode> m_stack; ///< internal stack
	uint32_t m_max_leaf_size;///< maximum leaf size
};

inline Scalar min(const Scalar x, const Scalar y)
{
	return x < y ? x : y;
}
inline Scalar max(const Scalar x, const Scalar y)
{
	return x > y ? x : y;
}

inline bool do_overlap(const BoundingBox<Scalar>& bbox1, const BoundingBox<Scalar>& bbox2)
{
	if (bbox1.max.x() < bbox2.min.x() || bbox2.max.x() < bbox1.min.x() || bbox1.max.y() < bbox2.min.y() || bbox2.max.y() < bbox1.min.y()
			|| bbox1.max.z() < bbox2.min.z() || bbox2.max.z() < bbox1.min.z())
		return false;
	else
		return true;
}

inline bool is_contained(const BoundingBox<Scalar>& bbox1, const BoundingBox<Scalar>& bbox2, const Scalar tol = 0.0f)
{
	if (bbox1.max.x() > bbox2.max.x() + tol || bbox1.min.x() < bbox2.min.x() - tol || bbox1.max.y() > bbox2.max.y() + tol || bbox1.min.y()
			< bbox2.min.y() - tol || bbox1.max.z() > bbox2.max.z() + tol || bbox1.min.z() < bbox2.min.z() - tol)
		return false;
	else
		return true;
}

inline BoundingBox<Scalar> intersection(const BoundingBox<Scalar>& bbox1, const BoundingBox<Scalar>& bbox2)
{
	BoundingBox<Scalar> bb;

	bb.min = Max(bbox1.min, bbox2.min);
	bb.max = Min(bbox1.max, bbox2.max);

	return bb;
}

///
/// A bbox-index pair vector
///
struct BBoxVector
{
	/// clear
	void clear()
	{
		m_bbox.clear();
		m_index.clear();
	}

	/// reserve
	void reserve(const uint32_t size)
	{
		m_bbox.reserve(size);
		m_index.reserve(size);
	}

	/// push back a new bbox
	void push_back(const BBoxType& bbox)
	{
		m_index.push_back((uint32_t) m_bbox.size());
		m_bbox.push_back(bbox);
	}
	/// push back a new bbox, with an accompanying index
	void push_back(const BBoxType& bbox, const uint32_t index)
	{
		m_index.push_back(index);
		m_bbox.push_back(bbox);
	}

	/// return vector size
	uint32_t size() const
	{
		return (uint32_t) m_bbox.size();
	}

	/// return i-th bbox
	const BBoxType& operator[](const uint32_t i) const
	{
		return m_bbox[i];
	}

	/// intersect the i-th bbox
	bool Intersect(const uint32_t i, const BBoxType& crop_bbox, BBoxType& result) const
	{
		if (do_overlap(m_bbox[i], crop_bbox) == false)
			return false;

		result = intersection(m_bbox[i], crop_bbox);
		return true;
	}

	/// return i-th index
	uint32_t index(const uint32_t i) const
	{
		return m_index[i];
	}

	void SetIndex(const uint32_t i, const uint32_t value)
	{
		m_index[i] = value;
	}

	const std::vector<uint32_t>& GetIndices() const
	{
		return m_index;
	}

	/// swap i-th and j-th element
	void swap(const uint32_t i, const uint32_t j)
	{
		std::swap(m_bbox[i], m_bbox[j]);
		std::swap(m_index[i], m_index[j]);
	}

	BBoxType* begin()
	{
		return m_bbox.size() ? (&m_bbox[0]) : NULL;
	}
	BBoxType* end()
	{
		return m_bbox.size() ? (&m_bbox[0]) + m_bbox.size() : NULL;
	}

private:
	std::vector<BBoxType> m_bbox;
	std::vector<uint32_t> m_index;
};

///
/// An SAH-based bvh builder for 3d bboxes
///
class SAH_BVHBuilder
{
public:
	typedef BVH Bvh_type;
	typedef BVHNode Bvh_node_type;

	/// constructor
	SAH_BVHBuilder() :
		m_max_leaf_size(4u), m_force_splitting(true), m_force_alignment(false)
	{
	}

	/// set bvh parameters
	void SetMaxLeafSize(const uint32_t max_leaf_size)
	{
		m_max_leaf_size = max_leaf_size;
	}

	/// set force splitting
	void SetForceSplitting(const bool flag)
	{
		m_force_splitting = flag;
	}

	/// set force 'max leaf size'-aligned splits
	void SetForceAlignment(const bool flag)
	{
		m_force_alignment = flag;
	}

	/// build
	///
	/// Iterator is supposed to dereference to a BBoxType
	///
	/// \param begin            first point
	/// \param end              last point
	/// \param bvh              output bvh
	template<typename Iterator>
	void build(Iterator begin, Iterator end, Bvh_type* bvh);

	/// remapped point index
	uint32_t index(const uint32_t i) const
	{
		return m_entities[m_indices[0][i]].m_index;
	}

private:
	class Predicate;
	class IndexSortPredicate;

	struct Entity
	{
		void setTagBit()
		{
			m_tag_bit = true;
		}
		void clearTagBit()
		{
			m_tag_bit = false;
		}
		bool hasTagBit() const
		{
			return m_tag_bit;
		}

		BBoxType m_bbox;
		uint32_t m_index;
		bool m_tag_bit;
	};

	struct EntityBBoxVector
	{
		/// constructor
		EntityBBoxVector(const std::vector<Entity>& entities, const std::vector<int>& indices) :
			m_entities(entities), m_indices(indices)
		{
		}

		/// return vector size
		uint32_t size() const
		{
			return uint32_t(m_entities.size());
		}

		/// return i-th bbox
		const BBoxType& operator[](const uint32_t i) const
		{
			return m_entities[m_indices[i]].m_bbox;
		}

		/// return i-th index
		uint32_t index(const uint32_t i) const
		{
			return m_entities[m_indices[i]].m_index;
		}

	private:
		const std::vector<Entity>& m_entities;
		const std::vector<int>& m_indices;
	};

	struct Node
	{
		BBoxType m_bbox;
		uint32_t m_begin;
		uint32_t m_end;
		uint32_t m_node;
		uint32_t m_depth;
	};
	typedef std::stack<Node> Node_stack;

	Scalar area(const BBoxType& bbox);

	void compute_bbox(const uint32_t begin, const uint32_t end, BBoxType& bbox);

	bool find_best_split(const Node& node, uint32_t& pivot, BBoxType& l_bb, BBoxType& r_bb);

	struct Bvh_partitioner;

	uint32_t m_max_leaf_size;
	bool m_force_splitting;
	bool m_force_alignment;
	std::vector<int> m_indices[3]; // indices into entities, one vector per dimension
	std::vector<Entity> m_entities;
	std::vector<BBoxType> m_bboxes;
	std::vector<int> m_tmpInt; // temp array used during build
};

//template <uint32_t DIM>
void build_parents(const BVH* bvh, std::vector<uint32_t>& parents);

/// performs a sanity check on a bvh, fixing it if it's in an incosistent state.
/// returns false if no error was found, true otherwise.
template<typename BBoxVectorT>
bool fix(BVH* bvh, const BBoxVectorT& bbox_vector, const Scalar tol = 1.0e-6f);

}


#include "BVH_inline.hh"

#endif
