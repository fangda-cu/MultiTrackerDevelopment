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
    const uint32_t m_max_leaf_size;///< maximum leaf size
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
    if (bbox1.max.x() < bbox2.min.x() || bbox2.max.x() < bbox1.min.x() || bbox1.max.y() < bbox2.min.y() || bbox2.max.y()
            < bbox1.min.y() || bbox1.max.z() < bbox2.min.z() || bbox2.max.z() < bbox1.min.z())
        return false;
    else
        return true;
}

inline bool is_contained(const BoundingBox<Scalar>& bbox1, const BoundingBox<Scalar>& bbox2, const Scalar tol = 0.0f)
{
    if (bbox1.max.x() > bbox2.max.x() + tol || bbox1.min.x() < bbox2.min.x() - tol || bbox1.max.y() > bbox2.max.y() + tol
            || bbox1.min.y() < bbox2.min.y() - tol || bbox1.max.z() > bbox2.max.z() + tol || bbox1.min.z() < bbox2.min.z()
            - tol)
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

typedef BoundingBox<Scalar> BBoxType;

inline bool is_left(const BBoxType& bbox, const uint32_t axis, const float pivot)
{
    switch (axis)
    {
    case 0:
        return bbox.min.x() + bbox.max.x() < pivot * 2.0f;
    case 1:
        return bbox.min.y() + bbox.max.y() < pivot * 2.0f;
    case 2:
        return bbox.min.z() + bbox.max.z() < pivot * 2.0f;
    }
    return false; // Never reached, this is just to suppress a warning.
}

inline void insert(BBoxType& bbox, const BBoxType& bbox2)
{
    bbox.min = Min(bbox.min, bbox2.min);
    bbox.max = Max(bbox.max, bbox2.max);

}

inline BBoxType merge(const BBoxType& bbox1, const BBoxType& bbox2)
{
    BBoxType bb;
    bb.min = Min(bbox1.min, bbox2.min);
    bb.max = Max(bbox1.max, bbox2.max);
    return bb;
}

template<typename BBoxFunctorT>
uint32_t partition(BBoxFunctorT& bboxes, const uint32_t begin, const uint32_t end, const uint32_t axis, const float pivot)
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

template<typename BBoxFunctorT>
BBoxType compute_bbox(BBoxFunctorT& bboxes, const uint32_t begin, const uint32_t end)
{
    BBoxType node_bbox;
    for (uint32_t i = begin; i < end; i++)
        insert(node_bbox, bboxes[i]);

    return node_bbox;
}

// split the kd_bbox a few times along its longest direction, until
// the next split would fall inside the node_bbox.
inline BBoxType MiddleBVHBuilder::presplit(const BBoxType& node_bbox, const BBoxType& kd_bbox)
{
    int tests[3];
    for (uint32_t i = 0; i < 3; i++)
        tests[i] = 8;

    BBoxType out_bbox = kd_bbox;
    while (tests[0] && tests[1] && tests[2])
    {
        const Point<Scalar> edge = out_bbox.max - out_bbox.min;
        const uint32_t split_dim = uint32_t(std::max_element(&edge[0], &edge[0] + 3) - &edge[0]);
        const float split_plane = (out_bbox.min[split_dim] + out_bbox.max[split_dim]) * 0.5f;

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

// build a bvh
//
// The BBoxVectorT template argument must obey the following interface:
//
//     BBox BBoxVectorT::operator[] (const uint32_t i);
//     void BBoxVectorT::swap(uint32_t i, uint32_t j);
//
// \param bboxes       a BBoxVectorT, representing the bboxes to be sorted
// \param bvh          the output bvh
//
template<typename BBoxFunctorT>
void MiddleBVHBuilder::build(BBoxFunctorT& bboxes, BVH_Type* bvh)
{
    //
    // The following code is essentially a 3d quicksort, which sorts all
    // the input bboxes in place, and stores all the splits as internal nodes
    // of the output bvh.
    //
    const uint32_t n = bboxes.size();

    // build an initial node containing all input bboxes and push it into the processing stack
    {
        BBoxType node_bbox = compute_bbox(bboxes, 0u, n);

        m_stack.push(StackNode(0u, // node index
                0u, // begin index
                n, // end index
                0u, // node depth
                node_bbox)); // node bbox
    }

    std::vector<BVHNode_Type>& nodes = bvh->GetNodeVector();
    nodes.resize(1u);

    while (m_stack.empty() == false)
    {
        // extract the next node to be processed from the stack
        StackNode node = m_stack.top();
        m_stack.pop();

        const BBoxType node_bbox = compute_bbox(bboxes, node.m_begin, node.m_end);

        // check if we can stop splitting and make a leaf
        if (node.m_end - node.m_begin <= m_max_leaf_size)
        {
            // make a leaf
            nodes[node.m_node_index] = BVHNode_Type(node_bbox, node.m_begin, node.m_end);
            continue;
        }
        else
        {
            // partition along the longest axis, in the middle
            BBoxType kd_bbox = presplit(node_bbox, node.m_kd_bbox);

            /*const*/
            Vector_Type edge = kd_bbox.max - kd_bbox.min;
            /*const*/
            uint32_t split_dim = uint32_t(std::max_element(&edge[0], &edge[0] + 3) - &edge[0]);

            /*const*/
            float split_plane = (kd_bbox.max[split_dim] + kd_bbox.min[split_dim]) * 0.5f;

            /*const*/
            uint32_t split_index = partition(bboxes, node.m_begin, node.m_end, split_dim, split_plane);

            // check if the split was unsuccessful
            if (split_index == node.m_begin || split_index == node.m_end)
            {
                // try a middle split of the real bbox
                edge = node_bbox.max - node_bbox.min;
                split_dim = uint32_t(std::max_element(&edge[0], &edge[0] + 3) - &edge[0]);
                split_plane = (node_bbox.max[split_dim] + node_bbox.min[split_dim]) * 0.5f;
                split_index = partition(bboxes, node.m_begin, node.m_end, split_dim, split_plane);
            }

            // check if the split was unsuccessful
            if (split_index == node.m_begin || split_index == node.m_end)
            {
                // try a mean split along the longest dimension
                edge = node_bbox.max - node_bbox.min;
                split_dim = uint32_t(std::max_element(&edge[0], &edge[0] + 3) - &edge[0]);

                float mean = 0.0f;
                for (uint32_t i = node.m_begin; i < node.m_end; i++)
                    mean += (bboxes[i].min[split_dim] + bboxes[i].max[split_dim]) * 0.5f;

                split_plane = mean / float(node.m_begin - node.m_end);
                split_index = partition(bboxes, node.m_begin, node.m_end, split_dim, split_plane);
            }

            // check if the split was unsuccessful
            if (split_index == node.m_begin || split_index == node.m_end)
            {
                // make a leaf
                nodes[node.m_node_index] = BVHNode_Type(node_bbox, node.m_begin, node.m_end);
                continue;
            }

            // split the node
            const uint32_t child_index = uint32_t(nodes.size());
            nodes.resize(nodes.size() + 2u);

            // store the newly defined inner node
            nodes[node.m_node_index] = BVHNode_Type(node_bbox, child_index);

            // push right child into the processing stack
            BBoxType right_bbox = node.m_kd_bbox;
            right_bbox.min[split_dim] = split_plane;
            StackNode right_node(child_index + 1u, split_index, node.m_end, node.m_depth + 1u, right_bbox);
            m_stack.push(right_node);

            // push left child into the processing stack
            BBoxType left_bbox = node.m_kd_bbox;
            left_bbox.max[split_dim] = split_plane;
            StackNode left_node(child_index, node.m_begin, split_index, node.m_depth + 1u, left_bbox);
            m_stack.push(left_node);
        }
    }
}

//  template <uint32_t DIM>
inline void build_parents(const BVH* bvh, const uint32_t node_index, std::vector<uint32_t>& parents)
{
    const BVHNode* nodes = bvh->GetNodes();
    const BVHNode* node = nodes + node_index;

    if (node->IsLeaf() == false)
    {
        const uint32_t child_index = node->ChildIndex();
        parents[child_index + 0u] = node_index;
        parents[child_index + 1u] = node_index;
        build_parents(bvh, child_index + 0u, parents);
        build_parents(bvh, child_index + 1u, parents);
    }
}
//  template <uint32_t DIM>
inline void build_parents(const BVH* bvh, std::vector<uint32_t>& parents)
{
    const uint32_t nnodes = uint32_t(bvh->GetNodeVector().size());

    parents.resize(nnodes);

    build_parents(bvh, 0u, parents);

    parents[0] = uint32_t(-1);
}

inline bool fix(BVH* bvh, const float tol, const uint32_t node_index, const std::vector<BBoxType>& bbox_vector)
{
    BVHNode* nodes = bvh->GetNodes();
    BVHNode* node = nodes + node_index;

    if (node->IsLeaf())
    {
        const uint32_t leaf_begin = node->LeafBegin();
        const uint32_t leaf_end = node->LeafEnd();

        BBoxType bbox;
        for (uint32_t leaf_index = leaf_begin; leaf_index < leaf_end; ++leaf_index)
            bbox = merge(bbox, bbox_vector[leaf_index]);

        if (is_contained(bbox, node->BBox(), tol) == false)
        {
            node->SetBBox(bbox);
            return true;
        }
        return false;
    }
    else
    {
        const uint32_t child_index = node->ChildIndex();

        const bool result = fix(bvh, tol, child_index + 0u, bbox_vector) | fix(bvh, tol, child_index + 1u, bbox_vector);

        const BBoxType bbox = merge(nodes[child_index + 0].BBox(), nodes[child_index + 1].BBox());

        if (is_contained(bbox, node->BBox(), tol) == false)
        {
            node->SetBBox(bbox);
            return true;
        }
        return result;
    }
}

inline bool fix(BVH* bvh, const std::vector<BBoxType>& bbox_vector, const float tol)
{
    return fix(bvh, tol, 0u, bbox_vector);
}

}

#endif
