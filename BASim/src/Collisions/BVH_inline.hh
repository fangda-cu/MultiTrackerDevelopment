// Base on PantaRay's BVH_inline.hh

namespace BASim
{

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

