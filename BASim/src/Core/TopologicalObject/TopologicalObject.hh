/**
 * \file TopologicalObject.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/12/2009
 */

#ifndef TOPOLOGICALOBJECT_HH
#define TOPOLOGICALOBJECT_HH

namespace BASim {

/** An object that represents a collection of vertices and edges with
    associated connectivity information. */
class TopologicalObject : public ObjectBase
{
public:

  typedef VertexHandle<TopologicalObject>               vertex_handle;
  typedef VertexIterator<TopologicalObject>             vertex_iter;
  typedef const VertexIterator<TopologicalObject>       const_vertex_iter;
  typedef VertexTopology<TopologicalObject>             vertex_topology;

  typedef EdgeHandle<TopologicalObject>                 edge_handle;
  typedef EdgeIterator<TopologicalObject>               edge_iter;
  typedef const EdgeIterator<TopologicalObject>         const_edge_iter;
  typedef EdgeTopology<TopologicalObject>               edge_topology;

  typedef VertexEdgeIterator<TopologicalObject>         VertexEdgeIter;
  typedef const VertexEdgeIterator<TopologicalObject>   ConstVertexEdgeIter;

  typedef EdgeVertexIterator<TopologicalObject>         EdgeVertexIter;
  typedef const EdgeVertexIterator<TopologicalObject>   ConstEdgeVertexIter;

  typedef VertexVertexIterator<TopologicalObject>       VertexVertexIter;
  typedef const VertexVertexIterator<TopologicalObject> ConstVertexVertexIter;

  TopologicalObject();

  virtual ~TopologicalObject() {}

  /** Number of vertices */
  int nv() const { return property(m_nv); }
  /** Number of edges */
  int ne() const { return property(m_ne); }

  /** Add a vertex to the object */
  vertex_handle addVertex();
  /** Add a directed edge that connects the two given vertices */
  edge_handle addEdge(const vertex_handle& v0, const vertex_handle& v1);

  /** Deletes the given vertex and associated edges */
  void deleteVertex(const vertex_handle& vertex);
  /** Deletes the given edge */
  void deleteEdge(const edge_handle& edge);

  /** Returns a handle to the vertex at the head of the edge */
  vertex_handle fromVertex(const edge_handle& eh) const;
  /** Returns a handle to the vertex at the tail of the edge */
  vertex_handle toVertex(const edge_handle& eh) const;

  /** \name Iterators */

  //@{

  /** Returns an iterator to the first vertex */
  vertex_iter vertices_begin();
  /** Const version */
  const_vertex_iter vertices_begin() const;

  /** Returns an iterator to the past-the-end vertex */
  vertex_iter vertices_end();
  /** Const version */
  const_vertex_iter vertices_end() const;

  /** Returns an iterator to the first edge */
  edge_iter edges_begin();
  /** Const version */
  const_edge_iter edges_begin() const;

  /** Returns an iterator to the past-the-end edge */
  edge_iter edges_end();
  /** Const version */
  const_edge_iter edges_end() const;

  //@}

  /** \name Circulators */

  //@{

  /** Circulator over vertices adjacent to an edge */
  VertexEdgeIter ve_iter(const vertex_handle& vh);
  /** Const version */
  ConstVertexEdgeIter ve_iter(const vertex_handle& vh) const;

  /** Circulator over edges adjacent to a vertex */
  EdgeVertexIter ev_iter(const edge_handle& eh);
  /** Const version */
  ConstEdgeVertexIter ev_iter(const edge_handle& eh) const;

  /** Circulator over vertices in the one-ring of the given vertex */
  VertexVertexIter vv_iter(const vertex_handle& vh);
  /** Const version */
  ConstVertexVertexIter vv_iter(const vertex_handle& vh) const;

  //@}

  /** \name Topology */

  //@{

  /** Returns topology associated to a vertex */
  vertex_topology& getVertexTopology(const vertex_handle& vh);
  /** Const version */
  const vertex_topology& getVertexTopology(const vertex_handle& vh) const;

  /** Returns topology associated to an edge */
  edge_topology& getEdgeTopology(const edge_handle& eh);
  /** Const version */
  const edge_topology& getEdgeTopology(const edge_handle& eh) const;

  //@}

  /** \name Vertex properties */
  //@{
  BA_CREATE_PROPERTY(VPropHandle, m_vertexProps, nv());
  //@}

  /** \name Edge properties */
  //@{
  BA_CREATE_PROPERTY(EPropHandle, m_edgeProps, ne());
  //@}

  BA_INHERIT_BASE(ObjectBase);

protected:

  PropertyContainer m_vertexProps; ///< Vertex property container
  PropertyContainer m_edgeProps;   ///< Edge property container

  ObjPropHandle<int> m_nv; ///< number of vertices
  ObjPropHandle<int> m_ne; ///< number of edges

  VPropHandle<vertex_topology> m_vertTop; ///< Topology of the vertices
  EPropHandle<edge_topology> m_edgeTop;   ///< Topology of the edges
};

#include "TopologicalObject.inl"

} // namespace BASim

#endif // TOPOLOGICALOBJECT_HH
