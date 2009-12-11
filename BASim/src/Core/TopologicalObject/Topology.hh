/**
 * \file Topology.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/13/2009
 */

#ifndef TOPOLOGY_HH
#define TOPOLOGY_HH

#include "BASim/src/Core/Definitions.hh"

namespace BASim {

/** Class that holds the topology (connectivity information) of a
    vertex */
template <class T>
class VertexTopology
{
public:

  typedef typename T::vertex_handle vertex_handle;
  typedef typename T::edge_handle   edge_handle;

  VertexTopology() {}

  VertexTopology(const VertexTopology& vt)
    : m_edges(vt.m_edges)
  {}

  VertexTopology& operator= (const VertexTopology& vt)
  {
    m_edges = vt.m_edges;
    return *this;
  }

  void addEdge(const edge_handle& eh)
  {
    m_edges.push_back(eh);
  }

  edge_handle& operator[] (size_t j)
  {
    assert(j < m_edges.size());
    return m_edges[j];
  }

  const edge_handle& operator[] (size_t j) const
  {
    assert(j < m_edges.size());
    return m_edges[j];
  }

  /** Valence of the vertex */
  size_t size() const
  {
    return m_edges.size();
  }

protected:

  std::vector<edge_handle> m_edges; ///< Edges that the vertex is part of
};

/** Class that holds the topology (connectivity information) of an
    edge */
template <class T>
class EdgeTopology
{
public:

  typedef typename T::vertex_handle vertex_handle;
  typedef typename T::edge_handle   edge_handle;

  EdgeTopology()
    : m_verts(2)
  {}

  EdgeTopology(const EdgeTopology& et)
    : m_verts(et.m_verts)
  {}

  EdgeTopology& operator= (const EdgeTopology& et)
  {
    m_verts = et.m_verts;
    return *this;
  }

  const vertex_handle& getFromVertex() const
  {
    return m_verts[0];
  }

  const vertex_handle& getToVertex() const
  {
    return m_verts[1];
  }

  void setFromVertex(const vertex_handle& vh)
  {
    m_verts[0] = vh;
  }

  void setToVertex(const vertex_handle& vh)
  {
    m_verts[1] = vh;
  }

  vertex_handle& operator[] (size_t i)
  {
    assert(i < m_verts.size());
    return m_verts[i];
  }

  const vertex_handle& operator[] (size_t i) const
  {
    assert(i < m_verts.size());
    return m_verts[i];
  }

  /** Number of vertices making up the edge (always 2) */
  size_t size() const
  {
    return m_verts.size();
  }

protected:

  std::vector<vertex_handle> m_verts; ///< Vertices adjacent to the edge
};

} // namespace BASim

#endif // TOPOLOGY_HH
