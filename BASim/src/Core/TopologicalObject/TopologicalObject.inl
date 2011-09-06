/**
 * \file TopologicalObject.inl
 *
 * \author miklos@cs.columbia.edu
 * \date 09/12/2009
 * \author smith@cs.columbia.edu
 * \date 03/14/2010
 */



inline TopologicalObject::vertex_iter 
TopologicalObject::vertices_begin()
{
  vertex_iter first_position(this, VertexHandle(-1));
  ++first_position;
  return first_position;
}

inline TopologicalObject::const_vertex_iter
TopologicalObject::vertices_begin() const
{
  vertex_iter first_position(this, VertexHandle(-1));
  ++first_position;
  return first_position;
}

inline TopologicalObject::vertex_iter 
TopologicalObject::vertices_end()
{
  return vertex_iter(this, VertexHandle(numVertexSlots()));
}

inline TopologicalObject::const_vertex_iter
TopologicalObject::vertices_end() const
{
  return vertex_iter(this, VertexHandle(numVertexSlots()));
}

inline TopologicalObject::edge_iter 
TopologicalObject::edges_begin()
{
  edge_iter iter(this, EdgeHandle(-1));
  ++iter;
  return iter;
}

inline TopologicalObject::const_edge_iter 
TopologicalObject::edges_begin() const
{
  edge_iter iter(this, EdgeHandle(-1));
  ++iter;
  return iter;
}


inline TopologicalObject::edge_iter 
TopologicalObject::edges_end()
{
  return edge_iter(this, EdgeHandle(numEdgeSlots()));
}

inline TopologicalObject::const_edge_iter 
TopologicalObject::edges_end() const
{
  return edge_iter(this, EdgeHandle(numEdgeSlots()));
}

inline TopologicalObject::face_iter 
TopologicalObject::faces_begin()
{
  FaceIterator fit(this, FaceHandle(-1));
  ++fit;
  return fit;
}

inline TopologicalObject::const_face_iter 
TopologicalObject::faces_begin() const
{
  FaceIterator fit(this, FaceHandle(-1));
  ++fit;
  return fit;
}

inline TopologicalObject::face_iter 
TopologicalObject::faces_end()
{
  return face_iter(this, FaceHandle(numFaceSlots()));
}

inline TopologicalObject::const_face_iter 
TopologicalObject::faces_end() const
{
  return face_iter(this, FaceHandle(numFaceSlots()));
}

inline TopologicalObject::tet_iter 
TopologicalObject::tets_begin()
{
  TetIterator tit(this, TetHandle(-1));
  ++tit;
  return tit;
}

inline TopologicalObject::const_tet_iter 
TopologicalObject::tets_begin() const
{
  TetIterator tit(this, TetHandle(-1));
  ++tit;
  return tit;
}

inline TopologicalObject::tet_iter 
TopologicalObject::tets_end()
{
  return tet_iter(this, tet_handle(numTetSlots()));
}

inline TopologicalObject::const_tet_iter 
TopologicalObject::tets_end() const
{
  return tet_iter(this, tet_handle(numTetSlots()));
}

inline VertexEdgeIterator
TopologicalObject::ve_iter(const VertexHandle& vh)
{
  return VertexEdgeIterator(this, vh);
}

inline const VertexEdgeIterator
TopologicalObject::ve_iter(const VertexHandle& vh) const
{
  return VertexEdgeIterator(this, vh);
}

inline EdgeVertexIterator
TopologicalObject::ev_iter(const EdgeHandle& eh)
{
  return EdgeVertexIterator(this, eh);
}

inline const EdgeVertexIterator
TopologicalObject::ev_iter(const EdgeHandle& eh) const
{
  return EdgeVertexIterator(this, eh);
}

inline VertexVertexIterator
TopologicalObject::vv_iter(const VertexHandle& vh)
{
  return VertexVertexIterator(this, vh);
}

inline const VertexVertexIterator
TopologicalObject::vv_iter(const VertexHandle& vh) const
{
  return VertexVertexIterator(this, vh);
}



inline FaceVertexIterator
TopologicalObject::fv_iter(const FaceHandle& fh)
{
  return FaceVertexIterator(this, fh);
}

inline const FaceVertexIterator
TopologicalObject::fv_iter(const FaceHandle& fh) const
{
  return FaceVertexIterator(this, fh);
}

inline FaceEdgeIterator
TopologicalObject::fe_iter(const FaceHandle& fh)
{
  return FaceEdgeIterator(this, fh);
}

inline const FaceEdgeIterator
TopologicalObject::fe_iter(const FaceHandle& fh) const
{
  return FaceEdgeIterator(this, fh);
}

inline VertexFaceIterator
TopologicalObject::vf_iter(const VertexHandle& vh)
{
  return VertexFaceIterator(this, vh);
}

inline const VertexFaceIterator
TopologicalObject::vf_iter(const VertexHandle& vh) const
{
  return VertexFaceIterator(this, vh);
}

inline EdgeFaceIterator
TopologicalObject::ef_iter(const EdgeHandle& eh)
{
  return EdgeFaceIterator(this, eh);
}

inline const EdgeFaceIterator
TopologicalObject::ef_iter(const EdgeHandle& eh) const
{
  return EdgeFaceIterator(this, eh);
}

