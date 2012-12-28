/**
 * \file TopologicalObject.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/12/2009
 * \author smith@cs.columbia.edu
 * \date 03/14/2010
 */

#ifndef TOPOLOGICALOBJECT_HH
#define TOPOLOGICALOBJECT_HH

#include "BASim/src/Core/ObjectBase.hh"
#include "BASim/src/Core/TopologicalObject/TopObjHandles.hh"
#include "BASim/src/Core/TopologicalObject/TopObjIterators.hh"
#include "BASim/src/Core/TopologicalObject/Topology.hh"
#include "BASim/src/Core/TopologicalObject/IncidenceMatrix.hh"
#include "BASim/src/Core/Definitions.hh"
#include <queue>
namespace BASim {

class TopObjPropertyBase;

/** An object that represents a collection of vertices, edges, faces and tets
   with associated connectivity information.  */
class TopologicalObject : public ObjectBase
{
public:

  //Allow iterators access to topological object internals
  friend class VertexIterator;  friend class EdgeIterator;  
  friend class FaceIterator;    friend class TetIterator;
  
  friend class VertexVertexIterator; friend class VertexEdgeIterator; friend class VertexFaceIterator;
  friend class EdgeVertexIterator; friend class EdgeFaceIterator;
  friend class FaceVertexIterator; friend class FaceEdgeIterator; friend class FaceTetIterator;
  friend class TetFaceIterator; friend class TetVertexIterator;
  
  //allow properties to access object internals (for registering/unregistering themselves) 
  template<class T>
  friend class VertexProperty;
  
  template<class T>
  friend class EdgeProperty;

  template<class T>
  friend class FaceProperty;

  template<class T>
  friend class TetProperty;


   //** Simplex handles and iterators - these typedefs are somewhat redundant now
  typedef VertexHandle                      vertex_handle;
  typedef VertexIterator                    vertex_iter;
  typedef const VertexIterator              const_vertex_iter;

  typedef EdgeHandle                        edge_handle;
  typedef EdgeIterator                      edge_iter;
  typedef const EdgeIterator                const_edge_iter;

  typedef FaceHandle                        face_handle;
  typedef FaceIterator                      face_iter;
  typedef const FaceIterator                const_face_iter;

  typedef TetHandle                         tet_handle;
  typedef TetIterator                       tet_iter;
  typedef const TetIterator                 const_tet_iter;


  TopologicalObject();
  virtual ~TopologicalObject() {}

  /** Number of vertices */
  int nv() const;
  /** Number of edges */
  int ne() const;
  /** Number of faces */
  int nf() const;
  /** Number of tetrahedra */
  int nt() const;

  int getRelativeOrientation(const TetHandle& th, const FaceHandle& fh) const;
  int getRelativeOrientation(const FaceHandle& fh, const EdgeHandle& eh) const;
  int getRelativeOrientation(const EdgeHandle& eh, const VertexHandle& vh) const;

  /** Incidence counts for different simplices */

  unsigned int vertexIncidentEdges(const VertexHandle& v) const { return m_VE.getNumEntriesInRow(v.idx()); }
  unsigned int edgeIncidentFaces(const EdgeHandle& e) const { return m_EF.getNumEntriesInRow(e.idx()); }
  unsigned int faceIncidentTets(const FaceHandle& f) const { return m_FT.getNumEntriesInRow(f.idx()); }

  /** Add a vertex to the object */
  VertexHandle addVertex();
  /** Add an edge that connects the two given vertices, directed from v0 to v1 */
  EdgeHandle addEdge(const VertexHandle& v0, const VertexHandle& v1);
  /** Add a face using three edges, with edge directions determined by order of edge parameters*/
  FaceHandle addFace(const EdgeHandle& e0, const EdgeHandle& e1, const EdgeHandle& e2);
  /** Add a tet using four faces. Set face directions to be consistent with the first face, with last parameter dictating whether it is flipped.*/
  TetHandle addTet(const FaceHandle& f0, const FaceHandle& f1,
                    const FaceHandle& f2, const FaceHandle& f3,
                    bool flip_face0 = false);
  
  /** Add a triangular face to the object, using the given vertices. Note, this will be slower than the above methods.*/ 
  FaceHandle addFace(const VertexHandle& v0, const VertexHandle& v1, const VertexHandle& v2);

  /** Deletes the given vertex */
  bool deleteVertex(const VertexHandle& vertex);
  /** Deletes the given edge, and if recurse==true, any resulting isolated vertices composing it*/
  bool deleteEdge(const EdgeHandle& edge, bool recurse);
  /** Deletes the given face, and if recurse==true, any resulting isolated lower dimensional simplices composing it*/
  bool deleteFace(const FaceHandle& face, bool recurse);
  /** Deletes the given tet, and if recurse==true, any resulting isolated lower dimensional simplices composing it*/
  bool deleteTet(const tet_handle& tet, bool recurse);

  //** Return true if the vertex handle points to a valid vertex. */
  bool vertexExists(const VertexHandle& vertex) const;
  //** Return true if the edge handle points to a valid edge. */
  bool edgeExists(const EdgeHandle& edge) const;
  //** Return true if the face handle points to a face vertex. */
  bool faceExists(const FaceHandle& face) const;
  //** Return true if the tet handle points to a tet vertex. */
  bool tetExists(const tet_handle& tet) const;
  
  //*** Return true if the VH is incident to a boundary edge. */
  bool isBoundary(const VertexHandle & vh) const;
  //*** Return true if the EH is a boundary edge. */
  bool isBoundary(const EdgeHandle & eh) const;

  /** Returns a handle to the vertex at the tail of the edge */
  VertexHandle fromVertex(const EdgeHandle& eh) const;
  /** Returns a handle to the vertex at the head of the edge */
  VertexHandle toVertex(const EdgeHandle& eh) const;

  /** Given an edge of a face, get the next edge in order */
  EdgeHandle nextEdge(const FaceHandle& face, const EdgeHandle& curEdge) const;
  /** Given an edge of a face, get the prev edge in order */
  EdgeHandle prevEdge(const FaceHandle& face, const EdgeHandle& curEdge) const;

  /** Collapse an edge, deleting the given adjacent vertex */
  VertexHandle collapseEdge(const EdgeHandle& eh, const VertexHandle& vertToRemove, std::vector<EdgeHandle>& deletedEdges);

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

  /** Returns an iterator to the first face */
  face_iter faces_begin();
  /** Const version */
  const_face_iter faces_begin() const;

  /** Returns an iterator to the past-the-end face */
  face_iter faces_end();
  /** Const version */
  const_face_iter faces_end() const;

  /** Returns an iterator to the first tet */
  tet_iter tets_begin();
  /** Const version */
  const_tet_iter tets_begin() const;

  /** Returns an iterator to the past-the-end face */
  tet_iter tets_end();
  /** Const version */
  const_tet_iter tets_end() const;


  //@}

  /** \name Iterators over (co-)boundaries and adjacencies */

  //@{

  //VertexXXXIterators first...

  /** Iterator over vertices in the one-ring of the given vertex */
  VertexVertexIterator vv_iter(const VertexHandle& vh);
  /** Const version */
  const VertexVertexIterator vv_iter(const VertexHandle& vh) const;

  /** Iterator over edges adjacent to a vertex */
  VertexEdgeIterator ve_iter(const VertexHandle& vh);
  /** Const version */
  const VertexEdgeIterator ve_iter(const VertexHandle& vh) const;

  /** Iterator over faces adjacent to a vertex */
  VertexFaceIterator vf_iter(const VertexHandle& vh);
  /** Const version */
  const VertexFaceIterator vf_iter(const VertexHandle& vh) const;


  //EdgeXXXIterators next...

  /** Iterator over vertices adjacent to an edge */
  EdgeVertexIterator ev_iter(const EdgeHandle& eh);
  /** Const version */
  const EdgeVertexIterator ev_iter(const EdgeHandle& eh) const;

  /** Iterator over faces adjacent to a vertex */
  EdgeFaceIterator ef_iter(const EdgeHandle& eh);
  /** Const version */
  const EdgeFaceIterator ef_iter(const EdgeHandle& eh) const;


  //FaceXXXIterators next...
  
  /** Iterator over edges adjacent to a face */
  FaceTetIterator ft_iter(const FaceHandle& th);
  /** Const version */
  const FaceTetIterator ft_iter(const FaceHandle& th) const;


  /** Iterator over edges adjacent to a face */
  FaceEdgeIterator fe_iter(const FaceHandle& fh);
  /** Const version */
  const FaceEdgeIterator fe_iter(const FaceHandle& fh) const;

  /** Iterator over vertices adjacent to a face */
  FaceVertexIterator fv_iter(const FaceHandle& fh);
  /** Const version */
  const FaceVertexIterator fv_iter(const FaceHandle& fh) const;


  //TetXXXIterators next...
  
  /** Iterator over vertices adjacent to a face */
  TetFaceIterator tf_iter(const TetHandle& th);
  /** Const version */
  const TetFaceIterator tf_iter(const TetHandle& th) const;


  /** Iterator over vertices adjacent to a face */
  TetVertexIterator tv_iter(const TetHandle& th);
  /** Const version */
  const TetVertexIterator tv_iter(const TetHandle& th) const;


  //@}

  static void serializeStructure(std::ofstream& of, const TopologicalObject& obj);
  static void loadStructure(std::ifstream& ifs, TopologicalObject& obj);
  

  /** \name Vertex properties */
  //@{
  BA_CREATE_PROPERTY(VPropHandle, m_vertexProps, numVertexSlots());
  //@}

  /** \name Edge properties */
  //@{
  BA_CREATE_PROPERTY(EPropHandle, m_edgeProps, numEdgeSlots());
  //@}

  /** \name Face properties */
  //@{
  BA_CREATE_PROPERTY(FPropHandle, m_faceProps, numFaceSlots());
  //@}

  /** \name Tet properties */
  //@{
  BA_CREATE_PROPERTY(TPropHandle, m_tetProps, numTetSlots());
  //@}

  BA_INHERIT_BASE(ObjectBase);

protected:

   //Functions for registering/unregistering properties associated to simplex elements
  void registerVertexProperty(TopObjPropertyBase* prop) { m_vertPropsNew.push_back(prop); }
  void removeVertexProperty(TopObjPropertyBase* prop) { 
    std::vector<TopObjPropertyBase*>::iterator it = std::find(m_vertPropsNew.begin(), m_vertPropsNew.end(), prop);
    if(it != m_vertPropsNew.end()) 
      m_vertPropsNew.erase(it); 
  }
  void registerEdgeProperty(TopObjPropertyBase* prop) { m_edgePropsNew.push_back(prop); }
  void removeEdgeProperty(TopObjPropertyBase* prop) { 
    std::vector<TopObjPropertyBase*>::iterator it = std::find(m_edgePropsNew.begin(), m_edgePropsNew.end(), prop);
    if(it != m_edgePropsNew.end()) 
      m_edgePropsNew.erase(it); 
  }
  void registerFaceProperty(TopObjPropertyBase* prop) { m_facePropsNew.push_back(prop); }
  void removeFaceProperty(TopObjPropertyBase* prop) { 
    std::vector<TopObjPropertyBase*>::iterator it = std::find(m_facePropsNew.begin(), m_facePropsNew.end(), prop);
    if(it != m_facePropsNew.end()) 
      m_facePropsNew.erase(it); 
  }
  void registerTetProperty(TopObjPropertyBase* prop) { m_tetPropsNew.push_back(prop); }
  void removeTetProperty(TopObjPropertyBase* prop) { 
    std::vector<TopObjPropertyBase*>::iterator it = std::find(m_tetPropsNew.begin(), m_tetPropsNew.end(), prop);
    if(it != m_tetPropsNew.end()) 
      m_tetPropsNew.erase(it); 
  }

  //The number of spaces currently allocated for each simplex type. (Note this is different
  //from the number of active simplices of each type.)
  unsigned int numVertexSlots() const {return m_V.size();}
  unsigned int numEdgeSlots() const {return m_EV.getNumRows();}
  unsigned int numFaceSlots() const {return m_FE.getNumRows();}
  unsigned int numTetSlots() const {return m_TF.getNumRows();}

  //New lists of simplex properties. These are just pointers so that memory can be managed by the SimplexMesh for adding/deleting.
  //But the data lives wherever it has been created.
  std::vector<TopObjPropertyBase*> m_vertPropsNew;
  std::vector<TopObjPropertyBase*> m_edgePropsNew;
  std::vector<TopObjPropertyBase*> m_facePropsNew;
  std::vector<TopObjPropertyBase*> m_tetPropsNew;
  
  //Given two faces, return the shared edge if one exists
  EdgeHandle getSharedEdge(const FaceHandle& f0, const FaceHandle& f1) const;

  //Old simplex property containers
  PropertyContainer m_vertexProps; ///< Vertex property container
  PropertyContainer m_edgeProps;   ///< Edge property container
  PropertyContainer m_faceProps;   ///< Face property container
  PropertyContainer m_tetProps;   ///< Face property container

  int m_nv, m_ne, m_nf, m_nt;

  //New Topological Data Structure (Incidence graph, a la "Building Your Own DEC At Home", Elcott & Schroder 2005
  
  //Fundamental mesh data
  IncidenceMatrix m_TF;  ///< tet-to-face relations
  IncidenceMatrix m_FE;  ///< face-to-edge relations
  IncidenceMatrix m_EV;  ///< edge-to-vert relations
  std::vector<bool> m_V; ///< vertex existence, to support isolated vertices

  //Transposes, needed for efficient deletion/traversal
  IncidenceMatrix m_FT; ///< face-to-tet relations
  IncidenceMatrix m_EF; ///< edge-to-face relations
  IncidenceMatrix m_VE; ///< vert-to-edge relations

  //Pools of empty rows/columns in the above matrices, to efficiently add
  //data in previously deleted slots.
  std::deque<unsigned int> m_deadVerts, m_deadEdges, m_deadFaces, m_deadTets;
  
  //flags to indicate whether derived data is currently valid, or needs recomputing
  mutable bool m_validVF;
  mutable bool m_validTV, m_validVT;
  mutable bool m_validTE, m_validET;

  //Auxiliary, derived data, cached for convenient iteration. 
  //These are only computed if/when someone needs the associated iterator.
  mutable IncidenceMatrix m_nbrsVF;
  mutable IncidenceMatrix m_nbrsTV, m_nbrsVT;
  mutable IncidenceMatrix m_nbrsTE, m_nbrsET;

  void compute_nbrsVF() const;
  /*void compute_nbrsTV();
  void compute_nbrsVT();
  void compute_nbrsTE();
  void compute_nbrsET();*/

 
};

#include "TopologicalObject.inl"

} // namespace BASim

#endif // TOPOLOGICALOBJECT_HH
