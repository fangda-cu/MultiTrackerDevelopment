#include "BASim/src/Core/TopologicalObject/TopologicalObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjProperty.hh"
#include <utility>
#include <set>
#include "BASim/src/IO/SerializationUtils.hh"

namespace BASim {

/** Number of vertices */
int TopologicalObject::nv() const { return m_nv;}
/** Number of edges */
int TopologicalObject::ne() const { return m_ne;}
/** Number of faces */
int TopologicalObject::nf() const { return m_nf;}
/** Number of tetrahedra */
int TopologicalObject::nt() const { return m_nt;}

TopologicalObject::TopologicalObject()
{
  
  m_nv = 0;
  m_ne = 0;
  m_nf = 0;
  m_nt = 0;

  m_validTE = false;
  m_validET = false;
  m_validVF = false;
  m_validTV = false;
  m_validVT = false;
}

bool TopologicalObject::vertexExists(const VertexHandle& vertex) const {
  if(vertex.idx() < 0 || vertex.idx() >= (int)m_V.size())
    return false;
  else
    return m_V[vertex.idx()];
}

bool TopologicalObject::edgeExists(const EdgeHandle& edge) const {
  if(edge.idx() < 0 || edge.idx() >= (int)m_EV.getNumRows())
    return false;
  else
    return m_EV.getNumEntriesInRow(edge.idx()) > 0;

}

bool TopologicalObject::faceExists(const FaceHandle& face) const {
  if(face.idx() < 0 || face.idx() >= (int)m_FE.getNumRows())
    return false;
  else
    return m_FE.getNumEntriesInRow(face.idx()) > 0;
}

bool TopologicalObject::tetExists(const TetHandle& tet) const {
  if(tet.idx() < 0 || tet.idx() >= (int)m_TF.getNumRows())
    return false;
  else
    return m_TF.getNumEntriesInRow(tet.idx()) > 0;
}

VertexHandle TopologicalObject::addVertex()
{

  int new_index;
  if(m_deadVerts.size() == 0) {
    //create new vertex in incidence matrices
    m_EV.addCols(1);
    m_VE.addRows(1);

    //add the vertex to the data/property array
    m_V.push_back(true);
    m_vertexProps.resize(m_V.size());
    for(unsigned int i = 0; i < m_vertPropsNew.size(); ++i) m_vertPropsNew[i]->resize(m_V.size());

    new_index = m_V.size()-1;
  }
  else {
    new_index = m_deadVerts.back();
    m_deadVerts.pop_back();

    assert(!m_V[new_index]);

    m_V[new_index] = true;
  }

  m_nv += 1;

  m_validVF = false;

  return VertexHandle(new_index);
}


EdgeHandle
TopologicalObject::addEdge(const VertexHandle& v0, const VertexHandle& v1)
{
  assert(vertexExists(v0) && vertexExists(v1));
  assert(v0.idx() != v1.idx()); //prevent self-edges - I don't think they are useful.

  //check if an edge matching this description exists!
  unsigned int loop_end = m_VE.getNumEntriesInRow(v0.idx());
  for(unsigned int i = 0; i < loop_end; ++i) {
    int edgeID = m_VE.getColByIndex(v0.idx(), i);
    unsigned int loop_end2 = m_EV.getNumEntriesInRow(edgeID);
    for(unsigned int j = 0; j < loop_end2; ++j){
      int vertID = m_EV.getColByIndex(edgeID, j);
      assert(vertID != v1.idx()); //this creates a duplicate edge
    }
  }

  //get the next free edge index, or add space
  int new_index;
  if(m_deadEdges.size() == 0) {
    //make room for the new edge
    m_FE.addCols(1);
    m_EF.addRows(1);

    m_EV.addRows(1);
    m_VE.addCols(1);

    new_index = m_EV.getNumRows()-1;

    //add property slots for the new edge
    m_edgeProps.resize(m_EV.getNumRows());
    for(unsigned int i = 0; i < m_edgePropsNew.size(); ++i) m_edgePropsNew[i]->resize(m_EV.getNumRows());

    assert(m_EV.getNumRows() == m_VE.getNumCols());
    assert(m_EV.getNumCols() == m_VE.getNumRows());
    assert(m_FE.getNumRows() == m_EF.getNumCols());
    assert(m_FE.getNumCols() == m_EF.getNumRows());
  }
  else {
    //grab the first dead edge off the pile
    new_index = m_deadEdges.back();
    m_deadEdges.pop_back();
  }

  //build new edge connectivity
  //Choose ordering within the row explicitly
  m_EV.setByIndex(new_index, 0, v0.idx(), -1);
  m_EV.setByIndex(new_index, 1, v1.idx(), +1);

  m_VE.set(v0.idx(), new_index, -1);
  m_VE.set(v1.idx(), new_index, 1);

  //adjust edge count
  m_ne += 1;

  
  return EdgeHandle(new_index);
}

FaceHandle
TopologicalObject::addFace(const EdgeHandle& e0, 
                           const EdgeHandle& e1, 
                           const EdgeHandle& e2)
{
  assert(edgeExists(e0) && edgeExists(e1) && edgeExists(e2));
  
  //TODO Add asserts to check that all edges use the same 3 vertices, twice each.
   
  //check if a face matching this description exists!
  //grab one edge, and check all of its faces.
  //This is probably expensive/unnecessary to do all the time so might want to comment it.
  
  /* //Commented out for efficiency.
  for(unsigned int i = 0; i < m_EF.getNumEntriesInRow(e0.idx()); ++i) {
     int faceInd = m_EF.getColByIndex(e0.idx(), i);
     for(unsigned int j = 0; j < m_FE.getNumEntriesInRow(faceInd); ++j) {
        int edgeInd = m_FE.getColByIndex(faceInd, j);
        assert(edgeInd != e1.idx() && edgeInd != e2.idx());
     }
  }
  for(unsigned int i = 0; i < m_EF.getNumEntriesInRow(e1.idx()); ++i) {
     int faceInd = m_EF.getColByIndex(e1.idx(), i);
     for(unsigned int j = 0; j < m_FE.getNumEntriesInRow(faceInd); ++j) {
        int edgeInd = m_FE.getColByIndex(faceInd, j);
        assert(edgeInd != e0.idx() && edgeInd != e2.idx());
     }
  }
  for(unsigned int i = 0; i < m_EF.getNumEntriesInRow(e2.idx()); ++i) {
     int faceInd = m_EF.getColByIndex(e2.idx(), i);
     for(unsigned int j = 0; j < m_FE.getNumEntriesInRow(faceInd); ++j) {
        int edgeInd = m_FE.getColByIndex(faceInd, j);
        assert(edgeInd != e0.idx() && edgeInd != e1.idx());
     }
  }
  */

  //get the next free face or add one
  int new_index;
  if(m_deadFaces.size() == 0) {
    //create space for new face
    m_TF.addCols(1);
    m_FT.addRows(1);

    m_FE.addRows(1);
    m_EF.addCols(1);

    new_index = m_FE.getNumRows()-1;

    //allocate space for properties associated to the edge
    m_faceProps.resize(m_FE.getNumRows());
    for(unsigned int i = 0; i < m_facePropsNew.size(); ++i) m_facePropsNew[i]->resize(m_FE.getNumRows());

    assert(m_FE.getNumRows() == m_EF.getNumCols());
    assert(m_FE.getNumCols() == m_EF.getNumRows());
    assert(m_FT.getNumRows() == m_TF.getNumCols());
    assert(m_FT.getNumCols() == m_TF.getNumRows());
  }
  else {
    //grab the next empty face of the pile
    new_index = m_deadFaces.back();
    m_deadFaces.pop_back();
  }

  //Signs are chosen to follow the ordering of edges provided as input,
  //so we flip to ensure edge vertices connect properly.

  //If the head of the first edge doesn't match either of the second edge's vertices, we must flip it,
  //since we want it oriented towards the second edge.
  bool flip0 = (toVertex(e0) != fromVertex(e1) && toVertex(e0) != toVertex(e1));

  //Now determine the shared vertex between edges 0 and 1.
  //Then flip edge1 if the shared_vertex is at the head; it should be at the tail.
  VertexHandle shared_vert0_hnd = flip0? fromVertex(e0) : toVertex(e0);
  bool flip1 = (shared_vert0_hnd != fromVertex(e1));

  //Determine shared vertex between edge 1 and 2.
  //Then flip edge2 if the shared_vertex is at the head.
  VertexHandle shared_vert1_hnd = flip1 ? fromVertex(e1) : toVertex(e1);
  bool flip2 = (shared_vert1_hnd != fromVertex(e2));

  //build face connectivity
  //get a unique ordering by arbitrarily choosing smallest index to go first
  int smallest = std::min(std::min(e0.idx(), e1.idx()), e2.idx());

  //add'em in requested ordering
  m_FE.setByIndex(new_index, 0, e0.idx(), flip0?-1:1);
  m_FE.setByIndex(new_index, 1, e1.idx(), flip1?-1:1);
  m_FE.setByIndex(new_index, 2, e2.idx(), flip2?-1:1);

  //cycle to get the smallest one first, for consistency
  while(m_FE.getColByIndex((unsigned int)new_index, (unsigned int)0) != (unsigned int)smallest)
     m_FE.cycleRow(new_index);

  m_EF.set(e0.idx(), new_index, flip0?-1:1);
  m_EF.set(e1.idx(), new_index, flip1?-1:1);
  m_EF.set(e2.idx(), new_index, flip2?-1:1);

  m_nf += 1;
  
  //invalidate the relevant cached neighbour data
  m_validVF = false;
  
  return FaceHandle(new_index);
}

TetHandle
TopologicalObject::addTet(const FaceHandle& f0,
                          const FaceHandle& f1,
                          const FaceHandle& f2,
                          const FaceHandle& f3, bool flip_face0)
{

  assert(faceExists(f0) && faceExists(f1) && faceExists(f2) && faceExists(f3));

  //TODO Add asserts that ensure these faces are all touching and share the same 6 edges

  //get the next free tet or add one
  int new_index;
  if(m_deadTets.size() == 0) {
    //create new tet
    m_TF.addRows(1);
    m_FT.addCols(1);

    new_index = m_TF.getNumRows()-1;

    //allocate space for the properties
    m_tetProps.resize(m_TF.getNumRows());
    for(unsigned int i = 0; i < m_tetPropsNew.size(); ++i) m_tetPropsNew[i]->resize(m_TF.getNumRows());

    assert(m_FT.getNumRows() == m_TF.getNumCols());
    assert(m_FT.getNumCols() == m_TF.getNumRows());
  }
  else {
    //grab the next unused tet
    new_index = m_deadTets.back();
    m_deadTets.pop_back();
  }

  //Need to figure out signs to be consistent with the choice of the first face
  
  //Determine the shared edge between two adjacent faces. Flip the 2nd so its direction is opposed (signs differ) after the flips.
  EdgeHandle shared_edge0 = getSharedEdge(f0, f1);
  bool flip1 = flip_face0 ? m_FE.get(f0.idx(), shared_edge0.idx()) == m_FE.get(f1.idx(), shared_edge0.idx()) : 
                            m_FE.get(f0.idx(), shared_edge0.idx()) != m_FE.get(f1.idx(), shared_edge0.idx());
  
  EdgeHandle shared_edge1 = getSharedEdge(f0, f2);
  bool flip2 = flip_face0 ? m_FE.get(f0.idx(), shared_edge1.idx()) == m_FE.get(f2.idx(), shared_edge1.idx()) : 
                            m_FE.get(f0.idx(), shared_edge1.idx()) != m_FE.get(f2.idx(), shared_edge1.idx());

  EdgeHandle shared_edge2 = getSharedEdge(f0, f3);
  bool flip3 = flip_face0 ? m_FE.get(f0.idx(), shared_edge2.idx()) == m_FE.get(f3.idx(), shared_edge2.idx()) : 
                            m_FE.get(f0.idx(), shared_edge2.idx()) != m_FE.get(f3.idx(), shared_edge2.idx());

  //build tet connectivity
  m_TF.set(new_index, f0.idx(), flip_face0?1:-1);
  m_TF.set(new_index, f1.idx(), flip1?1:-1);
  m_TF.set(new_index, f2.idx(), flip2?1:-1);
  m_TF.set(new_index, f3.idx(), flip3?1:-1);

  m_FT.set(f0.idx(), new_index, flip_face0?1:-1);
  m_FT.set(f1.idx(), new_index, flip1?1:-1);
  m_FT.set(f2.idx(), new_index, flip2?1:-1);
  m_FT.set(f3.idx(), new_index, flip3?1:-1);

  m_nt += 1;

  //invalidate the relevant cached neighbour data
  m_validTV = false;
  m_validVT = false;
  m_validET = false;
  m_validTE = false;

  return tet_handle(new_index);
}

bool TopologicalObject::deleteVertex(const VertexHandle& vertex)
{
  assert(vertexExists(vertex));

  //test for safety - don't perform the delete if the simplex is not orphaned, or we get inconsistency.
  if(m_VE.getNumEntriesInRow(vertex.idx()) != 0)
    return false;

  //set the vertex to inactive
  m_V[vertex.idx()] = false;
  m_deadVerts.push_front(vertex.idx());

  //adjust the vertex count
  m_nv -= 1;
  
  m_validVF = false;

  return true;
}

bool TopologicalObject::deleteEdge(const EdgeHandle& edge, bool recurse)
{
  assert(edgeExists(edge));

  //test for safety - don't perform delete if not orphaned, or it creates inconsistency.
  if(m_EF.getNumEntriesInRow(edge.idx()) != 0)
    return false;

  //determine the corresponding vertices
  unsigned int loop_end =  m_EV.getNumEntriesInRow(edge.idx());
  for(unsigned int i = 0; i < loop_end; ++i) {

    //delete the edge entry in the transpose
    int col = m_EV.getColByIndex(edge.idx(),i);
    m_VE.zero(col, edge.idx()); 

    //delete the composing vertices if desired
    if(recurse)
      deleteVertex(VertexHandle(col));       
  }

  //...and delete the row
  m_EV.zeroRow(edge.idx());
  m_deadEdges.push_front(edge.idx());

  //adjust the edge count
  m_ne -= 1;

  return true;
}

bool TopologicalObject::deleteFace(const FaceHandle& face, bool recurse)
{
  assert(faceExists(face));

  //test for safety - don't perform delete if not orphaned, or we'll have inconsistencies.
  if(m_FT.getNumEntriesInRow(face.idx()) != 0)
    return false;

  //determine the corresponding edges
  unsigned int loop_end = m_FE.getNumEntriesInRow(face.idx());
  for(unsigned int i = 0; i < loop_end; ++i) {

    //remove face entry from the transpose
    int col = m_FE.getColByIndex(face.idx(), i);
    m_EF.zero(col, face.idx());

    //delete the composing edges
    if(recurse)
      deleteEdge(EdgeHandle(col), recurse);
  }

  //...and delete the row
  m_FE.zeroRow(face.idx());
  m_deadFaces.push_front(face.idx());

  //adjust the face count
  m_nf -= 1;

  //invalidate cached relationships
  m_validVF = false;

  return true;
}

bool TopologicalObject::deleteTet(const tet_handle& tet, bool recurse) {
  assert(tetExists(tet));

  //There are no higher dimensional simplices in 3D, so deleting the tet cannot introduce inconsistencies
  //as it can in the lower cases.

  //determine the corresponding faces
  unsigned int loop_end = m_TF.getNumEntriesInRow(tet.idx());
  for(unsigned int i = 0; i < loop_end; ++i) {

    //clear the tet entries in the transpose
    int col = m_TF.getColByIndex(tet.idx(), i);
    m_FT.zero(col, tet.idx());

    //delete composing faces if desired
    if(recurse)
      deleteFace(FaceHandle(col), recurse);
  }

  //...and delete the row
  m_TF.zeroRow(tet.idx());
  m_deadTets.push_front(tet.idx());

  //adjust the tet count
  m_nt -= 1;
  
  //invalidate cached relationships
  m_validTV = false;
  m_validVT = false;
  m_validTE = false;
  m_validET = false;

  return true;

}


//This alternate method for adding a face requires 
//searching for the edges containing the desired vertices, and creating them if they don't exist.
//Likely to be slower, albeit more convenient.
FaceHandle
TopologicalObject::addFace(const VertexHandle& v0, const VertexHandle& v1, const VertexHandle& v2)
{
  assert(vertexExists(v0) && vertexExists(v1) && vertexExists(v2));

  //Find the edges we need, or create them.

  EdgeHandle e01(-1), e02(-1);
  for(VertexEdgeIterator ve_iter(this, v0); ve_iter; ++ve_iter) {
    EdgeHandle e_hnd = *ve_iter;

    for(VertexEdgeIterator ve_iter2(this, v1); ve_iter2; ++ve_iter2)
      if(e_hnd == *ve_iter2)
        e01 = e_hnd;

    for(VertexEdgeIterator ve_iter3(this, v2); ve_iter3; ++ve_iter3)
      if(e_hnd == *ve_iter3)
        e02 = e_hnd;
  }

  EdgeHandle e12(-1);
  for(VertexEdgeIterator ve_iter(this, v1); ve_iter; ++ve_iter) {
    EdgeHandle e_hnd = *ve_iter;

    for(VertexEdgeIterator ve_iter2(this, v2); ve_iter2; ++ve_iter2)
      if(e_hnd == *ve_iter2)
        e12 = e_hnd;
  }


  //careful about making the ordering match the desired input ordering of the vertices
  if(!e01.isValid())
    e01 = addEdge(v0,v1);
  if(!e12.isValid())
    e12 = addEdge(v1,v2);
  if(!e02.isValid())
    e02 = addEdge(v2,v0);

  return addFace(e01, e12, e02);
}


EdgeHandle TopologicalObject::getSharedEdge(const FaceHandle& f0, const FaceHandle& f1) const {
  assert(faceExists(f0) && faceExists(f1));
  assert(m_FE.getNumEntriesInRow(f0.idx()) == 3);
  assert(m_FE.getNumEntriesInRow(f1.idx()) == 3);

  //Iterate over all pairs, stop when we hit the match.
  for(unsigned int ind0 = 0; ind0 < 3; ++ind0) for(unsigned ind1 = 0; ind1 < 3; ++ind1) {
    int col0 = m_FE.getColByIndex(f0.idx(), ind0);
    int col1 = m_FE.getColByIndex(f1.idx(), ind1);
    if(col0 == col1)
      return EdgeHandle(col0);
  }
  
  return EdgeHandle(-1);
}

VertexHandle
TopologicalObject::fromVertex(const EdgeHandle& eh) const
{
  assert(m_EV.getNumEntriesInRow(eh.idx()) == 2);

  return VertexHandle(m_EV.getColByIndex(eh.idx(), 0));
  //return m_EV.getValueByIndex(eh.idx(), 0) == -1? 
  //  VertexHandle(m_EV.getColByIndex(eh.idx(), 0)) : VertexHandle(m_EV.getColByIndex(eh.idx(), 1));

}

VertexHandle
TopologicalObject::toVertex(const EdgeHandle& eh) const
{
   assert(m_EV.getNumEntriesInRow(eh.idx()) == 2);

   return VertexHandle(m_EV.getColByIndex(eh.idx(), 1));
   //return m_EV.getValueByIndex(eh.idx(), 0) == 1? 
   //   VertexHandle(m_EV.getColByIndex(eh.idx(), 0)) : VertexHandle(m_EV.getColByIndex(eh.idx(), 1));

}

int TopologicalObject::getRelativeOrientation(const TetHandle& th, const FaceHandle& fh) const {
  if(th.idx() < 0 || th.idx() >= (int)m_TF.getNumRows() || fh.idx() < 0 || fh.idx() >= (int)m_TF.getNumCols())
     return 0;
  return m_TF.get(th.idx(), fh.idx());
}

int TopologicalObject::getRelativeOrientation(const FaceHandle& fh, const EdgeHandle& eh) const {
  if(fh.idx() < 0 || fh.idx() >= (int)m_FE.getNumRows() || eh.idx() < 0 || eh.idx() >= (int)m_FE.getNumCols())
     return 0;
  return m_FE.get(fh.idx(), eh.idx());
}

int TopologicalObject::getRelativeOrientation(const EdgeHandle& eh, const VertexHandle& vh) const {
  if(eh.idx() < 0 || eh.idx() >= (int)m_EV.getNumRows() || vh.idx() < 0 || vh.idx() >= (int)m_EV.getNumCols())
     return 0;
  return m_EV.get(eh.idx(), vh.idx());
}
bool TopologicalObject::isBoundary(const VertexHandle & vh) const {
    for ( VertexEdgeIterator veit = ve_iter(vh); veit; ++veit){
        if( isBoundary(*veit) ){ return true; }
    }
    return false;
}
bool TopologicalObject::isBoundary(const EdgeHandle & eh) const{
    //TODO: extend this for non-manifold meshes
    return (edgeIncidentFaces(eh) == 1);
}



EdgeHandle TopologicalObject::nextEdge(const FaceHandle& face, const EdgeHandle& curEdge) const {
   //take advantage of known fixed ordering of indices
   int faceIdx = face.idx();
   int colIdx = curEdge.idx();
   int col0 = m_FE.getColByIndex(faceIdx, 0);
   int col1 = m_FE.getColByIndex(faceIdx, 1);
   if(col0 == colIdx)
      return EdgeHandle(col1);
   else if(col1 == colIdx)
      return EdgeHandle(m_FE.getColByIndex(faceIdx, 2));
   else 
      return EdgeHandle(col0);
}

EdgeHandle TopologicalObject::prevEdge(const FaceHandle& face, const EdgeHandle& curEdge) const {
   int faceIdx = face.idx();
   int colIdx = curEdge.idx();
   int col0 = m_FE.getColByIndex(faceIdx, 0);
   int col2 = m_FE.getColByIndex(faceIdx, 2);
   if(col0 == colIdx)
      return EdgeHandle(col2);
   else if(col2 == colIdx)
      return EdgeHandle(m_FE.getColByIndex(faceIdx, 1));
   else 
      return EdgeHandle(col0);
}

void TopologicalObject::compute_nbrsVF() const {
  if(m_nbrsVF.getNumRows() < numVertexSlots())
    m_nbrsVF.addRows(numVertexSlots() - m_nbrsVF.getNumRows());
  if(m_nbrsVF.getNumCols() < numFaceSlots())
    m_nbrsVF.addCols(numFaceSlots() - m_nbrsVF.getNumCols());
 
  m_nbrsVF.zeroAll();
  
  //consider each vertex
  for(unsigned int vert_idx = 0; vert_idx < m_VE.getNumRows(); ++vert_idx) {

    //consider each edge of the face
     unsigned int loop_end = m_VE.getNumEntriesInRow(vert_idx);
    for(unsigned int e = 0; e < loop_end; ++e) {
      unsigned int edge_idx = m_VE.getColByIndex(vert_idx, e);

      //consider each face of the edge
      unsigned int loop_end2 = m_EF.getNumEntriesInRow(edge_idx);
      for(unsigned int f = 0; f < loop_end2; ++f) {
        unsigned int face_idx = m_EF.getColByIndex(edge_idx, f);

        //set the vertex to true. (This will be visited ~2x for each face, but that's okay.)
        m_nbrsVF.set(vert_idx, face_idx, 1);
      }
    }
  }

  m_validVF = true;

}

VertexHandle TopologicalObject::collapseEdge(const EdgeHandle& eh, const VertexHandle& vertToRemove, std::vector<EdgeHandle>& deletedEdges) {
  
  VertexHandle fromV = fromVertex(eh);
  VertexHandle toV = toVertex(eh);
  assert(fromV == vertToRemove || toV == vertToRemove); //make sure the vertex selected for deletion is actually used by the edge
  
  //grab the vertex that we are keeping
  int vertToKeep = fromV == vertToRemove? toV.idx() : fromV.idx();

  //determine adjacent faces to the collapsing edge
  int faceCount = m_EF.getNumEntriesInRow(eh.idx());
  
  //look at faces on left and right of collapsing edge
  std::vector<int> facesToDelete(faceCount);
  for(int f = 0; f < faceCount; ++f) {
    unsigned int face_idx = m_EF.getColByIndex(eh.idx(), f);
    facesToDelete[f] = face_idx;
  }
  
  //delete the faces and then edge, leaving a hole to be stitched
  for(unsigned int i = 0; i < facesToDelete.size(); ++i) {
    bool success = deleteFace(FaceHandle(facesToDelete[i]), false);
    assert(success);
  }
  bool success = deleteEdge(eh, false);
  assert(success);
  
  
  //determine all existing edges using the vertex being eliminated
  std::vector< std::pair<int,int> > edgeIndices;
  unsigned int loop_end = m_VE.getNumEntriesInRow(vertToRemove.idx());
  for(unsigned int e = 0; e < loop_end; ++e) {
    unsigned int edgeInd = m_VE.getColByIndex(vertToRemove.idx(), e);
    int sign = m_VE.getValueByIndex(vertToRemove.idx(), e);
    edgeIndices.push_back(std::make_pair(edgeInd,sign));
  }
  
  //relabel all the edges' to-be-deleted endpoints to the vertex being kept.
  //doing it "in place" like this rather than using safer atomic add/deletes
  //ensures that the original data on the modified edges gets maintained.
  for(unsigned int i = 0; i < edgeIndices.size(); ++i) {
    unsigned int edgeInd = edgeIndices[i].first;
    int vertSign = edgeIndices[i].second;
    
    m_VE.zero(vertToRemove.idx(), edgeInd);
    m_VE.set(vertToKeep, edgeInd, vertSign);

    //replace the data in the correct slot (0 == from, 1 == to).
    m_EV.setByIndex(edgeInd, vertSign==-1? 0 : 1, vertToKeep, vertSign);
    
  }

  
  //now we have some edges that are duplicates, possibly pointing in opposite directions
  //identify duplicate edges for deletion.
  std::vector< std::pair<int,int> > duplicateEdges; //list of (edge,edge) pairs that are duplicates
  std::map<int,int> vertEdgeMap;//for each vertex in the set of nbrs, the first edge we hit that uses it.
  unsigned int loop_end2 = m_VE.getNumEntriesInRow(vertToKeep);
  for(unsigned int e = 0; e < loop_end2; ++e) {
    int edgeInd = m_VE.getColByIndex(vertToKeep,e);
    int fromV = fromVertex(EdgeHandle(edgeInd)).idx();
    int toV = toVertex(EdgeHandle(edgeInd)).idx();
    int otherVert = fromV == vertToKeep? toV : fromV;
    
    //check if an edge with this end vertex has been seen yet
    //if not, add to the set; if so, log it as a duplicate
    std::map<int,int>::iterator iter = vertEdgeMap.find(otherVert);
    if(iter != vertEdgeMap.end())
      duplicateEdges.push_back(std::make_pair(edgeInd,(*iter).second));
    else
      vertEdgeMap[otherVert] = edgeInd;
  }

  
  //now replace the duplicate edges with their partner everywhere they're used (relabelling again)
  for(unsigned int i = 0; i < duplicateEdges.size(); ++i) {
    EdgeHandle e0(duplicateEdges[i].first);
    EdgeHandle e1(duplicateEdges[i].second);
    
    //if the edges pointed in opposite directions, their directions in faces need to be swapped during the relabelling
    int flipSign = fromVertex(e0) != fromVertex(e1) ? -1 : 1;
    
    //let's choose to remove e1 arbitrarily.

    //collect all faces that use this edge
    std::vector< std::pair<int,int> > faceIndices;
    unsigned int loop_end3 = m_EF.getNumEntriesInRow(e1.idx());
    for(unsigned int f = 0; f < loop_end3; ++f) {
      unsigned int faceInd = m_EF.getColByIndex(e1.idx(), f);
      int sign = m_EF.getValueByIndex(e1.idx(), f);
      faceIndices.push_back(std::make_pair(faceInd,sign));
    }

    //relabel all the faces' to-be-deleted edge to the other duplicate edge being kept.
    //doing it "in place" like this rather than using safer atomic add/deletes
    //ensures that the original data on the retained edges gets maintained.
    for(unsigned int i = 0; i < faceIndices.size(); ++i) {
      unsigned int faceInd = faceIndices[i].first;
      int edgeSign = faceIndices[i].second;

      int newSign = flipSign*edgeSign;
      m_EF.zero(e1.idx(), faceInd);
      m_EF.set(e0.idx(), faceInd, newSign);
      
      //determine row index of this one
      int ind = m_FE.getIndexByCol(faceInd, e1.idx());
      m_FE.setByIndex(faceInd, ind, e0.idx(), newSign);

      //cycle the row to make sure the smallest column index is still stored first 
      //(this helps ensure the face is uniquely defined)
      while(m_FE.getColByIndex(faceInd, 0) > m_FE.getColByIndex(faceInd, 1) || m_FE.getColByIndex(faceInd, 0) > m_FE.getColByIndex(faceInd, 2))
         m_FE.cycleRow(faceInd);
    
    }

    //finally, delete the orphaned edge
    deletedEdges.push_back(e1);
    success = deleteEdge(e1, false);
    assert(success);
  }
  
  //invalidate the relevant cached neighbour data
  m_validVF = false;

  success = deleteVertex(vertToRemove);
  assert(success);
  

  return VertexHandle(vertToKeep);
}

void TopologicalObject::serializeStructure(std::ofstream& of, const TopologicalObject& obj) {
   assert(of.is_open());

   //core data
   IncidenceMatrix::serialize(of, obj.m_VE);
   IncidenceMatrix::serialize(of, obj.m_EF);
   IncidenceMatrix::serialize(of, obj.m_FT);
   serializeVectorBool(of, obj.m_V);

   //this should probably be reconstructed from the above instead of saved
   IncidenceMatrix::serialize(of, obj.m_EV);
   IncidenceMatrix::serialize(of, obj.m_FE);
   IncidenceMatrix::serialize(of, obj.m_TF);
   
  std::vector<unsigned int> dead;
  dead.assign(obj.m_deadVerts.begin(), obj.m_deadVerts.end());
  serializeVectorUint(of, dead);
  dead.assign(obj.m_deadEdges.begin(), obj.m_deadEdges.end());
  serializeVectorUint(of, dead);
  dead.assign(obj.m_deadFaces.begin(), obj.m_deadFaces.end());
  serializeVectorUint(of, dead);
  dead.assign(obj.m_deadTets.begin(), obj.m_deadTets.end());
  serializeVectorUint(of, dead);
   
}

void TopologicalObject::loadStructure(std::ifstream& ifs, TopologicalObject& obj) {
   assert(ifs.is_open());

   //core data
   IncidenceMatrix::load(ifs, obj.m_VE);
   IncidenceMatrix::load(ifs, obj.m_EF);
   IncidenceMatrix::load(ifs, obj.m_FT);
   loadVectorBool(ifs, obj.m_V);

   //this should probably be reconstructed from the above instead of saved
   IncidenceMatrix::load(ifs, obj.m_EV);
   IncidenceMatrix::load(ifs, obj.m_FE);
   IncidenceMatrix::load(ifs, obj.m_TF);

  std::vector<unsigned int> dead;
  loadVectorUint(ifs, dead);
  obj.m_deadVerts.assign(dead.begin(), dead.end());
  loadVectorUint(ifs, dead);
  obj.m_deadEdges.assign(dead.begin(), dead.end());
  loadVectorUint(ifs, dead);
  obj.m_deadFaces.assign(dead.begin(), dead.end());
  loadVectorUint(ifs, dead);
  obj.m_deadTets.assign(dead.begin(), dead.end());

}

} //namespace BASim
