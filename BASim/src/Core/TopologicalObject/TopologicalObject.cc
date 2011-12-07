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
  
  return VertexHandle(new_index);
}


EdgeHandle
TopologicalObject::addEdge(const VertexHandle& v0, const VertexHandle& v1)
{
  assert(vertexExists(v0) && vertexExists(v1));
  assert(v0.idx() != v1.idx()); //prevent self-edges - I don't think they are useful.

  //check if an edge matching this description exists!
  for(unsigned int i = 0; i < m_VE.getNumEntriesInRow(v0.idx()); ++i) {
    int edgeID = m_VE.getColByIndex(v0.idx(), i);
    for(unsigned int j = 0; j < m_EV.getNumEntriesInRow(edgeID); ++j){
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
  m_EV.set(new_index, v0.idx(), -1);
  m_EV.set(new_index, v1.idx(), 1);

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
  m_FE.set(new_index, e0.idx(), flip0?-1:1);
  m_FE.set(new_index, e1.idx(), flip1?-1:1);
  m_FE.set(new_index, e2.idx(), flip2?-1:1);

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
                            m_FE.get(f0.idx(), shared_edge2.idx()) != m_FE.get(f3.idx(), shared_edge1.idx());

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
  m_deadVerts.push_back(vertex.idx());

  //adjust the vertex count
  m_nv -= 1;

  return true;
}

bool TopologicalObject::deleteEdge(const EdgeHandle& edge, bool recurse)
{
  assert(edgeExists(edge));

  //test for safety - don't perform delete if not orphaned, or it creates inconsistency.
  if(m_EF.getNumEntriesInRow(edge.idx()) != 0)
    return false;

  //determine the corresponding vertices
  for(unsigned int i = 0; i < m_EV.getNumEntriesInRow(edge.idx()); ++i) {

    //delete the edge entry in the transpose
    int col = m_EV.getColByIndex(edge.idx(),i);
    m_VE.zero(col, edge.idx()); 

    //delete the composing vertices if desired
    if(recurse)
      deleteVertex(VertexHandle(col));       
  }

  //...and delete the row
  m_EV.zeroRow(edge.idx());
  m_deadEdges.push_back(edge.idx());

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
  for(unsigned int i = 0; i < m_FE.getNumEntriesInRow(face.idx()); ++i) {

    //remove face entry from the transpose
    int col = m_FE.getColByIndex(face.idx(), i);
    m_EF.zero(col, face.idx());

    //delete the composing edges
    if(recurse)
      deleteEdge(EdgeHandle(col), recurse);
  }

  //...and delete the row
  m_FE.zeroRow(face.idx());
  m_deadFaces.push_back(face.idx());

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
  for(unsigned int i = 0; i < m_TF.getNumEntriesInRow(tet.idx()); ++i) {

    //clear the tet entries in the transpose
    int col = m_TF.getColByIndex(tet.idx(), i);
    m_FT.zero(col, tet.idx());

    //delete composing faces if desired
    if(recurse)
      deleteFace(FaceHandle(col), recurse);
  }

  //...and delete the row
  m_TF.zeroRow(tet.idx());
  m_deadTets.push_back(tet.idx());

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

  return m_EV.getValueByIndex(eh.idx(), 0) == -1? 
    VertexHandle(m_EV.getColByIndex(eh.idx(), 0)) : VertexHandle(m_EV.getColByIndex(eh.idx(), 1));

}

VertexHandle
TopologicalObject::toVertex(const EdgeHandle& eh) const
{
  assert(m_EV.getNumEntriesInRow(eh.idx()) == 2);

  return m_EV.getValueByIndex(eh.idx(), 0) == 1? 
    VertexHandle(m_EV.getColByIndex(eh.idx(), 0)) : VertexHandle(m_EV.getColByIndex(eh.idx(), 1));

}

EdgeHandle TopologicalObject::nextEdge(const FaceHandle& face, const EdgeHandle& curEdge, int direction) const {
  assert(faceExists(face));
  assert(edgeExists(curEdge));
  assert(m_FE.getNumEntriesInRow(face.idx()) == 3);
  assert(direction == 1 || direction == -1);

  //Determine the edge sign (and hence direction wrt the face).
  int curSign = m_FE.get(face.idx(), curEdge.idx());
  
  assert(curSign != 0);

  //Grab the end vertex, from the face's perspective.
  VertexHandle sharedVertex = direction*curSign == 1 ? toVertex(curEdge) : fromVertex(curEdge);

  //Check the face's edges for the one that also shares this vertex, as its initial vertex, wrt the face ordering.
  for(int i = 0; i < 3; ++i) {
    EdgeHandle e_hnd(m_FE.getColByIndex(face.idx(), i));
    if(e_hnd == curEdge) continue;
    int sign = m_FE.getValueByIndex(face.idx(), i);

    //If the 'from' vertex matches the target shared vertex, return this edge.
    VertexHandle fromVert = direction*sign == 1? fromVertex(e_hnd) : toVertex(e_hnd);
    if(fromVert == sharedVertex) {
      return e_hnd;
    }
  }

  assert(false); //failed to find an edge that shares the face

  return EdgeHandle(-1);
}


EdgeHandle TopologicalObject::nextEdge(const FaceHandle& face, const EdgeHandle& curEdge) const {
  return nextEdge(face, curEdge, 1);
}

EdgeHandle TopologicalObject::prevEdge(const FaceHandle& face, const EdgeHandle& curEdge) const {
  return nextEdge(face, curEdge, -1);
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
    for(unsigned int e = 0; e < m_VE.getNumEntriesInRow(vert_idx); ++e) {
      unsigned int edge_idx = m_VE.getColByIndex(vert_idx, e);

      //consider each face of the edge
      for(unsigned int f = 0; f < m_EF.getNumEntriesInRow(edge_idx); ++f) {
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
  
    //visit the face's other edges
    std::set<int> neighbourEdges;
    for(unsigned int i = 0; i < m_FE.getNumEntriesInRow(face_idx); ++i) {
      int edge_idx = m_FE.getColByIndex(face_idx, i);
      if(edge_idx != eh.idx()) {

        //walk over the other faces that this edge belongs to, checking for a shared edge
        for(unsigned int j = 0; j < m_EF.getNumEntriesInRow(edge_idx); ++j) {
          int curFace = m_EF.getColByIndex(edge_idx, j);
          
          //if we see a shared edge, then the collapsing edge merges two faces
          //and that's unacceptable.
          if(curFace != face_idx) {
            for(unsigned int k = 0; k < m_FE.getNumEntriesInRow(curFace); ++k) {
              int curEdge = m_FE.getColByIndex(curFace, k);
              
              if(neighbourEdges.find(curEdge) != neighbourEdges.end()) {
                return VertexHandle(-1); //there's a shared face/edge, don't collapse
              }
              else {
                neighbourEdges.insert(curEdge);
              }
            }
          }
        }
      }
    }

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
  for(unsigned int e = 0; e < m_VE.getNumEntriesInRow(vertToRemove.idx()); ++e) {
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
    m_EV.zero(edgeInd, vertToRemove.idx());
    
    m_VE.set(vertToKeep, edgeInd, vertSign);
    m_EV.set(edgeInd, vertToKeep, vertSign);
  }

  
  //now we have some edges that are duplicates, possibly pointing in opposite directions
  //identify duplicate edges for deletion.
  std::vector< std::pair<int,int> > duplicateEdges; //list of (edge,edge) pairs that are duplicates
  std::map<int,int> vertEdgeMap;//for each vertex in the set of nbrs, the first edge we hit that uses it.
  for(unsigned int e = 0; e < m_VE.getNumEntriesInRow(vertToKeep); ++e) {
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
    for(unsigned int f = 0; f < m_EF.getNumEntriesInRow(e1.idx()); ++f) {
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
      m_FE.zero(faceInd, e1.idx());

      m_EF.set(e0.idx(), faceInd, newSign);
      m_FE.set(faceInd, e0.idx(), newSign);
    }

    //finally, delete the orphaned edge
    deletedEdges.push_back(e1);
    bool success = deleteEdge(e1, false);
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
   
   serializeVectorUint(of, obj.m_deadVerts);
   serializeVectorUint(of, obj.m_deadEdges);
   serializeVectorUint(of, obj.m_deadFaces);
   serializeVectorUint(of, obj.m_deadTets);
   
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

   loadVectorUint(ifs, obj.m_deadVerts);
   loadVectorUint(ifs, obj.m_deadEdges);
   loadVectorUint(ifs, obj.m_deadFaces);
   loadVectorUint(ifs, obj.m_deadTets);

}

} //namespace BASim
