#include "TopObjUtil.hh"
#include <set>

namespace BASim {

VertexHandle getSharedVertex(const TopologicalObject& obj, const EdgeHandle& e0, const EdgeHandle& e1) {
  VertexHandle v0 = obj.fromVertex(e0), v1 = obj.toVertex(e0),
    v2 = obj.fromVertex(e1), v3 = obj.toVertex(e1);

  if(v0 == v2 || v0 == v3) return v0;
  else if(v1 == v2 || v1 == v3) return v1;
  else return VertexHandle(-1);
}

VertexHandle getEdgesOtherVertex(const TopologicalObject& obj, const EdgeHandle &eh, const VertexHandle& vh) {
  VertexHandle v0 = obj.fromVertex(eh), v1 = obj.toVertex(eh);
  assert(v0 == vh || v1 == vh);
  return v0 == vh? v1 : v0;
}

FaceHandle getEdgeOtherFace(const TopologicalObject& obj, const EdgeHandle& eh, const FaceHandle& fh) {
  EdgeFaceIterator efit = obj.ef_iter(eh);
  for(; efit; ++efit) {
    FaceHandle cur = *efit;
    if(cur != fh)
      return cur;
  }
  return FaceHandle(-1);
}

bool getEdgeFacePair(const TopologicalObject& obj, const EdgeHandle& eh, FaceHandle& f0, FaceHandle& f1) {

  EdgeFaceIterator ef_it = obj.ef_iter(eh);
  if(!ef_it) return false; //if there is no face at all
  f0 = *ef_it;
  ++ef_it;
  if(!ef_it) return false; //there is no 2nd face, we also can't flip
  f1 = *ef_it;
  ++ef_it;

  assert(f0 != f1); //we should never get the same face
  assert(!ef_it); //this should probably not be used in non-manifold situations
  
  return true;
}

VertexHandle splitEdge(TopologicalObject& obj, const EdgeHandle& splitEdge, std::vector<FaceHandle>& newFaces) {
  EdgeFaceIterator ef_iter = obj.ef_iter(splitEdge);

  newFaces.clear();

  //get the edge's vertices
  VertexHandle from_vh, to_vh;
  from_vh = obj.fromVertex(splitEdge);
  to_vh = obj.toVertex(splitEdge);

  //add a new midpoint vertex
  VertexHandle newVert = obj.addVertex();

  //add the two new edges that are part of the original split edge
  EdgeHandle e_0 = obj.addEdge(from_vh, newVert);
  EdgeHandle e_1 = obj.addEdge(to_vh, newVert);

  //now iterate over the existing faces, splitting them in two appropriately
  std::vector<FaceHandle> facesToDelete;
  EdgeFaceIterator ef_it = obj.ef_iter(splitEdge);
  for(;ef_it; ++ef_it) {
    const FaceHandle& fh = *ef_it;

    //store this for deletion later.
    facesToDelete.push_back(fh);

    //find the other vertex in the first face
    FaceVertexIterator fv_it = obj.fv_iter(fh);
    while((*fv_it == from_vh) || (*fv_it == to_vh)) ++fv_it;
    VertexHandle other_vh = *fv_it;

    //create the new edge that splits this face
    EdgeHandle e_faceSplit = obj.addEdge(other_vh, newVert);

    //build the two new faces
    FaceEdgeIterator fe_it = obj.fe_iter(fh);
    for(;fe_it; ++fe_it) {
      EdgeHandle cur = *fe_it;

      //for each edge that isn't the splitEdge...      
      if(cur == splitEdge) continue;

      //accumulate a list of the new edges for the associated new face.
      //it is done in this manner to ensure new orientation is consistent with the old
      std::vector<EdgeHandle> edgeList;

      //determine which half of the splitEdge to use for this new face
      EdgeHandle halfEdge = obj.fromVertex(cur) == from_vh || obj.toVertex(cur) == from_vh? e_0 : e_1;

      FaceEdgeIterator fe_it2 = obj.fe_iter(fh);
      for(;fe_it2; ++fe_it2) {
        EdgeHandle cur2 = *fe_it2;
        if(cur2 == cur) edgeList.push_back(cur); //the current edge
        else if(cur2 == splitEdge) edgeList.push_back(halfEdge); //half of the original split edge
        else edgeList.push_back(e_faceSplit); //the new edge that splits the face
      }
      FaceHandle newFace = obj.addFace(edgeList[0], edgeList[1], edgeList[2]);
      newFaces.push_back(newFace);
    }

  }

  //Delete the previous faces and edge. Done as a post-process so as not to mess up the iterators.
  for(unsigned int i = 0; i < facesToDelete.size(); ++i)
    obj.deleteFace(facesToDelete[i], false);
  obj.deleteEdge(splitEdge, false);
  
  //Return a handle to the vertex we created
  return newVert;
}

EdgeHandle flipEdge(TopologicalObject& obj, const EdgeHandle& eh) {
  assert(obj.edgeExists(eh));

  VertexHandle from_vh, to_vh;
  from_vh = obj.fromVertex(eh);
  to_vh = obj.toVertex(eh);

  assert(from_vh != to_vh);

  //Just use the first two faces we hit.
  EdgeFaceIterator ef_it = obj.ef_iter(eh);
  if(!ef_it) return EdgeHandle(-1); //if there is no face at all, can't flip
  const FaceHandle& fh = *ef_it;
  ++ef_it;
  if(!ef_it) return EdgeHandle(-1); //there is no 2nd face, we also can't flip
  const FaceHandle& fh2 = *ef_it;
  ++ef_it;
  assert(fh != fh2); //we should never hit the same face
  assert(!ef_it); //this edge should have no more faces. A non-manifold edge flip doesn't make sense.

  //find the 3rd vertex in the first face
  FaceVertexIterator fv_it = obj.fv_iter(fh);
  while((*fv_it == from_vh) || (*fv_it == to_vh)) ++fv_it;
  VertexHandle f1_vh = *fv_it;

  //find the 3rd vertex in the second face
  FaceVertexIterator fv_it2 = obj.fv_iter(fh2);
  while((*fv_it2 == from_vh) || (*fv_it2 == to_vh)) ++fv_it2;
  VertexHandle f2_vh = *fv_it2;
  
  assert(f1_vh != f2_vh);
  
  assert(from_vh != f1_vh);
  assert(from_vh != f2_vh);
  assert(to_vh != f1_vh);
  assert(to_vh != f2_vh);

  //check for an edge already matching this description... and don't do the flip.
  EdgeHandle edge = findEdge(obj, f1_vh, f2_vh);
  if(edge.isValid()) return EdgeHandle(-1);

  EdgeHandle newEdge = obj.addEdge(f1_vh, f2_vh);

  //grab all the current edges in proper order, starting from the shared face
  EdgeHandle e0 = obj.nextEdge(fh, eh);
  EdgeHandle e1 = obj.nextEdge(fh, e0);

  EdgeHandle e2 = obj.nextEdge(fh2, eh);
  EdgeHandle e3 = obj.nextEdge(fh2, e2);

  assert(e0 != e1);
  assert(e0 != e2);
  assert(e0 != e3);
  assert(e0 != eh);
  assert(e0 != newEdge);
  
  assert(e1 != e2);
  assert(e1 != e3);
  assert(e1 != eh);
  assert(e1 != newEdge);

  assert(e2 != e3);
  assert(e2 != eh);
  assert(e2 != newEdge);

  assert(e3 != eh);
  assert(e3 != newEdge);

  assert(eh != newEdge);

  //flip the edges of the second face so it matches the first face 
  //(if the original faces didn't sync, it doesn't matter, since there's no way to flip and maintain consistency)
  VertexHandle sharedV = getSharedVertex(obj, e1, e2);
  if(!sharedV.isValid())
    std::swap(e2,e3);

  //add in the new faces
  obj.addFace(e1, e2, newEdge);
  obj.addFace(e3, e0, newEdge);

  //Delete the old patch
  bool success = obj.deleteFace(fh, false);
  assert(success);
  success = obj.deleteFace(fh2, false);
  assert(success);
  success = obj.deleteEdge(eh, false);
  assert(success);
  
  return newEdge;
}

//check if a vertex is on the boundary of a simplex mesh.
//We define this by the fact that one of its' edges doesn't have
//two faces attached.
bool isVertexOnBoundary(TopologicalObject& obj, VertexHandle& v) {
  VertexEdgeIterator ve_iter = obj.ve_iter(v);
  for(;ve_iter; ++ve_iter) {
    FaceHandle f0, f1;
    bool success = getEdgeFacePair(obj, *ve_iter, f0, f1);
    if(!success)
      return true;
  }
  return false;
}


bool getEdgeOppositeVertices(const TopologicalObject& obj, const EdgeHandle& eh, VertexHandle& v0, VertexHandle& v1) {
  FaceHandle f0,f1;
  bool success = getEdgeFacePair(obj, eh, f0, f1);
  
  if(!success) //there is only one other vertex
    return false;

  VertexHandle edgeV0 = obj.fromVertex(eh);
  VertexHandle edgeV1 = obj.toVertex(eh);
  
  for(FaceVertexIterator fv_it = obj.fv_iter(f0); fv_it; ++fv_it) {
    VertexHandle cur = *fv_it;
    if(cur != edgeV0 && cur != edgeV1)
      v0 = cur;
  }
  for(FaceVertexIterator fv_it = obj.fv_iter(f1); fv_it; ++fv_it) {
    VertexHandle cur = *fv_it;
    if(cur != edgeV0 && cur != edgeV1)
      v1 = cur;
  }
  return true;
}

bool getFaceThirdVertex(const TopologicalObject& obj, const FaceHandle& fh, const EdgeHandle&eh, VertexHandle& vertex) {
   //assumes edge eh if one of the edges of face fh
   for(FaceVertexIterator fvit = obj.fv_iter(fh); fvit; ++fvit) {
      if(*fvit  != obj.fromVertex(eh) && *fvit != obj.toVertex(eh)) {
         vertex = *fvit;
         return true;
      }
   }
   return false;
}

void sanityCheckTopology(TopologicalObject& obj) {
  
  for(EdgeIterator e_it = obj.edges_begin(); e_it != obj.edges_end(); ++e_it) {
    assert(obj.fromVertex(*e_it) != obj.toVertex(*e_it));
    assert(obj.edgeExists(*e_it));
    assert(obj.vertexExists(obj.fromVertex(*e_it)));
    assert(obj.vertexExists(obj.toVertex(*e_it)));

    FaceHandle f0, f1;
    bool success = getEdgeFacePair(obj, *e_it, f0, f1);
    if(success) {
      //count the number of unique edges
      std::set<EdgeHandle> edges;
      for(FaceEdgeIterator fe_it = obj.fe_iter(f0); fe_it; ++fe_it) edges.insert(*fe_it);
      for(FaceEdgeIterator fe_it = obj.fe_iter(f1); fe_it; ++fe_it) edges.insert(*fe_it);
      assert(edges.size() == 5);

      VertexHandle v0, v1;
      assert(f0 != f1);
      bool success2 = getEdgeOppositeVertices(obj, *e_it, v0, v1);
      if(success2) {
        assert(v0 != v1);
      }
    }

  }

  for(FaceIterator f_it = obj.faces_begin(); f_it != obj.faces_end(); ++f_it) {
    std::vector<VertexHandle> verts;
    for(FaceVertexIterator fv_it = obj.fv_iter(*f_it); fv_it; ++fv_it) {
      verts.push_back(*fv_it);

    }
    assert(verts.size() == 3);
    assert(verts[0] != verts[1] && verts[0] != verts[2]);
    
    std::vector<EdgeHandle> edges;
    for(FaceEdgeIterator fe_it = obj.fe_iter(*f_it); fe_it; ++fe_it) {
      edges.push_back(*fe_it);
    }
    assert(edges.size() == 3);
    assert(edges[0] != edges[1] && edges[0] != edges[2]);

    //check for other faces with the same edges
    for(FaceEdgeIterator fe_it = obj.fe_iter(*f_it); fe_it; ++fe_it) {
      EdgeHandle edge = *fe_it;
      FaceHandle f0, f1;
      bool success = getEdgeFacePair(obj, edge, f0, f1);
      if(!success) continue;
      assert(f0 == *f_it || f1 == *f_it);
      assert(f0 != f1);
      VertexHandle v0, v1;
      success = getEdgeOppositeVertices(obj, edge, v0, v1);
      if(!success) continue;
      assert(v0 != v1);
    }

  }

}

EdgeHandle findEdge( const TopologicalObject& obj, const VertexHandle& v0, const VertexHandle& v1 ) {
  for(VertexEdgeIterator ve_it = obj.ve_iter(v0); ve_it; ++ve_it) {
    EdgeHandle eh = *ve_it;
    if(obj.fromVertex(eh) == v1 || obj.toVertex(eh) == v1)
      return eh;
  }
  return EdgeHandle(-1);
}

}