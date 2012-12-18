#ifndef TOPOBJUTIl_H
#define TOPOBJUTIl_H

//Utility operations on TopologicalObjects, eg. messing with connectivity.
//These use the iterators and add/delete[Simplex] functions, rather than messing with internals.
//This is perhaps less optimal, but likely far simpler than writing the special case code.
#include "BASim/src/Core/TopologicalObject/TopObjHandles.hh"
#include "BASim/src/Core/TopologicalObject/TopologicalObject.hh"

namespace BASim {

//Given an edge, return the two wings (adjacent faces) of that edge.
bool getEdgeFacePair(const TopologicalObject& obj, const EdgeHandle& eh, FaceHandle& f0, FaceHandle& f1);

//Given two edges, find a shared vertex if one exists.
VertexHandle getSharedVertex(const TopologicalObject& obj, const EdgeHandle& e0, const EdgeHandle& e1);

//Given two faces, find a shared edge if one exists.
EdgeHandle getSharedEdge(const TopologicalObject& obj, const FaceHandle& f0, const FaceHandle& f1);

//Given an edge and a vertex which is on the edge, return the edge's other vertex
VertexHandle getEdgesOtherVertex(const TopologicalObject& obj, const EdgeHandle &eh, const VertexHandle& vh);

//Split an edge and insert a new vertex in between, subdividing all the faces sharing the edge.
//Should be non-manifold friendly. Orientation is maintained.
VertexHandle splitEdge(TopologicalObject& obj, const EdgeHandle& h, std::vector<FaceHandle>& newFaces);

//Takes an edge with two adjacent faces comprising a quad, and replaces the edge with the other diagonal of the quad.
//Handle non-manifold scenarios.
EdgeHandle flipEdge(TopologicalObject& obj, const EdgeHandle& eh, 
                    const FaceHandle& f1, const FaceHandle& f2, 
                    FaceHandle& f1New, FaceHandle& f2New);

bool getEdgeOppositeVertices(const TopologicalObject& obj, const EdgeHandle& eh, VertexHandle& v0, VertexHandle& v1);

bool getFaceThirdVertex(const TopologicalObject& obj, const FaceHandle& fh, const EdgeHandle&eh, VertexHandle& vertex);
bool getTetFourthVertex(const TopologicalObject& obj, const TetHandle & th, const FaceHandle&fh, VertexHandle& vertex);
  
FaceHandle getVertexOppositeFaceInTet(const TopologicalObject & obj, const TetHandle & th, const VertexHandle & vh);
EdgeHandle getVertexOppositeEdgeInFace(const TopologicalObject & obj, const FaceHandle & fh, const VertexHandle & vh);
  
bool isVertexOnBoundary(TopologicalObject& obj, VertexHandle& v);

void sanityCheckTopology(TopologicalObject& obj);

FaceHandle getEdgeOtherFace(const TopologicalObject& obj, const EdgeHandle& eh, const FaceHandle& fh);
TetHandle  getFaceOtherTet (const TopologicalObject& obj, const FaceHandle& fh, const TetHandle& th);

EdgeHandle findEdge( const TopologicalObject& obj, const VertexHandle& v0, const VertexHandle& v1 );
//bool faceExists(const EdgeHandle& e0, const EdgeHandle& e1, const EdgeHandle& e2)
//bool faceExists(const VertexHandle& v0, const VertexHandle& v1, const VertexHandle& v2)
void tearInteriorEdge(TopologicalObject& obj,const EdgeHandle& e, const VertexHandle &va, const VertexHandle & vb,
        std::vector<VertexHandle> &newVerts,
        std::vector<FaceHandle> &newFaces,
        std::vector<FaceHandle> &facesToDelete,
        std::vector<EdgeHandle> &edgesToDelete);
void tearEdge(TopologicalObject& obj,const EdgeHandle& e, const VertexHandle &va, const VertexHandle & vb,
        VertexHandle & newVerta, VertexHandle & newVertb, std::vector<FaceHandle> &newFaces,
        std::vector<FaceHandle> &facesToDelete, std::vector<EdgeHandle> &edgesToDelete);
void tearVertexAlong(TopologicalObject& obj,const EdgeHandle& e, const VertexHandle &va,
        VertexHandle & newVert, std::vector<FaceHandle> &newFaces,
        std::vector<FaceHandle> &facesToDelete, std::vector<EdgeHandle> &edgesToDelete);
void addPrevSide(TopologicalObject & obj, const FaceHandle &f, const EdgeHandle &e,
        const VertexHandle& pivot, const VertexHandle &newVert, std::vector<EdgeHandle> & oldEdges,
        std::vector<FaceHandle> & oldFaces, std::vector<FaceHandle> &newFaces);
void addNextSide(TopologicalObject & obj, const FaceHandle &f, const EdgeHandle &e,
        const VertexHandle& pivot, const VertexHandle &newVert, std::vector<EdgeHandle> & oldEdges,
        std::vector<FaceHandle> & oldFaces, std::vector<FaceHandle> &newFaces);

bool isFaceMatch(const TopologicalObject& obj, const FaceHandle& fh, const VertexHandle& v0, const VertexHandle& v1, const VertexHandle& v2);
bool faceContainsVertex(const TopologicalObject & obj, const FaceHandle & fh, const VertexHandle & vh);
bool faceContainsEdge  (const TopologicalObject & obj, const FaceHandle & fh, const EdgeHandle & eh);
bool tetContainsVertex (const TopologicalObject & obj, const TetHandle & th,  const VertexHandle & vh);

FaceHandle findFace( const TopologicalObject& obj, const VertexHandle& v0, const VertexHandle& v1, const VertexHandle& v2 );
FaceHandle findFace( const TopologicalObject& obj, const EdgeHandle& e0,   const EdgeHandle& e1,   const EdgeHandle& e2 );
TetHandle  findTet ( const TopologicalObject& obj, const VertexHandle& v0, const VertexHandle& v1, const VertexHandle& v2, const VertexHandle& v3 );
TetHandle  findTet ( const TopologicalObject& obj, const FaceHandle& f0,   const FaceHandle& f1,   const FaceHandle& f2,   const FaceHandle& f3 );

bool isVertexManifold(const TopologicalObject & obj, const VertexHandle & v);
bool isEdgeManifold(const TopologicalObject & obj, const EdgeHandle & e);
bool isFaceManifold(const TopologicalObject & obj, const FaceHandle & f);
  
} //namespace BASim 

#endif //TOPOBJUTIl_H
