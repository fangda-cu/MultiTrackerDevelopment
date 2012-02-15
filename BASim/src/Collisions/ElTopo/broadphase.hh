// ---------------------------------------------------------
//
//  broadphase.h
//  Tyson Brochu 2008
//  
//  Interface for abstract broad phase collision detector class.  
//  Used to avoid performing collision detection between all 
//  primitives. Abstract so we can try different strategies.
//
// ---------------------------------------------------------

#ifndef ELTOPO_BROADPHASE_H
#define ELTOPO_BROADPHASE_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------
#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Core/TopologicalObject/TopologicalObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjProperty.hh"

#include "vec.hh"

namespace ElTopoCode {
// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

// ---------------------------------------------------------
//  Interface declarations
// ---------------------------------------------------------

// --------------------------------------------------------
///
/// Abstract broad phase collision detector
///
// --------------------------------------------------------

class BroadPhase
{   
public:
   
   virtual ~BroadPhase() 
   {}
   
   /// Rebuild the broad phase using current vertex positions
   ///
   virtual void update_broad_phase_static( const BASim::TopologicalObject& m_obj, const BASim::VertexProperty<BASim::Vec3d>& vertices, double proximity_eps ) = 0;

   /// Rebuild the broad phase using current and predicted vertex positions
   ///
   //virtual void update_broad_phase_continuous( const BASim::TopologicalObject& m_obj, const BASim::VertexProperty<Vec3d>& vertices  ) = 0;

   virtual void add_vertex( unsigned int index, const Vec3d& aabb_low, const Vec3d& aabb_high ) = 0; 
   virtual void add_edge( unsigned int index, const Vec3d& aabb_low, const Vec3d& aabb_high ) = 0; 
   virtual void add_triangle( unsigned int index, const Vec3d& aabb_low, const Vec3d& aabb_high ) = 0; 

   virtual void update_vertex( unsigned int index, const Vec3d& aabb_low, const Vec3d& aabb_high ) = 0; 
   virtual void update_edge( unsigned int index, const Vec3d& aabb_low, const Vec3d& aabb_high ) = 0; 
   virtual void update_triangle( unsigned int index, const Vec3d& aabb_low, const Vec3d& aabb_high ) = 0; 
   
   virtual void remove_vertex( unsigned int index ) = 0; 
   virtual void remove_edge( unsigned int index ) = 0; 
   virtual void remove_triangle( unsigned int index ) = 0; 
   
   /// Get the set of vertices whose bounding volumes overlap the specified bounding volume
   ///
   virtual void get_potential_vertex_collisions( const Vec3d& query_low, const Vec3d& query_high, std::vector<unsigned int>& overlapping_vertices ) = 0;

   /// Get the set of edges whose bounding volumes overlap the specified bounding volume
   ///
   virtual void get_potential_edge_collisions( const Vec3d& query_low, const Vec3d& query_high, std::vector<unsigned int>& overlapping_triangles ) = 0;
   
   /// Get the set of triangles whose bounding volumes overlap the specified bounding volume
   ///
   virtual void get_potential_triangle_collisions( const Vec3d& query_low, const Vec3d& query_high, std::vector<unsigned int>& overlapping_triangles ) = 0;
   
};

}

#endif

