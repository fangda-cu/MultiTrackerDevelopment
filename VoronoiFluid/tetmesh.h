
#ifndef TETMESH_H
#define TETMESH_H

#include "geometryutils.h"
#include <vec.h>

#include "CGAL_wrapper.h"


class TetMesh
{
   
public:
   
   //
   // simplices
   //
   
   // tet vertices are in positive-volume orientation order
   // tets are unordered in the vector
   std::vector<ElTopo::Vec4ui> tets;

   // tris and edges have their vertices sorted and then are stored in lexicographical order in their vectors
   std::vector<ElTopo::Vec3ui> tris;
   std::vector<ElTopo::Vec2ui> edges;
   
   std::vector<ElTopo::Vec3f> vertices;
   
   //
   // incidence relationships
   //
   
   std::vector<ElTopo::Vec4ui> tet_to_tri_map; 
   std::vector<ElTopo::Vec6ui> tet_to_edge_map;
   std::vector<ElTopo::Vec3ui> tri_to_edge_map;   
   std::vector<bool> closed_tet_neighbourhood;     // whether the tets around this edge form a closed neighbourhood
   std::vector< std::vector<unsigned int> > edge_to_tet_map;   
   std::vector< std::vector<unsigned int> > vert_to_edge_map;
   std::vector< std::vector<unsigned int> > vert_to_tri_map;
   std::vector< std::vector<unsigned int> > vert_to_tet_map;

   
   //
   // cached geometry
   //
   
   std::vector<ElTopo::Vec3f> tet_circumcentres;
   
   std::vector<float> tet_volumes;
   std::vector<float> edge_lengths;
   
   std::vector<ElTopo::Vec3f> tet_edge_vectors;
   std::vector<ElTopo::Vec3f> tet_edge_midpoints;
   
   std::vector<float> voronoi_face_areas;
   std::vector<ElTopo::Vec3f> voronoi_face_centroids;

   // acceleration grid
   
   ElTopo::Array3< std::vector<int> > acceleration_grid;   
   float accel_dx;
   ElTopo::Vec3f accel_origin;
   unsigned int accel_ni, accel_nj, accel_nk;
   
   //CGAL tet structure for optimal point location in the Delaunay mesh
   Triangulation cgal_T;
   Cell_handle last_cell; //a hint for where to start from - let's see if this is helpful
   //
   // initialization functions
   //

   void add_tri( unsigned int a, unsigned int b, unsigned int c );
   void add_edge( unsigned int a, unsigned int b );
   
   void initialize( const std::vector<ElTopo::Vec4ui>& new_tets, const std::vector<ElTopo::Vec3f>& new_vertices, const Triangulation & T );
   
   void build_simplices( const std::vector<ElTopo::Vec4ui>& new_tets, const std::vector<ElTopo::Vec3f>& new_vertices );
   
   void sort_edge_to_tet_map();
   void build_incidences();
      
   void build_cached_geometry();
      
   void build_acceleration_grid();
   
   bool verify();
   
   //
   // query functions
   //
   
   inline int get_containing_tet( const ElTopo::Vec3f& point );
   inline int get_containing_voronoi( const ElTopo::Vec3f& point );

   void get_closed_tet_neighbourhood( unsigned int edge_index, std::vector<unsigned int>& sorted_incident_tets );
   
   float compute_voronoi_face_area( unsigned int edge_index );
   ElTopo::Vec3f compute_voronoi_face_centroid( unsigned int edge_index );
   
   //TODO If we can convert the next function to use CGAL's acceleration structure instead
   //we can eliminate the other acceleration grid.
   void get_overlapping_tets( const ElTopo::Vec3f& aabb_low, const ElTopo::Vec3f& aabb_high, std::vector<int>& query_results );

};


// ---------------------------------------------------------
// Inline functions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Return the index of the tet containing the specified point, or -1 if the point is not in any tet.
///
// ---------------------------------------------------------

extern unsigned int g_num_pit_tests;
extern unsigned int g_num_pit_hits;

inline int TetMesh::get_containing_tet( const ElTopo::Vec3f& point )
{
  
   Point p(point[0], point[1], point[2]);
   Cell_handle c = cgal_T.locate(p, last_cell);
   last_cell = c;

   return c->info();
  
  
}


// ---------------------------------------------------------
///
/// Return the index of the tet containing the specified point, or -1 if the point is not in any tet.
///
// ---------------------------------------------------------

inline int TetMesh::get_containing_voronoi( const ElTopo::Vec3f& point )
{
   

   Point p(point[0], point[1], point[2]);
   Vertex_handle vh = cgal_T.nearest_vertex(p, last_cell);
   last_cell = vh->cell();
   return vh->info();
  

}


// ---------------------------------------------------------
//
// Return true if the triangles have the same orientation
//
// ---------------------------------------------------------

inline bool same_oriented_triangle( const ElTopo::Vec3ui& tri_a, const ElTopo::Vec3ui& tri_b )
{
   if ( tri_a[0] == tri_b[0] && tri_a[1] == tri_b[1] && tri_a[2] == tri_b[2] ) { return true; }
   if ( tri_a[0] == tri_b[1] && tri_a[1] == tri_b[2] && tri_a[2] == tri_b[0] ) { return true; }
   if ( tri_a[0] == tri_b[2] && tri_a[1] == tri_b[0] && tri_a[2] == tri_b[1] ) { return true; }   
   return false;
}

// ---------------------------------------------------------
//
// Return true if the triangle shares the same vertices as a face on the tet, and has the same orientation.
//
// ---------------------------------------------------------

inline bool same_orientation( const ElTopo::Vec4ui& tet, const ElTopo::Vec3ui& tri )
{
   
   ElTopo::Vec3ui face( tet[0], tet[1], tet[2] );
   if ( same_oriented_triangle( face, tri ) ) { return true; }

   face = ElTopo::Vec3ui( tet[1], tet[2], tet[3] );
   if ( same_oriented_triangle( face, tri ) ) { return true; }

   face = ElTopo::Vec3ui( tet[2], tet[3], tet[0] );
   if ( same_oriented_triangle( face, tri ) ) { return true; }

   face = ElTopo::Vec3ui( tet[0], tet[1], tet[3] );
   if ( same_oriented_triangle( face, tri ) ) { return true; }

   return false;
}

// ---------------------------------------------------------
//
//
//
// ---------------------------------------------------------

inline bool tets_are_adjacent( const ElTopo::Vec4ui& tet_a, const ElTopo::Vec4ui& tet_b )
{
   unsigned int num_shared_vertices = 0;
   
   num_shared_vertices += ( tet_a[0] == tet_b[0] ) ? 1 : 0;
   num_shared_vertices += ( tet_a[0] == tet_b[1] ) ? 1 : 0;
   num_shared_vertices += ( tet_a[0] == tet_b[2] ) ? 1 : 0;
   num_shared_vertices += ( tet_a[0] == tet_b[3] ) ? 1 : 0;
   
   num_shared_vertices += ( tet_a[1] == tet_b[0] ) ? 1 : 0;
   num_shared_vertices += ( tet_a[1] == tet_b[1] ) ? 1 : 0;
   num_shared_vertices += ( tet_a[1] == tet_b[2] ) ? 1 : 0;
   num_shared_vertices += ( tet_a[1] == tet_b[3] ) ? 1 : 0;

   num_shared_vertices += ( tet_a[2] == tet_b[0] ) ? 1 : 0;
   num_shared_vertices += ( tet_a[2] == tet_b[1] ) ? 1 : 0;
   num_shared_vertices += ( tet_a[2] == tet_b[2] ) ? 1 : 0;
   num_shared_vertices += ( tet_a[2] == tet_b[3] ) ? 1 : 0;
   
   num_shared_vertices += ( tet_a[3] == tet_b[0] ) ? 1 : 0;
   num_shared_vertices += ( tet_a[3] == tet_b[1] ) ? 1 : 0;
   num_shared_vertices += ( tet_a[3] == tet_b[2] ) ? 1 : 0;
   num_shared_vertices += ( tet_a[3] == tet_b[3] ) ? 1 : 0;
   
   assert( num_shared_vertices >= 2 );
   
   return ( num_shared_vertices == 3 );   
}

#endif



