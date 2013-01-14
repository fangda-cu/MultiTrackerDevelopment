
#include "tetmesh.h"
#include "geometryutils.h"
#include "mat.h"


unsigned int g_num_pit_tests;
unsigned int g_num_pit_hits;

using namespace ElTopo;
// ---------------------------------------------------------

static void sort_three( unsigned int& a, unsigned int& b, unsigned int& c )
{
   if ( b < a ) { swap( a, b ); }
   if ( c < b ) { swap( b, c ); }
   if ( b < a ) { swap( a, b ); }   
}

// ---------------------------------------------------------


static bool edge_less( const Vec2st& edge_a, const Vec2st& edge_b )
{
   if ( edge_a[0] < edge_b[0] ) { return true; }
   if ( edge_a[0] > edge_b[0] ) { return false; }
   if ( edge_a[1] < edge_b[1] ) { return true; }
   assert( edge_b[1] < edge_a[1] || edge_a == edge_b );   
   return false;
}

// ---------------------------------------------------------

static bool tri_less( const Vec3st& tri_a, const Vec3st& tri_b )
{
   if ( tri_a[0] < tri_b[0] ) { return true; }
   if ( tri_a[0] > tri_b[0] ) { return false; }
   // a[0] == b[0]
   if ( tri_a[1] < tri_b[1] ) { return true; }
   if ( tri_a[1] > tri_b[1] ) { return false; }
   // a[1] == b[1]
   if ( tri_a[2] < tri_b[2] ) { return true; }
   if ( tri_a[2] > tri_b[2] ) { return false; }
   // a[2] == b[2]
   return false;
}


// ---------------------------------------------------------

static float convex_polygon_area( const std::vector<Vec3f>& polygon_vertices )
{
   Vec3f barycentre(0,0,0);
   for ( unsigned int i = 0; i < polygon_vertices.size(); ++i )
   {
      barycentre += polygon_vertices[i];
   }
   barycentre /= (float) polygon_vertices.size();
   
   float area = 0.0f;
   for ( unsigned int i = 0; i < polygon_vertices.size(); ++i )
   {
      unsigned int next = (i+1) % polygon_vertices.size();
      area += 0.5f * mag( cross( polygon_vertices[i] - barycentre, polygon_vertices[next] - barycentre ) );
   }
   
   return area;
   
}

// ---------------------------------------------------------

static bool tet_contains_edge( const Vec4st& tet, const Vec2st& edge )
{
   if (   ( tet[0] == edge[0] || tet[1] == edge[0] || tet[2] == edge[0] || tet[3] == edge[0] ) 
       && ( tet[0] == edge[1] || tet[1] == edge[1] || tet[2] == edge[1] || tet[3] == edge[1] ) )
   {
      return true;
   }
   return false;
}

// ---------------------------------------------------------

void TetMesh::add_tri( unsigned int a, unsigned int b, unsigned int c )
{
   assert( a != b && b != c && a != c );
   sort_three( a, b, c );
   assert( a < b && b < c );
   tris.push_back( Vec3st(a,b,c) );
}

// ---------------------------------------------------------

void TetMesh::add_edge( unsigned int a, unsigned int b )
{
   assert( a != b );
   if ( a > b ) { swap( a, b ); }
   assert( a < b );
   edges.push_back( Vec2st(a, b) );
}

// ---------------------------------------------------------

void TetMesh::initialize( const std::vector<Vec4st>& new_tets, const std::vector<Vec3f>& new_vertices, const Triangulation& T )
{

   cgal_T = T;

   std::cout << "building simplices..." << std::endl;
   
   build_simplices( new_tets, new_vertices );
   
   std::cout << "building incidence maps..." << std::endl;
   
   build_incidences();
   
   std::cout << "building cached geometry information..." << std::endl;
   
   build_cached_geometry();
   
   std::cout << "building acceleration grid..." << std::endl;
   
   build_acceleration_grid();
   
   //std::cout << "verifying..." << std::endl;
   //verify();
   
}


// ---------------------------------------------------------

void TetMesh::build_simplices( const std::vector<Vec4st>& new_tets, const std::vector<Vec3f>& new_vertices )
{
   
   //
   // verts
   //
   
   vertices = new_vertices;
   
   //
   // tets
   // 
   
   // add tets, in positive-volume orienation
   tets.clear();
   tets.reserve( new_tets.size() );
   for ( unsigned int i = 0; i < new_tets.size(); ++i )
   {
      unsigned int a = new_tets[i][0], b = new_tets[i][1], c = new_tets[i][2], d = new_tets[i][3];
      
      if ( a == b || a == c || a == d || b == c || b == d || c == d ) { continue; }
      
      if ( tet_signed_volume( vertices[a], vertices[b], vertices[c], vertices[d] ) > 0.0f )
      {
         if( tet_signed_volume( vertices[a], vertices[b], vertices[d], vertices[c] ) > 0.0f )
         {
            std::cout << "volABDC: " << tet_signed_volume( vertices[a], vertices[b], vertices[d], vertices[c] ) << std::endl;
         }
         
         tets.push_back( Vec4st(a,b,c,d) );
      }
      else
      {
         tets.push_back( Vec4st(a,b,d,c) );
      }
   }
   
   //
   // tris
   // 
   
   tris.clear();
   tris.reserve( 4 * tets.size() );
   
   // sort each triangle, but don't sort list of tris, and allow duplicates for now
   for ( unsigned int i = 0; i < tets.size(); ++i )
   {
      add_tri( tets[i][0], tets[i][1], tets[i][2] );
      add_tri( tets[i][1], tets[i][2], tets[i][3] );
      add_tri( tets[i][2], tets[i][3], tets[i][0] );
      add_tri( tets[i][0], tets[i][1], tets[i][3] );
   }
   
   // now sort and remove duplicates
   std::sort( tris.begin(), tris.end(), tri_less );
   tris.erase( std::unique( tris.begin(), tris.end() ), tris.end() );

   //
   // edges
   // 
   
   edges.clear();
   edges.reserve( 6 * tets.size() );
   // sort each edge, but don't sort list of edges, and allow duplicates for now
   for ( unsigned int i = 0; i < tets.size(); ++i )
   {
      add_edge( tets[i][0], tets[i][1] );
      add_edge( tets[i][1], tets[i][2] );
      add_edge( tets[i][2], tets[i][3] );
      add_edge( tets[i][3], tets[i][0] );
      add_edge( tets[i][0], tets[i][2] );
      add_edge( tets[i][1], tets[i][3] );
   }
   
   // now sort and remove duplicates
   std::sort( edges.begin(), edges.end(), edge_less );
   edges.erase( std::unique( edges.begin(), edges.end() ), edges.end() );
    
}

// ---------------------------------------------------------
// put each edge_to_tet list in order so that adjacent tets are next to each other in the list

void TetMesh::sort_edge_to_tet_map()
{
 
   closed_tet_neighbourhood.clear();
   closed_tet_neighbourhood.resize( edges.size(), false );

   for ( unsigned int edge_index = 0; edge_index < edges.size(); ++edge_index )
   {
      
      if ( edge_to_tet_map[edge_index].size() < 3 ) { continue; }     
      
      std::vector<unsigned int> incident_tets = edge_to_tet_map[edge_index];        // copy & destroy this...
      std::vector<unsigned int> sorted_incident_tets;                               // ...while building this
      
      sorted_incident_tets.push_back( incident_tets.front() );
      incident_tets.erase( incident_tets.begin() );

      bool neighbourhood_is_closed = true;
      
      while( !incident_tets.empty() )
      {
         
         // find a tet adjacent to the last tet in sorted_incident_tets
         
         bool found = false;
         const Vec4st& tet_b = tets[sorted_incident_tets.back()];
         
         for ( unsigned int i = 0; i < incident_tets.size(); ++i )
         {
            const Vec4st& tet_a = tets[incident_tets[i]]; 
            if ( tets_are_adjacent( tet_a, tet_b ) )
            {
               sorted_incident_tets.push_back( incident_tets[i] );
               incident_tets.erase( incident_tets.begin() + i );
               found = true;
               break;
            }
         }
         
         if ( !found )
         {
            // The incident tets don't form a closed ring.  This edge is on a boundary.
            neighbourhood_is_closed = false;
            break;
         }
         
      }
      
      if( neighbourhood_is_closed && !tets_are_adjacent( tets[sorted_incident_tets.back()], tets[sorted_incident_tets[0]] ) )
      {
         // The incident tets don't form a closed ring.  This edge is on a boundary.  The corresponding voronoi face has infinite area.
         neighbourhood_is_closed = false;
      }

      if ( !neighbourhood_is_closed )
      {
         // leave old map intact
         continue;
      }
      
      assert( sorted_incident_tets.size() == edge_to_tet_map[edge_index].size() );
      
      edge_to_tet_map[edge_index] = sorted_incident_tets;
      closed_tet_neighbourhood[edge_index] = true;
   }   
   
}


// ---------------------------------------------------------

void TetMesh::build_incidences()
{
   // verts -> edges
   vert_to_edge_map.clear();
   vert_to_edge_map.resize( vertices.size() );
   for ( unsigned int i = 0; i < edges.size(); ++i )
   {
      const Vec2st& e = edges[i];
      vert_to_edge_map[e[0]].push_back(i);
      vert_to_edge_map[e[1]].push_back(i);
   }

   // verts -> tris
   vert_to_tri_map.clear();
   vert_to_tri_map.resize( vertices.size() );
   for ( unsigned int i = 0; i < tris.size(); ++i )
   {
      const Vec3st& t = tris[i];
      vert_to_tri_map[t[0]].push_back(i);
      vert_to_tri_map[t[1]].push_back(i);
      vert_to_tri_map[t[2]].push_back(i);      
   }
   
   // verts -> tets
   vert_to_tet_map.clear();
   vert_to_tet_map.resize( vertices.size() );
   for ( unsigned int i = 0; i < tets.size(); ++i )
   {
      const Vec4st& t = tets[i];
      vert_to_tet_map[t[0]].push_back(i);
      vert_to_tet_map[t[1]].push_back(i);
      vert_to_tet_map[t[2]].push_back(i);      
      vert_to_tet_map[t[3]].push_back(i);            
   }
   
   // tets -> tris
   tet_to_tri_map.clear();
   tet_to_tri_map.resize( tets.size() );
   
   for ( unsigned int i = 0; i < tets.size(); ++i )
   {
      const Vec4st& t = tets[i];
      Vec3st tet_faces[4] = { Vec3st( t[0], t[1], t[2] ),
                              Vec3st( t[0], t[1], t[3] ),
                              Vec3st( t[0], t[2], t[3] ),
                              Vec3st( t[1], t[2], t[3] ) };
      
      for ( unsigned int j = 0; j < 4; ++j )
      {
         // make a triangle
         unsigned int a = tet_faces[j][0], b = tet_faces[j][1], c = tet_faces[j][2];
         sort_three( a, b, c );
         
         // find the triangle
         const std::vector<unsigned int>& incident_tris = vert_to_tri_map[a];
         unsigned int triangle_index = (unsigned int) ~0;
         for ( unsigned int t = 0; t < incident_tris.size(); ++t )
         {
            if ( tris[incident_tris[t]] == Vec3st(a,b,c) )
            {
               triangle_index = incident_tris[t];
            }
         }
         assert( triangle_index != (unsigned int) ~0 );
         
         // fill out the tet->tri map
         tet_to_tri_map[i][j] = triangle_index;
      }
   }
   
   // tets <-> edges
   tet_to_edge_map.clear();
   tet_to_edge_map.resize(tets.size());
   edge_to_tet_map.clear();
   edge_to_tet_map.resize(edges.size());
   for ( unsigned int i = 0; i < tets.size(); ++i )
   {
      const Vec4st& t = tets[i];
      Vec2st tet_edges[6] = { Vec2st(t[0], t[1]),
                              Vec2st(t[1], t[2]),
                              Vec2st(t[2], t[3]),
                              Vec2st(t[3], t[0]),
                              Vec2st(t[0], t[2]),
                              Vec2st(t[1], t[3]) };
   
      for ( unsigned int j = 0; j < 6; ++j )
      {
         if ( tet_edges[j][0] > tet_edges[j][1] ) { swap( tet_edges[j][0], tet_edges[j][1] ); }
         
         // find the edge
         unsigned int a = tet_edges[j][0];
         unsigned int b = tet_edges[j][1];
         
         const std::vector<unsigned int>& incident_edges = vert_to_edge_map[a];
         unsigned int edge_index = (unsigned int) ~0;
         for ( unsigned int e = 0; e < incident_edges.size(); ++e )
         {
            if ( edges[incident_edges[e]][1] == b )
            {
               edge_index = incident_edges[e];
            }
         }
         assert( edge_index != (unsigned int) ~0 );
         
         // fill out the tet->edge map
         tet_to_edge_map[i][j] = edge_index;
         
         // fill out the edge->tet map
         edge_to_tet_map[edge_index].push_back( i );
      }
   }
   
   // re-arrange incident tet lists so that adjacent tets are next to each other when the tet neighbourhood is closed
   sort_edge_to_tet_map();
   
   // tris -> edges
   tri_to_edge_map.clear();
   tri_to_edge_map.resize( tris.size(), Vec3st(~0, ~0, ~0) );
   for ( unsigned int i = 0; i < tris.size(); ++i )
   {
      Vec2st tri_edges[3] = { Vec2st( tris[i][0], tris[i][1] ),
                              Vec2st( tris[i][1], tris[i][2] ),
                              Vec2st( tris[i][0], tris[i][2] ) };
      
      for ( unsigned int j = 0; j < 3; ++j )
      {
         // find the edge
         unsigned int a = tri_edges[j][0];
         unsigned int b = tri_edges[j][1];
         
         const std::vector<unsigned int>& incident_edges = vert_to_edge_map[a];
         unsigned int edge_index = (unsigned int) ~0;
         for ( unsigned int e = 0; e < incident_edges.size(); ++e )
         {
            if ( edges[incident_edges[e]][1] == b )
            {
               edge_index = incident_edges[e];
            }
         }
         assert( edge_index != (unsigned int) ~0 );
         
         // fill out the tri->edge map
         tri_to_edge_map[i][j] = edge_index;
         
      }
   }
   
}


// ---------------------------------------------------------
//
// Returns all tets incident on an edge, in order s.t. adjecent tets are next to each other in the list.  If
// the set of incident tets doesn't form a closed ring, return an empty vector.
//
// ---------------------------------------------------------

void TetMesh::get_closed_tet_neighbourhood( unsigned int edge_index, std::vector<unsigned int>& sorted_incident_tets )
{   
   
   if ( !closed_tet_neighbourhood[edge_index] ) 
   { 
      sorted_incident_tets.clear();
   }
   else
   {
      sorted_incident_tets = edge_to_tet_map[edge_index];
   }
   
}


// ---------------------------------------------------------

float TetMesh::compute_voronoi_face_area( unsigned int edge_index )
{
   std::vector<unsigned int> sorted_incident_tets;
   get_closed_tet_neighbourhood( edge_index, sorted_incident_tets );
   
   if ( sorted_incident_tets.empty() ) { return 0.0f; }
   
   std::vector<Vec3f> voronoi_face_polygon;
   for ( unsigned int i = 0; i < sorted_incident_tets.size(); ++i )
   {
      voronoi_face_polygon.push_back( tet_circumcentres[ sorted_incident_tets[i] ] );
   }
   
   return convex_polygon_area( voronoi_face_polygon );
   
}

// ---------------------------------------------------------

Vec3f TetMesh::compute_voronoi_face_centroid( unsigned int edge_index ) 
{
   Vec3f barycentre(0,0,0);
   for ( unsigned int i = 0; i < edge_to_tet_map[edge_index].size(); ++i )
   {
      barycentre += tet_circumcentres[ edge_to_tet_map[edge_index][i] ];
   }
   return barycentre / (float)edge_to_tet_map[edge_index].size();
   
}

// ---------------------------------------------------------

void TetMesh::build_cached_geometry()
{
   
   // tet circumcentres
   
   tet_circumcentres.clear();
   tet_circumcentres.reserve( tets.size() );
   for ( unsigned int i = 0; i < tets.size(); ++i )
   {
      tet_circumcentres.push_back( tet_circumcentre( vertices[tets[i][0]], 
                                                     vertices[tets[i][1]], 
                                                     vertices[tets[i][2]], 
                                                     vertices[tets[i][3]] ) );
   }
   
   // tet volumes
   
   tet_volumes.resize( tets.size() );
   for ( unsigned int i = 0; i < tets.size(); ++i )
   {
      const Vec3f& a = vertices[ tets[i][0] ];
      const Vec3f& b = vertices[ tets[i][1] ];
      const Vec3f& c = vertices[ tets[i][2] ];
      const Vec3f& d = vertices[ tets[i][3] ];      
      tet_volumes[i] = tet_signed_volume( a, b, c, d );
      assert( tet_volumes[i] > 0.0f );
   }
   
   // edge lengths and vectors
   
   edge_lengths.resize( edges.size() );
   tet_edge_vectors.resize( edges.size() );
   tet_edge_midpoints.resize( edges.size() );
   for ( unsigned int i = 0; i < edges.size(); ++i )
   {
      tet_edge_vectors[i] = vertices[edges[i][1]] - vertices[edges[i][0]];
      edge_lengths[i] = mag( tet_edge_vectors[i] );
      assert ( edge_lengths[i] > 0 );
      tet_edge_vectors[i] /= edge_lengths[i];  
      assert( fabs( mag(tet_edge_vectors[i]) - 1.0f ) < 1e-6 );
      tet_edge_midpoints[i] = 0.5f * (vertices[edges[i][0]] + vertices[edges[i][1]]);
   }
   
   // compute voronoi face areas
   
   voronoi_face_areas.resize( edges.size() );
   voronoi_face_centroids.resize( edges.size() );
   for ( unsigned int i = 0; i < edges.size(); ++i )
   {   
      voronoi_face_areas[i] = compute_voronoi_face_area(i);
      voronoi_face_centroids[i] = compute_voronoi_face_centroid(i);
   }
   
}


// ---------------------------------------------------------

void TetMesh::build_acceleration_grid()
{
   // get average edge length
   
   float average_edge_length = 0.0f;
   for ( unsigned int i = 0; i < edge_lengths.size(); ++i )
   {
      average_edge_length += edge_lengths[i];
   }
   average_edge_length /= (float) ( edge_lengths.size() );
   
   // get bounds
   
   Vec3f min_x(1e30f), max_x(-1e30f);
   for ( unsigned int i = 0; i < vertices.size(); ++i )
   {
      min_x = min_union( min_x, vertices[i] );
      max_x = max_union( max_x, vertices[i] );
   }
   min_x -= Vec3f(average_edge_length);
   max_x += Vec3f(average_edge_length);
   
   assert( max_x[0] > min_x[0] );
   assert( max_x[1] > min_x[1] );
   assert( max_x[2] > min_x[2] );   
   
   Vec3st n_ijk = (Vec3st)((max_x - min_x) / average_edge_length);

   // build and populate the grid

   accel_dx = average_edge_length;
   accel_origin = min_x;
   
   accel_ni = n_ijk[0];
   accel_nj = n_ijk[1];
   accel_nk = n_ijk[2];;
   
   acceleration_grid.clear();
   acceleration_grid.resize( accel_ni, accel_nj, accel_nk );   
   
   for ( unsigned int t = 0; t < tets.size(); ++t )
   {
      const Vec3f& a = vertices[ tets[t][0] ];
      const Vec3f& b = vertices[ tets[t][1] ];
      const Vec3f& c = vertices[ tets[t][2] ];
      const Vec3f& d = vertices[ tets[t][3] ];
      
      Vec3f aabb_low = min_union( a, b );
      aabb_low = min_union( aabb_low, c );
      aabb_low = min_union( aabb_low, d );

      Vec3f aabb_high = max_union( a, b );
      aabb_high = max_union( aabb_high, c );
      aabb_high = max_union( aabb_high, d );

      //use the bound box to determine a range of cells to look at
      int i_lower = (int)((aabb_low[0] - accel_origin[0])/accel_dx-1);
      int j_lower = (int)((aabb_low[1] - accel_origin[1])/accel_dx-1);
      int k_lower = (int)((aabb_low[2] - accel_origin[2])/accel_dx-1);
      
      int i_upper = (int)((aabb_high[0] - accel_origin[0])/accel_dx+1);
      int j_upper = (int)((aabb_high[1] - accel_origin[1])/accel_dx+1);
      int k_upper = (int)((aabb_high[2] - accel_origin[2])/accel_dx+1);

      for(int i = i_lower; i < i_upper; ++i) 
      {
         for(int j = j_lower; j < j_upper; ++j) 
         {
            for(int k = k_lower; k < k_upper; ++k) 
            {
               acceleration_grid(i,j,k).push_back(t);
            }
         }
      }
      
   }
   printf("Completed building the regular acceleration grid.\n");

}


// ---------------------------------------------------------

bool TetMesh::verify()
{
   
   // verify positive volume
   for ( unsigned int i = 0; i < tets.size(); ++i )
   {
      assert( tet_signed_volume( vertices[tets[i][0]], vertices[tets[i][1]], vertices[tets[i][2]], vertices[tets[i][3]] ) > 0.0f );
   }
   
   // verify tris are sorted
   for ( unsigned int i = 0; i < tris.size()-1; ++i )
   {
      assert( tri_less( tris[i], tris[i+1] ) );
   }
   
   // verify edges are sorted   
   for ( unsigned int i = 0; i < edges.size() - 1; ++i )
   {
      assert( edge_less( edges[i], edges[i+1] ) );
   }
   
   // verts -> edges
   for ( unsigned int i = 0; i < edges.size(); ++i )
   {
      unsigned int a = edges[i][0];
      bool found = false;
      for ( unsigned int j = 0; j < vert_to_edge_map[a].size(); ++j )
      {
         if ( vert_to_edge_map[a][j] == i ) { found = true; }
      }
      assert( found );

      unsigned int b = edges[i][1];
      found = false;
      for ( unsigned int j = 0; j < vert_to_edge_map[b].size(); ++j )
      {
         if ( vert_to_edge_map[b][j] == i ) { found = true; }
      }
      assert( found );
   }
      
   // verts -> tris
   for ( unsigned int i = 0; i < tris.size(); ++i )
   {
      for ( unsigned int j = 0; j < 3; ++j )
      {
         unsigned int v = tris[i][j];
         bool found = false;
         for ( unsigned int k = 0; k < vert_to_tri_map[v].size(); ++k )
         {
            if ( vert_to_tri_map[v][k] == i ) { found = true; }
         }
         assert( found );
      }
   }
   
   // verts -> tets
   for ( unsigned int i = 0; i < tets.size(); ++i )
   {
      for ( unsigned int j = 0; j < 4; ++j )
      {
         unsigned int v = tets[i][j];
         assert( v < vertices.size() );
         bool found = false;
         for ( unsigned int k = 0; k < vert_to_tet_map[v].size(); ++k )
         {
            if ( vert_to_tet_map[v][k] == i ) { found = true; }
         }
         assert( found );
      }
   }
   
   // tets <-> tris
   for ( unsigned int i = 0; i < tets.size(); ++i )
   {
      
      float tet_volume = 0.0f;
      
      for ( unsigned int j = 0; j < 4; ++j )
      {
         unsigned int triangle_index = tet_to_tri_map[i][j];
         const Vec3st& tri = tris[ triangle_index ];
         assert( tri[0] == tets[i][0] || tri[0] == tets[i][1] || tri[0] == tets[i][2] || tri[0] == tets[i][3] );
         assert( tri[1] == tets[i][0] || tri[1] == tets[i][1] || tri[1] == tets[i][2] || tri[1] == tets[i][3] );
         assert( tri[2] == tets[i][0] || tri[2] == tets[i][1] || tri[2] == tets[i][2] || tri[2] == tets[i][3] );
         //assert( tri_to_tet_map[triangle_index][0] == i || tri_to_tet_map[triangle_index][1] == i );           
         
         float sign = same_orientation( tets[i], tri ) ? 1.0f : -1.0f;
         if ( j % 2 == 0 ) { sign *= -1.0f; }
         
         tet_volume += sign * 1.0f / 6.0f * triple( vertices[tri[0]], vertices[tri[1]], vertices[tri[2]] );
      }      
      
      if( tet_volume <= 0.0f )
      {
         std::cout << "tet_volume: " << tet_volume << std::endl;
         std::cout << "signed volume: " << tet_signed_volume( vertices[tets[i][0]], vertices[tets[i][1]], vertices[tets[i][2]], vertices[tets[i][3]] ) << std::endl;
         std::cout << "edge lengths: " << std::endl;
         std::cout << dist( vertices[tets[i][0]], vertices[tets[i][1]] ) << std::endl;
         std::cout << dist( vertices[tets[i][0]], vertices[tets[i][2]] ) << std::endl;
         std::cout << dist( vertices[tets[i][0]], vertices[tets[i][3]] ) << std::endl;
         std::cout << dist( vertices[tets[i][1]], vertices[tets[i][2]] ) << std::endl;
         std::cout << dist( vertices[tets[i][1]], vertices[tets[i][3]] ) << std::endl;
         std::cout << dist( vertices[tets[i][2]], vertices[tets[i][3]] ) << std::endl;
      }
      
   }
   
   // tets <-> edges
   for ( unsigned int i = 0; i < tets.size(); ++i )
   {
      for ( unsigned int j = 0; j < 6; ++j )
      {
         unsigned int edge_index = tet_to_edge_map[i][j];
         const Vec2st& edge = edges[ edge_index ];
         assert( edge[0] == tets[i][0] || edge[0] == tets[i][1] || edge[0] == tets[i][2] || edge[0] == tets[i][3] );
         assert( edge[1] == tets[i][0] || edge[1] == tets[i][1] || edge[1] == tets[i][2] || edge[1] == tets[i][3] );
         
         const std::vector<unsigned int>& incident_tets = edge_to_tet_map[edge_index];
         bool found = false;
         for ( unsigned int k = 0; k < incident_tets.size(); ++k )
         {
            if ( incident_tets[k] == i )
            {
               found = true;
               break;
            }   
         }
         assert(found);
      }      
   }
   
   for ( unsigned int i = 0; i < edges.size(); ++i )
   {
      if ( closed_tet_neighbourhood[i] )
      {
         for ( unsigned int j = 0; j < edge_to_tet_map[i].size() - 1; ++j )
         {
            assert( tets_are_adjacent( tets[edge_to_tet_map[i][j]], tets[edge_to_tet_map[i][j+1]] ) );
         }
      }
   }
   
   std::cout << "TetMesh verified" << std::endl;
   
   return true;
   
}



// ---------------------------------------------------------


void TetMesh::get_overlapping_tets( const Vec3f& aabb_low, const Vec3f& aabb_high, std::vector<int>& query_results )
{
   Vec3i cell_low = Vec3i(( aabb_low - accel_origin ) / accel_dx);
   Vec3i cell_high = Vec3i(( aabb_high - accel_origin ) / accel_dx);
   
   cell_low[0] = min( cell_low[0], (int)accel_ni - 1 );
   cell_low[1] = min( cell_low[1], (int)accel_nj - 1 );
   cell_low[2] = min( cell_low[2], (int)accel_nk - 1 );
   cell_low[0] = max( cell_low[0], 0 );
   cell_low[1] = max( cell_low[1], 0 );
   cell_low[2] = max( cell_low[2], 0 );

   cell_high[0] = min( cell_high[0], (int)accel_ni - 1 );
   cell_high[1] = min( cell_high[1], (int)accel_nj - 1 );
   cell_high[2] = min( cell_high[2], (int)accel_nk - 1 );
   cell_high[0] = max( cell_high[0], 0 );
   cell_high[1] = max( cell_high[1], 0 );
   cell_high[2] = max( cell_high[2], 0 );
   
   for ( int i = cell_low[0]; i <= cell_high[0]; ++i )
   {
      for ( int j = cell_low[1]; j <= cell_high[1]; ++j )
      {
         for ( int k = cell_low[2]; k <= cell_high[2]; ++k )
         {
            const std::vector<int>& cell = acceleration_grid(i,j,k);
            for ( unsigned int t = 0; t < cell.size(); ++t )
            {
               query_results.push_back( cell[t] );
            }
         }
      }
   }
   
}



