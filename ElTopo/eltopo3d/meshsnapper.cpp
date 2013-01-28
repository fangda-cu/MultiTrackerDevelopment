// ---------------------------------------------------------
//
//  meshsnapper.cpp
//  Christopher Batty 2013
//  
//  Functions supporting the "vertex snap" operation: merging nearby vertices.
//  The opposite of a pinch operation.
//
// ---------------------------------------------------------

#include <meshsnapper.h>

#include <broadphase.h>
#include <collisionpipeline.h>
#include <collisionqueries.h>
#include <nondestructivetrimesh.h>
#include <runstats.h>
#include <subdivisionscheme.h>
#include <surftrack.h>
#include <trianglequality.h>

// ---------------------------------------------------------
//  Extern globals
// ---------------------------------------------------------

namespace ElTopo {

extern RunStats g_stats;

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

// --------------------------------------------------------
///
/// Mesh Snapper constructor.  Takes a SurfTrack object.
///
// --------------------------------------------------------

MeshSnapper::MeshSnapper( SurfTrack& surf ) :
m_surf( surf )
{}


// --------------------------------------------------------
///
/// Get all triangles which are incident on either edge end vertex.
///
// --------------------------------------------------------

void MeshSnapper::get_moving_triangles(size_t source_vertex, 
   size_t destination_vertex, 
   std::vector<size_t>& moving_triangles )
{

   moving_triangles.clear();

   for ( size_t i = 0; i < m_surf.m_mesh.m_vertex_to_triangle_map[source_vertex].size(); ++i )
   {
      moving_triangles.push_back( m_surf.m_mesh.m_vertex_to_triangle_map[source_vertex][i] );
   }
   for ( size_t i = 0; i < m_surf.m_mesh.m_vertex_to_triangle_map[destination_vertex].size(); ++i )
   {
      moving_triangles.push_back( m_surf.m_mesh.m_vertex_to_triangle_map[destination_vertex][i] );
   }

}


// --------------------------------------------------------
///
/// Get all edges which are incident on either vertex.
///
// --------------------------------------------------------

void MeshSnapper::get_moving_edges( size_t source_vertex, 
   size_t destination_vertex, 
   std::vector<size_t>& moving_edges )
{

   moving_edges = m_surf.m_mesh.m_vertex_to_edge_map[ source_vertex ];
   moving_edges.insert( moving_edges.end(), m_surf.m_mesh.m_vertex_to_edge_map[ destination_vertex ].begin(), m_surf.m_mesh.m_vertex_to_edge_map[ destination_vertex ].end() );

}


// --------------------------------------------------------
///
/// Check the "pseudo motion" introduced by a collapsing edge for collision
///
// --------------------------------------------------------

bool MeshSnapper::snap_pseudo_motion_introduces_collision( size_t source_vertex, 
   size_t destination_vertex, 
   const Vec3d& )
{
   assert( m_surf.m_collision_safety );

   // Get the set of triangles which move because of this motion
   std::vector<size_t> moving_triangles;
   get_moving_triangles( source_vertex, destination_vertex, moving_triangles );

   // And the set of edges
   std::vector<size_t> moving_edges;
   get_moving_edges( source_vertex, destination_vertex, moving_edges );


   // Check for collisions, holding everything static except for the source and destination vertices

   CollisionCandidateSet collision_candidates;

   // triangle-point candidates
   for ( size_t i = 0; i < moving_triangles.size(); ++i )
   { 
      m_surf.m_collision_pipeline->add_triangle_candidates( moving_triangles[i], true, true, collision_candidates );
   }      

   // point-triangle candidates
   m_surf.m_collision_pipeline->add_point_candidates( source_vertex, true, true, collision_candidates );
   m_surf.m_collision_pipeline->add_point_candidates( destination_vertex, true, true, collision_candidates );

   // edge-edge candidates
   for ( size_t i = 0; i < moving_edges.size(); ++i )
   {
      m_surf.m_collision_pipeline->add_edge_candidates( moving_edges[i], true, true, collision_candidates );
   }

   // Prune collision candidates containing both the source and destination vertex (they will trivially be collisions )

   for ( size_t i = 0; i < collision_candidates.size(); ++i )
   {
      const Vec3st& candidate = collision_candidates[i];
      bool should_delete = false;

      if ( candidate[2] == 1 )
      {
         // edge-edge
         const Vec2st& e0 = m_surf.m_mesh.m_edges[ candidate[0] ];
         const Vec2st& e1 = m_surf.m_mesh.m_edges[ candidate[1] ];

         if ( e0[0] == source_vertex || e0[1] == source_vertex || e1[0] == source_vertex || e1[1] == source_vertex )
         {
            if ( e0[0] == destination_vertex || e0[1] == destination_vertex || e1[0] == destination_vertex || e1[1] == destination_vertex )
            {
               should_delete = true;
            }
         }

      }
      else
      {
         // point-triangle
         size_t t = candidate[0];
         const Vec3st& tri = m_surf.m_mesh.get_triangle(t);
         size_t v = candidate[1];

         if ( v == source_vertex || tri[0] == source_vertex || tri[1] == source_vertex || tri[2] == source_vertex )
         {
            if ( v == destination_vertex || tri[0] == destination_vertex || tri[1] == destination_vertex || tri[2] == destination_vertex )
            {
               should_delete = true;
            }
         }
      }

      if ( should_delete )
      {
         collision_candidates.erase( collision_candidates.begin() + i );
         --i;
      }

   }

   Collision collision;
   if ( m_surf.m_collision_pipeline->any_collision( collision_candidates, collision ) )
   {
      return true;
   }

   return false;

}


/*
// --------------------------------------------------------
///
/// Determine if the edge collapse operation would invert the normal of any incident triangles.
///
// --------------------------------------------------------

bool MeshSnapper::collapse_edge_introduces_normal_inversion( size_t source_vertex, 
size_t destination_vertex, 
size_t edge_index, 
const Vec3d& vertex_new_position )
{

// Get the set of triangles which are going to be deleted
std::vector< size_t >& triangles_incident_to_edge = m_surf.m_mesh.m_edge_to_triangle_map[edge_index];   

// Get the set of triangles which move because of this motion
std::vector<size_t> moving_triangles;
for ( size_t i = 0; i < m_surf.m_mesh.m_vertex_to_triangle_map[source_vertex].size(); ++i )
{
moving_triangles.push_back( m_surf.m_mesh.m_vertex_to_triangle_map[source_vertex][i] );
}
for ( size_t i = 0; i < m_surf.m_mesh.m_vertex_to_triangle_map[destination_vertex].size(); ++i )
{
moving_triangles.push_back( m_surf.m_mesh.m_vertex_to_triangle_map[destination_vertex][i] );
}

//
// check for normal inversion
//

for ( size_t i = 0; i < moving_triangles.size(); ++i )
{ 

// Disregard triangles which will end up being deleted - those triangles incident to the collapsing edge.
bool triangle_will_be_deleted = false;
for ( size_t j = 0; j < triangles_incident_to_edge.size(); ++j )
{
if ( moving_triangles[i] == triangles_incident_to_edge[j] )
{
triangle_will_be_deleted = true;
break;
}
}

if ( triangle_will_be_deleted ) { continue; }

const Vec3st& current_triangle = m_surf.m_mesh.get_triangle( moving_triangles[i] );
Vec3d old_normal = m_surf.get_triangle_normal( current_triangle );

Vec3d new_normal;

////////////////////////////////////////////////////////////
// FD 20121126
//
//  the new triangle always has the same orientation as the
//  old triangle so no change is needed
//
////////////////////////////////////////////////////////////

double new_area;
if ( current_triangle[0] == source_vertex || current_triangle[0] == destination_vertex )
{ 
new_normal = triangle_normal( vertex_new_position, m_surf.get_position(current_triangle[1]), m_surf.get_position(current_triangle[2]) ); 
new_area = triangle_area( vertex_new_position, m_surf.get_position(current_triangle[1]), m_surf.get_position(current_triangle[2]) ); 
}
else if ( current_triangle[1] == source_vertex || current_triangle[1] == destination_vertex ) 
{ 
new_normal = triangle_normal( m_surf.get_position(current_triangle[0]), vertex_new_position, m_surf.get_position(current_triangle[2]) ); 
new_area = triangle_area( m_surf.get_position(current_triangle[0]), vertex_new_position, m_surf.get_position(current_triangle[2]) ); 
}
else 
{ 
assert( current_triangle[2] == source_vertex || current_triangle[2] == destination_vertex ); 
new_normal = triangle_normal( m_surf.get_position(current_triangle[0]), m_surf.get_position(current_triangle[1]), vertex_new_position );
new_area = triangle_area( m_surf.get_position(current_triangle[0]), m_surf.get_position(current_triangle[1]), vertex_new_position );
}      

if ( dot( new_normal, old_normal ) < 1e-5 )
{
if ( m_surf.m_verbose ) { std::cout << "collapse edge introduces normal inversion" << std::endl; }

g_stats.add_to_int( "MeshSnapper:collapse_normal_inversion", 1 );

return true;
} 

if ( new_area < m_surf.m_min_triangle_area )
{
if ( m_surf.m_verbose ) { std::cout << "collapse edge introduces tiny triangle area" << std::endl; }

g_stats.add_to_int( "MeshSnapper:collapse_degenerate_triangle", 1 );

return true;
} 

}

return false;

}
*/

/*
// --------------------------------------------------------
///
/// Determine whether collapsing an edge will introduce an unacceptable change in volume.
///
// --------------------------------------------------------

bool MeshSnapper::collapse_edge_introduces_volume_change( size_t source_vertex, 
size_t edge_index, 
const Vec3d& vertex_new_position )
{
//
// If any incident triangle has a tiny area, collapse the edge without regard to volume change
//

const std::vector<size_t>& inc_tris = m_surf.m_mesh.m_edge_to_triangle_map[edge_index];

for ( size_t i = 0; i < inc_tris.size(); ++i )
{
if ( m_surf.get_triangle_area( inc_tris[i] ) < m_surf.m_min_triangle_area )
{
return false;
}
}

//
// Check volume change
//

const std::vector< size_t >& triangles_incident_to_vertex = m_surf.m_mesh.m_vertex_to_triangle_map[source_vertex];
double volume_change = 0;

for ( size_t i = 0; i < triangles_incident_to_vertex.size(); ++i )
{
const Vec3st& inc_tri = m_surf.m_mesh.get_triangle( triangles_incident_to_vertex[i] );
volume_change += signed_volume( vertex_new_position, m_surf.get_position(inc_tri[0]), m_surf.get_position(inc_tri[1]), m_surf.get_position(inc_tri[2]) );
}

if ( std::fabs(volume_change) > m_surf.m_max_volume_change )
{
if ( m_surf.m_verbose ) { std::cout << "collapse edge introduces volume change"  << std::endl; }
return true;
}

return false;

}
*/

/*
// ---------------------------------------------------------
///
/// Returns true if the edge collapse would introduce a triangle with a min or max angle outside of the specified min or max.
///
// ---------------------------------------------------------

bool MeshSnapper::collapse_edge_introduces_bad_angle(size_t source_vertex, 
size_t destination_vertex, 
const Vec3d& vertex_new_position )
{

std::vector<size_t> moving_triangles;
get_moving_triangles( source_vertex, destination_vertex,  moving_triangles );

for ( size_t i = 0; i < moving_triangles.size(); ++i )
{
const Vec3st& tri = m_surf.m_mesh.get_triangle( moving_triangles[i] );

///////////////////////////////////////////////////////////////////////
// FD 20121218
//
//  This is a bug in El Topo? This check always includes the two
//  triangles that are incident to the edge being collapsed (both
//  source_vertex and destination_vertex are in the triangle), and
//  for these triangles min_triangle_angle may return nan or a very
//  small (virtually zero) value. Random.
//
int count = 0;

Vec3d a = m_surf.get_position( tri[0] );

if ( tri[0] == source_vertex || tri[0] == destination_vertex )
{
a = vertex_new_position;
count++;
}

Vec3d b = m_surf.get_position( tri[1] );

if ( tri[1] == source_vertex || tri[1] == destination_vertex )
{
b = vertex_new_position;
count++;
}

Vec3d c = m_surf.get_position( tri[2] );

if ( tri[2] == source_vertex || tri[2] == destination_vertex )
{
c = vertex_new_position;
count++;
}

if (count == 2)
continue;

///////////////////////////////////////////////////////////////////////

double min_angle = min_triangle_angle( a, b, c );

if ( rad2deg(min_angle) < m_surf.m_min_triangle_angle )
{
return true;
}

double max_angle = max_triangle_angle( a, b, c );

if ( rad2deg(max_angle) > m_surf.m_max_triangle_angle )
{
return true;
}

}

return false;

}
*/

// --------------------------------------------------------
///
/// Delete an edge by moving its source vertex to its destination vertex
///
// --------------------------------------------------------

bool MeshSnapper::snap_vertex_pair( size_t vertex_to_keep, size_t vertex_to_delete)
{

   assert(m_surf.m_allow_non_manifold);

   bool keep_vert_is_boundary = m_surf.m_mesh.m_is_boundary_vertex[vertex_to_keep];
   bool del_vert_is_boundary = m_surf.m_mesh.m_is_boundary_vertex[vertex_to_delete];

   //either we're allowing remeshing of boundary edges, or this edge is not on the boundary.
   //TODO Figure out how to deal with boundaries here.

   ///////////////////////////////////////////////////////////////////////
   // FD 20130102
   //
   // do not snap if the two vertices are moving apart

   Vec3d edge_vec = m_surf.get_position(vertex_to_keep) - m_surf.get_position(vertex_to_delete);
   Vec3d rel_vel = m_surf.get_remesh_velocity(vertex_to_keep) - m_surf.get_remesh_velocity(vertex_to_delete);

   if (dot(rel_vel, edge_vec) > 0)
   {
      if (m_surf.m_verbose)
         std::cout << "The vertices are moving apart. No need to snap." << std::endl;
      return false;
   }

   ///////////////////////////////////////////////////////////////////////

   ///////////////////////////////////////////////////////////////////////
   // FD 20121229
   //
   // this should be allowed, in order to simulate pinching of films

   //  // *don't* collapse when both verts are on the boundary and the edge is not, since this would create a nonmanifold/singular vertex
   //  //assert( ((keep_vert_is_boundary && del_vert_is_boundary) && !edge_is_boundary) == false);   // these two asserts are equivalent
   //  assert(edge_is_boundary || !(keep_vert_is_boundary && del_vert_is_boundary));

   ///////////////////////////////////////////////////////////////////////

   if ( m_surf.m_verbose ) { std::cout << "Snapping vertices.  Doomed vertex: " << vertex_to_delete << " --- Vertex to keep: " << vertex_to_keep << std::endl; }
  

   // --------------

   {
      /*
      const std::vector< size_t >& r_triangles_incident_to_edge = m_surf.m_mesh.m_edge_to_triangle_map[edge];

      // Do not collapse edge on a degenerate tet or degenerate triangle
      for ( size_t i=0; i < r_triangles_incident_to_edge.size(); ++i )    // FD 20121126: Why can such deg tet/tri exist at all? trim_non_manifold() is called after every collapse_edge().
      {
         const Vec3st& triangle_i = m_surf.m_mesh.get_triangle( r_triangles_incident_to_edge[i] );

         if ( triangle_i[0] == triangle_i[1] || triangle_i[1] == triangle_i[2] || triangle_i[2] == triangle_i[0] )
         {
            if ( m_surf.m_verbose ) { std::cout << "duplicate vertices on triangle" << std::endl; }
            return false;
         }

         for ( size_t j=i+1; j < r_triangles_incident_to_edge.size(); ++j )
         {
            const Vec3st& triangle_j = m_surf.m_mesh.get_triangle( r_triangles_incident_to_edge[j] );

            if ( NonDestructiveTriMesh::triangle_has_these_verts( triangle_i, triangle_j ) )            
            {
               if ( m_surf.m_verbose ) { std::cout << "two triangles share vertices" << std::endl; }
               g_stats.add_to_int( "MeshSnapper:collapse_degen_tet", 1 );
               return false;
            }
         }
      }
      */
   }


   // --------------
   // decide on new vertex position
   // for now, we'll just take the average of the two endpoints.

   Vec3d vertex_new_position = 0.5*(m_surf.get_position(vertex_to_delete) + m_surf.get_position(vertex_to_keep));

   if ( m_surf.m_verbose ) { std::cout << "Collapsing edge.  Doomed vertex: " << vertex_to_delete << " --- Vertex to keep: " << vertex_to_keep << std::endl; }



   // --------------

   // Check vertex pseudo motion for collisions and volume change

   if ( mag ( m_surf.get_position(vertex_to_delete) - m_surf.get_position(vertex_to_keep) ) > 0 )
   {

      // Change source vertex predicted position to superimpose onto destination vertex
      m_surf.set_newposition( vertex_to_keep, vertex_new_position );
      m_surf.set_newposition( vertex_to_delete, vertex_new_position );

      /*
      bool volume_change = collapse_edge_introduces_volume_change( vertex_to_delete, edge, vertex_new_position );

      if ( volume_change )
      {
         // Restore saved positions which were changed by the function we just called.
         m_surf.set_newposition( vertex_to_keep, m_surf.get_position(vertex_to_keep) );
         m_surf.set_newposition( vertex_to_delete, m_surf.get_position(vertex_to_delete) );

         g_stats.add_to_int( "MeshSnapper:collapse_volume_change", 1 );

         if ( m_surf.m_verbose ) { std::cout << "collapse_volume_change" << std::endl; }
         return false;
      }

      bool normal_inversion = collapse_edge_introduces_normal_inversion(  vertex_to_delete, vertex_to_keep, edge, vertex_new_position );

      if ( normal_inversion )
      {
         // Restore saved positions which were changed by the function we just called.
         m_surf.set_newposition( vertex_to_keep, m_surf.get_position(vertex_to_keep) );
         m_surf.set_newposition( vertex_to_delete, m_surf.get_position(vertex_to_delete) );

         if ( m_surf.m_verbose ) { std::cout << "normal_inversion" << std::endl; }
         return false;
      }

      bool bad_angle = collapse_edge_introduces_bad_angle( vertex_to_delete, vertex_to_keep, vertex_new_position );

      if ( bad_angle )
      {
         // Restore saved positions which were changed by the function we just called.
         m_surf.set_newposition( vertex_to_keep, m_surf.get_position(vertex_to_keep) );
         m_surf.set_newposition( vertex_to_delete, m_surf.get_position(vertex_to_delete) );

         if ( m_surf.m_verbose ) { std::cout << "bad_angle" << std::endl; }

         g_stats.add_to_int( "MeshSnapper:collapse_bad_angle", 1 );
         return false;

      }
      */

      bool collision = false;

      if ( m_surf.m_collision_safety )
      {
         collision = snap_pseudo_motion_introduces_collision( vertex_to_delete, vertex_to_keep, vertex_new_position );
      }

      if ( collision ) 
      { 
         if ( m_surf.m_verbose ) { std::cout << "collision" << std::endl; }
         g_stats.add_to_int( "MeshSnapper:collapse_collisions", 1 ); 
      }

      // Restore saved positions which were changed by the function we just called.
      m_surf.set_newposition( vertex_to_keep, m_surf.get_position(vertex_to_keep) );
      m_surf.set_newposition( vertex_to_delete, m_surf.get_position(vertex_to_delete) );

      if ( collision )
      {
         // edge collapse would introduce collision or change volume too much or invert triangle normals
         return false;
      }
   }

   // --------------

   // start building history data
   MeshUpdateEvent collapse(MeshUpdateEvent::SNAP);
   collapse.m_vert_position = vertex_new_position;
   collapse.m_v0 = vertex_to_keep;
   collapse.m_v1 = vertex_to_delete;

   // move the vertex we decided to keep

   m_surf.set_position( vertex_to_keep, vertex_new_position );
   m_surf.set_newposition( vertex_to_keep, vertex_new_position );

   ///////////////////////////////////////////////////////////////////////
   // FD 20121229
   //
   // update the vertex constraint label
   //m_surf.m_mesh.set_vertex_constraint_label(vertex_to_keep, new_vert_constraint_label);

   ///////////////////////////////////////////////////////////////////////

  
   // Find anything pointing to the doomed vertex and change it

   // copy the list of triangles, don't take a reference to it
   std::vector< size_t > triangles_incident_to_vertex = m_surf.m_mesh.m_vertex_to_triangle_map[vertex_to_delete];    
   std::cout << "Vertices: " << vertex_to_keep << " and " << vertex_to_delete << std::endl;
   for ( size_t i=0; i < triangles_incident_to_vertex.size(); ++i )
   {
      //don't bother copying over dead tris... cuz why would we, eh?
      if(m_surf.m_mesh.triangle_is_deleted(triangles_incident_to_vertex[i])) continue;
      
      

      Vec3st new_triangle = m_surf.m_mesh.get_triangle( triangles_incident_to_vertex[i] );
      
      std::cout << "Original triangle: " << new_triangle << std::endl;

      if ( new_triangle[0] == vertex_to_delete )   { new_triangle[0] = vertex_to_keep; }
      if ( new_triangle[1] == vertex_to_delete )   { new_triangle[1] = vertex_to_keep; }
      if ( new_triangle[2] == vertex_to_delete )   { new_triangle[2] = vertex_to_keep; }

      std::cout << "Edited triangle: " << new_triangle << std::endl;
      if ( m_surf.m_verbose ) { std::cout << "adding updated triangle: " << new_triangle << std::endl; }

      size_t new_triangle_index = m_surf.add_triangle( new_triangle );
      collapse.m_created_tris.push_back( new_triangle_index );
      collapse.m_created_tri_data.push_back(new_triangle);

      ////////////////////////////////////////////////////////////
      // FD 20121126
      //
      // the old label carries over to the new triangle.
      // no need to test for orientation because the new triangle
      // generation code above does not change orientation.
      //
      Vec2i label = m_surf.m_mesh.get_triangle_label(triangles_incident_to_vertex[i]);

      m_surf.m_mesh.set_triangle_label(new_triangle_index, label);
      collapse.m_created_tri_labels.push_back(label);

      ////////////////////////////////////////////////////////////

      m_surf.m_dirty_triangles.push_back( new_triangle_index );
   }

   for ( size_t i=0; i < triangles_incident_to_vertex.size(); ++i )
   {  
      if ( m_surf.m_verbose )
      {
         std::cout << "removing vertex-incident triangle: " << m_surf.m_mesh.get_triangle( triangles_incident_to_vertex[i] ) << std::endl;
      }

      m_surf.remove_triangle( triangles_incident_to_vertex[i] );
      collapse.m_deleted_tris.push_back(triangles_incident_to_vertex[i]);
   }

   // Delete vertex
   assert( m_surf.m_mesh.m_vertex_to_triangle_map[vertex_to_delete].empty() );
   m_surf.remove_vertex( vertex_to_delete );
   collapse.m_deleted_verts.push_back( vertex_to_delete );

   m_surf.m_mesh.update_is_boundary_vertex( vertex_to_keep );

   // Store the history
   m_surf.m_mesh_change_history.push_back(collapse);
   return true;
}

// --------------------------------------------------------
///
/// Determine if the edge should be allowed to collapse
///
// --------------------------------------------------------

bool MeshSnapper::vert_pair_is_snappable( size_t vert0, size_t vert1, double& current_length )
{

   // skip deleted vertices
   if(m_surf.m_mesh.vertex_is_deleted(vert0) || m_surf.m_mesh.vertex_is_deleted(vert1) ) 
      return false;

   //shouldn't be calling this on a duplicate vert
   assert(vert0 != vert1);

   //Check if there is an edge connecting the two verts - if so, we can't do the snap
   size_t edge_id = m_surf.m_mesh.get_edge_index(vert0, vert1);
   if(edge_id != m_surf.m_mesh.m_edges.size())
      return false;

   
   //  if ( m_surf.edge_is_any_solid(edge_index) ) { return false; }

   //skip boundary edges if we're not remeshing those
   //if(!m_remesh_boundaries && m_surf.m_mesh.m_is_boundary_edge[edge_index]) { return false; }

   ///////////////////////////////////////////////////////////////////////
   // FD 20121229
   //
   // this should be allowed, in order to simulate pinching of films

   //  //this would introduce non-manifold ("singular") boundary vertex (this is an internal edge joining two boundary vertices)
   //  if ( m_surf.m_mesh.m_edge_to_triangle_map[edge_index].size() == 2 && 
   //      m_surf.m_mesh.m_is_boundary_vertex[m_surf.m_mesh.m_edges[edge_index][0]] &&
   //      m_surf.m_mesh.m_is_boundary_vertex[m_surf.m_mesh.m_edges[edge_index][1]] ) 
   //  {
   //    return false;
   //  }

   ///////////////////////////////////////////////////////////////////////

   current_length = mag(m_surf.get_position(vert0) - m_surf.get_position(vert1));
   
   return current_length < m_surf.m_merge_proximity_epsilon;  
   

}


// --------------------------------------------------------
///
/// Collapse all short edges
///
// --------------------------------------------------------

bool MeshSnapper::snap_pass()
{

   if ( m_surf.m_verbose )
   {
      std::cout << "\n\n\n---------------------- MeshSnapper: collapsing ----------------------" << std::endl;
      std::cout << "m_merge_proximity_epsilon: " << m_surf.m_merge_proximity_epsilon;
     
   }

   bool snap_occurred = false;

   assert( m_surf.m_dirty_triangles.size() == 0 );

   std::vector<SortablePair> sortable_pairs_to_try;


   //
   // get set of vertex pairs to snap!
   //

   for( size_t i = 0; i < m_surf.get_num_vertices(); i++ )
   {    
      if(m_surf.m_mesh.vertex_is_deleted(i)) continue;

      for(size_t j = i+1; j < m_surf.get_num_vertices(); ++j)  {
         
         if(m_surf.m_mesh.vertex_is_deleted(j)) continue;

         double cur_len;
         if(vert_pair_is_snappable(i,j,cur_len)) {
            sortable_pairs_to_try.push_back( SortablePair( i, j, cur_len ) );
         }
      }
   }

   //
   // sort in ascending order by length (collapse shortest edges first)
   //

   std::sort( sortable_pairs_to_try.begin(), sortable_pairs_to_try.end() );

   //if ( m_surf.m_verbose )
   //{
      std::cout << sortable_pairs_to_try.size() << " candidate pairs sorted" << std::endl;
   //}

   //
   // attempt to collapse each edge in the sorted list
   //

   for ( size_t si = 0; si < sortable_pairs_to_try.size(); ++si )
   {
      size_t v0 = sortable_pairs_to_try[si].m_vertex0;
      size_t v1 = sortable_pairs_to_try[si].m_vertex1;

      assert( v0 < m_surf.get_num_vertices() );
      assert( v1 < m_surf.get_num_vertices() );

      double dummy;
      if(vert_pair_is_snappable(v0, v1, dummy)) {
         std::cout << "Attempting vertex snapping.\n";
         bool result = snap_vertex_pair( v0, v1 );

         if ( result )
         { 
            std::cout << "Snapped vertices!\n";
            // clean up degenerate triangles and tets
            m_surf.trim_non_manifold( m_surf.m_dirty_triangles );            
         }
         else
            std::cout << "Snap failed!\n";

         snap_occurred |= result;
      }
   }

   return snap_occurred;

}

}
