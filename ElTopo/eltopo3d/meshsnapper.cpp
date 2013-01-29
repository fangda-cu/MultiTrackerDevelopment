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
#include <edgesplitter.h>
#include <facesplitter.h>

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
   m_surf( surf ), 
   m_edgesplitter(surf, false, false, 0), 
   m_facesplitter(surf)
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
/// Check the "pseudo motion" introduced by snapping vertices for collision
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
/// Determine whether the operation will introduce an unacceptable change in volume...
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



// --------------------------------------------------------
///
/// Snap a pair of vertices by moving the source vertex to the destination vertex
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

   /*
   Vec3d edge_vec = m_surf.get_position(vertex_to_keep) - m_surf.get_position(vertex_to_delete);
   Vec3d rel_vel = m_surf.get_remesh_velocity(vertex_to_keep) - m_surf.get_remesh_velocity(vertex_to_delete);

   if (dot(rel_vel, edge_vec) > 0)
   {
      if (m_surf.m_verbose)
         std::cout << "The vertices are moving apart. No need to snap." << std::endl;
      return false;
   }
   */

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
   // for now, we'll just take the average of the two points to be snapped.

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
   
   for ( size_t i=0; i < triangles_incident_to_vertex.size(); ++i )
   {
      //don't bother copying over dead tris. (Can this ever happen?)
      if(m_surf.m_mesh.triangle_is_deleted(triangles_incident_to_vertex[i])) continue;
      
      Vec3st new_triangle = m_surf.m_mesh.get_triangle( triangles_incident_to_vertex[i] );
   
      if ( new_triangle[0] == vertex_to_delete )   { new_triangle[0] = vertex_to_keep; }
      if ( new_triangle[1] == vertex_to_delete )   { new_triangle[1] = vertex_to_keep; }
      if ( new_triangle[2] == vertex_to_delete )   { new_triangle[2] = vertex_to_keep; }

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
/// Snap an edge-edge pair by splitting each, and snapping the vertices
///
// --------------------------------------------------------

bool MeshSnapper::snap_edge_pair( size_t edge0, size_t edge1)
{

   assert(m_surf.m_allow_non_manifold);

   //TODO Figure out how to deal with boundaries here.

   // TODO Check relative velocities
   // do not snap if the two edges are moving apart
   /*
   Vec3d edge_vec = m_surf.get_position(vertex_to_keep) - m_surf.get_position(vertex_to_delete);
   Vec3d rel_vel = m_surf.get_remesh_velocity(vertex_to_keep) - m_surf.get_remesh_velocity(vertex_to_delete);

   if (dot(rel_vel, edge_vec) > 0)
   {
      if (m_surf.m_verbose)
         std::cout << "The vertices are moving apart. No need to snap." << std::endl;
      return false;
   }
   */
   ///////////////////////////////////////////////////////////////////////


   //if ( m_surf.m_verbose ) { std::cout << "Snapping edges.  Doomed vertex: " << vertex_to_delete << " --- Vertex to keep: " << vertex_to_keep << std::endl; }
   double distance, s0, s2;
   Vec3d normal;
   const Vec2st& edge_data0 = m_surf.m_mesh.m_edges[edge0];
   const Vec2st& edge_data1 = m_surf.m_mesh.m_edges[edge1];

   const Vec3d& v0 = m_surf.get_position(edge_data0[0]);
   const Vec3d& v1 = m_surf.get_position(edge_data0[1]);
   const Vec3d& v2 = m_surf.get_position(edge_data1[0]);
   const Vec3d& v3 = m_surf.get_position(edge_data1[1]);
   
   check_edge_edge_proximity( v0, v1, v2, v3, 
      distance, s0, s2, normal );
   
   //compute the resulting closest points
   Vec3d midpoint0 = s0*v0 + (1-s0)*v1;
   Vec3d midpoint1 = s2*v2 + (1-s2)*v3;
   
   //Check if we're fairly close to an end vertex; if so just use the vertex directly for snapping.
   //otherwise, do a split to create a new point which will then be snapped.
   double snap_threshold = 0.05;
   
   size_t snapping_vert0;
   if(s0 < snap_threshold) {
      snapping_vert0 = edge_data0[0];
   }
   else if(s0 > 1 - snap_threshold) {
      snapping_vert0 = edge_data0[1];
   }
   else {
      size_t split_result;
      
      std::cout << "Attempting to edge split.\n";
      if(m_edgesplitter.edge_is_splittable(edge0) || !m_edgesplitter.split_edge(edge0, split_result, true, &midpoint0))
         return false;

      snapping_vert0 = split_result;
   }

   size_t snapping_vert1;
   if(s2 < snap_threshold) {
      snapping_vert1 = edge_data1[0];
   }
   else if(s2 > 1 - snap_threshold) {
      snapping_vert1 = edge_data1[1];
   }
   else {
      size_t split_result;

      std::cout << "Attempting to edge split.\n";
      if(m_edgesplitter.edge_is_splittable(edge1) || !m_edgesplitter.split_edge(edge1, split_result, true, &midpoint1))
         return false;

      snapping_vert1 = split_result;
   }
   
   std::cout << "Attempting to snap vertex pair.\n";

   bool success = vert_pair_is_snappable(snapping_vert0, snapping_vert1) && snap_vertex_pair(snapping_vert0, snapping_vert1);
   
   return success;
}

// --------------------------------------------------------
///
/// Snap an edge-edge pair by splitting each, and snapping the vertices
///
// --------------------------------------------------------

bool MeshSnapper::snap_face_vertex_pair( size_t face, size_t vertex)
{

   assert(m_surf.m_allow_non_manifold);

   //TODO Figure out how to deal with boundaries here.

   // TODO Check relative velocities
   // do not snap if the vertex and face are separating
   /*
   Vec3d edge_vec = m_surf.get_position(vertex_to_keep) - m_surf.get_position(vertex_to_delete);
   Vec3d rel_vel = m_surf.get_remesh_velocity(vertex_to_keep) - m_surf.get_remesh_velocity(vertex_to_delete);

   if (dot(rel_vel, edge_vec) > 0)
   {
      if (m_surf.m_verbose)
         std::cout << "The vertices are moving apart. No need to snap." << std::endl;
      return false;
   }
   */
   ///////////////////////////////////////////////////////////////////////


   double s0, s1, s2;
   Vec3d normal;
   Vec3st face_data = m_surf.m_mesh.m_tris[face];
   double dist;

   const Vec3d& v_pos = m_surf.get_position(vertex);
   const Vec3d& t0_pos = m_surf.get_position(face_data[0]);
   const Vec3d& t1_pos = m_surf.get_position(face_data[1]);
   const Vec3d& t2_pos = m_surf.get_position(face_data[2]);

   //determine the distance, and 
   //barycentric coords of the closest point
   check_point_triangle_proximity(v_pos, 
      t0_pos, t1_pos, t2_pos,
      dist, s0, s1, s2, normal );

   std::cout << "Barycentric coords: " << s0 << " " << s1 << " " << s2 << ".\n";

   //Depending on the barycentric coordinates, either snap to one of the face vertices,
   //split an edge and snap to it, or split the face and snap to it.

   double snap_threshold = 0.05;
   size_t snapping_vertex;
   if(s0 < snap_threshold) {
      if(s1 < snap_threshold) {
         snapping_vertex = face_data[2];
      }
      else if(s2 < snap_threshold) {
         snapping_vertex = face_data[1];
      }
      else {
       size_t result_vertex;
       size_t edge_to_split = m_surf.m_mesh.get_edge_index(face_data[1], face_data[2]);
       
       double edge_frac = s1 / (s1+s2);
       Vec3d split_point = edge_frac * t1_pos + (1-edge_frac) * t2_pos;
       
       std::cout << "Attempting to edge split.\n";
       if(!m_edgesplitter.edge_is_splittable(edge_to_split) || !m_edgesplitter.split_edge(edge_to_split, result_vertex, true, &split_point))
         return false;
       
       snapping_vertex = result_vertex;
      }
   }
   else if(s1 < snap_threshold) {
      if(s2 < snap_threshold) {
         snapping_vertex = face_data[0];
      }
      else {
         size_t result_vertex;
         size_t edge_to_split = m_surf.m_mesh.get_edge_index(face_data[0], face_data[2]);

         double edge_frac = s0 / (s0+s2);
         Vec3d split_point = edge_frac * t0_pos + (1-edge_frac) * t2_pos;

         std::cout << "Attempting to edge split.\n";
         if(!m_edgesplitter.edge_is_splittable(edge_to_split) || !m_edgesplitter.split_edge(edge_to_split, result_vertex, true, &split_point))
            return false;

         snapping_vertex = result_vertex;
      }
   }
   else if(s2 < snap_threshold) {
      size_t result_vertex;
      size_t edge_to_split = m_surf.m_mesh.get_edge_index(face_data[0], face_data[1]);
      
      double edge_frac = s0 / (s0+s1);
      Vec3d split_point = edge_frac * t0_pos + (1-edge_frac) * t1_pos;

      std::cout << "Attempting to edge split.\n";
      if(!m_edgesplitter.edge_is_splittable(edge_to_split) || !m_edgesplitter.split_edge(edge_to_split, result_vertex, true, &split_point))
         return false;
      
      snapping_vertex = result_vertex;
   }
   else{
      size_t result_vertex;
      Vec3d split_point = s0*t0_pos + s1*t1_pos + s2*t2_pos;
      
      //try to split the face, and if it fails, drop out.
      std::cout << "Attempting to split face.\n";
      if(!m_facesplitter.face_is_splittable(face) || !m_facesplitter.split_face(face, result_vertex, true, &split_point))
         return false;

      snapping_vertex = result_vertex;
   }

   std::cout << "Attempting to snap resulting vertex pair.\n";
   bool success = vert_pair_is_snappable(snapping_vertex, vertex) && snap_vertex_pair(snapping_vertex, vertex);
   
   return success;
}

// --------------------------------------------------------
///
/// Determine if the vertex pair should be allowed to snap
///
// --------------------------------------------------------

bool MeshSnapper::vert_pair_is_snappable( size_t vert0, size_t vert1 )
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

   //current_length = mag(m_surf.get_position(vert0) - m_surf.get_position(vert1));
   
   return true; //current_length < m_surf.m_merge_proximity_epsilon;  
   

}

// --------------------------------------------------------
///
/// Determine if the edge edge pair should be snapped
///
// --------------------------------------------------------

bool MeshSnapper::edge_pair_is_snappable( size_t edge0, size_t edge1, double& current_length )
{

   // skip deleted vertices
   if(m_surf.m_mesh.edge_is_deleted(edge0) || m_surf.m_mesh.edge_is_deleted(edge1) ) 
      return false;
   
   //always use the lower numbered edge, to avoid duplication
   if(edge0 >= edge1)
      return false;

   //edges shouldn't share a vertex
   const Vec2st& edge_data0 = m_surf.m_mesh.m_edges[edge0];
   const Vec2st& edge_data1 = m_surf.m_mesh.m_edges[edge1];
   if(edge_data0[0] == edge_data1[0] || edge_data0[0] == edge_data1[1] ||
      edge_data0[1] == edge_data1[0] || edge_data0[1] == edge_data1[1] ) 
      return false;

   //TODO extend to handle constraints, solids, and boundaries

   double s0, s2;
   Vec3d normal;

   check_edge_edge_proximity( m_surf.get_position(edge_data0[0]), 
      m_surf.get_position(edge_data0[1]), 
      m_surf.get_position(edge_data1[0]), 
      m_surf.get_position(edge_data1[1]), 
      current_length, s0, s2, normal );

   return current_length < m_surf.m_merge_proximity_epsilon;  
   
}


// --------------------------------------------------------
///
/// Determine if the face vertex pair should be snapped
///
// --------------------------------------------------------

bool MeshSnapper::face_vertex_pair_is_snappable( size_t face, size_t vertex, double& current_length )
{

   // skip deleted pairs
   if(m_surf.m_mesh.triangle_is_deleted(face) || m_surf.m_mesh.vertex_is_deleted(vertex) ) 
      return false;

   //face shouldn't contain the vertex
   const Vec3st& face_data = m_surf.m_mesh.m_tris[face];
   if(m_surf.m_mesh.triangle_contains_vertex(face_data, vertex))
      return false;
   
   //TODO extend to handle constraints, solids, and boundaries

   double s0, s1, s2;
   Vec3d normal;

   check_point_triangle_proximity( m_surf.get_position(vertex), 
      m_surf.get_position(face_data[0]), 
      m_surf.get_position(face_data[1]), 
      m_surf.get_position(face_data[2]), 
      current_length, s0, s1, s2, normal );

   return current_length < m_surf.m_merge_proximity_epsilon;  

}


bool MeshSnapper::snap_pass()
{

   if ( m_surf.m_verbose )
   {
      std::cout << "\n\n\n---------------------- MeshSnapper: collapsing ----------------------" << std::endl;
      std::cout << "m_merge_proximity_epsilon: " << m_surf.m_merge_proximity_epsilon;

   }

   bool snap_occurred = false;

   assert( m_surf.m_dirty_triangles.size() == 0 );

   std::vector<SortableProximity> sortable_pairs_to_try;


   //
   // get sets of geometry pairs to try snapping!
   //

   // first the face-vertex pairs
   for(size_t vertex = 0; vertex < m_surf.get_num_vertices(); ++vertex) {
      if(m_surf.m_mesh.vertex_is_deleted(vertex)) continue;

      Vec3d vmin, vmax;
      m_surf.vertex_static_bounds(vertex, vmin, vmax);
      vmin -= m_surf.m_merge_proximity_epsilon * Vec3d(1,1,1);
      vmax += m_surf.m_merge_proximity_epsilon * Vec3d(1,1,1);

      std::vector<size_t> overlapping_tris;
      m_surf.m_broad_phase->get_potential_triangle_collisions(vmin, vmax, false, true, overlapping_tris);
   
      for(size_t i = 0; i < overlapping_tris.size(); ++i) {
         size_t face = overlapping_tris[i];
         Vec3st tri_data = m_surf.m_mesh.m_tris[face];
         
         double len;
         if(face_vertex_pair_is_snappable(face, vertex, len))
         {
            SortableProximity prox(face, vertex, len, true);
            sortable_pairs_to_try.push_back(prox);
         }
      }
   }

   //now the edge-edge pairs
   for(size_t edge0 = 0; edge0 < m_surf.m_mesh.m_edges.size(); ++edge0) {
      if(m_surf.m_mesh.edge_is_deleted(edge0)) continue;

      Vec3d vmin, vmax;
      m_surf.edge_static_bounds(edge0, vmin, vmax);
      vmin -= m_surf.m_merge_proximity_epsilon * Vec3d(1,1,1);
      vmax += m_surf.m_merge_proximity_epsilon * Vec3d(1,1,1);

      std::vector<size_t> overlapping_edges;
      m_surf.m_broad_phase->get_potential_edge_collisions(vmin, vmax, false, true, overlapping_edges);

      const Vec2st& edge_data0 = m_surf.m_mesh.m_edges[edge0];

      for(size_t ind = 0; ind < overlapping_edges.size(); ++ind) {
         size_t edge1 = overlapping_edges[ind];
         
         double len;
         if (edge_pair_is_snappable(edge0, edge1, len))
         {
            SortableProximity prox(edge0, edge1, len, false);
            sortable_pairs_to_try.push_back(prox);
         }

      }
   }
   
   //
   // sort in ascending order by distance (prefer to merge nearby geometry first)
   //
   // TODO: can we update the local geometry after merges to find newly formed proximities
   // and then using a priority queue, make sure to test them immediately if they're closer?

   std::sort( sortable_pairs_to_try.begin(), sortable_pairs_to_try.end() );

   //if ( m_surf.m_verbose )
   //{
      std::cout << sortable_pairs_to_try.size() << " candidate pairs sorted" << std::endl;
   //}

   //
   // attempt to split and snap each pair in the sorted list
   //

   for ( size_t si = 0; si < sortable_pairs_to_try.size(); ++si )
   {
      std::cout << "\n\nConsidering pair " << si << " with proximity " << sortable_pairs_to_try[si].m_length << std::endl; 
      size_t ind0 = sortable_pairs_to_try[si].m_index0;
      size_t ind1 = sortable_pairs_to_try[si].m_index1;
      std::cout << "Data " << ind0 << " and " << ind1 << std::endl; 
      
      bool result = false;
      bool attempted = false;
      if(sortable_pairs_to_try[si].m_face_vert_proximity) {
         //perform face-vertex split-n-snap
         size_t face = ind0;
         size_t vertex = ind1;
         
         double cur_len;
         if(face_vertex_pair_is_snappable(face, vertex, cur_len)) {
            std::cout << "Attempting face-vertex snap.\n";
            result = snap_face_vertex_pair(face, vertex);
            attempted = true;
         }

      }
      else {
         //perform edge-edge split-n-snap
         size_t edge0 = ind0;
         size_t edge1 = ind1;

         double cur_len;
         if(edge_pair_is_snappable(edge0, edge1, cur_len)) {
            std::cout << "Attempting edge-edge snap.\n";
            result = snap_edge_pair(edge0, edge1);
            attempted = true;
         }
      }
      
      if ( result )
      { 
         std::cout << "Snapping succeeded. Cleaning up.\n\n\n";
         // clean up degenerate triangles and tets
         m_surf.trim_non_manifold( m_surf.m_dirty_triangles );            
      }
      else if(attempted) {
         std::cout << "Snapping failed.\n";
      }
      else {
         std::cout << "Pair no longer snappable.\n";
      }
      

      snap_occurred |= result;
      
   }

   return snap_occurred;

}


}
