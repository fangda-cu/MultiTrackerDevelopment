// ---------------------------------------------------------
//
//  edgecollapser.cpp
//  Tyson Brochu 2011
//  
//  Functions supporting the "edge collapse" operation: removing short edges from the mesh.
//
// ---------------------------------------------------------

#include <edgecollapser.h>

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
/// Edge collapser constructor.  Takes a SurfTrack object and curvature-adaptive parameters.
///
// --------------------------------------------------------

EdgeCollapser::EdgeCollapser( SurfTrack& surf, bool use_curvature, bool remesh_boundaries, double min_curvature_multiplier ) :
m_min_edge_length( UNINITIALIZED_DOUBLE ),
m_max_edge_length( UNINITIALIZED_DOUBLE ),
m_use_curvature( use_curvature ),
m_remesh_boundaries( remesh_boundaries ),
m_min_curvature_multiplier( min_curvature_multiplier ),
m_rank_region(-1),
m_surf( surf )
{}


// --------------------------------------------------------
///
/// Get all triangles which are incident on either edge end vertex.
///
// --------------------------------------------------------

void EdgeCollapser::get_moving_triangles(size_t source_vertex, 
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
/// Get all edges which are incident on either edge end vertex.
///
// --------------------------------------------------------

void EdgeCollapser::get_moving_edges( size_t source_vertex, 
                                     size_t destination_vertex, 
                                     size_t,
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

bool EdgeCollapser::collapse_edge_pseudo_motion_introduces_collision( size_t source_vertex, 
                                                                     size_t destination_vertex, 
                                                                     size_t edge_index, 
                                                                     const Vec3d& )
{
    assert( m_surf.m_collision_safety );
    
    // Get the set of triangles which move because of this motion
    std::vector<size_t> moving_triangles;
    get_moving_triangles( source_vertex, destination_vertex, moving_triangles );
    
    // And the set of edges
    std::vector<size_t> moving_edges;
    get_moving_edges( source_vertex, destination_vertex, edge_index, moving_edges );
    
    
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


// --------------------------------------------------------
///
/// Determine if the edge collapse operation would invert the normal of any incident triangles.
///
// --------------------------------------------------------

bool EdgeCollapser::collapse_edge_introduces_normal_inversion( size_t source_vertex, 
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
  
    double min_triangle_area = -1;
    for (size_t i = 0; i < moving_triangles.size(); i++)
    {
        Vec3st current_triangle = m_surf.m_mesh.get_triangle(moving_triangles[i]);
        double area = triangle_area(m_surf.get_position(current_triangle[0]), m_surf.get_position(current_triangle[1]), m_surf.get_position(current_triangle[2]));
        if (min_triangle_area < 0 || area < min_triangle_area)
            min_triangle_area = area;
    }
    assert(min_triangle_area > 0);
    min_triangle_area = std::min(min_triangle_area, m_surf.m_min_triangle_area);
    
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
            
            g_stats.add_to_int( "EdgeCollapser:collapse_normal_inversion", 1 );
            return true;
        } 
        
        if ( new_area < min_triangle_area )
        {
            if ( m_surf.m_verbose ) { std::cout << "collapse edge introduces tiny triangle area" << std::endl; }
            
            g_stats.add_to_int( "EdgeCollapser:collapse_degenerate_triangle", 1 );
            return true;
        } 
        
    }
    
    return false;
    
}


// --------------------------------------------------------
///
/// Determine whether collapsing an edge will introduce an unacceptable change in volume.
///
// --------------------------------------------------------

bool EdgeCollapser::collapse_edge_introduces_volume_change( size_t source_vertex, 
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


// ---------------------------------------------------------
///
/// Returns true if the edge collapse would introduce a triangle with a min or max angle outside of the specified min or max.
///
// ---------------------------------------------------------

bool EdgeCollapser::collapse_edge_introduces_bad_angle(size_t source_vertex, 
                                                       size_t destination_vertex, 
                                                       const Vec3d& vertex_new_position )
{
    
    std::vector<size_t> moving_triangles;
    get_moving_triangles( source_vertex, destination_vertex,  moving_triangles );
    
    double min_tri_angle = -1;
    double max_tri_angle = -1;
    for ( size_t i = 0; i < moving_triangles.size(); ++i )
    {
        const Vec3st& tri = m_surf.m_mesh.get_triangle( moving_triangles[i] );
        double mina = min_triangle_angle(m_surf.get_position(tri[0]), m_surf.get_position(tri[1]), m_surf.get_position(tri[2]));
        double maxa = max_triangle_angle(m_surf.get_position(tri[0]), m_surf.get_position(tri[1]), m_surf.get_position(tri[2]));
        assert(mina > 0);
        assert(maxa > mina);
        assert(M_PI > maxa);
        if (min_tri_angle < 0 || mina < min_tri_angle)
            min_tri_angle = mina;
        if (max_tri_angle < 0 || maxa > max_tri_angle)
            max_tri_angle = maxa;
    }
    assert(min_tri_angle > 0);
    assert(max_tri_angle > min_tri_angle);
    assert(M_PI > max_tri_angle);
    
    min_tri_angle = std::min(min_tri_angle, deg2rad(m_surf.m_min_triangle_angle));
    max_tri_angle = std::max(max_tri_angle, deg2rad(m_surf.m_max_triangle_angle));
    
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
        
        if ( min_angle < min_tri_angle )
        {
            return true;
        }
        
        double max_angle = max_triangle_angle( a, b, c );
        
        if ( max_angle > max_tri_angle )
        {
            return true;
        }
        
    }
    
    return false;
    
}





// --------------------------------------------------------
///
/// Choose the vertex to keep and delete, and the remaining vertices' position.
/// Return false if the edge turns out not to be collapsible
///

// --------------------------------------------------------
bool EdgeCollapser::get_new_vertex_position_dihedral(Vec3d& vertex_new_position, size_t& vertex_to_keep, size_t& vertex_to_delete, const size_t& edge, bool& new_vert_constraint_label) {


   // rank 1, 2, 3 = smooth, ridge, peak
   // if the vertex ranks don't match, keep the higher rank vertex

   int keep_incident_features = m_surf.vertex_feature_edge_count(vertex_to_keep);
   int delete_incident_features = m_surf.vertex_feature_edge_count(vertex_to_delete);

   //For now, map these to the existing rank scores
   unsigned int keep_rank = keep_incident_features < 2 ? 1 : (keep_incident_features >= 3 ? 3 : 2);
   unsigned int delete_rank = delete_incident_features < 2 ? 1 : (delete_incident_features >= 3 ? 3 : 2);
   
   bool edge_is_a_feature = m_surf.edge_is_feature(edge);
     
   double edge_length = mag(m_surf.get_position(vertex_to_keep) - m_surf.get_position(vertex_to_delete));
  
   bool keep_vert_is_boundary = m_surf.m_mesh.m_is_boundary_vertex[vertex_to_keep];
   bool del_vert_is_boundary = m_surf.m_mesh.m_is_boundary_vertex[vertex_to_delete];

   // boundary vertices have precedence
   if (keep_vert_is_boundary) keep_rank = 4;
   if (del_vert_is_boundary) delete_rank = 4;

   // constraint vertices have higher precedence
   bool keep_vert_is_constrained =   m_surf.m_mesh.get_vertex_constraint_label(vertex_to_keep);
   bool delete_vert_is_constrained = m_surf.m_mesh.get_vertex_constraint_label(vertex_to_delete);

   bool keep_vert_is_manifold = !m_surf.m_mesh.is_vertex_nonmanifold(vertex_to_keep);
   bool delete_vert_is_manifold = !m_surf.m_mesh.is_vertex_nonmanifold(vertex_to_delete);

   new_vert_constraint_label = false;
   if (keep_vert_is_constrained || delete_vert_is_constrained)
   {
      assert(m_surf.m_constrained_vertices_callback);
      new_vert_constraint_label = m_surf.m_constrained_vertices_callback->generate_collapsed_constraint_label(m_surf, vertex_to_keep, vertex_to_delete, m_surf.m_mesh.get_vertex_constraint_label(vertex_to_keep), m_surf.m_mesh.get_vertex_constraint_label(vertex_to_delete));
   }

   if (keep_vert_is_constrained)   keep_rank = 5;
   if (delete_vert_is_constrained) delete_rank = 5;

   // Handle different cases of constrained, boundary and interior vertices
   if (m_surf.m_allow_vertex_movement_during_collapse && !(keep_vert_is_boundary || del_vert_is_boundary) && !(keep_vert_is_constrained || delete_vert_is_constrained))
   {
      //Ranks dominate (i.e. use ranks to decide collapsing first, and if they match then use nonmanifoldness to decide).
      //-> This is particularly important for outward normal flow: it snaps the non-manifold curve back onto the 
      //feature curve produced at merge points. (Other scenarios might(?) work better with nonmanifoldness dominating; so I've left
      //the option in the code for now - Christopher Batty.)
      bool ranks_dominate = true; 
      if(ranks_dominate) {
         if ( keep_rank > delete_rank ) {
            vertex_new_position = m_surf.get_position(vertex_to_keep);
         }
         else if ( delete_rank > keep_rank ) {
            std::swap(vertex_to_keep, vertex_to_delete);
            vertex_new_position = m_surf.get_position(vertex_to_keep);
         }
         else
         {
            //same ranks, but one is non-manifold; may as well prefer to keep non-manifold points.
            if(!keep_vert_is_manifold && delete_vert_is_manifold) {
               vertex_new_position = m_surf.get_position(vertex_to_keep);
            }
            else if(keep_vert_is_manifold && !delete_vert_is_manifold) {
               std::swap(vertex_to_keep, vertex_to_delete);
               vertex_new_position = m_surf.get_position(vertex_to_keep);
            }
            else {
               // ranks are equal and manifoldness matches too
               m_surf.m_subdivision_scheme->generate_new_midpoint( edge, m_surf, vertex_new_position );
            }
         }
      }
      else {
         //Manifoldness dominates
         if(!keep_vert_is_manifold && delete_vert_is_manifold) {
            vertex_new_position = m_surf.get_position(vertex_to_keep);
         }
         else if(!delete_vert_is_manifold && keep_vert_is_manifold) {
            std::swap(vertex_to_keep, vertex_to_delete);
            vertex_new_position = m_surf.get_position(vertex_to_keep);
         }
         else {
            if ( keep_rank > delete_rank ) {
               vertex_new_position = m_surf.get_position(vertex_to_keep);
            }
            else if ( delete_rank > keep_rank ) {
               std::swap(vertex_to_keep, vertex_to_delete);
               vertex_new_position = m_surf.get_position(vertex_to_keep);
            }
            else {
               // ranks are equal and manifoldness matches too
               m_surf.m_subdivision_scheme->generate_new_midpoint( edge, m_surf, vertex_new_position );
            }
         }
      }

     if ((keep_rank >= 3 || !keep_vert_is_manifold) && (delete_rank >= 3 || !delete_vert_is_manifold) && edge_length >= m_t1_pull_apart_distance)
       return false;
   }
   else if (keep_vert_is_constrained || delete_vert_is_constrained)
   {
      assert(m_surf.m_constrained_vertices_callback);

      Vec3d newpos = (m_surf.get_position(vertex_to_keep) + m_surf.get_position(vertex_to_delete)) / 2;
      if (!m_surf.m_constrained_vertices_callback->generate_collapsed_position(m_surf, vertex_to_keep, vertex_to_delete, newpos))
      {
         // the callback decides this edge should not be collapsed
         if (m_surf.m_verbose)
            std::cout << "Constraint callback vetoed collapsing." << std::endl;
         return false;
      }
     
     if ((keep_rank >= 3 || !keep_vert_is_manifold || keep_vert_is_constrained) && (delete_rank >= 3 || !delete_vert_is_manifold || delete_vert_is_constrained) && edge_length >= m_t1_pull_apart_distance)
       return false;

      vertex_new_position = newpos;
   } else if (keep_vert_is_boundary || del_vert_is_boundary)
   {
      if (!keep_vert_is_boundary)
      {
         std::swap(keep_vert_is_boundary, del_vert_is_boundary);
         std::swap(vertex_to_keep, vertex_to_delete);
         std::swap(keep_rank, delete_rank);
      }

      vertex_new_position = m_surf.get_position(vertex_to_keep);
   } else  // m_surf.m_allow_vertex_movement_during_collapse == false
   {
      if (keep_rank < delete_rank)
      {
         std::swap(vertex_to_keep, vertex_to_delete);
         std::swap(keep_rank, delete_rank);
      }

      vertex_new_position = m_surf.get_position(vertex_to_keep);
   }

   return true;
}


// --------------------------------------------------------
///
/// Delete an edge by moving its source vertex to its destination vertex
///
// --------------------------------------------------------

bool EdgeCollapser::collapse_edge( size_t edge )
{
  size_t vertex_to_keep = m_surf.m_mesh.m_edges[edge][0];
  size_t vertex_to_delete = m_surf.m_mesh.m_edges[edge][1];
  
  bool keep_vert_is_boundary = m_surf.m_mesh.m_is_boundary_vertex[vertex_to_keep];
  bool del_vert_is_boundary = m_surf.m_mesh.m_is_boundary_vertex[vertex_to_delete];
  bool edge_is_boundary = m_surf.m_mesh.m_is_boundary_edge[edge];
  
  //either we're allowing remeshing of boundary edges, or this edge is not on the boundary.
  assert(m_remesh_boundaries || !edge_is_boundary);
  
  ///////////////////////////////////////////////////////////////////////
  // FD 20130102
  //
  // do not collapse the edge if the two vertices are moving apart
  
  Vec3d edge_vec = m_surf.get_position(vertex_to_keep) - m_surf.get_position(vertex_to_delete);
  Vec3d rel_vel = m_surf.get_remesh_velocity(vertex_to_keep) - m_surf.get_remesh_velocity(vertex_to_delete);
  double edge_len = mag(edge_vec);
  
  if ((dot(rel_vel, edge_vec) > 0 || edge_len >= m_t1_pull_apart_distance) && collapse_will_produce_irregular_junction(edge))
  {
    if (m_surf.m_verbose)
      std::cout << "The collapse will produce irregular junction, but the endpoints are moving apart. No need to collapse." << std::endl;
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
  
  if ( m_surf.m_verbose ) { std::cout << "Collapsing edge.  Doomed vertex: " << vertex_to_delete << " --- Vertex to keep: " << vertex_to_keep << std::endl; }

  // --------------

  // If we're disallowing topology changes, don't let an edge collapse form a degenerate tet

  if ( false == m_surf.m_allow_non_manifold )
  {

    bool would_be_non_manifold = false;

    // for each triangle that *would* be created, make sure that there isn't already a triangle with those 3 vertices

    for ( size_t i = 0; i < m_surf.m_mesh.m_vertex_to_triangle_map[vertex_to_delete].size(); ++i )
    {
      Vec3st new_triangle = m_surf.m_mesh.get_triangle( m_surf.m_mesh.m_vertex_to_triangle_map[vertex_to_delete][i] );
      if ( new_triangle[0] == vertex_to_delete )   { new_triangle[0] = vertex_to_keep; }
      if ( new_triangle[1] == vertex_to_delete )   { new_triangle[1] = vertex_to_keep; }
      if ( new_triangle[2] == vertex_to_delete )   { new_triangle[2] = vertex_to_keep; }

      for ( size_t j = 0; j < m_surf.m_mesh.m_vertex_to_triangle_map[vertex_to_keep].size(); ++j )
      {
        if ( NonDestructiveTriMesh::triangle_has_these_verts( m_surf.m_mesh.get_triangle( m_surf.m_mesh.m_vertex_to_triangle_map[vertex_to_keep][j] ), new_triangle ) )
        {
          if ( m_surf.m_verbose ) { std::cout << "would_be_non_manifold" << std::endl; }
          would_be_non_manifold = true;
          return false;
        }            
      }
    }

    assert ( !would_be_non_manifold );

    // Look for a vertex which is adjacent to both vertices on the edge, and which isn't on one of the incident triangles

    const std::vector< size_t >& triangles_incident_to_edge = m_surf.m_mesh.m_edge_to_triangle_map[edge];
    std::vector< size_t > third_vertices;

    for ( size_t i = 0; i < triangles_incident_to_edge.size(); ++i )
    {
      const Vec3st& inc_triangle = m_surf.m_mesh.get_triangle( triangles_incident_to_edge[i] );
      size_t opposite = m_surf.m_mesh.get_third_vertex( edge, inc_triangle );
      third_vertices.push_back( opposite );
    }

    std::vector<size_t> adj_vertices0, adj_vertices1;
    m_surf.m_mesh.get_adjacent_vertices( vertex_to_delete, adj_vertices0 );
    m_surf.m_mesh.get_adjacent_vertices( vertex_to_keep, adj_vertices1 );

    for ( size_t i = 0; i < adj_vertices0.size(); ++i )
    {
      for ( size_t j = 0; j < adj_vertices1.size(); ++j )
      {
        if ( adj_vertices0[i] == adj_vertices1[j] )
        {
          bool is_on_inc_triangle = false;
          for ( size_t k = 0; k < third_vertices.size(); ++k )
          {
            if ( adj_vertices0[i] == third_vertices[k] )
            {
              is_on_inc_triangle = true;
              break;
            }
          }

          if ( !is_on_inc_triangle )
          {
            // found a vertex adjacent to both edge vertices, which doesn't lie on the incident triangles

            //
            if ( m_surf.m_verbose )
            {
              std::cout << " --- Edge Collapser: found a vertex adjacent to both edge vertices, which doesn't lie on the incident triangles " << std::endl;
              std::cout << " --- Adjacent vertex: " << adj_vertices0[i] << ", incident triangles: ";
            }

            return false;
          }

        }
      }
    }


  }


  // --------------

  {
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
          g_stats.add_to_int( "EdgeCollapser:collapse_degen_tet", 1 );
          return false;
        }
      }
    }
  }


  // --------------
  // decide on new vertex position

  //Choose the vertex to keep and its new position.
  Vec3d vertex_new_position;
  bool new_vert_constraint_label;
  bool can_collapse = get_new_vertex_position_dihedral(vertex_new_position, vertex_to_keep, vertex_to_delete, edge, new_vert_constraint_label);
  if(!can_collapse)
     return false;

  if ( m_surf.m_verbose ) { std::cout << "Collapsing edge.  Doomed vertex: " << vertex_to_delete << " --- Vertex to keep: " << vertex_to_keep << std::endl; }

  // --------------

  // Check vertex pseudo motion for collisions and volume change

  if ( mag ( m_surf.get_position(m_surf.m_mesh.m_edges[edge][1]) - m_surf.get_position(m_surf.m_mesh.m_edges[edge][0]) ) > 0 )
  {

    // Change source vertex predicted position to superimpose onto destination vertex
    m_surf.set_newposition( vertex_to_keep, vertex_new_position );
    m_surf.set_newposition( vertex_to_delete, vertex_new_position );

    bool volume_change = collapse_edge_introduces_volume_change( vertex_to_delete, edge, vertex_new_position );

    if ( volume_change )
    {
      // Restore saved positions which were changed by the function we just called.
      m_surf.set_newposition( vertex_to_keep, m_surf.get_position(vertex_to_keep) );
      m_surf.set_newposition( vertex_to_delete, m_surf.get_position(vertex_to_delete) );

      g_stats.add_to_int( "EdgeCollapser:collapse_volume_change", 1 );

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

      g_stats.add_to_int( "EdgeCollapser:collapse_bad_angle", 1 );
      return false;

    }

    bool collision = false;

    if ( m_surf.m_collision_safety )
    {
      collision = collapse_edge_pseudo_motion_introduces_collision( vertex_to_delete, vertex_to_keep, edge, vertex_new_position );
    }

    if ( collision ) 
    { 
      if ( m_surf.m_verbose ) { std::cout << "collision" << std::endl; }
      g_stats.add_to_int( "EdgeCollapser:collapse_collisions", 1 ); 
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
  MeshUpdateEvent collapse(MeshUpdateEvent::EDGE_COLLAPSE);
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
  m_surf.m_mesh.set_vertex_constraint_label(vertex_to_keep, new_vert_constraint_label);
  
  ///////////////////////////////////////////////////////////////////////

  // Copy this vector, don't take a reference, as deleting will change the original
  std::vector< size_t > triangles_incident_to_edge = m_surf.m_mesh.m_edge_to_triangle_map[edge];

  // Delete triangles incident on the edge

  for ( size_t i=0; i < triangles_incident_to_edge.size(); ++i )
  {
    if ( m_surf.m_verbose )
    {
      std::cout << "removing edge-incident triangle: " << m_surf.m_mesh.get_triangle( triangles_incident_to_edge[i] ) << std::endl;
    }

    m_surf.remove_triangle( triangles_incident_to_edge[i] );
    collapse.m_deleted_tris.push_back(triangles_incident_to_edge[i] );
  }

  // Find anything pointing to the doomed vertex and change it

  // copy the list of triangles, don't take a reference to it
  std::vector< size_t > triangles_incident_to_vertex = m_surf.m_mesh.m_vertex_to_triangle_map[vertex_to_delete];    

  for ( size_t i=0; i < triangles_incident_to_vertex.size(); ++i )
  {
    for(size_t local_ind = 0; local_ind < triangles_incident_to_edge.size(); ++local_ind) {
      assert( triangles_incident_to_vertex[i] != triangles_incident_to_edge[local_ind] );
    }

    Vec3st new_triangle = m_surf.m_mesh.get_triangle( triangles_incident_to_vertex[i] );

    if ( new_triangle[0] == vertex_to_delete )   { new_triangle[0] = vertex_to_keep; }
    if ( new_triangle[1] == vertex_to_delete )   { new_triangle[1] = vertex_to_keep; }
    if ( new_triangle[2] == vertex_to_delete )   { new_triangle[2] = vertex_to_keep; }

    if ( m_surf.m_verbose ) { std::cout << "adding updated triangle: " << new_triangle << std::endl; }

    // the old label carries over to the new triangle.
    // no need to test for orientation because the new triangle
    // generation code above does not change orientation.
    //
    Vec2i label = m_surf.m_mesh.get_triangle_label(triangles_incident_to_vertex[i]);

    size_t new_triangle_index = m_surf.add_triangle( new_triangle, label );
    collapse.m_created_tris.push_back( new_triangle_index );
    collapse.m_created_tri_data.push_back(new_triangle);
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
  
  if (m_surf.m_mesheventcallback)
    m_surf.m_mesheventcallback->collapse(m_surf, edge);

  return true;
}

// --------------------------------------------------------
///
/// Determine if the edge should be allowed to collapse
///
// --------------------------------------------------------

bool EdgeCollapser::edge_is_collapsible( size_t edge_index, double& current_length )
{

  // skip deleted and solid edges
  if ( m_surf.m_mesh.edge_is_deleted(edge_index) ) { return false; }
//  if ( m_surf.edge_is_any_solid(edge_index) ) { return false; }
  
  //skip boundary edges if we're not remeshing those
  if(!m_remesh_boundaries && m_surf.m_mesh.m_is_boundary_edge[edge_index]) { return false; }

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
  
  current_length = m_surf.get_edge_length(edge_index);
  if ( m_use_curvature )
  {

    //collapse if we're below the lower limit
    if(current_length < m_min_edge_length) 
        return true;

    //don't collapse if we're near the upper limit, since it can produce edges above the limit
    if(current_length > m_max_edge_length*0.5)
        return false;

    //check all incident edges to see if any of them are super short, and if so, split this guy accordingly.
    //this enforces slow grading of the mesh.
    double min_nbr_len = current_length;
    size_t vertex_a = m_surf.m_mesh.m_edges[edge_index][0];
    size_t vertex_b = m_surf.m_mesh.m_edges[edge_index][1];
    for(size_t edge_id = 0; edge_id < m_surf.m_mesh.m_vertex_to_edge_map[vertex_a].size(); ++edge_id) {
        min_nbr_len = min(min_nbr_len, m_surf.get_edge_length(m_surf.m_mesh.m_vertex_to_edge_map[vertex_a][edge_id]));
    }
    for(size_t edge_id = 0; edge_id < m_surf.m_mesh.m_vertex_to_edge_map[vertex_b].size(); ++edge_id) {
        min_nbr_len = min(min_nbr_len, m_surf.get_edge_length(m_surf.m_mesh.m_vertex_to_edge_map[vertex_b][edge_id]));
    }

    if(current_length < 3*min_nbr_len )
        return false;

    double curvature_value = get_edge_curvature( m_surf, m_surf.m_mesh.m_edges[edge_index][0], m_surf.m_mesh.m_edges[edge_index][1] );
    
    //Assume we want to discretize any circle with at least ten segments.
    //Then give the curvature (i.e. inverse radius) the target edge length
    //here should be computed as... 
    int circlesegs = 32;
    double curvature_min_length = 2*M_PI / (double)circlesegs / max(curvature_value, 1e-8);

    //collapse if curvature dictates we should. 
    return current_length < curvature_min_length;
    
  }
  else
  {
    return current_length < m_min_edge_length;  
  }


}


// --------------------------------------------------------
///
/// Collapse all short edges
///
// --------------------------------------------------------

bool EdgeCollapser::collapse_pass()
{
    
    if ( m_surf.m_verbose )
    {
        std::cout << "\n\n\n---------------------- EdgeCollapser: collapsing ----------------------" << std::endl;
        std::cout << "m_min_edge_length: " << m_min_edge_length;
        std::cout << ", m_use_curvature: " << m_use_curvature;
        std::cout << ", m_min_curvature_multiplier: " << m_min_curvature_multiplier << std::endl;
        std::cout << "m_surf.m_collision_safety: " << m_surf.m_collision_safety << std::endl;
    }
    
    bool collapse_occurred = false;
    
    assert( m_surf.m_dirty_triangles.size() == 0 );
    
    std::vector<SortableEdge> sortable_edges_to_try;
    
    
    //
    // get set of edges to collapse
    //
    
    for( size_t i = 0; i < m_surf.m_mesh.m_edges.size(); i++ )
    {    
        double current_length;
        if(edge_is_collapsible(i, current_length)) 
          sortable_edges_to_try.push_back( SortableEdge( i, current_length ) );
    }
    
    //
    // sort in ascending order by length (collapse shortest edges first)
    //
    
    std::sort( sortable_edges_to_try.begin(), sortable_edges_to_try.end() );
    
    if ( m_surf.m_verbose )
    {
        std::cout << sortable_edges_to_try.size() << " candidate edges sorted" << std::endl;
        std::cout << "total edges: " << m_surf.m_mesh.m_edges.size() << std::endl;
    }
    
    //
    // attempt to collapse each edge in the sorted list
    //
    
    for ( size_t si = 0; si < sortable_edges_to_try.size(); ++si )
    {
        size_t e = sortable_edges_to_try[si].m_edge_index;
        
        assert( e < m_surf.m_mesh.m_edges.size() );
        
        double dummy;
        if(edge_is_collapsible(e, dummy)) {
          bool result = collapse_edge( e );

          if ( result )
          { 
            // clean up degenerate triangles and tets
            m_surf.trim_degeneracies( m_surf.m_dirty_triangles );                        
          }

          collapse_occurred |= result;
        }
    }
    
    return collapse_occurred;
    
}

bool EdgeCollapser::collapse_will_produce_irregular_junction(size_t edge)
{
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;
    
    size_t a = mesh.m_edges[edge][0];
    size_t b = mesh.m_edges[edge][1];
    
    std::vector<bool> regions;
    for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[a].size(); i++)
    {
        const Vec2i & l = mesh.get_triangle_label(mesh.m_vertex_to_triangle_map[a][i]);
        if (l[0] >= 0)
        {
            if (l[0] >= (int)regions.size()) regions.resize(l[0] + 1, false);
            regions[l[0]] = true;
        }
        if (l[1] >= 0)
        {
            if (l[1] >= (int)regions.size()) regions.resize(l[1] + 1, false);
            regions[l[1]] = true;
        }
    }
    for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[b].size(); i++)
    {
        const Vec2i & l = mesh.get_triangle_label(mesh.m_vertex_to_triangle_map[b][i]);
        if (l[0] >= 0)
        {
            if (l[0] >= (int)regions.size()) regions.resize(l[0] + 1, false);
            regions[l[0]] = true;
        }
        if (l[1] >= 0)
        {
            if (l[1] >= (int)regions.size()) regions.resize(l[1] + 1, false);
            regions[l[1]] = true;
        }
    }
    
    int nr = regions.size();
    bool * regiongraph = new bool[nr * nr];
    memset(regiongraph, 0, sizeof(bool) * nr * nr);
    for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[a].size(); i++)
    {
        const Vec2i & l = mesh.get_triangle_label(mesh.m_vertex_to_triangle_map[a][i]);
        if (l[0] >= 0 && l[1] >= 0)
            regiongraph[l[0] * nr + l[1]] = regiongraph[l[1] * nr + l[0]] = true;
    }
    for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[b].size(); i++)
    {
        const Vec2i & l = mesh.get_triangle_label(mesh.m_vertex_to_triangle_map[b][i]);
        if (l[0] >= 0 && l[1] >= 0)
            regiongraph[l[0] * nr + l[1]] = regiongraph[l[1] * nr + l[0]] = true;
    }

    bool irregular = false;
    for (int i = 0; i < nr; i++)
        for (int j = i + 1; j < nr; j++)
        {
            if (regions[i] && regions[j])
                if (!regiongraph[i * nr + j])
                    irregular = true;
        }
    
    delete []regiongraph;
    
    return irregular;
}
    
}

