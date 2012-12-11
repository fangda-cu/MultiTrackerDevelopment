// ---------------------------------------------------------
//
//  edgesplitter.cpp
//  Tyson Brochu 2011
//  
//  Functions supporting the "edge split" operation: subdividing an edge into two shorter edges.
//
// ---------------------------------------------------------

#include <edgesplitter.h>
#include <broadphase.h>
#include <collisionqueries.h>
#include <runstats.h>
#include <subdivisionscheme.h>
#include <surftrack.h>
#include <trianglequality.h>
#include <typeinfo>


// ---------------------------------------------------------
//  Extern globals
// ---------------------------------------------------------

namespace ElTopo {

extern RunStats g_stats;

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Constructor.  Active SurfTrack object must be supplied.
///
// ---------------------------------------------------------

EdgeSplitter::EdgeSplitter( SurfTrack& surf, bool use_curvature, bool remesh_boundaries, double max_curvature_multiplier ) :
m_max_edge_length( UNINITIALIZED_DOUBLE ),
m_min_edge_length( UNINITIALIZED_DOUBLE ),
m_use_curvature( use_curvature ),
m_remesh_boundaries( remesh_boundaries),
m_max_curvature_multiplier( max_curvature_multiplier ),
m_surf( surf )
{}


// --------------------------------------------------------
///
/// Check collisions between the edge [neighbour, new] and the given edge 
///
// --------------------------------------------------------

bool EdgeSplitter::split_edge_edge_collision( size_t neighbour_index, 
                                             const Vec3d& new_vertex_position, 
                                             const Vec3d& new_vertex_smooth_position, 
                                             const Vec2st& edge )
{
    
    size_t edge_vertex_0 = edge[0];
    size_t edge_vertex_1 = edge[1];
    size_t dummy_index = m_surf.get_num_vertices();
    
    if ( neighbour_index == edge_vertex_0 || neighbour_index == edge_vertex_1 )  { return false; }
    
    const std::vector<Vec3d>& x = m_surf.get_positions();
    
    double t_zero_distance; 
    check_edge_edge_proximity( new_vertex_position, 
                              x[ neighbour_index ], 
                              x[ edge_vertex_0 ], 
                              x[ edge_vertex_1 ],
                              t_zero_distance );
    
    if ( t_zero_distance < m_surf.m_improve_collision_epsilon )
    {
        return true;
    }
    
    if ( edge_vertex_1 < edge_vertex_0 ) { swap( edge_vertex_0, edge_vertex_1 ); }
    
    if ( segment_segment_collision(x[ neighbour_index ], x[ neighbour_index ], neighbour_index,
                                   new_vertex_position, new_vertex_smooth_position, dummy_index,
                                   x[ edge_vertex_0 ], x[ edge_vertex_0 ], edge_vertex_0,
                                   x[ edge_vertex_1 ], x[ edge_vertex_1 ], edge_vertex_1 ) )
        
    {      
        return true;
    }
    
    return false;
    
}


// ---------------------------------------------------------
///
/// Determine if the new vertex introduced by the edge split has a collision along its pseudo-trajectory.
///
// ---------------------------------------------------------

bool EdgeSplitter::split_triangle_vertex_collision( const Vec3st& triangle_indices, 
                                                   const Vec3d& new_vertex_position, 
                                                   const Vec3d& new_vertex_smooth_position, 
                                                   size_t overlapping_vert_index, 
                                                   const Vec3d& vert )
{
    
    if ( overlapping_vert_index == triangle_indices[0] || overlapping_vert_index == triangle_indices[1] || overlapping_vert_index == triangle_indices[2] )
    {
        return false;
    }
    
    Vec3st sorted_triangle = sort_triangle( triangle_indices );
    
    Vec3d tri_positions[3];
    Vec3d tri_smooth_positions[3];
    
    for ( unsigned int i = 0; i < 3; ++i )
    {
        if ( sorted_triangle[i] == m_surf.get_num_vertices() )
        {
            tri_positions[i] = new_vertex_position;
            tri_smooth_positions[i] = new_vertex_smooth_position;
        }
        else
        {
            tri_positions[i] = m_surf.get_position( sorted_triangle[i] );
            tri_smooth_positions[i] = m_surf.get_position( sorted_triangle[i] );
        }
    }
    
    
    // check distance at time t=0
    double t_zero_distance;
    check_point_triangle_proximity( vert, tri_positions[0], tri_positions[1], tri_positions[2], t_zero_distance );
    
    
    if ( t_zero_distance < m_surf.m_improve_collision_epsilon )
    {
        return true;
    }
    
    
    // now check continuous collision
    
    if ( point_triangle_collision( vert, vert, overlapping_vert_index,
                                  tri_positions[0], tri_smooth_positions[0], sorted_triangle[0],
                                  tri_positions[1], tri_smooth_positions[1], sorted_triangle[1],
                                  tri_positions[2], tri_smooth_positions[2], sorted_triangle[2] ) )
    {         
        return true;
    }
    
    return false;
    
    
}



// ---------------------------------------------------------
///
/// Determine if the pseudo-trajectory of the new vertex has a collision with the existing mesh.
///
// ---------------------------------------------------------
bool EdgeSplitter::split_edge_pseudo_motion_introduces_intersection( const Vec3d& new_vertex_position, 
  const Vec3d& new_vertex_smooth_position, 
  size_t edge,
  size_t vertex_a,
  size_t vertex_b,
  std::vector<size_t>& tris, 
  std::vector<size_t>& verts)
{

  NonDestructiveTriMesh& m_mesh = m_surf.m_mesh;

  if ( !m_surf.m_collision_safety)
  {
    return false;
  }

  // 
  // new point vs all triangles
  // 

  {

    Vec3d aabb_low, aabb_high;
    minmax( new_vertex_position, new_vertex_smooth_position, aabb_low, aabb_high );

    aabb_low -= m_surf.m_aabb_padding * Vec3d(1,1,1);
    aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);

    std::vector<size_t> overlapping_triangles;
    m_surf.m_broad_phase->get_potential_triangle_collisions( aabb_low, aabb_high, true, true, overlapping_triangles );

    for ( size_t i = 0; i < overlapping_triangles.size(); ++i )
    {

      // Exclude all incident triangles from the check
      bool overlap = false;
      for(size_t j = 0; j < tris.size(); ++j) {
        if(overlapping_triangles[i] == tris[j]) {
          overlap = true;
          break;
        }
      }
      
      if( overlap )  
        continue;
      
      size_t triangle_vertex_0 = m_mesh.get_triangle( overlapping_triangles[i] )[0];
      size_t triangle_vertex_1 = m_mesh.get_triangle( overlapping_triangles[i] )[1];
      size_t triangle_vertex_2 = m_mesh.get_triangle( overlapping_triangles[i] )[2];

      double t_zero_distance;

      check_point_triangle_proximity( new_vertex_position, 
        m_surf.get_position( triangle_vertex_0 ),
        m_surf.get_position( triangle_vertex_1 ),
        m_surf.get_position( triangle_vertex_2 ),
        t_zero_distance );

      size_t dummy_index = m_surf.get_num_vertices();

      if ( t_zero_distance < m_surf.m_improve_collision_epsilon )
      {
        return true;
      }

      Vec3st sorted_triangle = sort_triangle( Vec3st( triangle_vertex_0, triangle_vertex_1, triangle_vertex_2 ) );


      if ( point_triangle_collision(  new_vertex_position, new_vertex_smooth_position, dummy_index,
        m_surf.get_position( sorted_triangle[0] ), m_surf.get_position( sorted_triangle[0] ), sorted_triangle[0],
        m_surf.get_position( sorted_triangle[1] ), m_surf.get_position( sorted_triangle[1] ), sorted_triangle[1],
        m_surf.get_position( sorted_triangle[2] ), m_surf.get_position( sorted_triangle[2] ), sorted_triangle[2] ) )

      {
        return true;
      }
    }

  }

  //
  // new edges vs all edges
  //

  {

    Vec3d edge_aabb_low, edge_aabb_high;

    // do one big query into the broad phase for all new edges
    minmax( new_vertex_position, new_vertex_smooth_position, 
      m_surf.get_position( vertex_a ), m_surf.get_position( vertex_b ),
      edge_aabb_low, edge_aabb_high );
    for(size_t i = 0; i < verts.size(); ++i)
      update_minmax(m_surf.get_position(verts[i]), edge_aabb_low, edge_aabb_high);

    edge_aabb_low -= m_surf.m_aabb_padding * Vec3d(1,1,1);
    edge_aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);

    std::vector<size_t> overlapping_edges;
    m_surf.m_broad_phase->get_potential_edge_collisions( edge_aabb_low, edge_aabb_high, true, true, overlapping_edges );

    std::vector<size_t> vertex_neighbourhood;
    vertex_neighbourhood.push_back(vertex_a); vertex_neighbourhood.push_back(vertex_b);
    vertex_neighbourhood.insert(vertex_neighbourhood.end(), verts.begin(), verts.end());

    for ( size_t i = 0; i < overlapping_edges.size(); ++i )
    {

      if ( overlapping_edges[i] == edge ) { continue; }
      if ( m_mesh.m_edges[ overlapping_edges[i] ][0] == m_mesh.m_edges[ overlapping_edges[i] ][1] ) { continue; }

      for ( size_t v = 0; v < vertex_neighbourhood.size(); ++v )
      {
        bool collision = split_edge_edge_collision( vertex_neighbourhood[v], 
          new_vertex_position, 
          new_vertex_smooth_position, 
          m_mesh.m_edges[overlapping_edges[i]] );

        if ( collision ) { return true; }
      }
    }      
  }

  //
  // new triangles vs all points
  //

  {
    Vec3d triangle_aabb_low, triangle_aabb_high;

    // do one big query into the broad phase for all new triangles
    minmax( new_vertex_position, new_vertex_smooth_position, 
      m_surf.get_position( vertex_a ), m_surf.get_position( vertex_b ),
      triangle_aabb_low, triangle_aabb_high );
    for(size_t i = 0; i < verts.size(); ++i)
      update_minmax(m_surf.get_position(verts[i]), triangle_aabb_low, triangle_aabb_high);

    triangle_aabb_low -= m_surf.m_aabb_padding * Vec3d(1,1,1);
    triangle_aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);

    std::vector<size_t> overlapping_vertices;
    m_surf.m_broad_phase->get_potential_vertex_collisions( triangle_aabb_low, triangle_aabb_high, true, true, overlapping_vertices );

    size_t dummy_e = m_surf.get_num_vertices();

    std::vector< Vec3st > triangle_indices;

    //TODO Hopefully order is irrelevant here...
    for( size_t i = 0; i < verts.size(); ++i) {
      triangle_indices.push_back( Vec3st( vertex_a, dummy_e, verts[i] ) );    // triangle aec      
      triangle_indices.push_back( Vec3st( verts[i], dummy_e, vertex_b ) );    // triangle ceb      
    }

    for ( size_t i = 0; i < overlapping_vertices.size(); ++i )
    {
      if ( m_mesh.m_vertex_to_triangle_map[overlapping_vertices[i]].empty() ) 
      { 
        continue; 
      }

      size_t overlapping_vert_index = overlapping_vertices[i];
      const Vec3d& vert = m_surf.get_position(overlapping_vert_index);

      for ( size_t j = 0; j < triangle_indices.size(); ++j )
      {
        bool collision = split_triangle_vertex_collision( triangle_indices[j], 
          new_vertex_position, 
          new_vertex_smooth_position, 
          overlapping_vert_index, 
          vert );

        if ( collision ) 
        { 
          return true; 
        }

      }
    }

  }

  return false;

}


// --------------------------------------------------------
///
/// Split an edge, using subdivision_scheme to determine the new vertex location, if safe to do so.
///
// --------------------------------------------------------

bool EdgeSplitter::split_edge( size_t edge )
{   

  g_stats.add_to_int( "EdgeSplitter:edge_split_attempts", 1 );

  assert( edge_is_splittable(edge) );

  NonDestructiveTriMesh& mesh = m_surf.m_mesh;

  // --------------
  // Collect all the triangles around the edge

  std::vector<size_t> incident_tris = mesh.m_edge_to_triangle_map[edge];
  std::vector<double> tri_areas;
  for(size_t i = 0; i < incident_tris.size(); ++i) {
    tri_areas.push_back(m_surf.get_triangle_area(incident_tris[i]));
    
    // Splitting degenerate triangles causes problems
    if(tri_areas[i] < m_surf.m_min_triangle_area) {
      g_stats.add_to_int( "EdgeSplitter:split_edge_incident_to_tiny_triangle", 1 );
      return false;
    }
  }
  

  // --------------

  // convert each incident triangle abc into a pair of triangles aec, ebc

  size_t vertex_a = mesh.m_edges[edge][0];
  size_t vertex_b = mesh.m_edges[edge][1];
  

  // Collect the incident verts
  std::vector<size_t> other_verts;
  
  for(size_t i = 0; i < incident_tris.size(); ++i) {
    size_t cur_tri = incident_tris[i];
    Vec3st tri_data = mesh.get_triangle(cur_tri);
    other_verts.push_back(mesh.get_third_vertex(vertex_a, vertex_b, tri_data));
  }



  // --------------

  // get edge midpoint and the point on the smooth surface

  Vec3d new_vertex_position = 0.5 * ( m_surf.get_position( vertex_a ) + m_surf.get_position( vertex_b ) );
  Vec3d new_vertex_smooth_position;

  bool use_smooth_point = (incident_tris.size() == 2) || (incident_tris.size() == 1 && typeid(m_surf.m_subdivision_scheme) == typeid(ModifiedButterflyScheme));

  // generate the new midpoint according to the subdivision scheme, only if manifold and not on a boundary
  if(use_smooth_point)
    m_surf.m_subdivision_scheme->generate_new_midpoint( edge, m_surf, new_vertex_smooth_position );
  else
    new_vertex_smooth_position = new_vertex_position;

  // --------------

  // check if the generated point introduces an intersection
  if(use_smooth_point) {
    
    use_smooth_point = use_smooth_point && ! ( split_edge_pseudo_motion_introduces_intersection( new_vertex_position, 
      new_vertex_smooth_position, 
      edge, 
      vertex_a, 
      vertex_b, 
      incident_tris, 
      other_verts) );
    
    if ( !use_smooth_point ) { 
        g_stats.add_to_int( "EdgeSplitter:split_smooth_vertex_collisions", 1 ); }
  }
    
  
  // --------------

  // check normal inversion

  if ( use_smooth_point )
  {
    size_t tri0 = incident_tris[0];
    size_t tri1 = incident_tris[1];
    size_t vertex_c = other_verts[0];
    size_t vertex_d = other_verts[1];

    // ensure we're using the right triangle orientations (consistent with old splitting code)
    if ( !mesh.oriented( vertex_a, vertex_b, mesh.get_triangle(tri0) ) )
      swap(vertex_c, vertex_d);

    Vec3d tri0_normal = m_surf.get_triangle_normal( tri0 );
    Vec3d tri1_normal = m_surf.get_triangle_normal( tri1 );

    if ( dot( tri0_normal, tri1_normal ) >= 0.0 )
    {
      Vec3d new_normal = triangle_normal( m_surf.get_position(vertex_a), new_vertex_smooth_position, m_surf.get_position(vertex_c) );
      if ( dot( new_normal, tri0_normal ) < 0.0 || dot( new_normal, tri1_normal ) < 0.0 )
      {
        use_smooth_point = false;
      }
      new_normal = triangle_normal( m_surf.get_position(vertex_c), new_vertex_smooth_position, m_surf.get_position(vertex_b) );
      if ( dot( new_normal, tri0_normal ) < 0.0 || dot( new_normal, tri1_normal ) < 0.0 )
      {
        use_smooth_point = false;
      }         
      new_normal = triangle_normal( m_surf.get_position(vertex_d), m_surf.get_position(vertex_b), new_vertex_smooth_position );
      if ( dot( new_normal, tri0_normal ) < 0.0 || dot( new_normal, tri1_normal ) < 0.0 )
      {
        use_smooth_point = false;
      }         
      new_normal = triangle_normal( m_surf.get_position(vertex_d), new_vertex_smooth_position, m_surf.get_position(vertex_a) );
      if ( dot( new_normal, tri0_normal ) < 0.0 || dot( new_normal, tri1_normal ) < 0.0 )
      {
        use_smooth_point = false;
      }         
    }
  }

  // --------------

  // if the new point introduces an intersection, try using the edge midpoint
  

  if ( use_smooth_point == false )
  {

    if ( m_surf.m_verbose ) { std::cout << "not using smooth subdivision" << std::endl; }

    new_vertex_smooth_position = new_vertex_position;

    if ( split_edge_pseudo_motion_introduces_intersection( new_vertex_position, 
      new_vertex_smooth_position, 
      edge, 
      vertex_a, 
      vertex_b, 
      incident_tris, 
      other_verts ) )
    {

      g_stats.add_to_int( "EdgeSplitter:split_midpoint_collisions", 1 );

      if ( m_surf.m_verbose )  { std::cout << "Even mid-point subdivision introduces collision.  Backing out." << std::endl; }
      return false;
    }
  }
  else
  {
    if ( m_surf.m_verbose ) { std::cout << "using smooth subdivision" << std::endl; }
  }


  // --------------

  // Check angles on new triangles

  const Vec3d& va = m_surf.get_position( vertex_a );
  const Vec3d& vb = m_surf.get_position( vertex_b );
  std::vector<Vec3d> other_vert_pos;
  for(size_t i = 0; i < other_verts.size(); ++i)
    other_vert_pos.push_back(m_surf.get_position(other_verts[i]));

  double min_new_angle = 2*M_PI;
  for(size_t i = 0; i < other_vert_pos.size(); ++i) {
    min_new_angle = min( min_new_angle, min_triangle_angle( va, new_vertex_smooth_position, other_vert_pos[i] ) );
    min_new_angle = min( min_new_angle, min_triangle_angle( vb, new_vertex_smooth_position, other_vert_pos[i] ) );
  }
  
  if ( rad2deg(min_new_angle) < m_surf.m_min_triangle_angle )
  {
    g_stats.add_to_int( "EdgeSplitter:edge_split_small_angle", 1 );
    return false;
  }

  double max_current_angle = 0;
  for(size_t i = 0; i < incident_tris.size(); ++i) {
    max_current_angle = max( max_current_angle, max_triangle_angle( va, vb, other_vert_pos[i] ) );
  }

  double max_new_angle = 0;
  
  for(size_t i = 0; i < other_vert_pos.size(); ++i) {
    max_new_angle = max( max_new_angle, max_triangle_angle( va, new_vertex_smooth_position, other_vert_pos[i] ) );
    max_new_angle = max( max_new_angle, max_triangle_angle( vb, new_vertex_smooth_position, other_vert_pos[i] ) );
  }

  // if new angle is greater than the allowed angle, and doesn't 
  // improve the current max angle, prevent the split

  if ( rad2deg(max_new_angle) > m_surf.m_max_triangle_angle )
  {

    // if new triangle improves a large angle, allow it

    if ( rad2deg(max_new_angle) < rad2deg(max_current_angle) )
    {
      g_stats.add_to_int( "EdgeSplitter:edge_split_large_angle", 1 );      
      return false;
    }
  }

  // --------------
  
  // Do the actual splitting

  double new_vertex_mass = 0.5 * ( m_surf.m_masses[ vertex_a ] + m_surf.m_masses[ vertex_b ] );
  size_t vertex_e = m_surf.add_vertex( new_vertex_smooth_position, new_vertex_mass );

  // Add to change history
  m_surf.m_vertex_change_history.push_back( VertexUpdateEvent( VertexUpdateEvent::VERTEX_ADD, vertex_e, Vec2st( vertex_a, vertex_b) ) );
  
  if ( m_surf.m_verbose ) { std::cout << "new vertex: " << vertex_e << std::endl; }

  // Create new triangles with proper orientations (match their parents)
  std::vector<Vec3st> created_tri_data;
  for(size_t i = 0; i < other_verts.size(); ++i) {
    Vec3st newtri0, newtri1;
    if ( mesh.oriented( vertex_a, vertex_b, mesh.get_triangle( incident_tris[i] ) ) ) {
      newtri0 = Vec3st( vertex_a, vertex_e, other_verts[i] );
      newtri1 = Vec3st( other_verts[i], vertex_e, vertex_b );
    }
    else {
      newtri0 = Vec3st( vertex_a, other_verts[i], vertex_e );
      newtri1 = Vec3st( other_verts[i], vertex_b, vertex_e );
    }
    created_tri_data.push_back(newtri0);
    created_tri_data.push_back(newtri1);
  }

  // Delete the parent triangles
  for(size_t i = 0; i < incident_tris.size(); ++i) {
    m_surf.remove_triangle( incident_tris[i] );
  }

  // Now actually add the triangles to the mesh
  std::vector<size_t> created_tris;
  for(size_t i = 0; i < created_tri_data.size(); ++i) {
    //add the triangle
    size_t newtri0_id = m_surf.add_triangle( created_tri_data[i] );
    
    //record the data created
    created_tris.push_back(newtri0_id);
  }

  
  // Add to new history log
  MeshUpdateEvent split(MeshUpdateEvent::EDGE_SPLIT);
  split.m_v0 = vertex_a;
  split.m_v1 = vertex_b;
  split.m_vert_position = new_vertex_smooth_position;
  split.m_created_verts.push_back(vertex_e);
  for(size_t i = 0; i < incident_tris.size(); ++i)
    split.m_deleted_tris.push_back(incident_tris[i]);
  split.m_created_tris = created_tris;
  split.m_created_tri_data = created_tri_data;
  
  m_surf.m_mesh_change_history.push_back(split);

  return true;

}


// --------------------------------------------------------
///
/// Determine if an edge split is desirable
///
// --------------------------------------------------------


bool EdgeSplitter::edge_length_needs_split(size_t edge_index) {
    
    double edge_length = m_surf.get_edge_length(edge_index);
    size_t vertex_a = m_surf.m_mesh.m_edges[edge_index][0];
    size_t vertex_b = m_surf.m_mesh.m_edges[edge_index][1];

    if ( m_use_curvature )
    {
        //split if we're above the upper limit
        if(edge_length > m_max_edge_length)
            return true;

        //don't split if splitting would take us below the lower limit
        if(edge_length < 2*m_min_edge_length)
            return false;

        double curvature_value = get_edge_curvature( m_surf, vertex_a, vertex_b );
        int circlesegs = 16;
        double curvature_max_length = 2*M_PI / (double)circlesegs / max(curvature_value, 1e-8);
        
        //split if curvature dictates
        if(edge_length > curvature_max_length)
            return true;
        
        //check all incident edges to see if any of them are super short, and if so, split this guy accordingly.
        //this enforces slow grading of the mesh.
        double min_nbr_len = edge_length;
        for(size_t edge_id = 0; edge_id < m_surf.m_mesh.m_vertex_to_edge_map[vertex_a].size(); ++edge_id) {
            min_nbr_len = min(min_nbr_len, m_surf.get_edge_length(m_surf.m_mesh.m_vertex_to_edge_map[vertex_a][edge_id]));
        }
        for(size_t edge_id = 0; edge_id < m_surf.m_mesh.m_vertex_to_edge_map[vertex_b].size(); ++edge_id) {
            min_nbr_len = min(min_nbr_len, m_surf.get_edge_length(m_surf.m_mesh.m_vertex_to_edge_map[vertex_b][edge_id]));
        }
        
        if(edge_length > min_nbr_len * 3)
            return true;
              
    }
    else {
        return edge_length > m_max_edge_length;
    }
    
    return false;

}

// --------------------------------------------------------
///
/// Determine if edge should be allowed to be split
///
// --------------------------------------------------------

bool EdgeSplitter::edge_is_splittable( size_t edge_index )
{

  // skip deleted and solid edges
  if ( m_surf.m_mesh.edge_is_deleted(edge_index) ) { return false; }
  if ( m_surf.edge_is_all_solid(edge_index) ) { return false; }

  //if not remeshing boundary edges, skip those too
  if ( !m_remesh_boundaries && m_surf.m_mesh.m_is_boundary_edge[edge_index]) { return false; }

  return true;

}
// --------------------------------------------------------
///
/// Split edges opposite large angles
///
// --------------------------------------------------------

bool EdgeSplitter::large_angle_split_pass()
{

  NonDestructiveTriMesh& mesh = m_surf.m_mesh;

  bool split_occurred = false;

  for ( size_t e = 0; e < mesh.m_edges.size(); ++e )
  {

    if ( !edge_is_splittable(e) ) { continue; }

    // get edge end points
    const Vec2st& edge = m_surf.m_mesh.m_edges[e];      
    const Vec3d& edge_point0 = m_surf.get_position( edge[0] );
    const Vec3d& edge_point1 = m_surf.get_position( edge[1] );

    // get triangles incident to the edge
    std::vector<size_t> incident_tris = mesh.m_edge_to_triangle_map[e];
    for(size_t t = 0; t < incident_tris.size(); ++t) {
      const Vec3st& tri0 = mesh.get_triangle(incident_tris[t]);
    
      // get vertex opposite the edge for each triangle
      size_t opposite0 = mesh.get_third_vertex( e, tri0 );
    
      // compute the angle at each opposite vertex
      const Vec3d& opposite_point0 = m_surf.get_position(opposite0);
      
      double angle0 = rad2deg( acos( dot( normalized(edge_point0-opposite_point0), normalized(edge_point1-opposite_point0) ) ) );
    

      // if an angle is above the max threshold, split the edge

      if ( angle0 > m_surf.m_max_triangle_angle )
      {

        bool result = split_edge( e );

        if ( result )
        {
          g_stats.add_to_int( "EdgeSplitter:large_angle_split_success", 1 );
        }
        else
        {
          g_stats.add_to_int( "EdgeSplitter:large_angle_split_failed", 1 );
        }

        split_occurred |= result;
      }
    }
  }


  return split_occurred;

}


// --------------------------------------------------------
///
/// Split all long edges
///
// --------------------------------------------------------

bool EdgeSplitter::split_pass()
{
    
    if ( m_surf.m_verbose )
    {
        std::cout << "---------------------- Edge Splitter: splitting ----------------------" << std::endl;
    }
    
    assert( m_max_edge_length != UNINITIALIZED_DOUBLE );
    
    NonDestructiveTriMesh& mesh = m_surf.m_mesh;
    std::vector<SortableEdge> sortable_edges_to_try;
    
    for( size_t i = 0; i < mesh.m_edges.size(); i++ )
    {    
        if ( !edge_is_splittable(i) ) { continue; }
        
        bool should_split = edge_length_needs_split(i);
        if(should_split)
            sortable_edges_to_try.push_back( SortableEdge( i, m_surf.get_edge_length(i)) );
    }
    
    
    //
    // sort in ascending order, then iterate backwards to go from longest edge to shortest
    //
    
    // whether a split operation was successful in this pass
    
    bool split_occurred = false;
    
    std::sort( sortable_edges_to_try.begin(), sortable_edges_to_try.end() );
    
    std::vector<SortableEdge>::reverse_iterator iter = sortable_edges_to_try.rbegin();
    
    for ( ; iter != sortable_edges_to_try.rend(); ++iter )
    {
        size_t longest_edge = iter->m_edge_index;
        
        if ( !edge_is_splittable(longest_edge) ) { continue; }
        
        bool should_split = edge_length_needs_split(longest_edge);
        
        if(should_split) {
            bool result = split_edge(longest_edge);
            split_occurred |= result;
        }

    }
    
    
    // Now split to reduce large angles
    bool large_angle_split_occurred = large_angle_split_pass();
    
    return split_occurred || large_angle_split_occurred;
    
}

}
