// ---------------------------------------------------------
//
//  dynamicsurface.h
//  Tyson Brochu 2008
//  
//  A triangle mesh with associated vertex locations and  masses.  Query functions for getting geometry info.
//
//  The most important function is integrate(), which advances the mesh vertices from m_positions to m_newpositions, while 
//  performing collision detection and resolution.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_DYNAMICSURFACE_H
#define EL_TOPO_DYNAMICSURFACE_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <ccd_wrapper.h>
#include <nondestructivetrimesh.h>
#include <limits>

#define PBC_DOMAIN_SIZE_X    (1.0)
#define PBC_DOMAIN_SIZE_Y    (1.0)
#define PBC_DOMAIN_SIZE_Z    (1.0)

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

namespace ElTopo {

// Broad-phase collision detector.  Avoids performing collision detection between far-away primitives.
class BroadPhase;

// Class for encapsulating all collision detection and resolution functionality.
class CollisionPipeline;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// --------------------------------------------------------
///
/// Edge-triangle intersection info.
///
// --------------------------------------------------------

struct Intersection
{
    
    /// Constructor
    ///
    Intersection(size_t edge_index, 
                 size_t triangle_index ) :
    m_edge_index( edge_index ),
    m_triangle_index( triangle_index )
    {}
    
    /// The index of the edge intersecting the triangle
    ///
    size_t m_edge_index;
    
    /// Index of the triangle intersecting the edge
    ///
    size_t m_triangle_index;
    
};


// --------------------------------------------------------
///
/// A surface mesh.  Essentially consists of a NonDestructiveTriMesh object coupled with a set of vertex locations in 3D space.
///
// --------------------------------------------------------

class DynamicSurface
{
    
public:
    
    /// Create a DynamicSurface object from the specified vertices, triangles, and vertex masses.
    ///

    DynamicSurface( const std::vector<Vec3d>& vs, 
       const std::vector<Vec3st>& ts, 
       const std::vector<Vec2i>& labels,
       const std::vector<Vec3d>& masses,
       double in_proximity_epsilon = 1e-4,
       double in_friction_coefficient = 0.0,
       bool in_collision_safety = true,
       bool in_verbose = false );
    
    /// Destructor
    /// 
    virtual ~DynamicSurface(); 
    
private:
    
    /// Disallowed, do not implement
    ///
    DynamicSurface( const DynamicSurface& );
    
    /// Disallowed, do not implement
    ///
    DynamicSurface& operator=( const DynamicSurface& );
    
public:
    
    /// Advance from current state to a collision-free state as close as possible to predicted state.
    /// 
    virtual void integrate( double dt, double& actual_dt );

    //
    // Utility
    //

    /// Break up the triangle mesh into connected components, determine surface IDs for all vertices.
    ///
    void partition_surfaces( std::vector<size_t>& surface_ids, std::vector< std::vector< size_t> >& surfaces ) const;
    
    /// Compute the area of the specified triangle
    ///
    inline double get_triangle_area(size_t tri) const;
    
    /// Compute the area of the specified triangle
    ///    
    inline double get_triangle_area(const Vec3st& tri) const;
    
    /// Compute the area of the specified triangle
    ///
    inline double get_triangle_area(size_t v0, size_t v1, size_t v2) const;
    
    /// Get the smallest triangle in the mesh by area.
    ///
    inline double get_min_triangle_area( size_t& out_triangle_index ) const;
    
    /// Compute the vector normal to the specified triangle's plane
    ///
    inline Vec3d get_triangle_normal(size_t tri) const;
    
    /// Compute the vector normal to the specified triangle's plane
    ///    
    inline Vec3d get_triangle_normal(const Vec3st& tri) const;
    
    /// Compute the vector normal to the specified triangle's plane
    ///    
    inline Vec3d get_triangle_normal(size_t v0, size_t v1, size_t v2) const;
    
    /// Compute the vector normal to the specified triangle's plane
    ///
    inline Vec3d get_triangle_normal_by_region(size_t tri, int region) const;

    /// Compute the specified triangle's barycenter
    ///    
    inline Vec3d get_triangle_barycenter( size_t triangle_index ) const;
    
    /// Get an estimate for the surface normal at the specified vertex. Computed using an unweighted average of the normals of 
    /// incident triangles.
    ///
    inline Vec3d get_vertex_normal( size_t vertex ) const;
    
    /// Get an estimate for the surface normal at the specified vertex. Computed using angleweighted pseudonormal
    ///
    inline Vec3d get_vertex_normal_angleweighted( size_t vertex ) const;

    /// Get an estimate for the surface normal at the specified vertex. Computed using angleweighted pseudonormal
    ///
    inline Vec3d get_vertex_normal_angleweighted_by_label( size_t vertex, int label ) const;

    /// Compute all vertex normals, using an unweighted average of incident triangle normals.
    ///
    void get_all_vertex_normals( std::vector<Vec3d>& normals ) const;
    
    /// Get an estimate for the surface normal at the specified vertex. Computed using a weighted average as described in [Max 1999].
    ///
    inline Vec3d get_vertex_normal_max( size_t vertex_index ) const;
    
    /// Compute the edge length of the specified edge
    ///
    inline double get_edge_length( size_t edge_index ) const;
    
    /// Get the average edge length over all edges
    ///
    inline double get_average_edge_length() const;
    
    /// Get the average edge length, disregarding edges with vertices marked as solid
    ///
    inline double get_average_non_solid_edge_length() const;
    
    /// Determine if the vertex is on a solid surface (has infinite mass).
    ///
    inline bool vertex_is_solid( size_t vertex_index, int dof ) const;
    inline Vec3c vertex_is_solid_3 (size_t vertex_index ) const;
    inline bool vertex_is_all_solid( size_t vertex_index ) const;
    inline bool vertex_is_any_solid( size_t vertex_index ) const;
    
    /// Determine if the edge has any solid vertices
    ///
    inline bool edge_is_any_solid( size_t eedge_index ) const;

    /// Determine if the triangle has any solids vertices
    ///    
    inline bool triangle_is_any_solid( size_t triangle_index ) const;
    
    /// Determine if the edge has all solid vertices
    ///
    inline bool edge_is_all_solid( size_t eedge_index ) const;

    /// Determine if the triangle has all solid vertices
    ///    
    inline bool triangle_is_all_solid( size_t triangle_index ) const;

    /// Compute the total surface area defined by the mesh
    ///
    inline double get_surface_area( ) const;
    
    /// Compute the total surface area using predicted vertex locations (m_newpositions)
    ///
    inline double get_predicted_surface_area() const;
    
    /// Compute the volume enclosed by the surface
    ///
    inline double get_volume() const;
    
    /// Compute the volume enclosed by the surface, using predicted vertex locations (m_newpositions)
    ///
    inline double get_predicted_volume() const;
    
    /// Compute the distance from the given point to the surface.  Also return the index of the closest triangle.
    ///
    double distance_to_surface( const Vec3d& p, size_t& out_closest_triangle ) const;
    
    /// Determine the rank of the primary space at the given vertex (see Jiao07).
    /// Rank {1, 2, 3} == {smooth, ridge, peak}
    ///
    unsigned int vertex_primary_space_rank( size_t v, int region = -1) const;
    
    //unsigned int vertex_primary_space_rank_nonmanifold( size_t v ) const;

    unsigned int compute_rank_from_triangles(const std::vector<size_t>& tris) const;
    
    /// Get edge dihedral angle
    double get_largest_dihedral(size_t edge) const;
    double get_largest_dihedral(size_t edge, const std::vector<Vec3d>& cached_normals) const;

    /// Determine which region the point is inside by ray-casting and looking at the normal
    /// of the first intersection, and comparing that with the triangle's labeling
    int get_region_containing_point( const Vec3d& p);
    int test_region_via_ray_and_normal(const Vec3d& p, const Vec3d& ray_end);

    //
    // Broad phase collision detector
    //
    
    /// Delete and rebuild the broad phase object, using AABBs defined from m_positions.
    ///
    void rebuild_static_broad_phase( );
    
    /// Delete and rebuild the broad phase object, using AABBs defined from m_positions and m_newpositions.
    ///
    void rebuild_continuous_broad_phase( );
    
    /// Assume that the specified vertex has moved.  Update the broadphase entries of the vertex, its incident edges and incident 
    /// triangles.
    ///
    void update_static_broad_phase( size_t vertex_index );
    
    /// Assume that the specified vertex has moved.  Update the broadphase entries of the vertex, its incident edges and incident 
    /// triangles, using current and predicted positions.
    ///    
    void update_continuous_broad_phase( size_t vertex_index );  
    
    /// Get the padded AABB of the specified vertex using current positions.
    ///
    void vertex_static_bounds(size_t v, Vec3d &xmin, Vec3d &xmax) const;
    
    /// Get the padded AABB of the specified edge using current positions.
    ///
    void edge_static_bounds(size_t e, Vec3d &xmin, Vec3d &xmax) const;
    
    /// Get the padded AABB of the specified triangle using current positions.
    ///
    void triangle_static_bounds(size_t t, Vec3d &xmin, Vec3d &xmax) const; 
    
    /// Get the padded AABB of the specified vertex using current and predicted positions.
    ///
    void vertex_continuous_bounds(size_t v, Vec3d &xmin, Vec3d &xmax) const;
    
    /// Get the padded AABB of the specified edge using current and predicted positions.
    ///    
    void edge_continuous_bounds(size_t e, Vec3d &xmin, Vec3d &xmax) const;
    
    /// Get the padded AABB of the specified triangle using current and predicted positions.
    ///        
    void triangle_continuous_bounds(size_t t, Vec3d &xmin, Vec3d &xmax) const;
    
    /// Caution: slow!
    /// Check the consistency of the broad phase by comparing against the N^2 broadphase.
    ///
    void check_static_broad_phase_is_up_to_date() const;
    
    /// Caution: slow!
    /// Check the consistency of the broad phase by comparing against the N^2 broadphase.  Checks using current and predicted vertex 
    /// positions.
    void check_continuous_broad_phase_is_up_to_date() const;
    
    // 
    // Intersection detection 
    //
    
    /// Get the set of triangles intersected by the given segment. Also returns the barycentric coordinate along the segment of each
    /// intersection.
    ///
    void get_triangle_intersections(const Vec3d& segment_point_a, 
                                    const Vec3d& segment_point_b,
                                    std::vector<double>& hit_ss,
                                    std::vector<size_t>& hit_triangles,
                                    bool verbose = false) const;
    
    /// Count the number of triangles intersected by the given segment.
    ///
    size_t get_number_of_triangle_intersections( const Vec3d& segment_point_a, 
                                                const Vec3d& segment_point_b ) const;
    
    /// Using exact intersection testing, count the number of triangle intersected by the given segment.
    ///
    size_t get_number_of_triangle_intersections_exact( const Vec3d& segment_point_a, 
                                                      const Vec3d& segment_point_b ) const;
    
    /// Test the given triangle against all other triangles in the mesh for intersection.
    ///
    bool check_triangle_vs_all_triangles_for_intersection( size_t tri_index );
    
    /// Test the given triangle against all other triangles in the mesh for intersection.
    ///
    bool check_triangle_vs_all_triangles_for_intersection( const Vec3st& tri );
    
    /// Get all self-intersections in the surface
    ///
    void get_intersections(bool degeneracy_counts_as_intersection, 
                           bool use_new_positions, 
                           std::vector<Intersection>& intersections );
    
    /// Look for self-intersections, but stop when the first one is found
    ///
    void get_first_intersection( bool degeneracy_counts_as_intersection, 
                                bool use_new_positions, 
                                Intersection& intersections );
    
    /// Fire an assert if the mesh contains a self-intersection. Uses m_positions as the vertex locations.
    ///
    void assert_mesh_is_intersection_free( bool degeneracy_counts_as_intersection );              
    
    /// Using m_newpositions as the vertex locations, fire an assert if the mesh contains a self-intersection.
    ///
    void assert_predicted_mesh_is_intersection_free( bool degeneracy_counts_as_intersection ); 

    /// Returns the number of vertices in the mesh, including any vertices marked as deleted
    ///
    inline size_t get_num_vertices() const;
    
    /// Returns the current positions of a vertex.
    ///
    inline const Vec3d& get_position( size_t index ) const; // PBC: this returns the position of the base instance of vertex index, which is the same with the non-PBC implementation
    inline       Vec3d  get_position( size_t index, size_t ref_vert ) const; // PBC: this returns the position of one instance (out of the infinitely many instances) of vertex index that's closest to the position of the base instance of vertex ref_vert
    inline       Vec3d  get_position( size_t index, const Vec3d & ref_pos ) const; // PBC: this returns the position of one instance (out of the infinitely many instances) of vertex index that's closest to the position ref_pos
    inline       Vec3d  get_position( const Vec3d & pos, size_t ref_vert ) const; // PBC: this returns the position of one instance (out of the infinitely many instances) of position pos that's closest to the position of the base instance of vertex ref_vert
    inline       Vec3d  get_position( const Vec3d & pos, const Vec3d & ref_pos ) const; // PBC: this returns the position of one instance (out of the infinitely many instances) of position pos that's closest to the position ref_pos
  
    /// Returns the periodic domain offset of a position, an offset such that each component of (pos - offset) is in [0, 1)
    inline Vec3d get_domain_offset(const Vec3d & pos) const;

    /// Returns the set of all current vertex positions.
    ///
    inline const std::vector<Vec3d>& get_positions( ) const;
    
    /// Set the current position of an individual vertex.
    ///
    inline void set_position( size_t index, const Vec3d& x );   
    
    /// Set the current positions of all vertices in the mesh.
    ///
    inline void set_all_positions( const std::vector<Vec3d>& xs );
    
    /// Set the current positions of all vertices in the mesh, from a C-array of doubles
    ///
    inline void set_all_positions( size_t n, const double* xs );
    
    /// Copy predicted vertex positions into the current positions
    ///
    inline void set_positions_to_newpositions();

    /// Returns the predicted position of a vertex.
    ///
    inline const Vec3d& get_newposition( size_t index ) const; // PBC: this returns the position of the base instance of vertex index, which is the same with the non-PBC implementation
    inline       Vec3d  get_newposition( size_t index, size_t ref_vert ) const; // PBC: this returns the position of one instance (out of the infinitely many instances) of vertex index that's closest to the position of the base instance of vertex ref_vert
    inline       Vec3d  get_newposition( size_t index, const Vec3d & ref_pos ) const; // PBC: this returns the position of one instance (out of the infinitely many instances) of vertex index that's closest to the position ref_pos
    
    /// Returns the set of all predicted vertex positions.
    ///
    inline const std::vector<Vec3d>& get_newpositions( ) const;
    
    /// Set the predicted position of an individual vertex.
    ///
    inline void set_newposition( size_t index, const Vec3d& x );   
    
    /// Set the predicted positions of all vertices in the mesh.
    ///
    inline void set_all_newpositions( const std::vector<Vec3d>& xs );
    
    /// Set the predicted positions of all vertices in the mesh from a C-array.
    ///
    inline void set_all_newpositions( size_t n, const double* xs );
    
    ///////////////////////////////////////////////////////////////////////
    // FD 20130102
    
    // Vertex velocity, for remeshing only, not for collisions
    inline void set_all_remesh_velocities( const std::vector<Vec3d> & v );
  
    inline Vec3d get_remesh_velocity( size_t n );
    
    inline void set_remesh_velocity( size_t n, const Vec3d & v );
    


    ///////////////////////////////////////////////////////////////////////

    //
    // Data members
    //
    
    /// Elements closer than this have repulsion forces applied
    ///
    double m_proximity_epsilon;
    
    /// Dump lots of details to stdout
    ///
    bool m_verbose;
    
    /// Ensure that no mesh elements intersect, during mesh moving and mesh maintenance
    ///
    bool m_collision_safety;
    
    /// Vertex positions, predicted locations, velocities and masses
    ///
    std::vector<Vec3d> m_masses;
    
    /// The mesh graph
    ///
    NonDestructiveTriMesh m_mesh;
    
    /// collision acceleration structures
    ///
    BroadPhase* m_broad_phase;
    
    /// Encapsulates the collision detection functionality
    ///
    CollisionPipeline* m_collision_pipeline;
    
    /// Amount to pad AABBs by when doing broad-phase collision detection
    ///
    double m_aabb_padding;
    
    /// Angle threshold for declaring an edge to be a feature
    ///
    double m_feature_edge_angle_threshold;

public:
  
    friend class CollisionPipeline;
    friend class ImpactZoneSolver;
    friend class MeshSmoother;
    
    /// Current and predicted vertex positions
    ///
    std::vector<Vec3d> pm_positions, pm_newpositions;
    
    ///////////////////////////////////////////////////////////////////////
    // FD 20130102

    // Vertex velocity, for remeshing only, not for collision
    std::vector<Vec3d> pm_velocities;
    
    ///////////////////////////////////////////////////////////////////////

    /// Temporary velocities field
    ///
    std::vector<Vec3d> m_velocities;
    
};


// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------

// --------------------------------------------------------
///
/// Compute area of a triangle specified by three vertices
///
// --------------------------------------------------------

inline double triangle_area( const Vec3d& v0, const Vec3d &v1, const Vec3d &v2 )
{
    return 0.5 * mag( cross( v1 - v0, v2 - v0 ) );
}

// --------------------------------------------------------
///
/// Compute area of a triangle specified by a triangle index
///
// --------------------------------------------------------

inline double DynamicSurface::get_triangle_area(size_t tri) const
{
    const Vec3st &t = m_mesh.get_triangle( tri ); 
    return get_triangle_area(t[0], t[1], t[2]);
}

// --------------------------------------------------------
///
/// Compute area of a triangle specified by a triple of vertex indices
///
// --------------------------------------------------------

inline double DynamicSurface::get_triangle_area(const Vec3st& tri) const
{
    return get_triangle_area(tri[0], tri[1], tri[2]);
}

// --------------------------------------------------------
///
/// Compute area of a triangle specified by a three vertex indices
///
// --------------------------------------------------------

inline double DynamicSurface::get_triangle_area(size_t v0, size_t v1, size_t v2) const
{
    const Vec3d &p0 = get_position(v0);
    const Vec3d &p1 = get_position(v1, v0);
    const Vec3d &p2 = get_position(v2, v0);
    
    return 0.5 * mag(cross(p1-p0, p2-p0));
}

// --------------------------------------------------------
///
/// Compute the normal of a triangle specified by three vertices
///
// --------------------------------------------------------

inline Vec3d triangle_normal( const Vec3d& v0, const Vec3d &v1, const Vec3d &v2 )
{
    Vec3d u = v1 - v0;
    Vec3d v = v2 - v0;
    return normalized(cross(u, v));
}

// --------------------------------------------------------
///
/// Compute the normal of a triangle specified by a triangle index
///
// --------------------------------------------------------

inline Vec3d DynamicSurface::get_triangle_normal(size_t tri) const
{
    const Vec3st &t = m_mesh.get_triangle( tri ); 
    return get_triangle_normal(t[0], t[1], t[2]);
}

// --------------------------------------------------------
///
/// Compute the normal of a triangle specified by a triangle index, and the relevant region
///
// --------------------------------------------------------

inline Vec3d DynamicSurface::get_triangle_normal_by_region(size_t tri, int region) const
{
   const Vec3st &t = m_mesh.get_triangle( tri ); 
   Vec2i label = m_mesh.get_triangle_label(tri);
   if(label[0] != region && label[1] != region) return Vec3d(0,0,0);
   Vec3d normal = get_triangle_normal(tri);
   return label[1] == region? normal : -normal;
}
// --------------------------------------------------------
///
/// Compute the normal of a triangle specified by a triple of vertex indices
///
// --------------------------------------------------------

inline Vec3d DynamicSurface::get_triangle_normal(const Vec3st& tri) const
{
    return get_triangle_normal(tri[0], tri[1], tri[2]);
}

// --------------------------------------------------------
///
/// Compute the normal of a triangle specified by three vertex indices
///
// --------------------------------------------------------

inline Vec3d DynamicSurface::get_triangle_normal(size_t v0, size_t v1, size_t v2) const
{
    Vec3d start = get_position(v0);
    Vec3d u = get_position(v1, v0) - start;
    Vec3d v = get_position(v2, v0) - start;
    Vec3d res = cross(u, v);
    normalize(res);
    return res;
}

// --------------------------------------------------------

inline Vec3d DynamicSurface::get_triangle_barycenter( size_t triangle_index ) const
{
    const Vec3st& tri = m_mesh.get_triangle( triangle_index );
    return 1.0 / 3.0 * ( get_position( tri[0] ) + get_position( tri[1], tri[0] ) + get_position( tri[2], tri[0] ) );
}

// --------------------------------------------------------
///
/// Return the triangle with the the smallest area, and that area
///
// --------------------------------------------------------

inline double DynamicSurface::get_min_triangle_area( size_t& triangle_index ) const
{
    double min_area = BIG_DOUBLE;
    for ( size_t i = 0; i < m_mesh.num_triangles(); ++i )
    {
        if ( m_mesh.get_triangle(i)[0] == m_mesh.get_triangle(i)[1] )
        {
            continue;
        }
        
        double area = get_triangle_area(i);
        if ( area < min_area )
        {
            min_area = area;
            triangle_index = i;
        }
    }
    
    return min_area;
}

// --------------------------------------------------------
///
/// Compute surface normal at the specified vertex (unweighted average of incident triangle normals).
///
// --------------------------------------------------------

inline Vec3d DynamicSurface::get_vertex_normal( size_t vertex ) const
{
    Vec3d normal(0,0,0);
    for ( size_t i = 0; i < m_mesh.m_vertex_to_triangle_map[vertex].size(); ++i )
    {
        normal += get_triangle_normal( m_mesh.m_vertex_to_triangle_map[vertex][i] );
    }
    normal /= double(m_mesh.m_vertex_to_triangle_map[vertex].size());
    normal /= mag(normal);
    
    return normal;
}

// --------------------------------------------------------
///
/// Compute surface normal at the specified vertex (angle-weighted pseudo-normal) per region.
///
// --------------------------------------------------------

inline Vec3d DynamicSurface::get_vertex_normal_angleweighted_by_label( size_t vertex_index, int label ) const
{
   const std::vector<size_t>& inc_tris = m_mesh.m_vertex_to_triangle_map[vertex_index];

   Vec3d normal_sum(0,0,0);

   for ( unsigned int i = 0; i < inc_tris.size(); ++i )
   {
      const Vec3st & curr_tri = m_mesh.m_tris[inc_tris[i]];
      
      Vec2i tri_label = m_mesh.get_triangle_label(inc_tris[i]);
      if(tri_label[0] != label && tri_label[1] != label) continue; //no matching label, skip it.

      if ( curr_tri[0] == curr_tri[1] ) { continue; }

      Vec2st other_two;

      NonDestructiveTriMesh::index_in_triangle( curr_tri, vertex_index, other_two );

      unsigned int verti = curr_tri[other_two[0]];
      unsigned int vertnext = curr_tri[other_two[1]];

      Vec3d vi = get_position(verti) - get_position(vertex_index);
      Vec3d vnext = get_position(vertnext) - get_position(vertex_index);
      normalize(vi);
      normalize(vnext);

      double dotproduct = dot(vi,vnext);
      double angle = acos(dotproduct);
      Vec3d normal = cross(vi,vnext);
      normalize(normal);

      if(tri_label[0] == label) normal *= -1; //flip the normal in the case of label being on the reverse.

      normal_sum += angle*normal;
   }

   normalize(normal_sum);

   return normal_sum;
}

// --------------------------------------------------------
///
/// Compute surface normal at the specified vertex (angle-weighted pseudo-normal).
/// 
///
// --------------------------------------------------------

inline Vec3d DynamicSurface::get_vertex_normal_angleweighted( size_t vertex_index ) const
{
   const std::vector<size_t>& inc_tris = m_mesh.m_vertex_to_triangle_map[vertex_index];

   Vec3d normal_sum(0,0,0);

   for ( size_t i = 0; i < inc_tris.size(); ++i )
   {
      const Vec3st& curr_tri = m_mesh.m_tris[inc_tris[i]];

      if ( curr_tri[0] == curr_tri[1] ) { continue; } //testing for flaps?

      Vec2st other_two;

      NonDestructiveTriMesh::index_in_triangle( curr_tri, vertex_index, other_two );

      unsigned int verti = curr_tri[other_two[0]];
      unsigned int vertnext = curr_tri[other_two[1]];

      Vec3d vi = get_position(verti) - get_position(vertex_index);
      Vec3d vnext = get_position(vertnext) - get_position(vertex_index);
      normalize(vi);
      normalize(vnext);

      double dotproduct = dot(vi,vnext);
      double angle = acos(dotproduct);
      Vec3d normal = cross(vi,vnext);
      normalize(normal);

      normal_sum += angle*normal;
   }

   normalize(normal_sum);

   return normal_sum;
}

// --------------------------------------------------------
///
/// Compute surface normal at the specified vertex (weighted according to [Max 1999]).
///
// --------------------------------------------------------

inline Vec3d DynamicSurface::get_vertex_normal_max( size_t vertex_index ) const
{
    const std::vector<size_t>& inc_tris = m_mesh.m_vertex_to_triangle_map[vertex_index];
    
    Vec3d sum_cross_products(0,0,0);
    
    for ( size_t i = 0; i < inc_tris.size(); ++i )
    {
        const Vec3st& curr_tri = m_mesh.get_triangle( inc_tris[i] );
        
        if ( curr_tri[0] == curr_tri[1] ) { continue; }
        
        Vec2st other_two;
        NonDestructiveTriMesh::index_in_triangle( curr_tri, vertex_index, other_two );
        
        size_t verti = curr_tri[other_two[0]];
        size_t vertnext = curr_tri[other_two[1]];
        
        Vec3d vi = get_position(verti) - get_position(vertex_index);
        Vec3d vnext = get_position(vertnext) - get_position(vertex_index);
        
        sum_cross_products += cross( vi, vnext ) / ( mag2(vi)*mag2(vnext) );
    }
    
    sum_cross_products /= mag( sum_cross_products );
    
    return sum_cross_products;
}


// --------------------------------------------------------
///
/// Compute length of the specified edge
///
// --------------------------------------------------------

inline double DynamicSurface::get_edge_length( size_t edge_index ) const
{
    return mag( get_position( m_mesh.m_edges[edge_index][1], m_mesh.m_edges[edge_index][0] ) - get_position( m_mesh.m_edges[edge_index][0] ) );
}

// --------------------------------------------------------
///
/// Compute average length over all mesh edges
///
// --------------------------------------------------------

inline double DynamicSurface::get_average_edge_length() const
{
    double sum_lengths = 0;
    for ( size_t i = 0; i < m_mesh.m_edges.size(); ++i )
    {
        const Vec2st& e = m_mesh.m_edges[i]; 
        if ( e[0] == e[1] )  { continue; }
        sum_lengths += get_edge_length(i);
    }
    return sum_lengths / (double) m_mesh.m_edges.size();   
}

// --------------------------------------------------------
///
/// Compute average length over edges on non-solid meshes
///
// --------------------------------------------------------

inline double DynamicSurface::get_average_non_solid_edge_length() const
{
    double sum_lengths = 0;
    size_t counted_edges = 0;
    for ( size_t i = 0; i < m_mesh.m_edges.size(); ++i )
    {
        const Vec2st& e = m_mesh.m_edges[i]; 
        if ( e[0] == e[1] )  { continue; }
        if ( edge_is_all_solid(i) ) { continue; }
        sum_lengths += get_edge_length(i);
        ++counted_edges;
    }
    return sum_lengths / (double) counted_edges;   
}

// --------------------------------------------------------
///
/// Compute the surface area
///
// --------------------------------------------------------

inline double DynamicSurface::get_surface_area( ) const
{
    double area=0;
    const std::vector<Vec3st>& tris = m_mesh.get_triangles();
    for(size_t t=0; t < tris.size(); ++t )
    {
        if ( tris[t][0] ==  tris[t][1] ) { continue; }
        area += get_triangle_area(t);
    }
    return area;
}
// --------------------------------------------------------
///
/// Compute the surface area using predicted vertex locations
///
// --------------------------------------------------------

inline double DynamicSurface::get_predicted_surface_area( ) const
{
    double area=0;
    const std::vector<Vec3st>& tris = m_mesh.get_triangles();
    for(size_t t=0; t < tris.size(); ++t )
    {
        if ( tris[t][0] ==  tris[t][1] ) { continue; }
        const Vec3d &p0 = get_newposition(tris[t][0]);
        const Vec3d &p1 = get_newposition(tris[t][1], tris[t][0]);
        const Vec3d &p2 = get_newposition(tris[t][2], tris[t][0]);
        area += 0.5 * mag(cross(p1-p0, p2-p0));
    }
    return area;
}

// --------------------------------------------------------
///
/// Compute the volume enclosed by this surface
///
// --------------------------------------------------------

inline double DynamicSurface::get_volume( ) const
{
    static const double inv_six = 1.0/6.0;
    double volume=0;
    const std::vector<Vec3st>& tris = m_mesh.get_triangles();
    for(size_t t=0; t < tris.size(); ++t )
    {
        if ( tris[t][0] == tris[t][1] ) { continue; }
        const Vec3st& tri = tris[t];
        volume += inv_six * triple(get_position(tri[0]), get_position(tri[1], tri[0]), get_position(tri[2], tri[0]));
    }
    return volume;
}

// --------------------------------------------------------
///
/// Compute the volume using predicted vertex locations
///
// --------------------------------------------------------

inline double DynamicSurface::get_predicted_volume( ) const
{
    static const double inv_six = 1.0/6.0;
    double volume=0;
    const std::vector<Vec3st>& tris = m_mesh.get_triangles();
    for(size_t t=0; t < tris.size(); ++t )
    {
        if ( tris[t][0] ==  tris[t][1] ) { continue; }
        const Vec3st& tri = tris[t];
        volume += inv_six * triple(get_newposition(tri[0]), get_newposition(tri[1], tri[0]), get_newposition(tri[2], tri[0]));
    }
    return volume;
}

// --------------------------------------------------------
///
/// Return true if the specified vertex is solid (should be treated as having infinite mass).
///
// --------------------------------------------------------

inline bool DynamicSurface::vertex_is_solid( size_t v, int dof ) const
{
    assert( v < m_masses.size() );
    return ( m_masses[v][dof] == std::numeric_limits<double>::infinity() );
}

inline Vec3c DynamicSurface::vertex_is_solid_3( size_t v ) const
{
    assert( v < m_masses.size() );
    return (Vec3c(vertex_is_solid(v, 0), vertex_is_solid(v, 1), vertex_is_solid(v, 2)));
}

inline bool DynamicSurface::vertex_is_all_solid( size_t v ) const
{
    return (vertex_is_solid(v, 0) && vertex_is_solid(v, 1) && vertex_is_solid(v, 2));
}

inline bool DynamicSurface::vertex_is_any_solid( size_t v ) const
{
    return (vertex_is_solid(v, 0) || vertex_is_solid(v, 1) || vertex_is_solid(v, 2));
}



// --------------------------------------------------------
///
/// Return true if either end vertex of the specified edge is solid (should be treated as having infinite mass).
///
// --------------------------------------------------------

inline bool DynamicSurface::edge_is_any_solid( size_t e ) const
{
    const Vec2st& edge = m_mesh.m_edges[e];
    return ( vertex_is_any_solid(edge[0]) || vertex_is_any_solid(edge[1]) );
}

// --------------------------------------------------------
///
/// Return true if any corner vertex of the specified triangle is solid (should be treated as having infinite mass).
///
// --------------------------------------------------------

inline bool DynamicSurface::triangle_is_any_solid( size_t t ) const
{
    const Vec3st& tri = m_mesh.get_triangle(t);
    return ( vertex_is_any_solid(tri[0]) || vertex_is_any_solid(tri[1]) || vertex_is_any_solid(tri[2]) );
}

// --------------------------------------------------------
///
/// Return true if both end vertices of the specified edge is solid (should be treated as having infinite mass).
///
// --------------------------------------------------------

inline bool DynamicSurface::edge_is_all_solid( size_t e ) const
{
  const Vec2st& edge = m_mesh.m_edges[e];
  return ( vertex_is_all_solid(edge[0]) && vertex_is_all_solid(edge[1]) );
}

// --------------------------------------------------------
///
/// Return true if all corner vertices of the specified triangle is solid (should be treated as having infinite mass).
///
// --------------------------------------------------------

inline bool DynamicSurface::triangle_is_all_solid( size_t t ) const
{
  const Vec3st& tri = m_mesh.get_triangle(t);
  return ( vertex_is_all_solid(tri[0]) && vertex_is_all_solid(tri[1]) && vertex_is_all_solid(tri[2]) );
}

// ---------------------------------------------------------
///
/// Returns the number of vertices in the mesh, including any vertices marked as deleted
///
// ---------------------------------------------------------

inline size_t DynamicSurface::get_num_vertices() const
{
    return pm_positions.size();
}

// ---------------------------------------------------------
///
/// Returns the current positions of a vertex.
///
// ---------------------------------------------------------

inline const Vec3d& DynamicSurface::get_position( size_t index ) const
{
    assert( index < pm_positions.size() );
    return pm_positions[index];   
}

inline Vec3d DynamicSurface::get_position( size_t index, size_t ref_vert ) const
{
    return get_position(get_position(index), get_position(ref_vert));
}
 
inline Vec3d DynamicSurface::get_position( size_t index, const Vec3d & ref_pos ) const
{
    return get_position(get_position(index), ref_pos);
}
    
inline Vec3d DynamicSurface::get_position( const Vec3d & pos, size_t ref_vert ) const
{
    return get_position(pos, get_position(ref_vert));
}

inline Vec3d DynamicSurface::get_position( const Vec3d & pos, const Vec3d & ref_pos ) const
{
    return pos - get_domain_offset(pos + (Vec3d(PBC_DOMAIN_SIZE_X, PBC_DOMAIN_SIZE_Y, PBC_DOMAIN_SIZE_Z) * 0.5 - ref_pos));
}
  
    
// ---------------------------------------------------------
///
/// Returns the domain offset of a position
///
// ---------------------------------------------------------
inline Vec3d DynamicSurface::get_domain_offset(const Vec3d & pos) const
{
    // essentially dividing each component by the corresponding PBC_DOMAIN_SIZE_? and taking the floor
    Vec3d offset;
    offset[0] = std::floor(pos[0] / PBC_DOMAIN_SIZE_X) * PBC_DOMAIN_SIZE_X;
    offset[1] = std::floor(pos[1] / PBC_DOMAIN_SIZE_Y) * PBC_DOMAIN_SIZE_Y;
    offset[2] = std::floor(pos[2] / PBC_DOMAIN_SIZE_Z) * PBC_DOMAIN_SIZE_Z;
    
    return offset;
}
    
// ---------------------------------------------------------
///
/// Returns the set of all current vertex positions.
///
// ---------------------------------------------------------

inline const std::vector<Vec3d>& DynamicSurface::get_positions( ) const
{
    return pm_positions;
}

// ---------------------------------------------------------
///
/// Set the current position of an individual vertex.
///
// ---------------------------------------------------------

inline void DynamicSurface::set_position( size_t index, const Vec3d& x )
{
    assert( index < pm_positions.size() );
    pm_positions[index] = x - get_domain_offset(x);
    
    // update broad phase
    if ( m_collision_safety )
    {
        update_continuous_broad_phase( index );
    }
}

// ---------------------------------------------------------
///
/// Set the current positions of all vertices in the mesh.
///
// ---------------------------------------------------------

inline void DynamicSurface::set_all_positions( const std::vector<Vec3d>& xs )
{
    pm_positions = xs;
    pm_newpositions = xs;
    
    // update broad phase
    if ( m_collision_safety )
    {
        rebuild_continuous_broad_phase();
    }
}

// ---------------------------------------------------------
///
/// Set the current positions of all vertices in the mesh, from a C-array of doubles
///
// ---------------------------------------------------------

inline void DynamicSurface::set_all_positions( size_t n, const double* xs )
{
    pm_positions.resize(n);
    for ( size_t i = 0; i < n; ++i )
    {
        pm_positions[i][0] = xs[3*i+0];
        pm_positions[i][1] = xs[3*i+1];
        pm_positions[i][2] = xs[3*i+2];
    }
    
    pm_newpositions = pm_positions;
    
    // update broad phase
    if ( m_collision_safety )
    {
        rebuild_continuous_broad_phase();
    }
}

// ---------------------------------------------------------
///
/// Copy predicted vertex positions into the current positions
///
// ---------------------------------------------------------

inline void DynamicSurface::set_positions_to_newpositions()
{
    pm_positions = pm_newpositions;
    
    if ( m_collision_safety )
    {
        rebuild_continuous_broad_phase();
    }
}

// ---------------------------------------------------------
///
/// Returns the predicted position of a vertex.
///
// ---------------------------------------------------------

inline const Vec3d& DynamicSurface::get_newposition( size_t index ) const
{
    assert( index < pm_newpositions.size() );
    return pm_newpositions[index];   
}

inline Vec3d DynamicSurface::get_newposition( size_t index, size_t ref_vert ) const
{
    const Vec3d & ref_base = get_newposition(ref_vert);
    const Vec3d & voi_base = get_newposition(index);
    Vec3d voi_pos = voi_base;
    if      (voi_pos[0] < ref_base[0] - 0.5 * PBC_DOMAIN_SIZE_X) voi_pos[0] += PBC_DOMAIN_SIZE_X;
    else if (voi_pos[0] > ref_base[0] + 0.5 * PBC_DOMAIN_SIZE_X) voi_pos[0] -= PBC_DOMAIN_SIZE_X;
    if      (voi_pos[1] < ref_base[1] - 0.5 * PBC_DOMAIN_SIZE_Y) voi_pos[1] += PBC_DOMAIN_SIZE_Y;
    else if (voi_pos[1] > ref_base[1] + 0.5 * PBC_DOMAIN_SIZE_Y) voi_pos[1] -= PBC_DOMAIN_SIZE_Y;
    if      (voi_pos[2] < ref_base[2] - 0.5 * PBC_DOMAIN_SIZE_Z) voi_pos[2] += PBC_DOMAIN_SIZE_Z;
    else if (voi_pos[2] > ref_base[2] + 0.5 * PBC_DOMAIN_SIZE_Z) voi_pos[2] -= PBC_DOMAIN_SIZE_Z;
    
    return voi_pos;
}

inline Vec3d DynamicSurface::get_newposition( size_t index, const Vec3d & ref_pos ) const
{
    Vec3d ref_offset = get_domain_offset(ref_pos);
    const Vec3d & ref_base = ref_pos - ref_offset;
    const Vec3d & voi_base = get_newposition(index);
    Vec3d voi_pos = voi_base;
    if      (voi_pos[0] < ref_base[0] - 0.5 * PBC_DOMAIN_SIZE_X) voi_pos[0] += PBC_DOMAIN_SIZE_X;
    else if (voi_pos[0] > ref_base[0] + 0.5 * PBC_DOMAIN_SIZE_X) voi_pos[0] -= PBC_DOMAIN_SIZE_X;
    if      (voi_pos[1] < ref_base[1] - 0.5 * PBC_DOMAIN_SIZE_Y) voi_pos[1] += PBC_DOMAIN_SIZE_Y;
    else if (voi_pos[1] > ref_base[1] + 0.5 * PBC_DOMAIN_SIZE_Y) voi_pos[1] -= PBC_DOMAIN_SIZE_Y;
    if      (voi_pos[2] < ref_base[2] - 0.5 * PBC_DOMAIN_SIZE_Z) voi_pos[2] += PBC_DOMAIN_SIZE_Z;
    else if (voi_pos[2] > ref_base[2] + 0.5 * PBC_DOMAIN_SIZE_Z) voi_pos[2] -= PBC_DOMAIN_SIZE_Z;
    
    voi_pos += ref_base;
    return voi_pos;
}


// ---------------------------------------------------------
///
/// Set the predicted position of an individual vertex.
///
// ---------------------------------------------------------

inline void DynamicSurface::set_newposition( size_t index, const Vec3d& x )
{
    assert( index < pm_newpositions.size() );
    
    pm_newpositions[index] = x - get_domain_offset(x);
    
    // update broad phase
    if ( m_collision_safety )
    {
        update_continuous_broad_phase( index );
    }
    
}

// ---------------------------------------------------------
///
/// Set the predicted positions of all vertices in the mesh.
///
// ---------------------------------------------------------

inline void DynamicSurface::set_all_newpositions( const std::vector<Vec3d>& xs )
{
    pm_newpositions = xs;
    
    // update broad phase
    if ( m_collision_safety )
    {
        rebuild_continuous_broad_phase();
    }
}

// ---------------------------------------------------------
///
/// Set the predicted positions of all vertices in the mesh from a C-array.
///
// ---------------------------------------------------------

inline void DynamicSurface::set_all_newpositions( size_t n, const double* xs )
{
    pm_newpositions.resize(n);
    for ( size_t i = 0; i < n; ++i )
    {
        pm_newpositions[i][0] = xs[3*i+0];
        pm_newpositions[i][1] = xs[3*i+1];
        pm_newpositions[i][2] = xs[3*i+2];
    }
    
    // update broad phase
    if ( m_collision_safety )
    {
        rebuild_continuous_broad_phase();
    }
}

// ---------------------------------------------------------
///
/// Returns the set of all predicted vertex positions.
///
// ---------------------------------------------------------

inline const std::vector<Vec3d>& DynamicSurface::get_newpositions( ) const
{
    return pm_newpositions;
}

///////////////////////////////////////////////////////////////////////
// FD 20130102

// Vertex velocity, for remeshing only, not for collisions
inline void DynamicSurface::set_all_remesh_velocities( const std::vector<Vec3d> & v )
{
  pm_velocities = v;
}

inline Vec3d DynamicSurface::get_remesh_velocity( size_t n )
{
  assert( n < pm_velocities.size() );
  return pm_velocities[n];
}

inline void DynamicSurface::set_remesh_velocity( size_t n, const Vec3d & v )
{
    assert( n < pm_velocities.size() );
    pm_velocities[n] = v;
}

///////////////////////////////////////////////////////////////////////

// --------------------------------------------------------
///
/// Returns true if the specified edge is intersecting the specified triangle
///
// --------------------------------------------------------

inline bool check_edge_triangle_intersection_by_index( size_t edge_a, 
                                                      size_t edge_b, 
                                                      size_t triangle_a, 
                                                      size_t triangle_b, 
                                                      size_t triangle_c,
                                                      DynamicSurface * ds,
                                                      bool usenewpos,
                                                      bool verbose )
{
    if (    edge_a == triangle_a || edge_a == triangle_b || edge_a == triangle_c 
        || edge_b == triangle_a || edge_b == triangle_b || edge_b == triangle_c )
    {
        return false;
    }
    
    static const bool DEGEN_COUNTS_AS_INTERSECTION = true;
    
    if (usenewpos)
    {
        Vec3d ref = ds->get_position(edge_a);
        segment_triangle_intersection(ds->get_position(edge_a, ref), edge_a,
                                      ds->get_position(edge_b, ref), edge_b,
                                      ds->get_position(triangle_a, ref), triangle_a,
                                      ds->get_position(triangle_b, ref), triangle_b,
                                      ds->get_position(triangle_c, ref), triangle_c,
                                      DEGEN_COUNTS_AS_INTERSECTION, verbose);
    } else
    {
        Vec3d ref = ds->get_newposition(edge_a);
        segment_triangle_intersection(ds->get_newposition(edge_a, ref), edge_a,
                                      ds->get_newposition(edge_b, ref), edge_b,
                                      ds->get_newposition(triangle_a, ref), triangle_a,
                                      ds->get_newposition(triangle_b, ref), triangle_b,
                                      ds->get_newposition(triangle_c, ref), triangle_c,
                                      DEGEN_COUNTS_AS_INTERSECTION, verbose);
    }
    
}


// --------------------------------------------------------
///
/// Returns true if the an edge from one of the triangles intersects the other triangle.
/// NOTE: Using this routine will produce duplicate checks.  Better to use check_edge_triangle_intersection where possible.
///
// --------------------------------------------------------

inline bool check_triangle_triangle_intersection( Vec3st triangle_a, 
                                                 Vec3st triangle_b,
                                                 DynamicSurface * ds,
                                                 bool usenewpos )
{
    if ( triangle_a[0] == triangle_a[1] || triangle_b[0] == triangle_b[1] )    
    { 
        return false; 
    }
    
    if ( check_edge_triangle_intersection_by_index( triangle_a[0], triangle_a[1], 
                                                   triangle_b[0], triangle_b[1], triangle_b[2],
                                                   ds, usenewpos, false ) )
    {
        return true;
    }
    
    if ( check_edge_triangle_intersection_by_index( triangle_a[1], triangle_a[2], 
                                                   triangle_b[0], triangle_b[1], triangle_b[2], 
                                                   ds, usenewpos, false ) )
    {
        return true;
    }
    
    if ( check_edge_triangle_intersection_by_index( triangle_a[2], triangle_a[0], 
                                                   triangle_b[0], triangle_b[1], triangle_b[2], 
                                                   ds, usenewpos, false ) )
    {
        return true;
    }
    
    if ( check_edge_triangle_intersection_by_index( triangle_b[0], triangle_b[1], 
                                                   triangle_a[0], triangle_a[1], triangle_a[2], 
                                                   ds, usenewpos, false ) )
    {
        return true;
    }
    
    if ( check_edge_triangle_intersection_by_index( triangle_b[1], triangle_b[2], 
                                                   triangle_a[0], triangle_a[1], triangle_a[2], 
                                                   ds, usenewpos, false ) )
    {
        return true;
    }
    
    if ( check_edge_triangle_intersection_by_index( triangle_b[2], triangle_b[0], 
                                                   triangle_a[0], triangle_a[1], triangle_a[2], 
                                                   ds, usenewpos, false ) )
    {
		return true;
    }
    
    return false;
}

}

#endif


