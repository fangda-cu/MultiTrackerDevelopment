// ---------------------------------------------------------
//
//  surftrack.h
//  Tyson Brochu 2008
//  
//  The SurfTrack class: a dynamic mesh with topological changes and mesh maintenance operations.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_SURFTRACK_H
#define EL_TOPO_SURFTRACK_H

#include <dynamicsurface.h>
#include <edgecollapser.h>
#include <edgeflipper.h>
#include <edgesplitter.h>
#include <meshmerger.h>
#include <meshpincher.h>
#include <meshsmoother.h>
#include <meshcutter.h>
#include <t1transition.h>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

namespace ElTopo {

class SubdivisionScheme;
typedef std::vector<size_t> TriangleSet;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Structure for setting up a SurfTrack object with some initial parameters.  This is passed to the SurfTrack constructor.
///
// ---------------------------------------------------------

struct SurfTrackInitializationParameters
{
    
    ///  Constructor. Sets default values for parameters which are not likely to be specified.
    ///
    SurfTrackInitializationParameters();
    
    /// Elements closer than this are considered "near" (or proximate)
    ///
    double m_proximity_epsilon;
    
    /// Coefficient of friction to apply during collisions
    ///
    double m_friction_coefficient;
    
    /// Smallest triangle area to allow
    ///
    double m_min_triangle_area;
    
    /// Whether to enable T1 transition operations
    ///
    bool m_t1_transition_enabled;
    
    /// The velocity field callback
    ///
    T1Transition::VelocityFieldCallback * m_velocity_field_callback;
    
    /// Collision epsilon to use during mesh improvement operations (i.e. if any mesh elements are closer than this, the operation is 
    /// aborted).  NOTE: This should be greater than collision_epsilon, to prevent improvement operations from moving elements into 
    /// a collision configuration.
    ///
    double m_improve_collision_epsilon;
    
    /// Whether to set the min and max edge lengths as fractions of the initial average edge length
    ///
    bool m_use_fraction;
    
    // If use_fraction is true, the following three values are taken to be fractions of the average edge length of the new surface.
    // If use_fraction is false, these are absolute.
    
    /// Smallest edge length allowed
    ///
    double m_min_edge_length;
    
    /// Longest edge length allowed
    ///
    double m_max_edge_length; 
    
    /// Maximum change in volume allowed for one operation
    ///
    double m_max_volume_change;
    
    /// Smallest interior angle at a triangle vertex allowed
    ///
    double m_min_triangle_angle;
    
    /// Largest interior angle at a triangle vertex allowed
    ///
    double m_max_triangle_angle;   
    
    ///////////////////////////////////////////////////////////////////////
    // FD 20121229
    //
    /// Largest interior angle at a triangle vertex allowed before large-angle split pass kicks in
    ///
    double m_large_triangle_angle_to_split;
    
    ///////////////////////////////////////////////////////////////////////

    /// Whether to scale by curvature when computing edge lengths, in order to refine high-curvature regions
    ///
    bool m_use_curvature_when_splitting;

    /// Whether to scale by curvature when computing edge lengths, in order to coarsen low-curvature regions
    ///
    bool m_use_curvature_when_collapsing;
    
    /// The minimum curvature scaling allowed
    ///
    double m_min_curvature_multiplier;
    
    /// The maximum curvature scaling allowed
    ///
    double m_max_curvature_multiplier;
    
    /// boolean, whether to allow vertices to move during improvement
    int m_allow_vertex_movement_during_collapse;

    /// boolean, whether to allow vertices to move during improvement
    int m_perform_smoothing;
    
    /// Minimum edge length improvement in order to flip an edge
    //
    double m_edge_flip_min_length_change;
    
    /// Elements within this distance will trigger a merge attempt   
    ///
    double m_merge_proximity_epsilon;
    
    /// Type of subdivision to use when collapsing or splitting (butterfly, quadric error minimization, etc.)
    ///
    SubdivisionScheme *m_subdivision_scheme;   
    
    /// Whether to enforce collision-free surfaces (including during mesh maintenance operations)
    ///
    bool m_collision_safety;
    
    /// Whether to allow changes in topology
    ///
    bool m_allow_topology_changes;
    
    /// Whether to allow non-manifold (edges incident on more than two triangles)
    ///
    bool m_allow_non_manifold;
    
    /// Whether to allow mesh improvement
    ///
    bool m_perform_improvement;

    /// Whether to perform remeshing on mesh boundary edges (in the case of open surfaces, e.g. sheets)
    ///
    bool m_remesh_boundaries;

    /// Whether to be verbose in outputting data
    ///
    bool m_verbose;
    
};

// ---------------------------------------------------------
///
/// Used to build a list of edges sorted in order of increasing length.
/// 
// ---------------------------------------------------------

struct SortableEdge
{    
    /// Constructor
    ///
    SortableEdge( size_t ei, double el ) : 
    m_edge_index(ei), 
    m_edge_length(el) 
    {}
    
    /// Comparison operator for sorting
    ///
    bool operator<( const SortableEdge& other ) const
    {
        return (this->m_edge_length < other.m_edge_length);
    }
    
    /// The index of the edge
    ///
    size_t m_edge_index;
    
    /// The stored edge length
    ///
    double m_edge_length;

};


// ---------------------------------------------------------
///
/// Keeps track of a vertex removal or addition.  If it's an addition, it also points to the edge that was split to create it.
///
// ---------------------------------------------------------

struct VertexUpdateEvent
{
    /// Constructor
    ///
    VertexUpdateEvent(bool is_remove = false, 
                      size_t vertex_index = (size_t)~0, 
                      const Vec2st& split_edge = Vec2st((size_t)~0) ) :
    m_is_remove( is_remove ),
    m_vertex_index( vertex_index ),
    m_split_edge( split_edge )
    {}
    
    /// Tag for identifying a vertex removal
    ///
    static const bool VERTEX_REMOVE = true;
    
    /// Tag for identifying a vertex addition
    ///
    static const bool VERTEX_ADD = false;
    
    /// Whether this event is a vertex removal
    ///
    bool m_is_remove;
    
    /// The index of the vertex being added or removed
    ///
    size_t m_vertex_index;   
    
    /// If this is a vertex addition due to edge splitting, the edge that was split
    ///
    Vec2st m_split_edge;
    
};


// ---------------------------------------------------------
///
/// Keeps track of a triangle removal or addition. If addition, contains the three vertices that form the new triangle.
///
// ---------------------------------------------------------

struct TriangleUpdateEvent
{
    /// Constructor
    ///
    TriangleUpdateEvent(bool is_remove = false, 
                        size_t triangle_index = (size_t)~0, 
                        const Vec3st& triangle = Vec3st((size_t)~0) ) :
    m_is_remove( is_remove ),
    m_triangle_index( triangle_index ),
    m_tri( triangle )
    {}
    
    /// Tag for identifying a triangle removal
    ///
    static const bool TRIANGLE_REMOVE = true;
    
    /// Tag for identifying a triangle addition
    ///
    static const bool TRIANGLE_ADD = false;
    
    /// Whether this event is a triangle removal
    ///
    bool m_is_remove;
    
    /// The index of the triangle being added or removed
    ///
    size_t m_triangle_index;  
    
    /// If this is a triangle addition, the triangle added
    ///
    Vec3st m_tri;
    
};


// ---------------------------------------------------------
///
/// Keeps track of a triangle removal or addition. If addition, contains the three vertices that form the new triangle.
///
// ---------------------------------------------------------

struct MeshUpdateEvent
{
  enum EventType {
    EDGE_SPLIT,
    EDGE_FLIP,
    EDGE_COLLAPSE,
    EDGE_CUT,     //for fracturing/cutting
    FLAP_DELETE,  //remove non-manifold flap
    PINCH,        //separate singular vertex to allow topology change
    MERGE,        //zipper two edges together
      
    ////////////////////////////////////////////////////////////
    // FD 20130109
    //
    // T1 transition events
    
    EDGE_POP,
    VERTEX_POP,
    
    ////////////////////////////////////////////////////////////

  };

  /// Constructors
  ///
  MeshUpdateEvent(EventType eType):m_type(eType), m_deleted_tris(0), m_created_tris(0), m_created_tri_data(0)
  {}

  /// What type of mesh event this is
  ///
  EventType m_type;

  /// The start and end vertices of the edge
  size_t m_v0, m_v1;
  
  // Another identifying vertex, needed for internal cuts
  size_t m_v2; 
    
  // A fourth vertex needed only for zippering
  size_t m_v3;

  /// The index of the triangles involved. 
  ///
  std::vector<size_t> m_deleted_tris;
  std::vector<size_t> m_created_tris;

  /// The data of the new triangles (vertex indices)
  ///
  std::vector<Vec3st> m_created_tri_data;
  
  ////////////////////////////////////////////////////////////
  // FD 20121126
  //
  
  /// The label data of the new triangles
  ///
  std::vector<Vec2i> m_created_tri_labels;
  
  /// Dirty triangles whose labels have been changed
  ///
  std::vector<std::pair<size_t, Vec2i> > m_dirty_tris;

  ////////////////////////////////////////////////////////////

  /// The indices of the vertices involved. 
  ///
  std::vector<size_t> m_deleted_verts;
  std::vector<size_t> m_created_verts;
  
  /// The positions of the created vertices
  ///
  std::vector<Vec3d> m_created_vert_data;

  /// The location of the final vertex (for a split, collapse, or possibly smooth)
  ///
  Vec3d m_vert_position;


};

// ---------------------------------------------------------
///
/// A DynamicSurface with topological and mesh maintenance operations.
///
// ---------------------------------------------------------

class SurfTrack : public DynamicSurface
{
    
public:
    
    /// Create a SurfTrack object from a set of vertices and triangles using the specified parameters
    ///
    SurfTrack(const std::vector<Vec3d>& vs, 
              const std::vector<Vec3st>& ts,
              const std::vector<Vec2i>& labels,
              const std::vector<double>& masses,
              const SurfTrackInitializationParameters& initial_parameters );
    
    /// Destructor
    ///
    ~SurfTrack();
    
private:
    
    /// Disallow copying and assignment by declaring private
    ///
    SurfTrack( const SurfTrack& );
    
    /// Disallow copying and assignment by declaring private
    ///
    SurfTrack& operator=( const SurfTrack& );
    
    
public:
    

    //
    // Mesh bookkeeping
    //
    
    /// Add a triangle to the surface.  Update the underlying TriMesh and acceleration grid. 
    ///
    size_t add_triangle(const Vec3st& t);
    
    /// Remove a triangle from the surface.  Update the underlying TriMesh and acceleration grid. 
    ///
    void remove_triangle(size_t t);  
    
    /// Add a vertex to the surface.  Update the acceleration grid. 
    ///
    size_t add_vertex( const Vec3d& new_vertex_position, double new_vertex_mass );
    
    /// Remove a vertex from the surface.  Update the acceleration grid. 
    ///
    void remove_vertex(size_t v);
    
    /// Remove deleted vertices and triangles from the mesh data structures
    ///
    void defrag_mesh();

    /// Check for labels with -1 as their value, or the same label on both sides.
    /// 
    void assert_no_bad_labels();

    //
    // Main operations
    //
    
    /// Run mesh maintenance operations
    ///
    void improve_mesh( );
    
    /// Run edge-edge merging
    ///
    void topology_changes( );
    
    /// Run mesh cutting operations on a given set of edges
    ///
    void cut_mesh( const std::vector< std::pair<size_t, size_t> >& edges);


    //
    // Mesh cleanup
    //
    
    /// Check for and delete flaps and zero-area triangles among the given triangle indices, then separate singular vertices.
    ///
    void trim_non_manifold( std::vector<size_t>& triangle_indices );
    
    /// Check for and delete flaps and zero-area triangles among *all* triangles, then separate singular vertices.
    ///
    inline void trim_non_manifold();
    
    /// Fire an assert if any degenerate triangles or tets (flaps) are found.
    /// 
    void assert_no_degenerate_triangles();
    
    //
    // Member variables
    //
    
    /// Edge collapse operation object
    ///
    EdgeCollapser m_collapser;
    
    /// Edge split operation object
    ///
    EdgeSplitter m_splitter;
    
    /// Edge flip operation object
    ///
    EdgeFlipper m_flipper;
    
    /// NULL-space surface smoothing
    /// 
    MeshSmoother m_smoother;
    
    /// Surface merging object
    ///
    MeshMerger m_merger;
    
    /// Surface splitting operation object
    ///
    MeshPincher m_pincher;
    
    /// Surface cutting (tearing) operation object
    ///
    MeshCutter m_cutter;
    
    ////////////////////////////////////////////////////////////
    // FD 20130109
    
    /// T1 transition operation object
    ///
    T1Transition m_t1transition;
    
    /// An option to indicate whether T1 transition operations are enabled
    ///
    bool m_t1_transition_enabled;

    ////////////////////////////////////////////////////////////

    /// Collision epsilon to use during mesh improvement operations
    ///
    double m_improve_collision_epsilon;
    
    /// Minimum edge length improvement in order to flip an edge
    ///
    double m_edge_flip_min_length_change;
    
    /// Maximum volume change allowed when flipping or collapsing an edge
    ///
    double m_max_volume_change;
    
    /// Minimum edge length.  Edges shorter than this will be collapsed.
    ///
    double m_min_edge_length;   
    
    /// Maximum edge length.  Edges longer than this will be subdivided.
    ///
    double m_max_edge_length;   
    
    /// Elements within this distance will trigger a merge attempt
    ///
    double m_merge_proximity_epsilon;
    
    /// Try to prevent triangles with area less than this
    ///
    double m_min_triangle_area;
    
    /// Don't create triangles with angles less than this.  If angles less than this do exist, try to remove them.
    ///
    double m_min_triangle_angle;
    
    /// Don't create triangles with angles greater than this.  If angles greater than this do exist, try to remove them.
    ///
    double m_max_triangle_angle;
    
    ///////////////////////////////////////////////////////////////////////
    // FD 20121229
    //
    /// Split triangles with angles greater than this.
    ///
    double m_large_triangle_angle_to_split;
    
    ///////////////////////////////////////////////////////////////////////
    
    /// Interpolation scheme, determines edge midpoint location
    ///
    SubdivisionScheme *m_subdivision_scheme;
    
    /// If we allocate our own SubdivisionScheme object, we must delete it in this object's deconstructor.
    ///
    bool should_delete_subdivision_scheme_object;
    
    /// Triangles which are involved in connectivity changes which may introduce degeneracies
    ///
    std::vector<size_t> m_dirty_triangles;
    
    /// Whether to allow merging and separation
    ///
    bool m_allow_topology_changes;
    
    /// Whether to allow non-manifold (edges incident on more than two triangles)
    ///
    bool m_allow_non_manifold;
    
    /// Whether to perform adaptivity operations
    ///
    bool m_perform_improvement;

    /// Whether to perform remeshing on mesh boundary edges (in the case of open surfaces, e.g. sheets)
    ///
    bool m_remesh_boundaries;
    

    /// boolean, whether to allow vertices to move during collapses (i.e. use points other than the endpoints)
    int m_allow_vertex_movement_during_collapse;

    /// boolean, whether to do null space smoothing on vertex positions
    int m_perform_smoothing;
    
    
    ///////////////////////////////////////////////////////////////////////
    // FD 20121229
    
    class ConstrainedVerticesCallback
    {
    public:
        virtual bool generate_collapsed_position(SurfTrack & st, size_t v0, size_t v1, Vec3d & pos) = 0;
        
        virtual bool generate_splitted_position(SurfTrack & st, size_t v0, size_t v1, Vec3d & pos) = 0;
        
        virtual bool generate_collapsed_constraint_label(SurfTrack & st, size_t v0, size_t v1, bool label0, bool label1) = 0;
        
        virtual bool generate_splitted_constraint_label(SurfTrack & st, size_t v0, size_t v1, bool label0, bool label1) = 0;
        
        virtual bool generate_edge_popped_positions(SurfTrack & st, size_t oldv, const Vec2i & cut, Vec3d & pos_upper, Vec3d & pos_lower) = 0;
        
        virtual bool generate_vertex_popped_positions(SurfTrack & st, size_t oldv, int A, int B, Vec3d & pos_a, Vec3d & pos_b) = 0;
        
    };
    
    ConstrainedVerticesCallback * m_constrained_vertices_callback;

    ///////////////////////////////////////////////////////////////////////
        
    /// History of vertex removal or addition events
    ///
    std::vector<VertexUpdateEvent> m_vertex_change_history;
    
    /// History of triangle removal or addition events
    ///    
    std::vector<TriangleUpdateEvent> m_triangle_change_history;
    
    /// Map of triangle indices, mapping pre-defrag triangle indices to post-defrag indices
    ///
    std::vector<Vec2st> m_defragged_triangle_map;
    
    /// Map of vertex indices, mapping pre-defrag vertex indices to post-defrag indices
    ///
    std::vector<Vec2st> m_defragged_vertex_map;

    /// History of higher level mesh update events (split, flip, collapse, smooth)
    ///    
    std::vector<MeshUpdateEvent> m_mesh_change_history;

    
};

// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Search the entire mesh for non-manifold elements and remove them
/// NOTE: SHOULD USE THE VERSION THAT ACCEPTS A SET OF TRIANGLE INDICES INSTEAD.
///
// ---------------------------------------------------------

inline void SurfTrack::trim_non_manifold()
{
    
    std::vector<size_t> triangle_indices;
    triangle_indices.resize( m_mesh.num_triangles() );
    for ( size_t i = 0; i < triangle_indices.size(); ++i )
    {
        triangle_indices[i] = i;
    }
    
    trim_non_manifold( triangle_indices );
}

}

#endif

