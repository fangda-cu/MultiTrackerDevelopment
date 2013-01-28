// ---------------------------------------------------------
//
//  meshsnapper.h
//  Christopher Batty 2013
//  
//  Functions to handle snapping of meshes together when they get close.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_MESHSNAPPER_H
#define EL_TOPO_MESHSNAPPER_H

// ---------------------------------------------------------
//  Nested includes
// ---------------------------------------------------------

#include <cstddef>
#include <vector>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

namespace ElTopo {

class SurfTrack;
template<unsigned int N, class T> struct Vec;
typedef Vec<3,double> Vec3d;
typedef Vec<3,size_t> Vec3st;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Mesh snapper object.  Snaps together vertices that are close together.
///
// ---------------------------------------------------------

class MeshSnapper
{
    
public:
    
    /// Mesh snapper constructor.  Takes a SurfTrack object.
    ///
    MeshSnapper( SurfTrack& surf );

    /// Collapse all proximal vertices
    ///
    bool snap_pass();
    
  
private:
    
    friend class SurfTrack;
    
    /// The mesh this object operates on
    /// 
    SurfTrack& m_surf;

    /// Get all triangles which are incident on either involved vertex.
    ///
    void get_moving_triangles(size_t source_vertex, 
                              size_t destination_vertex, 
                              std::vector<size_t>& moving_triangles );
    
    
    /// Get all edges which are incident on either involved vertex.
    ///
    void get_moving_edges(size_t source_vertex, 
                          size_t destination_vertex, 
                          std::vector<size_t>& moving_edges );
    
    /// Check the "pseudo motion" introduced by snapping for collision
    ///
    bool snap_pseudo_motion_introduces_collision( size_t source_vertex, 
                                                          size_t destination_vertex, 
                                                          const Vec3d& vertex_new_position );
    
    /*
    /// Determine if the snap operation would invert the normal of any incident triangles.
    ///
    bool snap_edge_introduces_normal_inversion( size_t source_vertex, 
                                                   size_t destination_vertex, 
                                                   size_t edge_index, 
                                                   const Vec3d& vertex_new_position );
    
    /// Determine whether snapping will introduce an unacceptable change in volume.
    ///
    bool snap_edge_introduces_volume_change( size_t source_vertex, 
                                                size_t edge_index, 
                                                const Vec3d& vertex_new_position );   
    
    /// Returns true if the snap collapse would introduce a triangle with a min or max angle outside of the specified min or max.
    ///
    bool snap_edge_introduces_bad_angle( size_t source_vertex, 
                                            size_t destination_vertex, 
                                            const Vec3d& vertex_new_position );
    */

    /// Snap a vertex pair by moving both to their average point
    ///
    bool snap_vertex_pair( size_t vert0, size_t vert1);
    
    /// Determine if the edge should be allowed to collapse
    ///
    bool vert_pair_is_snappable( size_t vert0, size_t vert1, double& cur_length );
    
    
};

}

#endif //EL_TOPO_MESHSNAPPER_H


