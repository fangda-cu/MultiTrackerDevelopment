// ---------------------------------------------------------
//
//  t1transition.h
//  Fang Da 2013
//  
//  Functions handling T1 transitions (edge popping and vertex popping).
//
// ---------------------------------------------------------

#ifndef EL_TOPO_T1TRANSITION_H
#define EL_TOPO_T1TRANSITION_H

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
typedef Vec<2,size_t> Vec2st;
typedef Vec<3,size_t> Vec3st;

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// T1Transition class. Pull apart X-junction edges and X-junction vertices
///
// ---------------------------------------------------------

class T1Transition
{
    
public:
    
    /// Constructor
    ///
    T1Transition(SurfTrack & surf, bool remesh_boundaries);
    
    /// Perform edge popping (first step of T1 transition)
    ///
    bool pop_edges();
    
    /// Perform vertex popping (second step of T1 transition)
    ///
    bool pop_vertices();
    
    /// Whether or not to remesh the boundary
    ///
    bool m_remesh_boundaries;
    
private:
    
    /// The mesh this object operates on
    /// 
    SurfTrack & m_surf;   
    
//    /// Determine if the pseudo-trajectories of the new vertices have a collision with the existing mesh.
//    ///
//    bool t1_transition_pseudo_motion_introduces_intersection(const Vec3d& new_vertex_position, 
//                                                             const Vec3d& new_vertex_smooth_position, 
//                                                             size_t edge,
//                                                             size_t vertex_a,
//                                                             size_t vertex_b,
//                                                             std::vector<size_t>& tris, 
//                                                             std::vector<size_t>& verts);

};

}

#endif
