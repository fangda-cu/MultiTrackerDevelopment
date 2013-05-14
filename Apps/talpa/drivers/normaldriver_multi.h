// ---------------------------------------------------------
//
//  normaldriver.h
//  Tyson Brochu 2008
//
//  Mesh driver for motion in the normal direction using vertex normals.
//  NOTE: This was implemented only for comparison against FaceOff.  Not recommended for general use.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_NORMALDRIVER_MULTI_H
#define EL_TOPO_NORMALDRIVER_MULTI_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include "meshdriver.h"

#include <mat.h>
#include <vec.h> 
#include <vector>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------
namespace ElTopo {
class NonDestructiveTriMesh;
class SurfTrack;
}

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Motion in the normal direction using vertex normals 
/// (area-weighted average of incident triangle normals).
///
/// This driver starts with two spheres, moves in the positive
/// normal direction, then in the negative normal direction.
/// Error computations are done vs. an analytical solution.
///
// ---------------------------------------------------------

class NormalDriverMulti : public MeshDriver
{
public:
    
    NormalDriverMulti( const std::vector< std::vector<double> >& in_speeds,  int outer_surface ) : 
    speed_matrix(in_speeds), expanding_surface(outer_surface)
    {}
    
    void initialize( const ElTopo::SurfTrack& ) {}
    
    // Assign a velocity vector to each mesh vertex
    // 
    void set_predicted_vertex_positions( const ElTopo::SurfTrack& surf, std::vector<ElTopo::Vec3d>& predicted_positions, double current_t, double& adaptive_dt );
    
    /// Compute and output both L1 and L_inf errors
    ///
    void compute_error( const ElTopo::SurfTrack& surf, double current_t );
    
    /// dictate speeds based on pairwise labels
    std::vector< std::vector<double> > speed_matrix; 
    
    /// (outer) surface region to use for offsetting and null-space smoothing (i.e. for handling non-manifold cases)
    int expanding_surface;
};


// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------

#endif


