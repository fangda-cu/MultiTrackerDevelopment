// ---------------------------------------------------------
//
//  meancurvature.h
//  Tyson Brochu 2008
//
//  Mesh driver for motion in the normal direction scaled by mean curvature
//
// ---------------------------------------------------------

#ifndef EL_TOPO_MEANCURVATURE_H
#define EL_TOPO_MEANCURVATURE_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include "meshdriver.h"

#include <vec.h>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Motion in the normal direction scaled by mean curvature.
///
// ---------------------------------------------------------

class MeanCurvatureMultiDriver : public MeshDriver
{
public:
    /// Constructor specifying speed of flow and final signed distance field to compare against.
    /// 
    MeanCurvatureMultiDriver(double in_curvature_multiplier);
    
    /// Compute MC * normal at a vertex
    ///
    void add_triangle_contribution_of_mean_curvature_normal(size_t triangle_index, const ElTopo::SurfTrack & surf, ElTopo::Vec3d & vert0, ElTopo::Vec3d & vert1, ElTopo::Vec3d & vert2);
    
    /// Set velocities on each mesh vertex
    ///
    void set_predicted_vertex_positions(const ElTopo::SurfTrack & surf, std::vector<ElTopo::Vec3d> & predicted_positions, double current_t, double & adaptive_dt);
    
    /// Speed of motion
    ///
    double curvature_multiplier;
    
};

// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------


#endif
