// ---------------------------------------------------------
//
//  meancurvature.h
//  Tyson Brochu 2008
//
//  Mesh driver for motion in the normal direction scaled by mean curvature
//
// ---------------------------------------------------------

#ifndef EL_TOPO_MEANCURVATURE_MULTI_H
#define EL_TOPO_MEANCURVATURE_MULTI_H

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

class MeanCurvatureMultiDriver : public MeshDriver, ElTopo::SurfTrack::SolidVerticesCallback
{
public:
    /// Constructor specifying speed of flow and final signed distance field to compare against.
    /// 
    MeanCurvatureMultiDriver(double in_curvature_multiplier);
    
    void initialize(ElTopo::SurfTrack & st);
    
    /// Compute MC * normal at a vertex
    ///
    void add_triangle_contribution_of_mean_curvature_normal(size_t triangle_index, const ElTopo::SurfTrack & surf, ElTopo::Vec3d & vert0, ElTopo::Vec3d & vert1, ElTopo::Vec3d & vert2);
    
    /// Set velocities on each mesh vertex
    ///
    void set_predicted_vertex_positions(const ElTopo::SurfTrack & surf, std::vector<ElTopo::Vec3d> & predicted_positions, double current_t, double & adaptive_dt);
    
    /// Speed of motion
    ///
    double curvature_multiplier;
    
    // SolidVerticesCallback
    bool generate_collapsed_position(ElTopo::SurfTrack & st, size_t v0, size_t v1, ElTopo::Vec3d & pos);
    bool generate_split_position(ElTopo::SurfTrack & st, size_t v0, size_t v1, ElTopo::Vec3d & pos);
    ElTopo::Vec3c generate_collapsed_solid_label(ElTopo::SurfTrack & st, size_t v0, size_t v1, const ElTopo::Vec3c & label0, const ElTopo::Vec3c & label1);
    ElTopo::Vec3c generate_split_solid_label(ElTopo::SurfTrack & st, size_t v0, size_t v1, const ElTopo::Vec3c & label0, const ElTopo::Vec3c & label1);
    bool generate_edge_popped_positions(ElTopo::SurfTrack & st, size_t oldv, const ElTopo::Vec2i & cut, ElTopo::Vec3d & pos_upper, ElTopo::Vec3d & pos_lower);
    bool generate_vertex_popped_positions(ElTopo::SurfTrack & st, size_t oldv, int A, int B, ElTopo::Vec3d & pos_a, ElTopo::Vec3d & pos_b);
    bool solid_edge_is_feature(const ElTopo::SurfTrack & st, size_t e);

};

// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------


#endif
