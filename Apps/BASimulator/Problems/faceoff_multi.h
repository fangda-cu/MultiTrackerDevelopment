// ---------------------------------------------------------
//
//  faceoff_multi.h
//  Tyson Brochu 2008
//  Christopher Batty 2013 - extensions for multimaterial
//  Fang Da 2013 - modifications to incorporate into BASimulation
//
//  Mesh driver for motion in the normal direction using the faceoff method (entropy solution).
//
// ---------------------------------------------------------

#ifndef BASIM_EL_TOPO_FACEOFF_MULTI_H
#define BASIM_EL_TOPO_FACEOFF_MULTI_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <vector>
#include "ElTopo/common/mat.h"
#include "ElTopo/common/vec.h"
#include "ElTopo/eltopo3d/nondestructivetrimesh.h"
#include "ElTopo/eltopo3d/surftrack.h"

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

namespace ElTopo {
    class NonDestructiveTriMesh;
}

// ---------------------------------------------------------
//  Class definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Motion in the normal direction using Face Offsetting [Jiao 2007].
///
/// Specifically, this driver starts with two spheres, moves in the 
/// positive normal direction, then in the negative normal direction.
/// Error computations are done vs. an analytical solution to this problem.
///
// ---------------------------------------------------------

class FaceOffMultiDriver
{
public:
    
    /// Constructor 
    /// 
    FaceOffMultiDriver( const std::vector< std::vector<double> >& in_speeds, int outer_surf, bool nmf_stationary) : 
    speed_matrix(in_speeds) ,
    expanding_surface(outer_surf),
    nonmanifold_stationary(nmf_stationary),
    do_reverse(false)
    {}
    

    void set_reversing(double at_time) {
       do_reverse = true;
       reverse_time = at_time;
    }

    void set_solution_data(const ElTopo::Vec3d& in_sphere_a_centre, const ElTopo::Vec3d& in_sphere_b_centre, double in_max_radius, double in_interior_radius) {
       sphere_a_centre = in_sphere_a_centre;
       sphere_b_centre = in_sphere_b_centre;
       max_radius      = in_max_radius;
       interior_radius = in_interior_radius;
    }

    /// Get the quadric metric tensor at a vertex from the given incident triangles
    ///   
    void compute_quadric_metric_tensor( const std::vector<ElTopo::Vec3d>& triangle_normals, 
                                       const std::vector<double>& triangle_areas, 
                                       const std::vector<size_t>& incident_triangles,
                                       ElTopo::Mat33d& quadric_metric_tensor );
    
    /// Return intersection point between a set of planes in the least-squares sense
    ///
    void intersection_point( const std::vector<ElTopo::Vec3d>& triangle_normals, 
                            const std::vector<double>& triangle_plane_distances,
                            const std::vector<double>& triangle_areas, 
                            const std::vector<ElTopo::Vec2i>& triangle_labels, 
                            const std::vector<size_t>& incident_triangles,
                            ElTopo::Vec3d& out);
    
    ///  Analytic entropy solution to motion in the normal direction of two spheres
    ///
    double signed_distance_entropy( const ElTopo::Vec3d& pt,
                                   const ElTopo::Vec3d& sphere_a_centre,
                                   const ElTopo::Vec3d& sphere_b_centre,
                                   double sphere_max_radius,
                                   double sphere_interior_radius );
    
    /// Assign a velocity vector to each mesh vertex
    /// 
    void set_predicted_vertex_positions( const ElTopo::SurfTrack& surf, std::vector<ElTopo::Vec3d>& predicted_positions, double current_t, double& adaptive_dt );
    
    /// Compute the distance from vertices to analytic surface, weighted by associated vertex area
    /// 
    double compute_l1_error( const ElTopo::SurfTrack& surf );
    
    /// Compute the maximum distance from vertices to analytic surface
    /// 
    double compute_inf_error( const ElTopo::SurfTrack& surf );
    
    /// Compute and output both L1 and L_inf errors
    ///
    void compute_error( const ElTopo::SurfTrack& surf, double current_t );
        
    /// Speed of normal motion
    double speed;
    std::vector< std::vector<double> > speed_matrix; //dictate speeds based on pairwise labels
    
    /// a fixed time at which to switch directions for testing purposes
    double reverse_time;
    bool do_reverse;

    //
    // The following members are used to compute the analytic solution
    //
    
    /// Geometry to compare against
    ElTopo::Vec3d sphere_a_centre, sphere_b_centre;
    
    /// radius of spheres when motion switches from positive to negative
    double max_radius;         
    
    /// difference between maximum radius and final radius
    double interior_radius; 

    //
    // The following members are used to dictate
    //
    
    /// this flag dictates whether to use null-space smoothing on the interior branch of non-manifold junctions
    /// for shrinking, we do want to null-space smooth considering the interior branch, for expanding we don't.
    bool smooth_based_on_exterior;

    /// (outer) surface region to use for offsetting and null-space smoothing (i.e. for handling non-manifold cases)
    int expanding_surface;

    /// turn on if the non-manifold vertices are not to move (e.g. in curling sphere example)
    bool nonmanifold_stationary;

};


// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------

#endif //BASIM_EL_TOPO_FACEOFF_MULTI_H


