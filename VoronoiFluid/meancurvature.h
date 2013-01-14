// ---------------------------------------------------------
//
//  meancurvature.h
//  Tyson Brochu 2008
//
//  Mesh driver for motion in the normal direction scaled by mean curvature
//
// ---------------------------------------------------------

#ifndef MEANCURVATURE_H
#define MEANCURVATURE_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

//#include <meshdriver.h>
#include <vec.h>
#include "array3.h"

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------
namespace ElTopo{
class SurfTrack;
class DynamicSurface;
}

// ---------------------------------------------------------
//  Interface declarations
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Motion in the normal direction scaled by mean curvature.
///
// ---------------------------------------------------------


class MeanCurvatureDriver //: public MeshDriver
{
public:
      
   /// Constructor specifying speed of flow and final signed distance field to compare against.
   /// 
   MeanCurvatureDriver( double in_curvature_multiplier, 
                        const ElTopo::Array3d& in_final_signed_distance, 
                        const ElTopo::Vec3d& in_final_domain_low, 
                        double in_final_domain_dx );
   
   /// Compute the area of a triangle associated with the specified vertex.  (See [Meyer et al. 2002].)
   ///
   static double mixed_area( unsigned int vertex_index, unsigned int triangle_index, const ElTopo::DynamicSurface& surf );
   
   /// Compute mean curvature times normal at a vertex
   ///
   static void vertex_mean_curvature_normal( unsigned int vertex_index, const ElTopo::DynamicSurface& surf, ElTopo::Vec3d& out );
   
   /// Compute MC * normal at a vertex and return the sum of weights used (for computing the time step restriction)
   ///
   static void vertex_mean_curvature_normal( unsigned int vertex_index, const ElTopo::DynamicSurface& surf, ElTopo::Vec3d& out, double& weight_sum );
   
   /// Compute mean curvature times normal at a vertex, using per-triangle approach that supports non-manifoldness
   ///
   static void vertex_mean_curvature_normal_nonmanifold( unsigned int vertex_index, const ElTopo::DynamicSurface& surf, ElTopo::Vec3d& out );

   /// Set velocities on each mesh vertex
   ///
   void set_surface_velocity( const ElTopo::DynamicSurface& surf, std::vector<ElTopo::Vec3d>& out_velocity, double current_t, double& adaptive_dt );
   
   /// Compute the distance from vertices to analytic surface, weighted by associated vertex area
   /// 
   double compute_l1_error( const ElTopo::DynamicSurface& surf ) const;
   
   /// Compute the maximum distance from vertices to analytic surface
   /// 
   double compute_inf_error( const ElTopo::DynamicSurface& surf ) const;
   
   /// Compute and output both L1 and L_inf errors
   ///
   void compute_error( const ElTopo::DynamicSurface& surf, double current_t );
  

   /// Speed of motion
   ///
   double curvature_multiplier;
   
   //
   // For error computation: what the signed distance function "should" be at the end time
   //
   
   ElTopo::Array3d final_signed_distance;
   ElTopo::Vec3d final_domain_low;
   double final_domain_dx;
      
};


// ---------------------------------------------------------
//  Inline functions
// ---------------------------------------------------------


#endif
