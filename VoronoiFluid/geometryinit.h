// ---------------------------------------------------------
//
//  geometryinit.h
//  Tyson Brochu 2008
//
//  A set of geometry initialization functions.
//
// ---------------------------------------------------------

#ifndef GEOMETRYINIT_H
#define GEOMETRYINIT_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <vec.h>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------
namespace ElTopo{
class SurfTrack;
}

// ---------------------------------------------------------
//  Interface declarations
// ---------------------------------------------------------


/// Create a circular, open surface
///
void create_circle( std::vector<ElTopo::Vec3d>& verts, 
                    std::vector<ElTopo::Vec3st>& tris,
                    std::vector<double>& masses,
                    const ElTopo::Vec3d& centre,
                    double radius,
                    unsigned int nx );

/// Create a vertical or horizontal, triangulated sheet.
///
void create_sheet( std::vector<ElTopo::Vec3d>& verts, 
                   std::vector<ElTopo::Vec3st>& tris, 
                   const ElTopo::Vec3d& low_corner, 
                   const ElTopo::Vec3d& plane_normal, 
                   double width, 
                   unsigned int nx,
                   unsigned int ny );

/// Create a vertical or horizontal, triangulated sheet, curved like a curtain
///
void create_curved_sheet( std::vector<ElTopo::Vec3d>& verts, 
                          std::vector<ElTopo::Vec3st>& tris, 
                          const ElTopo::Vec3d& low_corner, 
                          const ElTopo::Vec3d& plane_normal, 
                          double width, 
                          unsigned int nx );

/// Create a rectangular cubic surface mesh
///
void create_cube( const ElTopo::Vec3d& cube_low,
                  const ElTopo::Vec3d& cube_high,
                  const ElTopo::Vec3st& resolution,
                  std::vector<ElTopo::Vec3d>& verts, 
                  std::vector<ElTopo::Vec3st>& tris );


/// 
///
void create_sphere( const ElTopo::Vec3d& sphere_centre,
                    double sphere_radius,
                    double dx,
                    std::vector<ElTopo::Vec3d>& verts, 
                    std::vector<ElTopo::Vec3st>& tris );



//
// Helpers
//


/// 
///
void project_to_exact_cube( std::vector<ElTopo::Vec3d>& verts, 
                            const ElTopo::Vec3d& cube_low, 
                            const ElTopo::Vec3d& cube_high );

/// Project mesh vertices onto an analytic sphere
///
void project_to_exact_sphere( std::vector<ElTopo::Vec3d>& verts, 
                              const ElTopo::Vec3d& sphere_centre, 
                              double radius );

/// Project mesh vertices onto two analytic spheres
///
void project_to_exact_2_spheres( ElTopo::SurfTrack* surf, 
                                 const ElTopo::Vec3d& sphere_a_centre, 
                                 const ElTopo::Vec3d& sphere_b_centre, 
                                 double max_radius, 
                                 double interior_radius );

/// Project mesh vertices onto analytic dumbbell
///
void project_to_exact_dumbbell( std::vector<ElTopo::Vec3d>& verts, 
                                const ElTopo::Vec3d& sphere_a_centre, 
                                const ElTopo::Vec3d sphere_b_centre, 
                                double sphere_radius, 
                                double handle_width );

/// Read a signed distance field from an ASCII file
///
void read_signed_distance( const char* filename, ElTopo::Array3d& signed_distance );

/// Create signed distance field for a capsule with the given geometry
///
void create_capsule_signed_distance( const ElTopo::Vec3d& capsule_end_a, 
                                     const ElTopo::Vec3d& capsule_end_b, 
                                     double capsule_radius,
                                     double dx,
                                     const ElTopo::Vec3d& domain_low,
                                     const ElTopo::Vec3d& domain_high,                                     
                                     ElTopo::Array3d& phi );

/// Create signed distance field for a cube with the given geometry
///
void create_cube_signed_distance( const ElTopo::Vec3d& cube_low, 
                                  const ElTopo::Vec3d& cube_high, 
                                  double dx,
                                  const ElTopo::Vec3d& domain_low,
                                  const ElTopo::Vec3d& domain_high,                                     
                                  ElTopo::Array3d& phi );

/// Create signed distance field for a sphere with the given geometry
///
void create_sphere_signed_distance( const ElTopo::Vec3d& sphere_centre, 
                                   double sphere_radius, 
                                   double dx,
                                   const ElTopo::Vec3d& domain_low,
                                   const ElTopo::Vec3d& domain_high,                                     
                                   ElTopo::Array3d& phi );

///  Analytic entropy solution to motion in the normal direction of two spheres
///
double signed_distance_entropy( const ElTopo::Vec3d& pt,
                               const ElTopo::Vec3d& sphere_a_centre, 
                               const ElTopo::Vec3d& sphere_b_centre, 
                               double sphere_max_radius,
                               double sphere_interior_radius );

/// Create signed distance field for two spheres
///
void create_two_sphere_signed_distance( const ElTopo::Vec3d& sphere_a_centre, 
                                       const ElTopo::Vec3d& sphere_b_centre, 
                                       double sphere_radius, 
                                       double dx,
                                       const ElTopo::Vec3d& domain_low,
                                       const ElTopo::Vec3d& domain_high,                                     
                                       ElTopo::Array3d& phi );

/// Analytic signed distance function for a dumbbell
///
double signed_distance_dumbbell( const ElTopo::Vec3d& pt,
                                const ElTopo::Vec3d& sphere_a_centre, 
                                const ElTopo::Vec3d& sphere_b_centre, 
                                double sphere_radius, 
                                double handle_width );

/// Create signed distance field for a dumbbell
///
void create_dumbbell_signed_distance( const ElTopo::Vec3d& sphere_a_centre, 
                                     const ElTopo::Vec3d& sphere_b_centre, 
                                     double sphere_radius, 
                                     double handle_width,
                                     double dx,
                                     const ElTopo::Vec3d& domain_low,
                                     const ElTopo::Vec3d& domain_high,                                     
                                     ElTopo::Array3d& phi );

#endif


