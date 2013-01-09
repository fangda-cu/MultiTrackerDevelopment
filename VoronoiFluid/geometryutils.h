
#ifndef GEOMETRYUTILS_H
#define GEOMETRYUTILS_H

#include <vec.h>

// ---------------------------------------------------------
///
/// Estimate the fraction of the polygon which lies above phi = 0.
///
// ---------------------------------------------------------

float compute_polygon_weight( const std::vector<ElTopo::Vec3f>& polygon_vertices, 
                              const std::vector<float>& polygon_phis,
                              const ElTopo::Vec3f& polygon_barycentre,
                              float phi_at_barycentre );

inline float tet_signed_volume( const ElTopo::Vec3f& a, const ElTopo::Vec3f& b, const ElTopo::Vec3f& c, const ElTopo::Vec3f& d );

ElTopo::Vec4f tet_barycentric_weights( const ElTopo::Vec3f& point, const ElTopo::Vec3f& a, const ElTopo::Vec3f& b, const ElTopo::Vec3f& c, const ElTopo::Vec3f& d );
ElTopo::Vec4f tet_barycentric_weights_careful( const ElTopo::Vec3f& point, const ElTopo::Vec3f& a, const ElTopo::Vec3f& b, const ElTopo::Vec3f& c, const ElTopo::Vec3f& d );

ElTopo::Vec3f tet_circumcentre( const ElTopo::Vec3f& a, const ElTopo::Vec3f& b, const ElTopo::Vec3f& c, const ElTopo::Vec3f& d );

bool point_in_tet( const ElTopo::Vec3f &p, const ElTopo::Vec3f &x1, const ElTopo::Vec3f &x2, const ElTopo::Vec3f &x3, const ElTopo::Vec3f &x4, float epsilon );


// ---------------------------------------------------------
///
/// Signed distance of a box (negative inside).
///
// ---------------------------------------------------------

float solid_box_phi( const ElTopo::Vec3f& box_centre, const ElTopo::Vec3f& box_extents, const ElTopo::Vec3f& eval_point );

ElTopo::Vec3f solid_box_gradient( const ElTopo::Vec3f& box_centre, const ElTopo::Vec3f& box_extents, const ElTopo::Vec3f& eval_point );



// ---------------------------------------------------------
// Inline functions
// ---------------------------------------------------------

inline float tet_signed_volume( const ElTopo::Vec3f& a, const ElTopo::Vec3f& b, const ElTopo::Vec3f& c, const ElTopo::Vec3f& d )
{
   static const float OVER_SIX = 1.0f / 6.0f;
   return triple(b-a, c-a, d-a) * OVER_SIX;
}

#endif
