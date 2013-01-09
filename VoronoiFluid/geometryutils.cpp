/*
 *  geometryutils.cpp
 *  voronoifluid3d_project
 *
 *  Created by tyson on 10/12/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "geometryutils.h"
#include <mat.h>

using namespace ElTopo;

// ---------------------------------------------------------
///
/// Helper for compute_positive_area.
/// Determine the phi == 0 point on a line segment, given values of phi on the segment's end point.
///
// ---------------------------------------------------------

static Vec3f crossing_point( const Vec3f& x0, const Vec3f& x1, float p0, float p1 )
{
   assert( ( p0 >= 0 && p1 < 0 ) || ( p0 < 0 && p1 >= 0 ) );
   
   float a=p1/(p1-p0);

   assert( a == a );
   assert( a >= 0 );
   assert( a <= 1 );
   
   return a*x0 + (1.0f-a)*x1;
}

// ---------------------------------------------------------
///
/// Helper for compute_polygon_weight.
/// Estimates the fraction of a triangle lying above phi == 0 from the phi values on its corners.
///
// ---------------------------------------------------------

static float compute_positive_area( const Vec3f& a, const Vec3f& b, const Vec3f& c, float phi_a, float phi_b, float phi_c )
{
   std::vector<Vec3f> clipped_polygon;
   
   if(phi_a<0)
   {
      if(phi_b<0)
      {
         if(phi_c<0) { return 0.0f; }
         else 
         {
            // bc, c, ca
            clipped_polygon.push_back( crossing_point( b, c, phi_b, phi_c ) );
            clipped_polygon.push_back( c );
            clipped_polygon.push_back( crossing_point( c, a, phi_c, phi_a ) );
         }
      }
      else
      { 
         // phi_b > 0
         if ( phi_c < 0) 
         { 
            // ab, b, bc
            clipped_polygon.push_back( crossing_point( a, b, phi_a, phi_b ) );
            clipped_polygon.push_back( b );
            clipped_polygon.push_back( crossing_point( b, c, phi_b, phi_c ) );
         }
         else 
         { 
            // ab b c ca
            clipped_polygon.push_back( crossing_point( a, b, phi_a, phi_b ) );
            clipped_polygon.push_back( b );               
            clipped_polygon.push_back( c );               
            clipped_polygon.push_back( crossing_point( c, a, phi_c, phi_a ) );               
         }
      }
   }
   else
   { 
      // phi_a > 0
      
      if( phi_b < 0)
      {
         if(phi_c < 0) 
         { 
            // a ab ca
            clipped_polygon.push_back( a );                 
            clipped_polygon.push_back( crossing_point( a, b, phi_a, phi_b ) );               
            clipped_polygon.push_back( crossing_point( c, a, phi_c, phi_a ) );               
         }
         else 
         { 
            // a ab bc c
            clipped_polygon.push_back( a );                 
            clipped_polygon.push_back( crossing_point( a, b, phi_a, phi_b ) );               
            clipped_polygon.push_back( crossing_point( b, c, phi_b, phi_c ) );               
            clipped_polygon.push_back( c );                 
         }
      }
      else
      { 
         // phi_b > 0
         if( phi_c < 0 ) 
         { 
            // a b bc ca
            clipped_polygon.push_back( a );                 
            clipped_polygon.push_back( b );   
            clipped_polygon.push_back( crossing_point( b, c, phi_b, phi_c ) );               
            clipped_polygon.push_back( crossing_point( c, a, phi_c, phi_a ) );               
         }
         else 
         { 
            // all positive
            return 0.5f * mag( cross( b - a, c - a ) ); 
         }
      }
   }
   
   if ( clipped_polygon.size() == 3 )
   {
      return 0.5f * mag( cross( clipped_polygon[1] - clipped_polygon[0], clipped_polygon[2] - clipped_polygon[0] ) ); 
   }
   
   if ( clipped_polygon.size() == 4 )
   {
      float area = 0.5f * mag( cross( clipped_polygon[0] - clipped_polygon[1], clipped_polygon[2] - clipped_polygon[1] ) ) 
                 + 0.5f * mag( cross( clipped_polygon[0] - clipped_polygon[3], clipped_polygon[2] - clipped_polygon[3] ) );
      
      return area;
   }
   
   assert( !"Clipped polygon doesn't have 3 or 4 points\n" );
   return 0;
}


// ---------------------------------------------------------
///
/// Estimate the fraction of the polygon which lies above phi = 0.  The algorithm triangulates the polygon using its
/// barycentre, then estimates the positive-phi fraction for each triangle by clipping it and computing the area of the
/// clipped polygon.  Whew.
///
// ---------------------------------------------------------

float compute_polygon_weight( const std::vector<Vec3f>& polygon_vertices, 
                             const std::vector<float>& polygon_phis,
                             const Vec3f& polygon_barycentre,
                             float phi_at_barycentre )
{
   float total_area = 0.0f;
   float positive_area = 0.0f;
   
   assert( polygon_vertices.size() > 2 );
   
   for ( unsigned int i = 0; i < polygon_vertices.size(); ++i )
   {
      unsigned int next = (i+1)%polygon_vertices.size();
      
      float triangle_area = 0.5f * mag( cross( polygon_vertices[i] - polygon_barycentre, polygon_vertices[next] - polygon_barycentre ) );
      
      float positive_triangle_area = compute_positive_area( polygon_barycentre, polygon_vertices[i], polygon_vertices[next],
                                                            phi_at_barycentre, polygon_phis[i], polygon_phis[next] );
      
      positive_triangle_area = clamp( positive_triangle_area, 0.0f, triangle_area );
            
      total_area += triangle_area;
      positive_area += positive_triangle_area;
   }

   if ( total_area == 0.0f )
   {
      // degenerate polygon
      
      // TODO: Review this!!  What's the correct thing to do here?
      
      if ( phi_at_barycentre > 0.0f ) 
      { 
         return 1.0f; 
      }
      else
      {
         return 0.0f;
      }
   }
   
   float positive_fraction = positive_area / total_area;
   positive_fraction = max( min( positive_fraction, 1.0f ), 0.0f );
   
   assert( positive_fraction == positive_fraction );
   
   return positive_fraction;
   
}


// ---------------------------------------------------------

Vec4f tet_barycentric_weights( const Vec3f& point, const Vec3f& a, const Vec3f& b, const Vec3f& c, const Vec3f& d )
{
   float volABCD = tet_signed_volume( a, b, c, d );
   
   Vec4f weights;
   
   weights[0] = tet_signed_volume( point, b, c, d ) / volABCD;
   weights[1] = tet_signed_volume( point, a, d, c ) / volABCD;
   weights[2] = tet_signed_volume( point, a, b, d ) / volABCD;
   weights[3] = 1.0f - weights[0] - weights[1] - weights[2];
   
   assert( weights[0] == weights[0] );
   assert( weights[1] == weights[1] );
   assert( weights[2] == weights[2] );
   assert( weights[3] == weights[3] );
   
   return weights;
}


Vec4f tet_barycentric_weights_careful( const Vec3f& point, const Vec3f& a, const Vec3f& b, const Vec3f& c, const Vec3f& d )
{
   float volABCD = tet_signed_volume( a, b, c, d );
   
   assert( volABCD > 0.0f );
   
   Vec4f weights;
   
   weights[0] = tet_signed_volume( point, b, c, d ) / volABCD;
   weights[1] = tet_signed_volume( point, a, d, c ) / volABCD;
   weights[2] = tet_signed_volume( point, a, b, d ) / volABCD;
   weights[3] = 1.0f - weights[0] - weights[1] - weights[2];

   weights[0] = clamp( weights[0], 0.0f, 1.0f );
   weights[1] = clamp( weights[1], 0.0f, 1.0f );
   weights[2] = clamp( weights[2], 0.0f, 1.0f );
   weights[3] = clamp( weights[3], 0.0f, 1.0f );
   
   float sum = weights[0] + weights[1] + weights[2] + weights[3];

   assert( sum > 0.0f );
   
   weights[0] /= sum;
   weights[1] /= sum;
   weights[2] /= sum;
   weights[3] /= sum;
   
   return weights;
}
// ---------------------------------------------------------
///
/// due to J R Shewchuk from Eppstein's geometry junkyard
///
// ---------------------------------------------------------


Vec3f tet_circumcentre( const Vec3f& a, const Vec3f& b, const Vec3f& c, const Vec3f& d )
{
   
   Mat33f M( b[0]-a[0], c[0]-a[0], d[0]-a[0], 
            b[1]-a[1], c[1]-a[1], d[1]-a[1], 
            b[2]-a[2], c[2]-a[2], d[2]-a[2] );
   
   float denom =  2.0f * determinant(M);
   
   Vec3f numer = mag2(d-a) * cross(b-a, c-a) + mag2(c-a) * cross(d-a, b-a) + mag2(b-a) * cross(c-a,d-a);
   
   return a + numer / denom;
}


// ---------------------------------------------------------
///
/// We'll probably want to sub in an exact arithmetic version at some point.
///
// ---------------------------------------------------------

bool point_in_tet( const Vec3f &p, const Vec3f &x1, const Vec3f &x2, const Vec3f &x3, const Vec3f &x4, float epsilon )
{
   
   // triangle 1 - x1 x2 x3
   float a = tet_signed_volume(p, x1, x2, x3);
   
   if (fabs(a) < epsilon)        // degenerate -- should check for point in triangle
   {        
      return true;
   }
   
   // triangle 2 - x2 x4 x3
   float b = tet_signed_volume(p, x2, x4, x3);
   
   if (fabs(b) < epsilon)         // degenerate -- should check for point in triangle
   {
      return true;
   }
   
   if ((a > epsilon) ^ (b > epsilon))
   {
      return false;
   }
   
   // triangle 3 - x1 x4 x2
   float c = tet_signed_volume(p, x1, x4, x2);
   
   if (fabs(c) < epsilon)     // degenerate -- should check for point in triangle
   {
      return true;
   }
   
   if ((a > epsilon) ^ (c > epsilon))
   {
      return false;
   }
   
   // triangle 4 - x1 x3 x4
   float d = tet_signed_volume(p, x1, x3, x4);
   
   if (fabs(d) < epsilon)     // degenerate -- should check for point in triangle
   { 
      return true;
   }
   
   if ((a > epsilon) ^ (d > epsilon))
   {
      return false;
   }
   
   return true;    // point is on the same side of all triangles
}


// ---------------------------------------------------------
///
/// Signed distance of a box (negative inside).
///
// ---------------------------------------------------------

float solid_box_phi( const Vec3f& box_centre, const Vec3f& box_extents, const Vec3f& eval_point )
{
   Vec3f offset = fabs(eval_point - box_centre);
   
   offset -= 0.5f * box_extents;
   
   if ( offset[0] < 0 && offset[1] < 0 && offset[2] < 0 )
   {
      return max(offset);
   }
   else 
   {
      float sum = 0;
      if(offset[0] > 0)
         sum += sqr(offset[0]);
      if(offset[1] > 0)
         sum += sqr(offset[1]);
      if(offset[2] > 0)
         sum += sqr(offset[2]);
      return sqrt(sum);
   }
   
}


Vec3f solid_box_gradient( const Vec3f& box_centre, const Vec3f& box_extents, const Vec3f& eval_point )
{

   float eval_sample = solid_box_phi( box_centre, box_extents, eval_point );

   Vec3f sample_point = eval_point + Vec3f( 1e-3f, 0, 0 );
   float sample = solid_box_phi( box_centre, box_extents, sample_point );
   
   float dx = sample - eval_sample;
   
   sample_point = eval_point + Vec3f( 0, 1e-3f, 0 );
   sample = solid_box_phi( box_centre, box_extents, sample_point );

   float dy = sample - eval_sample;

   sample_point = eval_point + Vec3f( 0, 0, 1e-3f );
   sample = solid_box_phi( box_centre, box_extents, sample_point );
   
   float dz = sample - eval_sample;

   return Vec3f(dx,dy,dz) / 1e-3f;
   
}



