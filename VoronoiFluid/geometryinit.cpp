// ---------------------------------------------------------
//
//  geometryinit.cpp
//  Tyson Brochu 2008
//
//  A set of geometry initialization functions.
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include "geometryinit.h"
#include "marching_tiles_hires.h"
#include "surftrack.h"
#include "array2.h"

// ---------------------------------------------------------
// Global externs
// ---------------------------------------------------------

// ---------------------------------------------------------
// Local constants, typedefs, macros
// ---------------------------------------------------------

// ---------------------------------------------------------
// Static and non-member function definitions
// ---------------------------------------------------------

using namespace ElTopo;

void create_circle( std::vector<Vec3d>& verts, 
                    std::vector<Vec3st>& tris, 
                    std::vector<double>& masses,
                    const Vec3d& centre,
                    double radius,
                    unsigned int nx )
{
   Vec3d low = centre - Vec3d( radius, 0.0, radius );
   
   double dx = 2.0 * radius / (double)nx;

   Array2<Vec3d> grid_points;
   grid_points.resize(nx, nx);
   for ( unsigned int i = 0; i < nx; ++i )
   {
      for ( unsigned int j = 0; j < nx; ++j )
      {
         grid_points(i,j) = low + dx * Vec3d( i, 0, j );
         verts.push_back( low + dx * Vec3d( i, 0, j ) );
      }
   }

   // cells
   for ( unsigned int i = 0; i < nx-1; ++i )
   {
      for ( unsigned int j = 0; j < nx-1; ++j )
      {
         unsigned int a = i   + nx*j;
         unsigned int b = i+1 + nx*j;
         unsigned int c = i+1 + nx*(j+1);
         unsigned int d = i   + nx*(j+1);

         
         if ( ( dist( verts[a], centre ) < radius ) &&
              ( dist( verts[b], centre ) < radius ) &&
              ( dist( verts[d], centre ) < radius ) )
         {
            tris.push_back( Vec3st( a, b, d ) );
         }

         if ( ( dist( verts[b], centre ) < radius ) &&
              ( dist( verts[c], centre ) < radius ) &&
              ( dist( verts[d], centre ) < radius ) )            
         {
            tris.push_back( Vec3st( b, c, d ) );
         }
      }
   }
   
   masses.clear();
   masses.resize( verts.size(), 1.0 );
   
   NonDestructiveTriMesh temp_tri_mesh;
   temp_tri_mesh.m_tris = tris;
   
   temp_tri_mesh.update_connectivity();
   
   for ( unsigned int i = 0; i < temp_tri_mesh.m_edges.size(); ++i )
   {
      if ( temp_tri_mesh.m_edge_to_triangle_map[i].size() == 1 )
      {
         unsigned int a = temp_tri_mesh.m_edges[i][0];
         unsigned int b = temp_tri_mesh.m_edges[i][1];
         
         masses[ a ] = 200.0;
         masses[ b ] = 200.0;
         
         // project to circle
         verts[ a ] += (radius - dist( verts[a], centre )) * normalized( verts[a] - centre );
         verts[ b ] += (radius - dist( verts[b], centre )) * normalized( verts[b] - centre );
      }
   }
   
}

// ---------------------------------------------------------
///
/// Create a vertical or horizontal, triangulated sheet.
///
// ---------------------------------------------------------

void create_sheet( std::vector<Vec3d>& verts, 
                   std::vector<Vec3st>& tris, 
                   const Vec3d& low_corner, 
                   const Vec3d& plane_normal, 
                   double dx, 
                   unsigned int nx, 
                   unsigned int ny )
{
   std::cout << "sheet: " << low_corner << ", " << plane_normal << ", " << dx << std::endl;
   std::cout << "resolution: " << nx << "x" << ny << std::endl;
   std::cout << "dimensions: " << dx*nx << "x" << dx*ny << std::endl;
   
   for(unsigned int i = 0; i < ny; i++)
   {
      for(unsigned int j = 0; j < nx; j++)
      {
         // COMPLETE HACK!
         
         if ( plane_normal[1] == 1.0 )
         {
            // plane normal is pointing in +y direction
            double x = dx*j;  //width*((double)j/nx); 
            double y = 0.0;
            double z = dx*i;  //width*((double)i/ny);              
            verts.push_back(Vec3d(x,y,z) + low_corner);
         }
         else
         {
            unsigned int num_curves = 4;
            double amplitude = dx*nx / 20.0;
            
            double x = dx*j; //width*((double)j/nx); 
            double y = dx*i; //width*((double)i/ny);    // make y correspond to rows

            double theta = 2.0 * M_PI * (double)j / (double)nx * (double)num_curves;
            double z = 1e-3*y + amplitude * sin( theta );  

            //double z = 1e-3*y;  
            
            verts.push_back(Vec3d(x,y,z) + low_corner);
         }
         
      }
   }
      
   for(unsigned int i = 0; i < ny-1; i++)
   {
      for(unsigned int j = 0; j < nx-1; j++)
      {
         unsigned int idx = i*(nx)+j;
         tris.push_back(Vec3st(idx, idx+(nx), idx+1));
         tris.push_back(Vec3st(idx+1, idx+(nx), idx+(nx)+1));
      }
   }
}


// ---------------------------------------------------------
///
///
///
// ---------------------------------------------------------

void merge_vertices( unsigned int a,   // delete
                     unsigned int b,   // keep
                     std::vector<Vec3st>& tris )
{
   for ( unsigned int i = 0; i < tris.size(); ++i )
   {
      if ( tris[i][0] == a ) { tris[i][0] = b; }
      if ( tris[i][1] == a ) { tris[i][1] = b; }
      if ( tris[i][2] == a ) { tris[i][2] = b; }
   }
}


// ---------------------------------------------------------
///
/// Create a cube surface mesh.
///
// ---------------------------------------------------------

void create_cube( const Vec3d& cube_low,
                  const Vec3d& cube_high,
                  const Vec3st& resolution,
                  std::vector<Vec3d>& verts, 
                  std::vector<Vec3st>& tris )
{

   //
   // Create each face seperately, then merge vertices.
   //
   
   unsigned int ni = resolution[0], nj = resolution[1], nk = resolution[2];
   
   double dx = (cube_high[0] - cube_low[0]) / (double) ni;
   double dy = (cube_high[1] - cube_low[1]) / (double) nj;
   double dz = (cube_high[2] - cube_low[2]) / (double) nk;
             
   unsigned int offset = verts.size();
   
   // face x = low
   for ( unsigned int j = 0; j <= nj; ++j ) for ( unsigned int k = 0; k <= nk; ++k )
   {
      verts.push_back( cube_low + Vec3d( 0, j*dy, k*dz ) );
      if ( j < nj && k < nk )
      {
         unsigned int idx = offset + j*(nk+1) + k;
         tris.push_back(Vec3st(idx+nk+1, idx, idx+1));
         tris.push_back(Vec3st(idx+nk+1, idx+1, idx+nk+2));
      }
   }
   
   // face x = high
   offset = verts.size();
   for ( unsigned int j = 0; j <= nj; ++j ) for ( unsigned int k = 0; k <= nk; ++k )
   {
      verts.push_back( cube_low + Vec3d( ni*dx, j*dy, k*dz ) );
      if ( j < nj && k < nk )
      {
         unsigned int idx = offset + j*(nk+1) + k;
         tris.push_back(Vec3st(idx,   idx+nk+1, idx+1));
         tris.push_back(Vec3st(idx+1, idx+nk+1, idx+nk+2));
      }
   }
   
   // face y = low
   offset = verts.size();
   for ( unsigned int i = 0; i <= ni; ++i ) for ( unsigned int k = 0; k <= nk; ++k )
   {
      verts.push_back( cube_low + Vec3d( i*dx, 0, k*dz ) );
      if ( i < ni && k < nk )
      {
         unsigned int idx = offset + i*(nk+1) + k;
         tris.push_back(Vec3st(idx,   idx+nk+1, idx+1));
         tris.push_back(Vec3st(idx+1, idx+nk+1, idx+nk+2));
      }
   }
   
   // face y = high
   offset = verts.size();
   for ( unsigned int i = 0; i <= ni; ++i ) for ( unsigned int k = 0; k <= nk; ++k )
   {
      verts.push_back( cube_low + Vec3d( i*dx, nj*dy, k*dz ) );
      if ( i < ni && k < nk )
      {
         unsigned int idx = offset + i*(nk+1) + k;
         tris.push_back(Vec3st(idx+nk+1, idx, idx+1));
         tris.push_back(Vec3st(idx+nk+1, idx+1, idx+nk+2));
      }
   }
   
   // face z = low
   offset = verts.size();
   for ( unsigned int i = 0; i <= ni; ++i ) for ( unsigned int j = 0; j <= nj; ++j )
   {
      verts.push_back( cube_low + Vec3d( i*dx, j*dy, 0 ) );
      if ( i < ni && j < nj )
      {
         unsigned int idx = offset + i*(nj+1) + j;
         tris.push_back(Vec3st(idx+nj+1, idx, idx+1));
         tris.push_back(Vec3st(idx+nj+1, idx+1, idx+nj+2));
      }
   }
   
   // face z = high
   offset = verts.size();
   for ( unsigned int i = 0; i <= ni; ++i ) for ( unsigned int j = 0; j <= nj; ++j )
   {
      verts.push_back( cube_low + Vec3d( i*dx, j*dy, nk*dz ) );
      if ( i < ni && j < nj )
      {
         unsigned int idx = offset + i*(nj+1) + j;
         tris.push_back(Vec3st(idx, idx+nj+1, idx+1));
         tris.push_back(Vec3st(idx+1, idx+nj+1, idx+nj+2));
      }
   }
   
   
   // now merge vertices (change the triangles pointing to the vertex to be deleted)
   
   for ( unsigned int i = 0; i < verts.size(); ++i )
   {
      for ( unsigned int j = i+1; j < verts.size(); ++j )
      {
         if ( dist( verts[i], verts[j] ) < 0.1 * dx )
         {
            merge_vertices( j, i, tris );
         }
      }      
   }

   
}


// ---------------------------------------------------------
///
/// Create a sphere surface
///
// ---------------------------------------------------------

void create_sphere( const Vec3d& sphere_centre,
                    double sphere_radius,
                    double dx,
                    std::vector<Vec3d>& verts, 
                    std::vector<Vec3st>& tris )
{

   const Vec3d domain_low = sphere_centre - Vec3d(sphere_radius + 3*dx);
   const Vec3d domain_high = sphere_centre + Vec3d(sphere_radius + 3*dx);
   Array3d phi;
   create_sphere_signed_distance( sphere_centre, sphere_radius, 1.5*dx, domain_low, domain_high, phi );  
   
   MarchingTilesHiRes marching_tiles( domain_low, 1.5*dx, phi );
   marching_tiles.contour();
   marching_tiles.improve_mesh();
   
   std::vector<Vec3d> new_verts( marching_tiles.x.size() );
   for ( unsigned int i = 0; i < new_verts.size(); ++i )
   {
      new_verts[i] = marching_tiles.x[i];
   }

   std::vector<Vec3st> new_tris( marching_tiles.tri.size() );
   for ( unsigned int i = 0; i < new_tris.size(); ++i )
   {
      new_tris[i] = marching_tiles.tri[i];
   }
   
   project_to_exact_sphere( new_verts, sphere_centre, sphere_radius );
   
   Vec3st offset( verts.size() );
   for ( unsigned int i = 0; i < new_verts.size(); ++i )
   {
      verts.push_back( new_verts[i] );
   }
   
   for ( unsigned int i = 0; i < new_tris.size(); ++i )
   {
      tris.push_back( new_tris[i] + offset );
   }
   
}


// ---------------------------------------------------------
///
/// Project mesh vertices onto an analytic cube
///
// ---------------------------------------------------------

void project_to_exact_cube( std::vector<Vec3d>& verts, 
                            const Vec3d& cube_low, 
                            const Vec3d& cube_high )
{

   for ( unsigned int i = 0; i < verts.size(); ++i )
   {
      Vec3d& v = verts[i];
      Vec3d p(0,0,0);
      double min_dist = 1e+30;
      
      if ( fabs( v[0] - cube_low[0] ) < min_dist ) { min_dist = fabs( v[0] - cube_low[0] ); p = v; p[0] = cube_low[0]; }
      if ( fabs( v[1] - cube_low[1] ) < min_dist ) { min_dist = fabs( v[1] - cube_low[1] ); p = v; p[1] = cube_low[1]; }
      if ( fabs( v[2] - cube_low[2] ) < min_dist ) { min_dist = fabs( v[2] - cube_low[2] ); p = v; p[2] = cube_low[2]; }
      if ( fabs( v[0] - cube_high[0] ) < min_dist ) { min_dist = fabs( v[0] - cube_high[0] ); p = v; p[0] = cube_high[0]; }
      if ( fabs( v[1] - cube_high[1] ) < min_dist ) { min_dist = fabs( v[1] - cube_high[1] ); p = v; p[1] = cube_high[1]; }
      if ( fabs( v[2] - cube_high[2] ) < min_dist ) { min_dist = fabs( v[2] - cube_high[2] ); p = v; p[2] = cube_high[2]; }
      
      v = p;
   }   
}

// ---------------------------------------------------------
///
/// Project mesh vertices onto an analytic sphere
///
// ---------------------------------------------------------

void project_to_exact_sphere( std::vector<Vec3d>& verts, const Vec3d& sphere_centre, double sphere_radius )
{
   for ( unsigned int i = 0; i < verts.size(); ++i )
   {
      Vec3d& v = verts[i];
      
      double dist = mag( v - sphere_centre ) - sphere_radius;
      
      v -= dist * ( v - sphere_centre ) / mag( v - sphere_centre );
   }
}


// ---------------------------------------------------------
///
/// Project mesh vertices onto two spheres
///
// ---------------------------------------------------------

void project_to_exact_2_spheres( SurfTrack* surf, 
                                 const Vec3d& sphere_a_centre, 
                                 const Vec3d& sphere_b_centre, 
                                 double max_radius, 
                                 double interior_radius )
{
   
   for ( unsigned int i = 0; i < surf->get_num_vertices(); ++i )
   {
      if ( surf->m_mesh.m_vertex_to_triangle_map[i].empty() )
      {
         continue;
      }
      
      Vec3d v = surf->get_position(i);
      
      double dist = signed_distance_entropy( v, sphere_a_centre, sphere_b_centre, max_radius, interior_radius );
      
      if ( mag( v - sphere_a_centre ) < mag ( v - sphere_b_centre ) )
      {
         v -= dist * ( v - sphere_a_centre ) / mag( v - sphere_a_centre );
      }
      else
      {
         v -= dist * ( v - sphere_b_centre ) / mag( v - sphere_b_centre );
      }
      surf->set_position(i, v);
      
   }
}


// ---------------------------------------------------------
///
/// Project mesh vertices onto an analytic dumbbell
///
// ---------------------------------------------------------

void project_to_exact_dumbbell( std::vector<Vec3d>& verts, 
                               const Vec3d& sphere_a_centre, 
                               const Vec3d sphere_b_centre, 
                               double sphere_radius, 
                               double handle_width )
{
   
   for ( unsigned int i = 0; i < verts.size(); ++i )
   {
      Vec3d& pt = verts[i];
      double dx = 1e-6;
      
      double x_pos = signed_distance_dumbbell( pt + Vec3d(dx,0,0), sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
      double x_neg = signed_distance_dumbbell( pt - Vec3d(dx,0,0), sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
      double y_pos = signed_distance_dumbbell( pt + Vec3d(0,dx,0), sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
      double y_neg = signed_distance_dumbbell( pt - Vec3d(0,dx,0), sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
      double z_pos = signed_distance_dumbbell( pt + Vec3d(0,0,dx), sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
      double z_neg = signed_distance_dumbbell( pt - Vec3d(0,0,dx), sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
      
      Vec3d n( x_pos - x_neg, y_pos - y_neg, z_pos - z_neg );
      normalize( n );
      
      double dist = signed_distance_dumbbell( pt, sphere_a_centre, sphere_b_centre, sphere_radius, handle_width );
      
      pt -= dist * n;
      
      assert( abs( signed_distance_dumbbell( pt, sphere_a_centre, sphere_b_centre, sphere_radius, handle_width ) < 1e-8 ) ); 
   }
   
}


// ---------------------------------------------------------
///
/// Read a signed distance field from an ASCII file
///
// ---------------------------------------------------------

void read_signed_distance( const char* filename, Array3d& signed_distance )
{
   FILE* file = fopen( filename, "r" );
   if ( file == NULL )
   {
      printf("failed to open file: %s\n", filename );
      assert(0);
      return;
   }
   
   int num_dim;
   fscanf( file, "%d ", &num_dim );
   
   if ( num_dim != 3 )
   {
      printf("num_dim: %d\n", num_dim );
      assert(0);
   }
   
   int dims[3];
   fscanf( file, "%d %d %d ", &(dims[0]), &(dims[1]), &(dims[2]) );
   
   unsigned int num_elements = dims[0]*dims[1]*dims[2];
   
   signed_distance.ni = dims[0];
   signed_distance.nj = dims[1];
   signed_distance.nk = dims[2];
   
   signed_distance.a.resize(num_elements);
   
   //fread( &buf, sizeof(double), num_elements, file );
   for ( unsigned int i = 0; i < num_elements; ++i )
   {
      float buf;
      fscanf( file, "%f ", &buf );
      signed_distance.a[i] = buf;
   }
   
   fclose(file);
   
   printf( "Read mat file\n" );
}



// ---------------------------------------------------------
///
/// Create signed distance field for a capsule with the given geometry
///
// ---------------------------------------------------------

void create_capsule_signed_distance( const Vec3d& capsule_end_a, 
                                     const Vec3d& capsule_end_b, 
                                     double capsule_radius,
                                     double dx,
                                     const Vec3d& domain_low,
                                     const Vec3d& domain_high,                                     
                                     Array3d& phi )
{ 
   phi.resize( (int) ceil( (domain_high[0]-domain_low[0]) / dx), (int) ceil( (domain_high[1]-domain_low[1]) / dx), (int) ceil( (domain_high[2]-domain_low[2]) / dx) );
   
   std::cout << "Generating signed distance function.  Grid resolution: " << phi.ni << " x " << phi.nj << " x " << phi.nk << std::endl;
   
   for ( int i = 0; i < phi.ni; ++i )
   {
      for ( int j = 0; j < phi.nj; ++j )
      {
         for ( int k = 0; k < phi.nk; ++k )
         {
            Vec3d pt = domain_low + dx * Vec3d(i,j,k);
            
            double distance;
            
            Vec3d central_segment(capsule_end_b - capsule_end_a);
            double m2=mag2(central_segment);
            
            // find parameter value of closest point on infinite line
            double s = dot(capsule_end_b - pt, central_segment)/m2;
            
            if ( s < 0.0 )
            {
               // dist = distance to the cylinder disc at end b
               
               distance = dist(pt, capsule_end_b);
               distance -= capsule_radius;
               
               
            } 
            else if ( s > 1.0 )
            {
               // dist = distance to the cylinder disc at end b
               distance = dist(pt, capsule_end_a );
               distance -= capsule_radius;
               
            }
            else
            {
               // dist = distance to the cylinder's central axis
               
               distance = dist(pt, s*capsule_end_a + (1-s)*capsule_end_b);
               distance -= capsule_radius;
            }
            
            phi(i,j,k) = distance;
            
         }
      }
   }
}


// ---------------------------------------------------------
///
/// Create signed distance field for a cube with the given geometry
///
// ---------------------------------------------------------

void create_cube_signed_distance( const Vec3d& cube_low, 
                                 const Vec3d& cube_high, 
                                 double dx,
                                 const Vec3d& domain_low,
                                 const Vec3d& domain_high,                                     
                                 Array3d& phi )
{
   phi.resize( (int) ceil( (domain_high[0]-domain_low[0]) / dx), (int) ceil( (domain_high[1]-domain_low[1]) / dx), (int) ceil( (domain_high[2]-domain_low[2]) / dx) );
   
   std::cout << "Generating signed distance function.  Grid resolution: " << phi.ni << " x " << phi.nj << " x " << phi.nk << std::endl;
   
   for ( int i = 0; i < phi.ni; ++i )
   {
      for ( int j = 0; j < phi.nj; ++j )
      {
         for ( int k = 0; k < phi.nk; ++k )
         {
            Vec3d pt = domain_low + dx * Vec3d(i,j,k);
            
            double dist_low_x = cube_low[0] - pt[0];
            double dist_high_x = pt[0] - cube_high[0];
            
            double dist_low_y = cube_low[1] - pt[1];
            double dist_high_y = pt[1] - cube_high[1];
            
            double dist_low_z = cube_low[2] - pt[2];
            double dist_high_z = pt[2] - cube_high[2];
            
            phi(i,j,k) = max( dist_low_x, dist_high_x, dist_low_y, dist_high_y, dist_low_z, dist_high_z );
            
         }
      }
   }
}


// ---------------------------------------------------------
///
/// Create signed distance field for a sphere with the given geometry
///
// ---------------------------------------------------------

void create_sphere_signed_distance( const Vec3d& sphere_centre, 
                                   double sphere_radius, 
                                   double dx,
                                   const Vec3d& domain_low,
                                   const Vec3d& domain_high,                                     
                                   Array3d& phi )
{
   phi.resize( (int) ceil( (domain_high[0]-domain_low[0]) / dx), (int) ceil( (domain_high[1]-domain_low[1]) / dx), (int) ceil( (domain_high[2]-domain_low[2]) / dx) );
   
   std::cout << "Generating signed distance function.  Grid resolution: " << phi.ni << " x " << phi.nj << " x " << phi.nk << std::endl;
   
   for ( int i = 0; i < phi.ni; ++i )
   {
      for ( int j = 0; j < phi.nj; ++j )
      {
         for ( int k = 0; k < phi.nk; ++k )
         {
            Vec3d pt = domain_low + dx * Vec3d(i,j,k);
            double dist = mag( pt - sphere_centre );
            phi(i,j,k) = dist - sphere_radius;
         }
      }
   }
}


// ---------------------------------------------------------
///
/// Analytic entropy solution to motion in the normal direction of two spheres
///
// ---------------------------------------------------------

double signed_distance_entropy(const Vec3d& x,
                               const Vec3d& a,
                               const Vec3d& b,
                               double max_radius,
                               double interior_radius)
{
   double d=dist(a, b);
   double alpha=dot(x-a, b-a)/d; // distance along vector from a to b of projection of x onto line
   
   if(alpha<=0){ // x is on the far side of a from b
      return dist(x, a) - max_radius + interior_radius; // regular sphere distance
      
   }else if(alpha>=d){ // x is on the far side of b from a
      return dist(x, b) - max_radius + interior_radius; // regular sphere distance
      
   }else if(alpha<=0.5*d){ // x is between a and b, but closer to a
      if(max_radius<=0.5*d){ // if the spheres at maximum radius never intersected, life is simple
         return dist(x, a) - max_radius + interior_radius; // regular sphere distance
      }
      double beta=std::sqrt(dist2(x, a) - sqr(alpha)); // distance between x and projection onto line through a and b
      double gamma=std::sqrt(sqr(max_radius) - sqr(0.5*d)); // radius of intersection curve between spheres at max_radiu
      if(beta/alpha>=gamma/(0.5*d)){ // if closest point is still on a's sphere (not in the intersection region)
         return dist(x, a) - max_radius + interior_radius; // regular sphere distance
      }else{
         // use the distance from closest point on the intersection circle of the two spheres
         return interior_radius - sqrt(sqr(gamma-beta) + sqr(0.5*d-alpha));
      }
      
   }else{ // x is between a and b, but closer to b
      if(max_radius<=0.5*d){ // if the spheres at maximum radius never intersected, life is simple
         return dist(x, b) - max_radius + interior_radius; // regular sphere distance
      }
      double beta=std::sqrt(dist2(x, a) - sqr(alpha)); // distance between x and projection onto line through a and b
      double gamma=std::sqrt(sqr(max_radius) - sqr(0.5*d)); // radius of intersection curve between spheres at max_radiu
      if(beta/(d-alpha)>=gamma/(0.5*d)){ // if closest point is still on b's sphere (not in the intersection region)
         return dist(x, b) - max_radius + interior_radius; // regular sphere distance
      }else{
         // use the distance from closest point on the intersection circle of the two spheres
         return interior_radius - sqrt(sqr(gamma-beta) + sqr(alpha-0.5*d));
         
      }
   }
}


// ---------------------------------------------------------
///
/// Create signed distance field for two spheres
///
// ---------------------------------------------------------

void create_two_sphere_signed_distance( const Vec3d& sphere_a_centre, 
                                       const Vec3d& sphere_b_centre, 
                                       double sphere_radius, 
                                       double dx,
                                       const Vec3d& domain_low,
                                       const Vec3d& domain_high,                                     
                                       Array3d& phi )
{  
   phi.resize( (int) ceil( (domain_high[0]-domain_low[0]) / dx), (int) ceil( (domain_high[1]-domain_low[1]) / dx), (int) ceil( (domain_high[2]-domain_low[2]) / dx) );
   
   std::cout << "Generating signed distance function.  Grid resolution: " << phi.ni << " x " << phi.nj << " x " << phi.nk << std::endl;
   
   for ( int i = 0; i < phi.ni; ++i )
   {
      for ( int j = 0; j < phi.nj; ++j )
      {
         for ( int k = 0; k < phi.nk; ++k )
         {
            Vec3d pt = domain_low + dx * Vec3d(i,j,k);
            phi(i,j,k) = signed_distance_entropy( pt, sphere_a_centre, sphere_b_centre, sphere_radius, 0.0 );
         }
      }
   }
   
}


// ---------------------------------------------------------
///
/// Analytic signed distance function for a dumbbell
///
// ---------------------------------------------------------

double signed_distance_dumbbell( const Vec3d& pt,
                                const Vec3d& sphere_a_centre, 
                                const Vec3d& sphere_b_centre, 
                                double sphere_radius, 
                                double handle_width )
{
   double dist_spheres = ( min( mag(pt - sphere_a_centre), mag(pt - sphere_b_centre) ) - sphere_radius );
   return min( dist_spheres, max( fabs(pt[0]) - sphere_b_centre[0], sqrt( pt[1]*pt[1] + pt[2]*pt[2] ) - handle_width ) );   
}


// ---------------------------------------------------------
///
/// Create signed distance field for a dumbbell
///
// ---------------------------------------------------------

void create_dumbbell_signed_distance( const Vec3d& sphere_a_centre, 
                                     const Vec3d& sphere_b_centre, 
                                     double sphere_radius, 
                                     double handle_width,
                                     double dx,
                                     const Vec3d& domain_low,
                                     const Vec3d& domain_high,                                     
                                     Array3d& phi )
{
   
   phi.resize( (int) ceil( (domain_high[0]-domain_low[0]) / dx), (int) ceil( (domain_high[1]-domain_low[1]) / dx), (int) ceil( (domain_high[2]-domain_low[2]) / dx) );
   
   std::cout << "Generating signed distance function.  Grid resolution: " << phi.ni << " x " << phi.nj << " x " << phi.nk << std::endl;
   
   for ( int i = 0; i < phi.ni; ++i )
   {
      for ( int j = 0; j < phi.nj; ++j )
      {
         for ( int k = 0; k < phi.nk; ++k )
         {
            Vec3d pt = domain_low + dx * Vec3d(i,j,k);
            phi(i,j,k) = signed_distance_dumbbell( pt, sphere_a_centre, sphere_b_centre, sphere_radius, handle_width ); 
            
         }
      }
   }
}


// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------


