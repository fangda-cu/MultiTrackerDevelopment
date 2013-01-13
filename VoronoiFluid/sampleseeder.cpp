
#include "sampleseeder.h"

#include <surftrack.h>
#include <wallclocktime.h>
#include <set>

using namespace ElTopo;

static const float COLLOCATION_THRESHOLD = 1e-4f;

// ---------------------------------------------------------
///
///
///
// ---------------------------------------------------------

static inline bool add_unique( const Vec3f& x, const std::vector<Vec3f>& check_against_points, std::vector<Vec3f>& add_to_points )
{
   
   bool already_exists = false;
   for ( unsigned int i = 0; i < check_against_points.size(); ++i )
   {
      if ( dist2( x, check_against_points[i] ) < COLLOCATION_THRESHOLD*COLLOCATION_THRESHOLD ) 
      {
         already_exists = true;
         break;
      }               
   }
   
   if ( !already_exists )
   {
      add_to_points.push_back(x);
   }

   return !already_exists;
    
}

// ---------------------------------------------------------
///
/// Create a set of BCC vertices.
///
// ---------------------------------------------------------

void SampleSeeder::generate_bcc_points( const Vec3f& domain_low, const Vec3f& domain_high, float dx, std::vector<Vec3f>& xs )
{
   Vec3st n = Vec3st( (domain_high - domain_low) / dx ) + Vec3st(1);

   std::vector<Vec3f> new_points;
   
   // regular grid: cell corners
   for ( unsigned int i = 0; i < n[0]; ++i )
   {
      for ( unsigned int j = 0; j < n[1]; ++j )
      {
         for ( unsigned int k = 0; k < n[2]; ++k )
         {
            add_unique( domain_low + dx * Vec3f(Vec3st(i,j,k)), xs, new_points );
         }
      }
   }
   
   // and cell centres
   for ( unsigned int i = 0; i < n[0]-1; ++i )
   {
      for ( unsigned int j = 0; j < n[1]-1; ++j )
      {
         for ( unsigned int k = 0; k < n[2]-1; ++k )
         {
            add_unique( domain_low + Vec3f(0.5f*dx) + dx * Vec3f(Vec3st(i,j,k)), xs, new_points );
         }
      }
   }
   
   xs.insert( xs.end(), new_points.begin(), new_points.end() );
   
}



double g_air_sample_rejection_threshold = 0.1;


// ---------------------------------------------------------
///
/// Create a set of points fitting the given surface.
///
// ---------------------------------------------------------

void SampleSeeder::generate_adaptive_points( const SurfTrack& surface, 
                                             double desired_dx, 
                                             std::vector<Vec3f>& samples )
{
   
   
   std::cout << "g_air_sample_rejection_threshold: " << g_air_sample_rejection_threshold << std::endl;
   
   
   double start_time = get_time_in_seconds();
   
   unsigned int num_rejections = 0;

   // for each vertex on the surface
   
   for ( unsigned int i = 0; i < surface.get_positions().size(); ++i )
   {
      if ( surface.m_mesh.m_vertex_to_edge_map[i].size() < 1 ) { continue; }

      extern bool seed_samples_at_solid;
      
      if ( !seed_samples_at_solid )
      {
         // only add a sample if the surface vertex is on a free surface, *or* is adjacent to a vertex on a free surface
         bool free_surface = ( surface.m_masses[i] < 1.5 );
         if ( !free_surface )
         {
            const std::vector<unsigned int>& incident_edges = surface.m_mesh.m_vertex_to_edge_map[i];
            for ( unsigned int j = 0; j < incident_edges.size(); ++j )
            {
               unsigned int neighbour = surface.m_mesh.m_edges[ incident_edges[j] ][0];
               if ( neighbour == i ) { neighbour = surface.m_mesh.m_edges[ incident_edges[j] ][1]; }
               if ( surface.m_masses[neighbour] < 1.5 )
               {
                  free_surface = true;
                  break;
               }
            }
         }
         
         if ( !free_surface ) { continue; }
      }
      
      std::vector<Vec3d> ray_dirs;

      const Vec3d& ray_origin = surface.get_position(i);
      
      //get the list of labels
      std::set<int> labels;
      for(size_t t = 0; t < surface.m_mesh.m_vertex_to_triangle_map[i].size(); ++t) {
         int tri_ID = surface.m_mesh.m_vertex_to_triangle_map[i][t];
         labels.insert(surface.m_mesh.get_triangle_label(tri_ID)[0]);
         labels.insert(surface.m_mesh.get_triangle_label(tri_ID)[1]);
      }
      
      //compute normal for each one
      for(std::set<int>::iterator it = labels.begin(); it != labels.end(); ++it) {
         int cur_label = *it;
         Vec3d normal = surface.get_vertex_normal_angleweighted_by_label(i, cur_label);
         ray_dirs.push_back(-normal);
      }

      //
      // fire a ray in the various normal direction
      //
      
      // ignore incident triangles
      const std::vector<unsigned int>& incident_triangles = surface.m_mesh.m_vertex_to_triangle_map[i];
      
      for(size_t i = 0; i < ray_dirs.size(); ++i)
      {      
         Vec3d normal = ray_dirs[i];
         
         if ( mag(normal) == 0.0 ) { continue; }

         const Vec3d ray_end = ray_origin + desired_dx * normal;

         std::vector<double> hit_ss;
         std::vector<unsigned int> hit_triangles; 
         
         surface.get_triangle_intersections( ray_origin, ray_end, hit_ss, hit_triangles );
         
         // get first ray intersection time (bounded by desired_dx)
         double min_hit = desired_dx;
         for ( unsigned int j = 0; j < hit_ss.size(); ++j )
         {
            // skip incident triangles
            bool hit_incident = false;
            for ( unsigned int t = 0; t < incident_triangles.size(); ++t )
            {
               if ( incident_triangles[t] == hit_triangles[j] )
               {
                  hit_incident = true;
               }
            }
            if ( hit_incident ) { continue; }
            
            min_hit = min( min_hit, desired_dx * hit_ss[j] );
            
         }
         
        
         Vec3f new_sample = Vec3f( ray_origin + 0.33 * min_hit * normal );
         //TODO Reject really close samples in free-surface situations
         //to better encourage merging? Or some other thresholding in
         //general situations?
        /* if( min_hit < g_air_sample_rejection_threshold * desired_dx ) 
         {
            ++num_rejections;
            continue;
         }*/

                  
         assert( new_sample[0] == new_sample[0] );
         assert( new_sample[1] == new_sample[1] );
         
         add_unique( new_sample, samples, samples );
      }
   }
   
   std::vector<Vec3f> surf_samples;
   
   //try adding samples at tri barycenters to break up regularity and avoid slivers.
   bool parity = false;
   for(unsigned int i = 0; i < surface.m_mesh.m_tris.size(); ++i) {
      if(surface.m_mesh.triangle_is_deleted(i)) 
         continue;
      Vec3d new_sample(0,0,0);
      Vec3st tri = surface.m_mesh.m_tris[i];
      for(int j = 0; j < 3; ++j)
         new_sample += surface.get_position(tri[j]);
      new_sample /= 3;
      
      //nudge it in or out slightly just to make the point's region less ambiguous
      new_sample += 0.01*desired_dx * surface.get_triangle_normal(i) * (parity?1:-1);
      parity=!parity;

      bool accepted = add_unique( Vec3f(new_sample), samples, samples );
   }


   std::cout << "number of samples rejected: " << num_rejections << std::endl;
   std::cout << "seed time: " << get_time_in_seconds() - start_time << std::endl;
   
}


