
#include "dualfluidsim3d.h"

#include <bfstream.h>
#include "dualpressure3d.h"
#include "geometryutils.h"
#include "iomesh.h"
#include "lapack_wrapper.h"
#include "mat.h"

#include <queue>
#include "sampleseeder.h"

#include "wallclocktime.h"
#include "CGAL_wrapper.h"



using namespace ElTopo;



// always handy
template<unsigned int M, unsigned int N, class T>
static bool least_squares( Mat<M,N,T>&matrix, Vec<M,T>&rhs ) {
   int info = 0;
   int rank;

   double rcond = 0.1; //treat as 0's any singular values (SV) that are less than rcond * largest SV, for purposes of rank computation
   LAPACK::solve_least_squares('N', M, N, 1, matrix.a, M, rhs.v, max(M,N), info, rcond, rank);
   
   //if(info != 0)
   //   printf("Info %d\n", info);
   //assert( info == 0 );
   
   //return failure if the data is rank-deficient. (we lack good velocity information some axes, so the solution will be bad).
   if(rank < 3) return false;

   return info == 0;
}


// ---------------------------------------------------------
///
///
///
// ---------------------------------------------------------

DualFluidSim3D::DualFluidSim3D( const std::vector<Vec3d>& surface_vertices, 
                                const std::vector<Vec3ui>& surface_triangles, 
                                const std::vector<double>& surface_vertex_masses,
                                const SurfTrackInitializationParameters& initial_parameters )
{
  mesh = new TetMesh(); 
  surface_tracker = new SurfTrack( surface_vertices, surface_triangles, surface_vertex_masses, initial_parameters );
  /*surface_tracker->improve_mesh();
  surface_tracker->improve_mesh();
  surface_tracker->improve_mesh();*/
}


DualFluidSim3D::~DualFluidSim3D()
{
   delete surface_tracker;
   delete mesh;
}


// ---------------------------------------------------------
///
/// Set up fluid sim.
///
// ---------------------------------------------------------

void DualFluidSim3D::initialize()
{
   compute_solids();   
   compute_liquid_phi();
   //extrapolate_liquid_phi_into_solid();
   
   densities.resize(3);
   densities[0] = 1.0;
   densities[1] = 1.0f;
   densities[2] = 1.0f;

   tet_edge_velocities.resize( mesh->edges.size(), 0.0f );
   tet_vertex_velocities.resize( mesh->vertices.size(), Vec3f(0.0f) );
   tet_vertex_velocity_is_valid.resize( mesh->vertices.size(), false );
   edge_is_valid.resize( mesh->edges.size(), false );
   should_remesh = false;
   surface_tension_coefficient = 0;//0.000001;
   
   allow_solid_overlap = false;
   
   total_remesh_time = 0;
   total_semilagrangian_time = 0;
   
   gravity = Vec3f(0,-1,0);
   
   interpolation_scheme = BARYCENTRIC;

}


// ---------------------------------------------------------
///
/// Evaluate the solid signed distance function at each tet circumcentre, then get the solid fraction at each tet edge.
/// (Dual: Evaluate the solid signed distance function at each voronoi vertex, then get the solid fraction on each voronoi face.)
///
// ---------------------------------------------------------

void DualFluidSim3D::compute_solids()
{
   solid_phi.resize( mesh->tets.size() );
      
   Vec3f box_centre = 0.5f * (solid_low + solid_high);
   Vec3f box_extents = solid_high - solid_low;
   
   for ( unsigned int i = 0; i < mesh->tets.size(); ++i )
   {
      solid_phi[i] = solid_box_phi( box_centre, box_extents, mesh->tet_circumcentres[i] );
   }

   
   //
   // compute solid weights per voronoi face
   //
   
   solid_weights.resize( mesh->edges.size(), 1.0f );
   
   for ( unsigned int edge_index = 0; edge_index < mesh->edges.size(); ++edge_index )
   {
      std::vector<unsigned int> sorted_incident_tets;
      mesh->get_closed_tet_neighbourhood( edge_index, sorted_incident_tets );
      if ( sorted_incident_tets.empty() ) 
      { 
         // boundary edge
         solid_weights[edge_index] = 0.0; 
         continue;
      }
      
      // OPTIMIZE: A potential optimization might be to determine if all phis are the 
      // same sign and avoid computing the fractional weight.  This will likely be
      // the most common case.
      
      std::vector<Vec3f> polygon_vertices;
      std::vector<float> polygon_phis;
      
      Vec3f polygon_barycentre(0);
      for ( unsigned int i = 0; i < sorted_incident_tets.size(); ++i )
      {
         unsigned int tet_index = sorted_incident_tets[i];
         polygon_vertices.push_back( mesh->tet_circumcentres[tet_index] );
         polygon_phis.push_back( solid_phi[tet_index] );
         polygon_barycentre += mesh->tet_circumcentres[tet_index];
      }
      
      polygon_barycentre /= (float) sorted_incident_tets.size();
      float phi_at_barycentre = solid_box_phi( box_centre, box_extents, polygon_barycentre );
      
      // this computes the fraction of the polygon which is positive phi...
      float positive_fraction = compute_polygon_weight( polygon_vertices, polygon_phis, polygon_barycentre, phi_at_barycentre );
      
      assert( positive_fraction >= 0.0f );
      assert( positive_fraction <= 1.0f );
      
      // ...so the fraction of the polygon which is *non-solid* is:
      solid_weights[edge_index] = (float)(1.0 - positive_fraction);
      
      if( solid_weights[edge_index] != solid_weights[edge_index] )
      {
         std::cout << "----------" << std::endl;
         std::cout << "polygon: " << std::endl;
         for ( unsigned int i = 0; i < polygon_vertices.size(); ++i ) { std::cout << polygon_vertices[i] << std::endl; }
         for ( unsigned int i = 0; i < sorted_incident_tets.size(); ++i )
         {
            std::cout << "tet " << sorted_incident_tets[i] << ": " << mesh->tets[ sorted_incident_tets[i] ] << std::endl;
            std::cout << "circumcentre: " << mesh->tet_circumcentres[ sorted_incident_tets[i] ] << std::endl;
         }
         assert(0);
      }
   }

}

//


// ---------------------------------------------------------
///
/// Use barycentric interpolation to get a 3d velocity vector from the tet vertices.
///
// ---------------------------------------------------------

Vec3f DualFluidSim3D::get_velocity_from_tet_vertices( const Vec3f& point )
{
   int tet_index = mesh->get_containing_tet( point );
   
   if ( tet_index < 0 ) 
   { 
      return Vec3f(0.0f); 
   }

   const Vec4ui& tet = mesh->tets[tet_index];
   const Vec3f& a = mesh->vertices[tet[0]];
   const Vec3f& b = mesh->vertices[tet[1]];
   const Vec3f& c = mesh->vertices[tet[2]];
   const Vec3f& d = mesh->vertices[tet[3]];
   
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
   
   Vec3f weighted_velocity = weights[0] * tet_vertex_velocities[tet[0]] 
                           + weights[1] * tet_vertex_velocities[tet[1]] 
                           + weights[2] * tet_vertex_velocities[tet[2]] 
                           + weights[3] * tet_vertex_velocities[tet[3]];
   
//   if ( mag(weighted_velocity) > 2.0 * characteristic_distance )
//   {
//      std::cout << "weights: " << weights[0] << ", " << weights[1] << ", " << weights[2] << ", " << weights[3] << std::endl;
//      std::cout << "tet: " << std::endl;
//      std::cout << a << std::endl;
//      std::cout << b << std::endl;
//      std::cout << c << std::endl;
//      std::cout << d << std::endl;
//      std::cout << "tet volume: " << std::endl;
//      std::cout << "velocities: " << volABCD << std::endl;
//      std::cout << tet_vertex_velocities[tet[0]] << std::endl;
//      std::cout << tet_vertex_velocities[tet[1]] << std::endl;
//      std::cout << tet_vertex_velocities[tet[2]] << std::endl;
//      std::cout << tet_vertex_velocities[tet[3]] << std::endl;
//   }
   
   return weighted_velocity;
   
}



// ---------------------------------------------------------
///
/// Use sharper barycentric interpolation to get a 3d velocity vector from the constructed sub-voronoi tetrahedra
///
// ---------------------------------------------------------

Vec3f DualFluidSim3D::get_sharper_barycentric_velocity( const Vec3f& point )
{
   int voronoi_index = mesh->get_containing_voronoi( point );
   
   if ( voronoi_index < 0 ) 
   { 
      //printf("Failed to find a voronoi cell.\n");
      return Vec3f(0.0f); 
   }
   
   //determine which sub-tet the point belongs to, and interpolate velocity from its points
   
   Vec3f voronoi_central_vertex = mesh->vertices[voronoi_index];

   //consider each voronoi face
   for(unsigned int i = 0; i < mesh->vert_to_edge_map[voronoi_index].size(); ++i) {
      
      int face_ind = mesh->vert_to_edge_map[voronoi_index][i];

      Vec3f face_central_vertex = mesh->voronoi_face_centroids[face_ind];
      
      std::vector<unsigned int> sorted_incident_tets;
      mesh->get_closed_tet_neighbourhood(face_ind, sorted_incident_tets);
      
      //within each face, consider each sub-tet produced by splitting the face into triangles via the centroid, and joining to the voronoi site
      unsigned int limit = sorted_incident_tets.size();
      for(unsigned int v0 = limit-1, v1 = 0; v1 < limit; v0 = v1++) { //cute trick to enumerate polygon edges, there's a JGT paper
         
         //determine the two remaining vertices of the tetrahedra
         int v0_ind = sorted_incident_tets[v0];
         int v1_ind = sorted_incident_tets[v1];
         Vec3f vert0 = mesh->tet_circumcentres[v0_ind];
         Vec3f vert1 = mesh->tet_circumcentres[v1_ind];
         
         //now we have four vertices forming the tet:
         //v0, v1, face_centre, voronoi_centre.
         //determine if our point is in this tet
         if(point_in_tet(point, vert0, vert1, face_central_vertex, voronoi_central_vertex, 1e-10f)) {
            
            //use barycentric interpolation to get the velocity from the tet
            float volume = tet_signed_volume(  vert0, vert1, face_central_vertex, voronoi_central_vertex );
            
            Vec3f result;

            if( fabs(volume) < 1e-7) {
               //std::cout << "Points: " << vert0 << " " << vert1 << " " << face_central_vertex << " " << voronoi_central_vertex << std::endl;
               //printf("Degenerate case: %f\n", volume);
               //Hopefully this only happens near boundaries
               return tet_vertex_velocities[voronoi_index];
            }

            //order appropriately, by swapping the last two entries, so that the volume is positive
            if(volume <= 0) {
               Vec4f weights = tet_barycentric_weights_careful(point, vert0, vert1, voronoi_central_vertex, face_central_vertex);

               //compute the result by barycentric interpolation
               result = weights[0] * voronoi_vertex_velocities[v0_ind] + 
                        weights[1] * voronoi_vertex_velocities[v1_ind] + 
                        weights[2] * tet_vertex_velocities[voronoi_index] +
                        weights[3] * full_voronoi_face_velocities[face_ind];
            
            }
            else {
               Vec4f weights = tet_barycentric_weights_careful(point, vert0, vert1, face_central_vertex, voronoi_central_vertex);

               //compute the result by barycentric interpolation
               result = weights[0] * voronoi_vertex_velocities[v0_ind] + 
                        weights[1] * voronoi_vertex_velocities[v1_ind] + 
                        weights[2] * full_voronoi_face_velocities[face_ind] + 
                        weights[3] * tet_vertex_velocities[voronoi_index];
            }
            return result;
         }
      }
   }
   
   //assert(false);
   //failed to determine a containing tet; just return the velocity of the voronoi site for lack of something better.
   
   return tet_vertex_velocities[voronoi_index];
   
}

// ---------------------------------------------------------
///
/// Use standard generalized barycentric interpolation to get a 3d velocity vector from the Voronoi vertices
///
// ---------------------------------------------------------

Vec3f DualFluidSim3D::get_generalized_barycentric_velocity( const Vec3f& point )
{
   int voronoi_index = mesh->get_containing_voronoi( point );
   
   if ( voronoi_index < 0 ) 
   { 
      //printf("Failed to find a voronoi cell.\n");
      return Vec3f(0.0f); 
   }
   
   Vec3f voronoi_central_vertex = mesh->vertices[voronoi_index];

   //consider each adjacent tetrahedron
   std::vector<float> weights(mesh->vert_to_tet_map[voronoi_index].size());

   for(unsigned int i = 0; i < mesh->vert_to_tet_map[voronoi_index].size(); ++i) {
      int tet_ID = mesh->vert_to_tet_map[voronoi_index][i];
      
      float cur_weight = 6 * mesh->tet_volumes[tet_ID];
      Vec3f tet_centre = mesh->tet_circumcentres[tet_ID];
      
      for(unsigned int j = 0; j < 4; ++j) {
         int tet_vert_ID = mesh->tets[tet_ID][j];
         
         //no need to process the point that corresponds to the Voronoi site
         if(tet_vert_ID == voronoi_index) 
            continue;
         
         Vec3f tet_point = mesh->vertices[tet_vert_ID];
            
         float denom_term = dot(tet_point - voronoi_central_vertex, tet_centre - point );
         if(denom_term < 1e-13)
            denom_term = 1e-13f;
         cur_weight /= denom_term;
      }
      weights[i] = cur_weight;
   }

   //determine the sum of the weights
   float weight_sum = 0;
   for(unsigned int i = 0; i < weights.size(); ++i) {
      weight_sum += weights[i];
   }

   //normalize them, and compute the final weight
   Vec3f result(0,0,0);
   for(unsigned int i = 0; i < weights.size(); ++i) {
      weights[i] /= weight_sum;
      int tet_ID = mesh->vert_to_tet_map[voronoi_index][i];
      result += voronoi_vertex_velocities[tet_ID] * weights[i];
   }

   return result;
   
}

// ---------------------------------------------------------
///
/// Use Whitney-element-style interpolation to get a velocity from tangential velocity components on tet edges
///
// ---------------------------------------------------------

Vec3f DualFluidSim3D::get_velocity_from_tet_edges( const Vec3f& point )
{
   int tet_index = mesh->get_containing_tet( point );
   
   if ( tet_index < 0 ) { return Vec3f(0.0f); }
   const Vec4ui& tet = mesh->tets[tet_index];
   Vec4f bary = tet_barycentric_weights( point, mesh->vertices[tet[0]], mesh->vertices[tet[1]], mesh->vertices[tet[2]], mesh->vertices[tet[3]] );
   
   const Vec6ui& tet_edges = mesh->tet_to_edge_map[tet_index];

   // for each vertex, which triangle is its opposite (in the tet_to_tri_map)
   unsigned int opposite_triangles[4] = { 3, 2, 1, 0 };
   
   Vec3f total_velocity(0);
   for ( unsigned int i = 0; i < 6; ++i )
   {
      unsigned int a = mesh->edges[tet_edges[i]][0];
      unsigned int b = mesh->edges[tet_edges[i]][1];                            
      
      // find the index of each vertex within the tet (0-3)
      
      unsigned int a_in_tet = ~0, b_in_tet = ~0;
      for ( unsigned int j = 0; j < 4; ++j )
      {
         if ( a == tet[j] ) { a_in_tet = j; }
         if ( b == tet[j] ) { b_in_tet = j; }         
      }

      assert( a_in_tet != (unsigned int) ~0 );
      assert( b_in_tet != (unsigned int) ~0 );
      
      float na = clamp( bary[a_in_tet], 0.0f, 1.0f );
      float nb = clamp( bary[b_in_tet], 0.0f, 1.0f );

      Vec3f grad_na;
      Vec3f grad_nb;
      
      // get the (out-facing) normal of the triangle opposite vertex a      
      
      {
         unsigned int tri_in_tet = opposite_triangles[a_in_tet];
         unsigned int opp_triangle_a = mesh->tet_to_tri_map[tet_index][tri_in_tet];
         const Vec3ui& tri = mesh->tris[opp_triangle_a];
         assert( tri[0] != a && tri[1] != a && tri[2] != a );         
         float sign = same_orientation( tet, tri ) ? -1.0f : 1.0f;
         if ( tri_in_tet % 2 == 0 ) { sign *= -1.0f; }
         Vec3f out_normal = sign * cross( mesh->vertices[tri[1]] - mesh->vertices[tri[0]], mesh->vertices[tri[2]] - mesh->vertices[tri[0]] );
         grad_na = out_normal;
      }
      
      // get the (out-facing) normal of the triangle opposite vertex b
      
      {
         unsigned int tri_in_tet = opposite_triangles[b_in_tet];
         unsigned int opp_triangle_b = mesh->tet_to_tri_map[tet_index][tri_in_tet];
         const Vec3ui& tri = mesh->tris[opp_triangle_b];
         assert( tri[0] != b && tri[1] != b && tri[2] != b );
         float sign = same_orientation( tet, tri ) ? -1.0f : 1.0f;
         if ( tri_in_tet % 2 == 0 ) { sign *= -1.0f; }
         Vec3f out_normal = sign * cross( mesh->vertices[tri[1]] - mesh->vertices[tri[0]], mesh->vertices[tri[2]] - mesh->vertices[tri[0]] );
         grad_nb = out_normal;
      }
      
      float edge_velocity_integral = tet_edge_velocities[tet_edges[i]] * mesh->edge_lengths[tet_edges[i]];
      total_velocity += (grad_nb * na - grad_na * nb) * edge_velocity_integral;
      
   }
   
   float vol = fabs( tet_signed_volume( mesh->vertices[tet[0]], mesh->vertices[tet[1]], mesh->vertices[tet[2]], mesh->vertices[tet[3]] ) );
   
   return total_velocity / (6.0f*vol);
   
}


// ---------------------------------------------------------
///
/// Construct a 3D velocity vector for each tet vertex from its incident edges.
/// (Dual: Construct a 3D velocity vector for each voronoi site from its incident voronoi faces.)
///
// ---------------------------------------------------------

void DualFluidSim3D::reset_vertex_velocities() {
   tet_vertex_velocities.clear();
   tet_vertex_velocities.resize( mesh->vertices.size() );

   tet_vertex_velocity_is_valid.clear();
   tet_vertex_velocity_is_valid.resize(mesh->vertices.size(), false);

   for(unsigned int i = 0; i < mesh->vertices.size(); ++i) 
   {
      tet_vertex_velocity_is_valid[i] = false;
      tet_vertex_velocities[i] = Vec3f(0,0,0);
   }

}

void DualFluidSim3D::tet_edge_to_vertex_velocities( )
{
   const int MAX_NBRS = 500;
   Mat<MAX_NBRS, 3, double> N;
   Vec<MAX_NBRS, double> z;
   
   
   reset_vertex_velocities();

   for(unsigned int i = 0; i < mesh->vertices.size(); ++i) 
   {
      
      //if(liquid_phi[i] >= 0)
      if(!is_liquid(i)) 
      {
         tet_vertex_velocities[i] = Vec3f(0,0,0);
         tet_vertex_velocity_is_valid[i] = false;
         continue;
      }
      
      const std::vector<unsigned int>& incident_edges = mesh->vert_to_edge_map[i];
      int num_incident_edges = incident_edges.size();
      
      if(num_incident_edges > MAX_NBRS) {
         printf("Too many incident edges: %d\n", num_incident_edges );
         assert(0);
      }
      
      // reset the matrix
      for(int j = 0; j < MAX_NBRS; ++j) 
      {
         N(j, 0) = 0;
         N(j, 1) = 0;
         N(j, 2) = 0;
         z[j] = 0;
      }
            
      int c = 0;
      
      for(int j = 0; j < num_incident_edges; ++j) 
      {
         float cur_edge_velocity = tet_edge_velocities[incident_edges[j]];
         Vec3f normal = mesh->tet_edge_vectors[incident_edges[j]];
         
         assert( fabs( mag(normal) - 1.0f ) < 1e-6 );
         
         if( solid_weights[incident_edges[j]] > 1e-7) 
         {
            N(c, 0) = normal[0];
            N(c, 1) = normal[1];
            N(c, 2) = normal[2];
            z[c] = cur_edge_velocity;
            ++c;
         }
      }
      
      if( c >= 3 && least_squares(N,z) ) 
      { 
         tet_vertex_velocity_is_valid[i] = true;
         tet_vertex_velocities[i] = Vec3f( (float)z[0], (float)z[1], (float)z[2] );
      }
      else 
      {
         tet_vertex_velocity_is_valid[i] = false;
         tet_vertex_velocities[i] = Vec3f(0,0,0);
      }
      
//      if ( mag(tet_vertex_velocities[i]) > 1e+7f )
//      {
//         std::cout << "huge reconstructed velocity" << std::endl;
//         std::cout << "normals: " << std::endl;
//         for(int j = 0; j < num_incident_edges; ++j) 
//         {
//            Vec3f normal = mesh.tet_edge_vectors[incident_edges[j]];
//            if( solid_weights[incident_edges[j]] > 1e-7) 
//            {
//               std::cout << normal << std::endl;
//             }
//         }
//
//         std::cout << "velocities: " << std::endl;
//         for(int j = 0; j < num_incident_edges; ++j) 
//         {
//            if( solid_weights[incident_edges[j]] > 1e-7) 
//            {
//               std::cout << tet_edge_velocities[incident_edges[j]] << std::endl;
//            }
//            else
//            {
//               std::cout << "solid: " << tet_edge_velocities[incident_edges[j]] << std::endl;
//            }
//         }
//         
//         std::cout << "N: " << std::endl;
//         for ( int j = 0; j < c; ++j )
//         {
//            std::cout << N(j,0) << " " << N(j,1) << " " << N(j,2) << std::endl;
//         }
//         std::cout << "z: " << std::endl;
//         for ( int j = 0; j < c; ++j )
//         {
//            std::cout << z[j] << std::endl;
//         }
//       
//         Vec3f box_centre = 0.5f * ( solid_low + solid_high );
//         Vec3f box_extents = solid_high - solid_low;
//
//         std::cout << "phi: " << solid_box_phi( box_centre, box_extents, mesh.vertices[i] ) << std::endl;
//         
//         assert(0);
//         
//      }
      
   }
   
   
   float max_velocity = -1.0f;
   
   for(unsigned int i = 0; i < mesh->vertices.size(); ++i) 
   {
      assert( tet_vertex_velocities[i][0] == tet_vertex_velocities[i][0] );
      assert( tet_vertex_velocities[i][1] == tet_vertex_velocities[i][1] );
      assert( tet_vertex_velocities[i][2] == tet_vertex_velocities[i][2] );
      
      max_velocity = max( max_velocity, mag( tet_vertex_velocities[i] ) );
   }
   
   std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% max vertex velocity: " << max_velocity << std::endl;
   
   
}


// ---------------------------------------------------------
///
/// Construct a 3D velocity vector for each tet vertex from its 4 incident edges.
/// (Dual: Construct a 3D velocity vector for each tet circumcentre from its 4 incident tet faces.)
///
// ---------------------------------------------------------

void DualFluidSim3D::tet_edge_to_circumcentre_velocities( )
{
   const int MAX_NBRS = 6;
   Mat<MAX_NBRS, 3, double> N;
   Vec<MAX_NBRS, double> z;
   
   voronoi_vertex_velocities.clear();
   voronoi_vertex_velocities.resize( mesh->tets.size() );
   
   voronoi_vertex_velocity_is_valid.clear();
   voronoi_vertex_velocity_is_valid.resize(mesh->tets.size(), false);
 
   
   for(unsigned int i = 0; i < mesh->tets.size(); ++i) 
   {
      //if any of the four vertices (voronoi sites/pressure samples) are liquid, then go ahead
      Vec4ui nbr_voronoi_cells = mesh->tets[i];
      int liquid_nbr_count = 0;
      for(int ind = 0; ind < 4; ++ind) {
         //liquid_nbr_count += (liquid_phi[nbr_voronoi_cells[ind]] < 0? 1: 0);
        liquid_nbr_count += (is_liquid(nbr_voronoi_cells[ind])?1:0);
      }
      
      //if no neighbouring voronoi cells are liquid, skip this voronoi vertex since it's not valid
      if(liquid_nbr_count == 0) 
      {
         voronoi_vertex_velocities[i] = Vec3f(0,0,0);
         voronoi_vertex_velocity_is_valid[i] = false;
         continue;
      }
      
      Vec6ui& incident_faces = mesh->tet_to_edge_map[i];
      int num_incident_faces = 6;
      
      float max_component = 0;

      // reset the matrix
      for(int j = 0; j < MAX_NBRS; ++j) 
      {
         N(j, 0) = 0;
         N(j, 1) = 0;
         N(j, 2) = 0;
         z[j] = 0;
      }
            
      int c = 0;
      
      for(int j = 0; j < num_incident_faces; ++j) 
      {
         int faceID = incident_faces[j];
         float cur_edge_velocity = tet_edge_velocities[faceID];
         Vec3f normal = mesh->tet_edge_vectors[faceID];
         
         //check if the edge is bordering on liquid, and if not, skip it
         Vec2ui nbrVertexIDs = mesh->edges[faceID];
         //if(liquid_phi[nbrVertexIDs[0]] >= 0 && liquid_phi[nbrVertexIDs[1]] >= 0)
         if(!is_liquid(nbrVertexIDs[0]) && !is_liquid(nbrVertexIDs[1]))
          continue;

         if( solid_weights[faceID] > 1e-7) 
         {
            N(c, 0) = normal[0];
            N(c, 1) = normal[1];
            N(c, 2) = normal[2];
            z[c] = cur_edge_velocity;
            max_component = max(max_component, (float)fabs(cur_edge_velocity));
            ++c;
         }
      }
      
      if( c >= 3 && least_squares(N,z) )  
      {
//         bool success = least_squares(N,z);
//         if( !success || z[0] != z[0] || z[1] != z[1] || z[2] != z[2] ) 
//         {
//           
//            std::cout << "N: " << std::endl;
//            for ( int j = 0; j < c; ++j )
//            {
//               std::cout << N(j,0) << " " << N(j,1) << " " << N(j,2) << std::endl;
//               
//               //Vec3f normal = mesh.tet_edge_vectors[incident_faces[j]];
//               //std::cout << "Vector: " << normal << " with velocity: " << tet_edge_velocities[incident_faces[j]] << std::endl;
//               
//            }
//            std::cout << "z: " << std::endl;
//            for ( int j = 0; j < c; ++j )
//            {
//               std::cout << z[j] << std::endl;
//            }
//            voronoi_vertex_velocity_is_valid[i] = false;
//            voronoi_vertex_velocities[i] = Vec3f(0,0,0);
//            exit(1);//continue;
//            
//         }
         Vec3f new_velocity = Vec3f( (float)z[0], (float)z[1], (float)z[2] );
         
         float max_expected = sqrt(3*sqr(max_component));
         if(mag(new_velocity) > 5*max_expected) {
            printf("\n\n\n***Warning: BAD VELOCITY***\nMax expected velocity magnitude: %f\nConstructed velocity: %f\n", max_expected, mag(new_velocity));
            std::cout << "Offending tetrahedron: " << i << std::endl;
            std::cout << "Final computed velocity: " << new_velocity << std::endl;
            std::cout << "Spatial position" << mesh->tet_circumcentres[i] << std::endl;
            for(int ver = 0; ver < 4; ++ver) {
               std::cout << "Vert pos: " << mesh->vertices[mesh->tets[i][ver]] << std::endl;
            }
            for(int j = 0; j < MAX_NBRS; ++j) {
               N(j, 0) = 0;
               N(j, 1) = 0;
               N(j, 2) = 0;
               z[j] = 0;
            }
                  
            int c = 0;
            
            for(int j = 0; j < num_incident_faces; ++j) 
            {
               int faceID = incident_faces[j];
               float cur_edge_velocity = tet_edge_velocities[faceID];
               Vec3f normal = mesh->tet_edge_vectors[faceID];
               
               //check if the edge is bordering on liquid, and if not, skip it
               Vec2ui nbrVertexIDs = mesh->edges[faceID];
               if(!is_liquid(nbrVertexIDs[0]) && !is_liquid(nbrVertexIDs[1])) {
                  std::cout << "Skipping because both ends are not liquid\n";
                  continue;
               }

               if( solid_weights[faceID] > 1e-7) 
               {
                  N(c, 0) = normal[0];
                  N(c, 1) = normal[1];
                  N(c, 2) = normal[2];
                  z[c] = cur_edge_velocity;
                  
                  ++c;
                  std::cout << "Normal: " << normal << " Velocity: " << cur_edge_velocity << std::endl;
               }
               else {
                  std::cout << "Skipping because the face is solid\n";
               }
            }         
         }
         
         voronoi_vertex_velocity_is_valid[i] = true;
         voronoi_vertex_velocities[i] = new_velocity;
      }
      else 
      {
         voronoi_vertex_velocity_is_valid[i] = false;
         voronoi_vertex_velocities[i] = Vec3f(0,0,0);
      }      
   }
   

   float max_velocity = -1.0f;
   
   for(unsigned int i = 0; i < mesh->tets.size(); ++i) 
   {
      assert( voronoi_vertex_velocities[i][0] == voronoi_vertex_velocities[i][0] );
      assert( voronoi_vertex_velocities[i][1] == voronoi_vertex_velocities[i][1] );
      assert( voronoi_vertex_velocities[i][2] == voronoi_vertex_velocities[i][2] );
      
      max_velocity = max( max_velocity, mag( voronoi_vertex_velocities[i] ) );
   }
   
   std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% max voronoi vertex velocity: " << max_velocity << std::endl;
   
   
}

// ---------------------------------------------------------
///
/// Assuming voronoi vertices (tet circumcentres) have velocities computed as above, 
/// construct full velocities on voronoi faces (tet edges)
///
// ---------------------------------------------------------

void DualFluidSim3D::construct_full_voronoi_face_velocities() {
   //voronoi faces already have an exact normal component
   //so we will just construct the tangential component
   //from the velocities on the voronoi face's vertices
   
   full_voronoi_face_velocities.resize(mesh->edges.size(), Vec3f(0,0,0) );
   full_voronoi_face_velocity_is_valid.resize(mesh->edges.size(), false );

   for(unsigned int i = 0; i < mesh->edges.size(); ++i) 
   {
      
      //if the edge contains no fluid region, skip it
      if ( solid_weights[i] <= 1e-7 ) { 
         full_voronoi_face_velocities[i] = Vec3f(0,0,0);
         full_voronoi_face_velocity_is_valid[i] = false;
         continue; 
      }
         
      //check if valid, ie. one of the two adjacent cells is liquid
      Vec2ui edge_nbrs = mesh->edges[i];
      //if(liquid_phi[edge_nbrs[0]] >= 0 && liquid_phi[edge_nbrs[1]] >= 0) {
      if(!is_liquid(edge_nbrs[0]) && !is_liquid(edge_nbrs[1])) {
         full_voronoi_face_velocities[i] = Vec3f(0,0,0);
         full_voronoi_face_velocity_is_valid[i] = false;
         continue;
      }

      //determine all the voronoi vertices (tet circumcentres) nearby
      //and find the average of their velocities
      Vec3f sum(0,0,0);
      std::vector<unsigned int> nbr_tets = mesh->edge_to_tet_map[i];
      unsigned int num_valid_voronoi_vertices = 0;
      for(unsigned int nbr = 0; nbr < nbr_tets.size(); ++nbr) 
      {
         int nbr_ind = nbr_tets[nbr];
         if ( voronoi_vertex_velocity_is_valid[nbr_ind] )
         {         
            sum += voronoi_vertex_velocities[nbr_ind];
            ++num_valid_voronoi_vertices;
         }
      }
      
      if ( num_valid_voronoi_vertices == 0 ) { continue; }
         
      Vec3f result = sum / (float)num_valid_voronoi_vertices; //just take the average, no real need to be "smarter" I don't think.

      
      //strip out just the tangential component, since we know the normal component already
      Vec3f face_normal = mesh->tet_edge_vectors[i];
      normalize(face_normal);
      Vec3f tangential_component = result - dot(result,face_normal)*face_normal;
      //construct result velocity as sum of known normal component, and constructed tangential component
      Vec3f full_velocity = tet_edge_velocities[i] * face_normal + tangential_component;
      
      full_voronoi_face_velocities[i] = full_velocity;//result;
      full_voronoi_face_velocity_is_valid[i] = true;
   }

}


// ---------------------------------------------------------
///
/// For every tet vertex inside the solid, project out the component normal to the solid wall.
///
// ---------------------------------------------------------

void DualFluidSim3D::clamp_solid_tet_vertex_velocities()
{
   Vec3f box_centre = 0.5f * ( solid_low + solid_high );
   Vec3f box_extents = solid_high - solid_low;
   
   for ( unsigned int i = 0; i < mesh->vertices.size(); ++i ) 
   {
      float solid_phi = solid_box_phi( box_centre, box_extents, mesh->vertices[i] );
      if ( solid_phi > 0.0f )
      {
         Vec3f normal = solid_box_gradient( box_centre, box_extents, mesh->vertices[i] );
         tet_vertex_velocities[i] -= normal * dot( normal, tet_vertex_velocities[i] );
      }
   }   
}

// ---------------------------------------------------------
///
/// For every tet vertex / voronoi vertex / voronoi face inside the solid, project out the component normal to the solid wall.
///
// ---------------------------------------------------------

void DualFluidSim3D::clamp_solid_all_velocities()
{
   Vec3f box_centre = 0.5f * ( solid_low + solid_high );
   Vec3f box_extents = solid_high - solid_low;
   
   for ( unsigned int i = 0; i < mesh->vertices.size(); ++i ) 
   {
      float solid_phi = solid_box_phi( box_centre, box_extents, mesh->vertices[i] );
      if ( solid_phi > 0.0f)
      {
         Vec3f normal = solid_box_gradient( box_centre, box_extents, mesh->vertices[i] );
         tet_vertex_velocities[i] -= normal * dot( normal, tet_vertex_velocities[i] );
      }
   }

   for ( unsigned int i = 0; i < mesh->tets.size(); ++i ) 
   {
      float solid_phi = solid_box_phi( box_centre, box_extents, mesh->tet_circumcentres[i] );
      if ( solid_phi > 0.0f )
      {
         Vec3f normal = solid_box_gradient( box_centre, box_extents, mesh->tet_circumcentres[i] );
         voronoi_vertex_velocities[i] -= normal * dot( normal, voronoi_vertex_velocities[i] );
      }
   }

   for ( unsigned int i = 0; i < mesh->edges.size(); ++i ) 
   {
      float solid_phi = solid_box_phi( box_centre, box_extents, mesh->voronoi_face_centroids[i] );
      if ( solid_phi > 0.0f )
      {
         Vec3f normal = solid_box_gradient( box_centre, box_extents, mesh->voronoi_face_centroids[i] );
         full_voronoi_face_velocities[i] -= normal * dot( normal, full_voronoi_face_velocities[i] );
      }
   }

}


// ---------------------------------------------------------
///
/// Take the velocities on the tet mesh vertices in the liquid, extrapolate out to the surrounding non-liquid vertices.
///
// ---------------------------------------------------------

void DualFluidSim3D::extrapolate_tet_vertex_velocities()
{
   
   //reset outside velocities to zero
   for ( unsigned int i = 0; i < mesh->vertices.size(); ++i ) 
   {
      if ( !tet_vertex_velocity_is_valid[i] )
      {
         tet_vertex_velocities[i] = Vec3f(0,0,0);
      }
   }
   
   std::vector<bool> new_valid = tet_vertex_velocity_is_valid;
   
   for( unsigned int pass = 0; pass < 5; ++pass ) 
   {
      for(unsigned int i = 0; i < mesh->vertices.size(); ++i) 
      {
         if ( tet_vertex_velocity_is_valid[i] ) { continue; }
         
         const std::vector<unsigned int>& incident_edges = mesh->vert_to_edge_map[i];
         
         float min_neighbour_distance = 1e30f;
         unsigned int closest_valid_neighbour = (unsigned int) ~0;
         
         for( unsigned int j = 0; j < incident_edges.size(); ++j) 
         {
            const Vec2ui& edge = mesh->edges[incident_edges[j]];
            if ( edge[0] == edge[1] ) { continue; }

            unsigned int neighbour = ( edge[0] == i ? edge[1] : edge[0] );
            
            if ( tet_vertex_velocity_is_valid[neighbour] ) 
            {
               if ( dist( mesh->vertices[i], mesh->vertices[neighbour] ) < min_neighbour_distance )
               {
                  min_neighbour_distance = dist( mesh->vertices[i], mesh->vertices[neighbour] );
                  closest_valid_neighbour = neighbour;
               }
            } 
         }
         
         if ( closest_valid_neighbour != (unsigned int) ~0 )
         {
            tet_vertex_velocities[i] = tet_vertex_velocities[closest_valid_neighbour];
            new_valid[i] = true;
         }
         
      }
      
      tet_vertex_velocity_is_valid = new_valid;      
   }
   
}


// ---------------------------------------------------------
///
/// Take the velocities on voronoi sites, vertices and faces in the liquid, 
/// and extrapolate out to the surrounding non-liquid data points.
///
// ---------------------------------------------------------

void DualFluidSim3D::extrapolate_all_velocities()
{
   
   //reset all outside velocities to zero
   for ( unsigned int i = 0; i < mesh->vertices.size(); ++i ) 
   {
      if ( !tet_vertex_velocity_is_valid[i] )
      {
         tet_vertex_velocities[i] = Vec3f(0,0,0);
      }
   }
   for ( unsigned int i = 0; i < mesh->tets.size(); ++i ) 
   {
      if ( !voronoi_vertex_velocity_is_valid[i] )
      {
         voronoi_vertex_velocities[i] = Vec3f(0,0,0);
      }
   }
   for ( unsigned int i = 0; i < mesh->edges.size(); ++i ) 
   {
      if ( !full_voronoi_face_velocity_is_valid[i] )
      {
         full_voronoi_face_velocities[i] = Vec3f(0,0,0);
      }
   }


   std::vector<bool> new_tet_vert_valid = tet_vertex_velocity_is_valid;
   std::vector<bool> new_vor_vert_valid = voronoi_vertex_velocity_is_valid;
   std::vector<bool> new_vor_face_valid = full_voronoi_face_velocity_is_valid;

   for( unsigned int pass = 0; pass < 12; ++pass ) 
   {
      
      //vor face + vor vert --> tet vert
      for(unsigned int i = 0; i < mesh->vertices.size(); ++i) 
      {
         if ( tet_vertex_velocity_is_valid[i] ) 
            continue;
         
         //determine incident voronoi faces and vertices

         const std::vector<unsigned int>& incident_edges = mesh->vert_to_edge_map[i];
         const std::vector<unsigned int>& incident_tets = mesh->vert_to_tet_map[i];

         Vec3f new_velocity_sum(0,0,0);
         int valid_nbrs = 0;

         for( unsigned int j = 0; j < incident_edges.size(); ++j) 
         {
            if(full_voronoi_face_velocity_is_valid[incident_edges[j]]) {
               new_velocity_sum += full_voronoi_face_velocities[incident_edges[j]];
               ++valid_nbrs;
            }
         }
         
         for(unsigned int j = 0; j < incident_tets.size(); ++j) {
            if(voronoi_vertex_velocity_is_valid[incident_tets[j]]) {
               new_velocity_sum += voronoi_vertex_velocities[incident_tets[j]];
               ++valid_nbrs;
            }
         }

         if ( valid_nbrs > 0 )
         {
            tet_vertex_velocities[i] = new_velocity_sum / (float)valid_nbrs;
            new_tet_vert_valid[i] = true;
         }
      }

      //vor face + tet vert --> voronoi vert
      for(unsigned int i = 0; i < mesh->tets.size(); ++i) {
         if(voronoi_vertex_velocity_is_valid[i]) 
            continue;
         
         Vec6ui incident_faces = mesh->tet_to_edge_map[i];
         Vec4ui incident_verts = mesh->tets[i];

         Vec3f new_velocity_sum(0,0,0);
         int valid_nbrs = 0;

         for(unsigned int j = 0; j < 6; ++j) {
            if(full_voronoi_face_velocity_is_valid[incident_faces[j]]) {
               new_velocity_sum += full_voronoi_face_velocities[incident_faces[j]];
               ++valid_nbrs;
            }
         }

         for(unsigned int j = 0; j < 4; ++j) {
            if(tet_vertex_velocity_is_valid[incident_verts[j]]) {
               new_velocity_sum += tet_vertex_velocities[incident_verts[j]];
               ++valid_nbrs;
            }

         }

         if ( valid_nbrs > 0 )
         {
            voronoi_vertex_velocities[i] = new_velocity_sum / (float)valid_nbrs;
            new_vor_vert_valid[i] = true;
         }
      }
      
      //tet vert + vor vert --> vor face
      for(unsigned int i = 0; i < mesh->edges.size(); ++i) {
         if(full_voronoi_face_velocity_is_valid[i])
            continue;
         
         Vec2ui incident_verts = mesh->edges[i];
         const std::vector<unsigned int>& incident_tets = mesh->edge_to_tet_map[i];
         
         Vec3f new_velocity_sum(0,0,0);
         int valid_nbrs = 0;

         for(unsigned int j = 0; j < 2; ++j) {
            if(tet_vertex_velocity_is_valid[incident_verts[j]]) {
               new_velocity_sum += tet_vertex_velocities[incident_verts[j]];
               ++valid_nbrs;
            }
         }

         for(unsigned int j = 0; j < incident_tets.size(); ++j) {
            if(voronoi_vertex_velocity_is_valid[incident_tets[j]]) {
               new_velocity_sum += voronoi_vertex_velocities[incident_tets[j]];
               ++valid_nbrs;
            }
         }

         if ( valid_nbrs > 0 )
         {
            full_voronoi_face_velocities[i] = new_velocity_sum / (float)valid_nbrs;
            new_vor_face_valid[i] = true;
         }
      }

      //copy over new validity information, and repeat
      tet_vertex_velocity_is_valid = new_tet_vert_valid;
      voronoi_vertex_velocity_is_valid = new_vor_vert_valid;
      full_voronoi_face_velocity_is_valid = new_vor_face_valid;

   }
   
}


// ---------------------------------------------------------
///
/// Compute the maximum timestep length.
///
// ---------------------------------------------------------

float DualFluidSim3D::get_cfl_limit()
{
   float max_vel = 0;
  /* for(unsigned int i = 0; i < tet_edge_velocities.size(); ++i)
   {
      max_vel = max(max_vel, tet_edge_velocities[i]);
   }*/
   
   for(unsigned int i = 0; i < tet_vertex_velocities.size(); ++i) {
      max_vel = max(max_vel, mag(tet_vertex_velocities[i]));
   }

   if ( interpolation_scheme == GENERALIZED_BARYCENTRIC || interpolation_scheme == IMPROVED_BARYCENTRIC )
   {
      for(unsigned int i = 0; i < voronoi_vertex_velocities.size(); ++i) {
         max_vel = max(max_vel, mag(voronoi_vertex_velocities[i]));
      }
      for(unsigned int i = 0; i < full_voronoi_face_velocities.size(); ++i) {
         max_vel = max(max_vel, mag(full_voronoi_face_velocities[i]));
      }
   }

   if ( surface_tension_coefficient != 0 )
   {
      // Not the real CFL condition for surface tension, but it doesn't seem to blow up...
      
      return 0.33f*characteristic_distance / max_vel;
   }
   else
   {
      return 0.33f*characteristic_distance / max_vel;
   }
}


// ---------------------------------------------------------
///
/// Extrapolate tangential edge velocities from inside the fluid to the surrounding air 
///
// ---------------------------------------------------------

void DualFluidSim3D::extrapolate_tet_edge_velocities()
{
   
   // mark valid edges (we know the velocities are right)
   // use the same creteria as the pressure solver
   
   edge_is_valid.clear( );
   edge_is_valid.resize( mesh->edges.size(), false );
   
   for(unsigned int i = 0; i < mesh->edges.size(); ++i) 
   {
      Vec2ui verts = mesh->edges[i];
      
      if( solid_weights[i] > 0 ) 
      {
         float phi0 = liquid_phi[verts[0]];
         float phi1 = liquid_phi[verts[1]];
         if( phi0 < 0 || phi1 < 0) 
         {
            edge_is_valid[i] = true;
         }
      }
   }
   
   // rebuild and extrapolate vertex velocities
   tet_edge_to_vertex_velocities();      
   extrapolate_tet_vertex_velocities();
   clamp_solid_tet_vertex_velocities();       

   unsigned int num_edges_extrapolated_to = 0;
   
   // look for invalid edges which have valid vertex velocities
   for(unsigned int i = 0; i < mesh->edges.size(); ++i) 
   {      
      if ( edge_is_valid[i] ) { continue; }
      
      Vec2ui verts = mesh->edges[i];
      
      // take the average of the velocities at the end points, if either one is valid
      
      unsigned int num_valid_end_vertices = 0;
      Vec3f sum_end_vertex_velocities(0,0,0);
      
      if ( tet_vertex_velocity_is_valid[verts[0]] )
      {
         ++num_valid_end_vertices;
         sum_end_vertex_velocities += tet_vertex_velocities[verts[0]];
      }
      
      if ( tet_vertex_velocity_is_valid[verts[1]] )
      {
         ++num_valid_end_vertices;
         sum_end_vertex_velocities += tet_vertex_velocities[verts[1]];
      }
      
      if ( num_valid_end_vertices > 0 )
      {
         edge_is_valid[i] = true;
         sum_end_vertex_velocities /= (float)num_valid_end_vertices;
         tet_edge_velocities[i] = dot( sum_end_vertex_velocities, mesh->tet_edge_vectors[i] );
         ++num_edges_extrapolated_to;
      }
      else
      {
         tet_edge_velocities[i] = 0.0f;
      }
   }
   
   std::cout << "num_edges_extrapolated_to: " << num_edges_extrapolated_to << std::endl;
   
}


// ---------------------------------------------------------
///
///
///
// ---------------------------------------------------------

void DualFluidSim3D::reconstruct_and_extrapolate_velocities()
{
   
   if ( interpolation_scheme == WHITNEY )
   {
      extrapolate_tet_edge_velocities();
   }
   else if ( interpolation_scheme == GENERALIZED_BARYCENTRIC || interpolation_scheme == IMPROVED_BARYCENTRIC )
   {
      std::cout << "---------------------- Voronoi Fluid Sim: Transferring velocities to high resolution vertices ----------------------" << std::endl;
      
      
      printf("...to tet circumcentres.\n");
      tet_edge_to_circumcentre_velocities();
      
      printf("...to tet edges.\n");
      construct_full_voronoi_face_velocities();
      
      printf("...to tet vertices.\n"); 
      reset_vertex_velocities(); //just reset these to empty, and extrapolation will fill them in from the neighbours.
      
      std::cout << "---------------------- Voronoi Fluid Sim: Extrapolating high resolution vertex velocities ----------------------" << std::endl;
      extrapolate_all_velocities();
      clamp_solid_all_velocities(); 
   }
   else
   {
      std::cout << "---------------------- Voronoi Fluid Sim: Transferring velocities to nodes ----------------------" << std::endl;
      tet_edge_to_vertex_velocities();      
      std::cout << "---------------------- Voronoi Fluid Sim: Extrapolating nodal velocities ----------------------" << std::endl;
      extrapolate_tet_vertex_velocities();
      clamp_solid_tet_vertex_velocities();     
   }
   
}

// ---------------------------------------------------------
///
///
///
// ---------------------------------------------------------

void DualFluidSim3D::correct_volume( )
{
   std::cout << "---------------------- Voronoi Fluid Sim: Correct volume ----------------------" << std::endl;
   
   double predicted_volume = surface_tracker->get_predicted_volume();
   double area = surface_tracker->get_predicted_surface_area();

   double dV = ( initial_volume - predicted_volume ) / area;
   
   for ( unsigned int i = 0; i < surface_tracker->get_num_vertices(); ++i )
   {
      if( surface_tracker->m_mesh.m_vertex_to_edge_map[i].size() == 0 ) { continue; }
      if( surface_tracker->m_masses[i] > 1.5 ) { continue; } //TODO: Use El Topo's actual constraint mechanism instead?
      
      Vec3d n = surface_tracker->get_vertex_normal( i );
      
      surface_tracker->set_newposition(i, surface_tracker->get_newposition(i) + 0.5 * dV * n);
      //surface_tracker->m_newpositions[i] += 0.5 * dV * n;
   }
   
}


// ---------------------------------------------------------
///
/// Just in case we get desperate.
///
// ---------------------------------------------------------

//void DualFluidSim3D::surface_laplacian_smoothing( double coefficient )
//{
//   std::vector<Vec3d> new_positions;
//   
//   for ( unsigned int i = 0; i < surface_tracker->m_newpositions.size(); ++i )
//   {
//      const std::vector<unsigned int>& incident_edges = surface_tracker->m_mesh.m_vtxedge[i];
//      Vec3d average_position(0,0,0);
//      unsigned int num_neighbours = 0;
//      for ( unsigned int j = 0; j < incident_edges.size(); ++j )
//      {
//         unsigned int neighbour = surface_tracker->m_mesh.m_edges[ incident_edges[j] ][0];
//         if ( neighbour == i ) { neighbour = surface_tracker->m_mesh.m_edges[ incident_edges[j] ][1]; }
//         average_position += surface_tracker->m_newpositions[neighbour];
//         ++num_neighbours;
//      }
//      average_position /= (double)num_neighbours;     
//      
//      new_positions[i] = surface_tracker->m_newpositions[i] + coefficient * ( average_position - surface_tracker->m_newpositions[i] );
//   }
//
//}


// ---------------------------------------------------------
///
/// Advance surface mesh / marker particles
///
// ---------------------------------------------------------

void DualFluidSim3D::advance_surface( float dt )
{
   
   VelocityFunctor3D* get_velocity;
   switch( interpolation_scheme )
   {
      case WHITNEY:
         get_velocity = new WhitneyEdgeVelocityFunctor( *this );
         break;
      case IMPROVED_BARYCENTRIC:
         get_velocity = new SharperBarycentricVelocityFunctor( *this);
         break;
      case GENERALIZED_BARYCENTRIC:
         get_velocity = new GeneralizedBarycentricVelocityFunctor( *this);
         break;
      case BARYCENTRIC:
         get_velocity = new BarycentricTetVelocityFunctor( *this );
         break;
      default:
         assert( !"Invalid interpolation scheme specified" );
   }
   
   // markers
   for ( unsigned int i = 0; i < markers.size(); ++i )
   {
      trace_rk2( markers[i], markers[i], dt, *get_velocity );
   }
   
   // El Topo: static operations
   
   std::vector<bool> vert_const_labels(surface_tracker->m_mesh.m_vertex_to_edge_map.size(), 0);
   surface_tracker->m_mesh.m_vertex_constraint_labels = vert_const_labels;
   std::vector<Vec3d> vert_vel(surface_tracker->get_num_vertices());
   for(unsigned int i = 0; i < vert_vel.size(); ++i) {
      vert_vel[i] = Vec3d((*get_velocity)(Vec3f(surface_tracker->get_position(i))));
   }
   surface_tracker->set_all_remesh_velocities(vert_vel);
   surface_tracker->assert_no_bad_labels();
   
   surface_tracker->improve_mesh();
   surface_tracker->assert_no_bad_labels();
   
   surface_tracker->topology_changes();
   surface_tracker->assert_no_bad_labels();
      
   // El Topo: surface advection

   Vec3f box_centre = 0.5f * ( solid_low + solid_high );
   Vec3f box_extents = solid_high - solid_low;
   std::vector<Vec3d> new_positions(surface_tracker->get_num_vertices());
   
   //surface_tracker->m_newpositions.resize( surface_tracker->get_num_vertices() );
   for ( unsigned int i = 0; i < surface_tracker->get_num_vertices(); ++i )
   {
      Vec3f new_position;
      trace_rk2( Vec3f(surface_tracker->get_position(i)), new_position, dt, *get_velocity );

      float vertex_solid_phi = solid_box_phi( box_centre, box_extents, new_position );
      
      // mark vertices against the solid wall
      if ( vertex_solid_phi > -1e-4 ) //TODO Make this somehow a tunable parameter or something...
      {
         surface_tracker->m_masses[i] = 2.0; //TODO use constraint mechanism instead?
      }
      
      if ( !allow_solid_overlap )
      {
         // snap to solid
         if ( surface_tracker->m_masses[i] > 1.0 ) //TODO use constraint mechanism instead?
         {
            new_position -= vertex_solid_phi * solid_box_gradient( box_centre, box_extents, new_position );
            static const float snap_distance = 1e-4f;
            Vec3f v(surface_tracker->get_position(i));
            if ( fabs( v[0] - solid_low[0] ) < snap_distance ) { new_position[0] = solid_low[0]; }
            if ( fabs( v[1] - solid_low[1] ) < snap_distance ) { new_position[1] = solid_low[1]; }
            if ( fabs( v[2] - solid_low[2] ) < snap_distance ) { new_position[2] = solid_low[2]; }
            if ( fabs( v[0] - solid_high[0] ) < snap_distance ) { new_position[0] = solid_high[0]; }
            if ( fabs( v[1] - solid_high[1] ) < snap_distance ) { new_position[1] = solid_high[1]; }
            if ( fabs( v[2] - solid_high[2] ) < snap_distance ) { new_position[2] = solid_high[2]; }      
         }
      }
      new_positions[i] = Vec3d(new_position);
      //surface_tracker->set_newposition(i, (Vec3d) new_position);
   }
   surface_tracker->set_all_newpositions(new_positions);

   
   if ( volume_correction )
   {
      //TODO Make volume correction work with multiphase mesh.
      //correct_volume( );
   }
   
   write_binary_file_with_newpositions( surface_tracker->m_mesh, 
                                        surface_tracker->get_positions(),
                                        surface_tracker->m_masses, 
                                        surface_tracker->get_newpositions(),
                                        dt,
                                        "/Users/tyson/scratch/pre-integration.bin" );
   double actual_dt;
   surface_tracker->integrate( dt, actual_dt );
   surface_tracker->assert_no_bad_labels();
   delete get_velocity;
}


// ---------------------------------------------------------
///
/// Generate a new tet mesh.  Get velocities for new tet mesh by using semi-Lagrangian advection looking into the old mesh.
///
// ---------------------------------------------------------

void DualFluidSim3D::remesh_and_advect_semilagrangian( float dt )
{
   
   double start_time = get_time_in_seconds();
   
   std::vector<Vec3f> input_xs;
   
   // Surface-adaptive points
   SampleSeeder::generate_adaptive_points( *surface_tracker, 0.5 * characteristic_distance, input_xs );
      
   
     
   // get delaunay mesh
   
   std::vector<Vec4ui> tets;
   std::vector<Vec3f> xs;

   
    SampleSeeder::generate_bcc_points( domain_low, domain_high, characteristic_distance, input_xs );
      
    std::cout << "num total pressure samples: " << input_xs.size() << std::endl;

    Triangulation cgal_T;
    //Use CGAL Delaunay mesher
    compute_delaunay_CGAL(input_xs, tets, cgal_T);
    xs = input_xs;
      
    //Use TetGen Delaunay mesher
    //TetGen::delaunay_mesh( input_xs, xs, tets );
   
   // create a new mesh
   TetMesh* new_mesh = new TetMesh;
   new_mesh->initialize( tets, xs, cgal_T );
   
   // transfer velocity from current mesh to new mesh
   // combine this step with semi-lagrangian advection to reduce interpolation error
   
   total_remesh_time += get_time_in_seconds() - start_time;
   start_time = get_time_in_seconds();

   std::cout << "---------------------- Voronoi Fluid Sim: Semi-lagrangian advection ----------------------" << std::endl;

   VelocityFunctor3D* get_velocity;
   switch( interpolation_scheme )
   {
      case WHITNEY:
         std::cout << "Whitney-style interpolation\n";
         get_velocity = new WhitneyEdgeVelocityFunctor( *this );
         break;
      case IMPROVED_BARYCENTRIC:
         std::cout << "Improved barycentric interpolation\n";
         get_velocity = new SharperBarycentricVelocityFunctor( *this);
         break;
      case GENERALIZED_BARYCENTRIC:
         std::cout << "Generalized barycentric interpolation\n";
         get_velocity = new GeneralizedBarycentricVelocityFunctor( *this);
         break;
      case BARYCENTRIC:
         std::cout << "Basic barycentric interpolation\n";
         get_velocity = new BarycentricTetVelocityFunctor( *this );
         break;
      default:
         assert( !"Invalid interpolation scheme specified" );
   }
   
   std::vector<float> new_edge_velocities( new_mesh->edges.size() );
   
   for ( unsigned int i = 0; i < new_mesh->tet_edge_midpoints.size(); ++i )
   {

      Vec3f previous_location;
      trace_rk2( new_mesh->tet_edge_midpoints[i], previous_location, -dt, *get_velocity );
      
      Vec3f previous_velocity = (*get_velocity)( previous_location );
      
      float tangential_component = dot( previous_velocity, new_mesh->tet_edge_vectors[i] );
      assert( tangential_component == tangential_component );
      new_edge_velocities[i] = tangential_component;
   }  

   total_semilagrangian_time += get_time_in_seconds() - start_time;
   start_time = get_time_in_seconds();

   // swap in new mesh with new velocities
   
   TetMesh* old_mesh = mesh;
   mesh = new_mesh;
   delete old_mesh;

   assert( new_edge_velocities.size() == mesh->tet_edge_midpoints.size() );
   tet_edge_velocities = new_edge_velocities;   
   
   compute_solids();

   total_remesh_time += get_time_in_seconds() - start_time;
   
}


// ---------------------------------------------------------
///
/// Semi-Lagrangian advection
///
// ---------------------------------------------------------

void DualFluidSim3D::semi_lagrangian_advection( float dt )
{
   double start_time = get_time_in_seconds();

   std::vector<float> new_edge_velocities( mesh->edges.size() );

   VelocityFunctor3D* get_velocity;
   switch( interpolation_scheme )
   {
      case WHITNEY:
         get_velocity = new WhitneyEdgeVelocityFunctor( *this );
         break;
      case IMPROVED_BARYCENTRIC:
         get_velocity = new SharperBarycentricVelocityFunctor( *this);
         break;
      case GENERALIZED_BARYCENTRIC:
         get_velocity = new GeneralizedBarycentricVelocityFunctor( *this);
         break;
      case BARYCENTRIC:
         get_velocity = new BarycentricTetVelocityFunctor( *this );
         break;
      default:
         assert( !"Invalid interpolation scheme specified" );
   }
      
   for ( unsigned int i = 0; i < mesh->tet_edge_midpoints.size(); ++i )
   {
      Vec3f previous_location;
      trace_rk2( mesh->tet_edge_midpoints[i], previous_location, -dt, *get_velocity );
      
      Vec3f previous_velocity = (*get_velocity)( previous_location );
      
      float tangential_component = dot( previous_velocity, mesh->tet_edge_vectors[i] );
      assert( tangential_component == tangential_component );
      new_edge_velocities[i] = tangential_component;
   }  
   
   tet_edge_velocities = new_edge_velocities;
   
   total_semilagrangian_time += get_time_in_seconds() - start_time;

}


// ---------------------------------------------------------


void DualFluidSim3D::add_thermal_buoyancy( float dt )
{
   static const float hot_particle_radius = 0.2f;
   
   // determine temperature from markers
   std::vector<float> temperature( mesh->edges.size(), 0.0f );
   for ( unsigned int i = 0; i < markers.size(); ++i )
   {
      int tet_index = mesh->get_containing_tet( markers[i] );
      if ( tet_index < 0 ) { continue; }
      
      const Vec6ui& es = mesh->tet_to_edge_map[tet_index];
      
      for ( unsigned int e = 0; e < 6; ++e )
      {
         unsigned int edge_index = es[e];
         const Vec3f& edge_midpoint = mesh->tet_edge_midpoints[edge_index];
         if ( dist( edge_midpoint, markers[i] ) < hot_particle_radius )
         {
            temperature[edge_index] += 1.0;
         }
      }      
   }
        
   // then add upward force in hot regions
   for ( unsigned int i = 0; i < mesh->edges.size(); ++i )
   {
      if ( temperature[i] > 0.0f )
      {
         tet_edge_velocities[i] += dt * dot( mesh->tet_edge_vectors[i], Vec3f(0, 3, 0) );
      }
   }
   
}

// ---------------------------------------------------------


void DualFluidSim3D::add_local_force( float dt )
{
   // localized force
   const Vec3f force_centre = 0.5f * (domain_low + domain_high) - Vec3f( 0.0f, 1.0f, 0.0f );
   const float force_radius = 0.25f;
   static const Vec3f force( 0.0f, 10.0f, 0.0f );
   
   for(unsigned int i = 0; i < mesh->edges.size(); ++i) 
   {
      if ( dist( mesh->tet_edge_midpoints[i], force_centre ) < force_radius )
      {
         float dv = dot( dt * force, mesh->tet_edge_vectors[i] );
         tet_edge_velocities[i] += dv;
      }
   }   
}

// ---------------------------------------------------------
///
/// Apply gravity, etc. to the tet edges / voronoi faces.
///
// ---------------------------------------------------------

void DualFluidSim3D::add_forces( float dt )
{

   if ( free_surface )
   {
      static const bool slosh = false;
      
      if ( slosh )
      {
         static unsigned int call = 40;
         
         if ( call < 80 )
         {
            gravity = Vec3f( 0.5f, -1.0f, 0.0f );
         }
         else
         {
            gravity = Vec3f( 0.0f, -1.0f, 0.0f );
         }
         
         gravity = normalized( gravity );
         
         ++call;         
      }

      for(unsigned int i = 0; i < mesh->edges.size(); ++i) 
      {
         float dv = dot( dt * gravity, mesh->tet_edge_vectors[i] );
         tet_edge_velocities[i] += dv;
      }
   }
   else
   {      
      add_thermal_buoyancy( dt ); 
   }

}



// ---------------------------------------------------------
///
/// Compute the liquid signed distance function from the surface.
/// Evaluate at tet vertices / voronoi sites.
///
// ---------------------------------------------------------

unsigned int num_distance_to_surface_calls;
unsigned int num_broadphase_queries;
unsigned int num_point_triangle_tests;

void DualFluidSim3D::compute_liquid_phi( )
{
      
   liquid_phi.resize( mesh->vertices.size() );
   if ( !free_surface )
   {
      for ( unsigned int i = 0; i < liquid_phi.size(); ++i )
      {
         liquid_phi[i] = -1;
      }
      
      return;
   }
   
   double redistance_start = get_time_in_seconds();
   
   surface_tracker->rebuild_static_broad_phase();
   
   // redistancing: propagate nearest-triangle information from the surface outward
 
   liquid_phi.clear();
   liquid_phi.resize( mesh->vertices.size(), 1e+30f );
   
   std::vector<int> closest_triangle( mesh->vertices.size(), -1 );
   
   std::vector<char> is_on_boundary_tet(mesh->vertices.size(), 0);

   std::queue<unsigned int> vertex_queue;
   
   // first compute nearest-triangle in a band around the surface
   
   std::vector<unsigned int> vertex_band;

   num_distance_to_surface_calls = 0;
   num_broadphase_queries = 0;
   num_point_triangle_tests = 0;
   
   for ( unsigned int i = 0; i < surface_tracker->m_mesh.m_tris.size(); ++i )
   {
      // get nearby tets for this triangle
      Vec3d aabb_low, aabb_high;
      surface_tracker->triangle_static_bounds( i, aabb_low, aabb_high );

      std::vector<int> nearby_tets;
      mesh->get_overlapping_tets( (Vec3f)aabb_low, (Vec3f)aabb_high, nearby_tets );
      
      for ( unsigned int t = 0; t < nearby_tets.size(); ++t )
      {

         for ( unsigned int v = 0; v < 4; ++v )
         {
            unsigned int tet_vertex = mesh->tets[ nearby_tets[t] ][v];
            
            is_on_boundary_tet[tet_vertex] = 1;

            if ( liquid_phi[tet_vertex] < 1e+30 ) { continue; }          
            unsigned int triangle;
            double dist = surface_tracker->distance_to_surface( Vec3d(mesh->vertices[tet_vertex]), triangle );
            
            if ( dist > 1e+10 ) { continue; }
            
            liquid_phi[tet_vertex] = (float)dist;
            closest_triangle[tet_vertex] = triangle;
            vertex_queue.push( tet_vertex );
            
         }
      }         
   }
      
   double narrow_band_end = get_time_in_seconds();

   std::cout << "num_distance_to_surface_calls: " << num_distance_to_surface_calls << std::endl;
   std::cout << "avg number of broadphase queries per call: " << num_broadphase_queries / (double) num_distance_to_surface_calls << std::endl;
   std::cout << "avg number of point-triangle distance tests per call: " << num_point_triangle_tests / (double) num_distance_to_surface_calls << std::endl;
   
   // now propagate 
   
   while( !vertex_queue.empty() )
   {
      // grab a point off the queue
      unsigned int curr_vertex = vertex_queue.front();
      vertex_queue.pop();
      
      const unsigned int curr_nearest_triangle = closest_triangle[curr_vertex];
      
      // for each of the point's neighbours, check if this point's nearest triangle is closer to what they already have
      
      const std::vector<unsigned int>& incident_edges = mesh->vert_to_edge_map[curr_vertex];
      
      for ( unsigned int i = 0; i < incident_edges.size(); ++i )
      {
         if ( mesh->edges[ incident_edges[i] ][0] == mesh->edges[ incident_edges[i] ][1] ) { continue; }
         
         unsigned int neighbour_vertex = mesh->edges[ incident_edges[i] ][0];
         if ( neighbour_vertex == curr_vertex ) { neighbour_vertex = mesh->edges[ incident_edges[i] ][1]; }
         assert( neighbour_vertex != curr_vertex );
      
         // check if the current vertex's nearest triangle is closer than what the neighbour has
         
         float new_dist = get_distance_to_surface_triangle( curr_nearest_triangle, mesh->vertices[neighbour_vertex] );
         
         if ( new_dist < liquid_phi[neighbour_vertex] )
         {
            liquid_phi[neighbour_vertex] = new_dist;
            closest_triangle[neighbour_vertex] = curr_nearest_triangle;
            vertex_queue.push( neighbour_vertex );
         }
      }
   }   
   
   double extrapolate_end = get_time_in_seconds();
   
   std::cout << "redistancing: time to compute narrow band phi: " << narrow_band_end - redistance_start << std::endl;
   std::cout << "redistancing: time to extrapolate phi: " << extrapolate_end - narrow_band_end << std::endl;
 
   //For points that are near the boundary (i.e. belong to a boundary tet)
   //we'll use raycasting on each one.
   region_IDs.assign(liquid_phi.size(), -1);
   int bvert_count = 0;
   for ( unsigned int i = 0; i < liquid_phi.size(); ++i )
   {
      assert( closest_triangle[i] >= 0 );
      assert( liquid_phi[i] < 1e30 );
      
      //use raycasting and triangle labels to determine containing region
      if(is_on_boundary_tet[i]) {
         bvert_count++;
         const Vec3f& x = mesh->vertices[i];   
         region_IDs[i] = surface_tracker->get_region_containing_point( Vec3d(x) );
      }

   }
   std::cout << "Boundary verts: " << bvert_count << " out of " << liquid_phi.size() << std::endl;

   //The remaining points away from the interface can be set by flood filling
   //from just a single check.
   for ( unsigned int i = 0; i < liquid_phi.size(); ++i )
   {
      
      //skip ones that are already processed
      if(region_IDs[i] != -1) continue;
     
      //use raycasting and triangle labels to determine containing region
      const Vec3f& x = mesh->vertices[i];
      int group_ID = surface_tracker->get_region_containing_point( Vec3d(x) );
      
      assert(group_ID != -1);

      std::queue<unsigned int> points_to_process;
      points_to_process.push(i);
      std::cout << "Processing interior region " << group_ID << std::endl;
      //flood fill neighboring non-boundary points, setting them to this region ID.
      while(points_to_process.size() != 0) {
         unsigned int cur_point = points_to_process.front();
         //std::cout << "\nS: " << cur_point <<  " to " << group_ID << std::endl;
         points_to_process.pop();
         
         if(region_IDs[cur_point] != -1) continue;

         region_IDs[cur_point] = group_ID;
         for(unsigned int j = 0; j < mesh->vert_to_edge_map[cur_point].size(); ++j) {
            int edge_ID = mesh->vert_to_edge_map[cur_point][j];
            Vec2st edge_data = mesh->edges[edge_ID];
            int other_vert = (edge_data[0] == cur_point ? edge_data[1] : edge_data[0]);
            
            if(region_IDs[other_vert] == -1) {
               points_to_process.push(other_vert);
            }
         }
      }

   }

   int baseline = 0;
   for(unsigned int i = 0; i < region_IDs.size(); ++i) {
      if(region_IDs[i] == 1)
         baseline++;
   }
   std::cout << "Preprocessed " << baseline << " out of " << baseline << std::endl;
 
   /*
   //Verify with explicit raycasting for all points.
   for ( unsigned int i = 0; i < liquid_phi.size(); ++i )
   {
      const Vec3f& x = mesh->vertices[i];
      int real_ID = surface_tracker->get_region_containing_point( Vec3d(x) );
      if(real_ID != region_IDs[i])
         std::cout << "Real vs. Fake mismatch.\n";
   }
   */

   double set_sign_end = get_time_in_seconds();
   
   std::cout << "redistancing: time to set region labels: " << set_sign_end - extrapolate_end << std::endl;
   

   assert(liquid_phi.size() == region_IDs.size());
   
}

// ---------------------------------------------------------
///
/// 
///
// ---------------------------------------------------------

void DualFluidSim3D::extrapolate_liquid_phi_into_solid( )
{

   Vec3f box_centre = 0.5f * (solid_low + solid_high);
   Vec3f box_extents = solid_high - solid_low;
   
   static const bool USE_INTERPOLATED_PHI = false;
   const unsigned int NUM_PASSES = USE_INTERPOLATED_PHI ? 5 : 1;

   for ( unsigned int pass = 0; pass < NUM_PASSES; ++pass )
   {
      std::vector<float> new_liquid_phi = liquid_phi;
      
      for( unsigned int i = 0; i < mesh->vertices.size(); ++i ) 
      {
         const Vec3f& phi_location = mesh->vertices[i];
         float solid_phi = solid_box_phi( box_centre, box_extents, phi_location );
         
         if ( solid_phi < 0.0 ) { continue; }
                           
         // voronoi cell is solid
         float sample_offset = solid_phi + 0.25f * characteristic_distance;
         Vec3f non_solid_sample = phi_location - sample_offset * solid_box_gradient( box_centre, box_extents, phi_location );
         
         if ( USE_INTERPOLATED_PHI )
         {
            int tet_index = mesh->get_containing_tet( non_solid_sample );
            if ( tet_index >= 0 )
            {
               Vec4f bary = tet_barycentric_weights( non_solid_sample, 
                                                     mesh->vertices[ mesh->tets[tet_index][0] ],
                                                     mesh->vertices[ mesh->tets[tet_index][1] ],
                                                     mesh->vertices[ mesh->tets[tet_index][2] ],
                                                     mesh->vertices[ mesh->tets[tet_index][3] ] );
               
               new_liquid_phi[i] = bary[0] * liquid_phi[mesh->tets[tet_index][0]] 
                                 + bary[1] * liquid_phi[mesh->tets[tet_index][1]] 
                                 + bary[2] * liquid_phi[mesh->tets[tet_index][2]] 
                                 + bary[3] * liquid_phi[mesh->tets[tet_index][3]] 
                                 - 0.2f * characteristic_distance;
               
            }
         }
         else
         {            
            // fire a ray, see if it goes through a "solid" triangle
            std::vector<double> hit_ss;
            std::vector<unsigned int> hit_triangles;
            surface_tracker->get_triangle_intersections( Vec3d(phi_location), Vec3d(non_solid_sample), hit_ss, hit_triangles );
            
            for ( unsigned int t = 0; t < hit_triangles.size(); ++t )
            {
               const Vec3ui& tri = surface_tracker->m_mesh.m_tris[ hit_triangles[t] ];
               if (   surface_tracker->m_masses[tri[0]] > 1.5 
                   && surface_tracker->m_masses[tri[1]] > 1.5 
                   && surface_tracker->m_masses[tri[2]] > 1.5 )
               {
                  new_liquid_phi[i] = -solid_phi; //-0.5f * characteristic_distance;
               }
            }   
         }         
      }
      
      liquid_phi = new_liquid_phi;
   }
   
}

// ---------------------------------------------------------
///
/// Project velocity field into divergence-free vector subspace.
///
// ---------------------------------------------------------

void DualFluidSim3D::solve_pressure()
{
   
   //Single-phase, free surface version
   //pressures = pressure_solve_voronoi( *mesh, *surface_tracker, surface_tension_coefficient, tet_edge_velocities, solid_weights, liquid_phi, wall_velocities);   
   
   //Multi-phase version (also supports free surfaces, when density of a region is 0.)
   pressures = pressure_solve_multi( *mesh, *surface_tracker, surface_tension_coefficient, tet_edge_velocities, solid_weights, liquid_phi, region_IDs, densities);   

}


// ---------------------------------------------------------


void DualFluidSim3D::test_interpolation( float dt ) 
{
   
   float accum_t = 0.0f;   
   
//   for(unsigned int i = 0; i < mesh.edges.size(); ++i) 
//   {
//      float dv = dot( Vec3f(0,1,0), mesh.tet_edge_vectors[i] );
//      tet_edge_velocities[i] = v;
//   }
   
   static const Vec3f gravity = Vec3f( 0.0f, -1.0f, 0.0f );

   for(unsigned int i = 0; i < mesh->edges.size(); ++i) 
   {
      float dv = dot( dt * gravity, mesh->tet_edge_vectors[i] );
      tet_edge_velocities[i] += dv;
   }

   
   while ( accum_t < dt - 1e-5 )
   {
      float sub_dt = dt - accum_t;
      sub_dt = min( sub_dt, get_cfl_limit() );
      
      std::cout << "---------------------- Voronoi Fluid Sim: ";
      std::cout << "sub dt: " << sub_dt << " ----------------------" << std::endl;

      for ( unsigned int i = 0; i < liquid_phi.size(); ++i )
      {
         liquid_phi[i] = -1;
      }
      
      // Transfer and extrapolate velocity
      reconstruct_and_extrapolate_velocities();
      
      // Advect surface
      std::cout << "---------------------- Voronoi Fluid Sim: Advancing fluid surface ----------------------" << std::endl;
      advance_surface( sub_dt );      
            
      // Remesh
      if ( should_remesh )
      {
         std::cout << "---------------------- Voronoi Fluid Sim: Remeshing and semi-lagrangian advection ----------------------" << std::endl;
         remesh_and_advect_semilagrangian( sub_dt );
      }
      else
      {
         std::cout << "---------------------- Voronoi Fluid Sim: Semi-lagrangian advection ----------------------" << std::endl;
         semi_lagrangian_advection( sub_dt );
      }
      
      accum_t += sub_dt;      
   }
   
   reconstruct_and_extrapolate_velocities();
   
}


// ---------------------------------------------------------


float DualFluidSim3D::max_velocity( )
{
   float max_vel = -1.0f;
   for ( unsigned int i = 0; i < tet_edge_velocities.size(); ++i )
   {
      max_vel = max( max_vel, (float)fabs( tet_edge_velocities[i] ) );
   }
   return max_vel;
}


extern unsigned int g_num_pit_tests;
extern unsigned int g_num_pit_hits;


// ---------------------------------------------------------
///
/// Advance the fluid sim by one time step
///
// ---------------------------------------------------------

void DualFluidSim3D::advance( float dt, unsigned int num_surface_substeps ) 
{
   assert( num_surface_substeps > 0 );
      
   float accum_t = 0.0f;
      
   double frame_start_time = get_time_in_seconds();
   
   double reconstruct_and_extrapolate_time = 0.0;
   double surface_tracking_time = 0.0;
   double add_force_time = 0.0;
   double redistancing_time  = 0.0;
   double pressure_solve_time = 0.0;
   
   static bool first_call = true;
   if ( first_call )
   {
      //If the velocities were initialized to some non-zero value, we need to do a pressure solve on
      //them because they may not be divergence-free, and we need them to be so the surface will move
      //in a logical way.

      // Transfer and extrapolate velocity
      std::cout << "---------------------- Voronoi Fluid Sim: First time step: Computing liquid phi ----------------------" << std::endl;
      compute_liquid_phi();
      
      if ( !allow_solid_overlap )
      {
         std::cout << "---------------------- Voronoi Fluid Sim: First time step: Extrapolating liquid phi ----------------------" << std::endl;
         extrapolate_liquid_phi_into_solid();
      }

      std::cout << "---------------------- Voronoi Fluid Sim: First time step: Pressure solve ----------------------" << std::endl;
      solve_pressure();


      std::cout << "---------------------- Voronoi Fluid Sim: First time step: Reconstructing and extrapolating velocities ----------------------" << std::endl;
      reconstruct_and_extrapolate_velocities();     
      
   }
   
   first_call = false;
   
   g_num_pit_tests = 0;
   g_num_pit_hits = 0;
      
   while ( accum_t < dt - 1e-5 )
   {
      float sub_dt = dt - accum_t;
      sub_dt = min( sub_dt, get_cfl_limit() );
   
      std::cout << std::endl << std::endl << std::endl << std::endl << std::endl;
      std::cout << "---------------------- Beginning Voronoi Fluid Sim iteration------------------ \n";
      
      std::cout << "Taking substep of length: " << sub_dt << " (" << 100 * (sub_dt / dt) << "% of a frame.)" << std::endl;
      std::cout << "Currently at sub-frame time: " << accum_t << " which is " << 100 * (accum_t / dt) << "% through the frame." << std::endl;

      std::cout << std::endl << std::endl << std::endl;

      // Advect surface
      std::cout << "---------------------- Voronoi Fluid Sim: Advancing fluid surface ----------------------" << std::endl;
      double start_time = get_time_in_seconds();
      
      std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% max velocity: " << max_velocity() << std::endl;
      

      double surface_dt = sub_dt / (double) num_surface_substeps;
      for ( unsigned int substep = 0; substep < num_surface_substeps; ++substep )
      {
         advance_surface( (float)surface_dt );      
      }
         
      double end_time = get_time_in_seconds();
      surface_tracking_time += ( end_time - start_time );
      // Remesh

      if ( should_remesh )
      {
         std::cout << "---------------------- Voronoi Fluid Sim: Remeshing and semi-lagrangian advection ----------------------" << std::endl;
         remesh_and_advect_semilagrangian( sub_dt );
      }
      else
      {
         std::cout << "---------------------- Voronoi Fluid Sim: Semi-lagrangian advection ----------------------" << std::endl;
         semi_lagrangian_advection( sub_dt );
      }
      
      std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% max velocity: " << max_velocity() << std::endl;

      std::cout << "---------------------- Voronoi Fluid Sim: Computing liquid phi ----------------------" << std::endl;
      start_time = get_time_in_seconds();
      compute_liquid_phi();
      end_time = get_time_in_seconds();
      redistancing_time += ( end_time - start_time );

      std::cout << "---------------------- Voronoi Fluid Sim: Adding forces ----------------------" << std::endl;
      start_time = get_time_in_seconds();
      add_forces( sub_dt );
      end_time = get_time_in_seconds();
      add_force_time += ( end_time - start_time );
      std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% max velocity: " << max_velocity() << std::endl;

      if ( !allow_solid_overlap )
      {
         std::cout << "---------------------- Voronoi Fluid Sim: Extrapolating liquid phi ----------------------" << std::endl;
         start_time = get_time_in_seconds();
         extrapolate_liquid_phi_into_solid();
         end_time = get_time_in_seconds();
         redistancing_time += ( end_time - start_time );
      }
            
      std::cout << "---------------------- Voronoi Fluid Sim: Pressure solve ----------------------" << std::endl;
      start_time = get_time_in_seconds();
      solve_pressure();
      end_time = get_time_in_seconds();
      pressure_solve_time += ( end_time - start_time );
      
      std::cout << " %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% max velocity: " << max_velocity() << std::endl;
      
      // Transfer and extrapolate velocity
      std::cout << "---------------------- Voronoi Fluid Sim: Reconstructing and extrapolating velocities ----------------------" << std::endl;
      start_time = get_time_in_seconds();
      reconstruct_and_extrapolate_velocities();
      end_time = get_time_in_seconds();
      reconstruct_and_extrapolate_time += ( end_time - start_time );
      
      accum_t += sub_dt;
      
   }
   
   double frame_end_time = get_time_in_seconds();
   double frame_time = frame_end_time - frame_start_time;
   std::cout << "---------------------- Voronoi Fluid Sim: Total frame time: " <<  frame_time << "  ----------------------" << std::endl;

   
   // track average execution time vs number of tets
   static unsigned int num_calls = 0;
   static unsigned int total_num_tets = 0;
   static double total_reconstruct_and_extrapolate_time = 0.0;
   static double total_surface_tracking_time = 0.0;
   static double total_add_force_time = 0.0;
   static double total_redistancing_time = 0.0;
   static double total_pressure_solve_time = 0.0;
   static double total_simulation_time = 0.0;
   
   ++num_calls;
   total_num_tets += mesh->tets.size();
   total_reconstruct_and_extrapolate_time += reconstruct_and_extrapolate_time;
   total_surface_tracking_time += surface_tracking_time;
   total_add_force_time += add_force_time;
   total_redistancing_time += redistancing_time;
   total_pressure_solve_time += pressure_solve_time;
   total_simulation_time += frame_time;
   
   unsigned int average_num_tets = total_num_tets / num_calls;
   
   std::cout << "====================== Voronoi Fluid Sim Timings ======================" << std::endl;
   std::cout << "   Average number of tetrahedra: " <<  average_num_tets << std::endl;
   std::cout << "   Average sim time: " <<  total_simulation_time / (double)num_calls << std::endl;
   std::cout << "   Average reconstruct_and_extrapolate_time: " <<  total_reconstruct_and_extrapolate_time / (double)num_calls << std::endl;
   std::cout << "   Average surface_tracking_time: " <<  total_surface_tracking_time / (double)num_calls << std::endl;
   std::cout << "   Average remeshing time: " << total_remesh_time / (double)num_calls << std::endl;
   std::cout << "   Average semi-Lagrangian advection time: " << total_semilagrangian_time / (double)num_calls << std::endl;
   std::cout << "   Average add_force_time: " <<  total_add_force_time / (double)num_calls << std::endl;
   std::cout << "   Average redistancing_time: " <<  total_redistancing_time / (double)num_calls << std::endl;
   std::cout << "   Average pressure_solve_time: " <<  total_pressure_solve_time / (double)num_calls << std::endl;
   std::cout << "=======================================================================" << std::endl;

   std::cout << "---------------------- Voronoi Fluid Sim: Number point-in-tet hits/tests: " <<  g_num_pit_hits << " / " << g_num_pit_tests << " ----------------------" << std::endl;   
   
   // Just for visualization.  Take this out to optimize.
   reconstruct_and_extrapolate_velocities();
   
}

