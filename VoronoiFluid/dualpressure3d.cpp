
#include "dualpressure3d.h"

#include "cgsolver.h"
#include "sparsematrix.h"
#include "tetmesh.h"
#include "meancurvature.h"
#include "dynamicsurface.h"
#include <fstream>
#include "util.h"

using namespace ElTopo;

void get_bary_coords(const Vec3d point, const Vec3d& x1, const Vec3d& x2, const Vec3d& x3, Vec3d& coords) {
   //Swiped from Robert's collision queries code.

   // do it the QR way for added robustness
   Vec3d x13=x1-x3;
   double r00=mag(x13)+1e-30;
   x13/=r00;
   Vec3d x23=x2-x3;
   double r01=dot(x23,x13);
   x23-=r01*x13;
   double r11=mag(x23)+1e-30;
   x23/=r11;
   Vec3d x03=point-x3;
   double s2=dot(x23,x03)/r11;
   double s1=(dot(x13,x03)-r01*s2)/r00;
   double s3=1-s1-s2;

   coords = Vec3d(s1,s2,s3);
}

float get_surface_curvature(TetMesh& mesh, DynamicSurface& surface, std::vector<double>& vertex_curvatures, int i, int neighbour_index) {
   
   //cast a ray from the interior pressure sample to the exterior one, to determine the crossing point
   Vec3d ray_origin = (Vec3d)mesh.vertices[i];
   Vec3d ray_end = (Vec3d)mesh.vertices[neighbour_index];
   std::vector<double> hit_ss;
   std::vector<size_t> hit_triangles; 
   surface.get_triangle_intersections( ray_origin, ray_end, hit_ss, hit_triangles );

   //estimate the mean curvature at the crossing point, and use it to dictate the surface tension force.
   double mean_curvature = 0;

   if(hit_triangles.size() > 0) {
      Vec3st tri = surface.m_mesh.m_tris[hit_triangles[0]];
      Vec3d cross_point = lerp(ray_origin, ray_end, hit_ss[0]);
      Vec3d v0 = surface.get_position(tri[0]);
      Vec3d v1 = surface.get_position(tri[1]);
      Vec3d v2 = surface.get_position(tri[2]);
      Vec3d coords;
      get_bary_coords(cross_point, v0,v1,v2, coords);

      mean_curvature = coords[0] * vertex_curvatures[tri[0]] + coords[1] * vertex_curvatures[tri[1]] + coords[2] * vertex_curvatures[tri[2]];
      
   }
   
   return (float)mean_curvature;
}

std::vector<double> pressure_solve_multi( TetMesh& mesh, 
   DynamicSurface& surface,
   float surface_tension_coeff,
   std::vector<float>& face_velocities,             // tangential to tet edges / normal to voronoi faces
   const std::vector<float>& solid_weights,         // on the voronoi faces (face fractions for solid boundary conditions)
   const std::vector<float>& liquid_phi,            // on the voronoi sites (distance field, treated as unsigned)
   const std::vector<int>& regions,                 // on the voronoi sites (region ID's, to distinguish where interfaces are)
   const std::vector<float>& densities)             // one per region ID  (zero implies a free surface region)
{

   // Verify some dimensions for good measure.
   assert(regions.size() == liquid_phi.size());
   assert(face_velocities.size() == solid_weights.size());

   std::vector<double> vertex_curvatures; //Per vertex surface mesh curvature estimates
   std::vector<double> face_curvatures;   //Saved curvatures at each interface crossing point, to avoid rayshooting twice.

   // Pre-compute vertex mean curvatures on the surface mesh
   if(surface_tension_coeff > 0) {
   
      vertex_curvatures.resize(surface.get_num_vertices());
      for(unsigned int i = 0; i < surface.get_num_vertices(); ++i) {
         Vec3d surf_normal = surface.get_vertex_normal_angleweighted(i);

         Vec3d curvatureNormal;
         MeanCurvatureDriver::vertex_mean_curvature_normal(i, surface, curvatureNormal);

         vertex_curvatures[i] = -(float)mag(curvatureNormal) * (dot(surf_normal,curvatureNormal) > 0?1:-1);
      }

      // Smooth the curvatures a bit
      for(int pass = 0; pass < 0; ++pass) {

         std::vector<double> temp_curvatures(surface.get_num_vertices());

         for(unsigned int i = 0; i < surface.get_num_vertices(); ++i) {
            double sum = 0;
            for(unsigned int edge = 0; edge < surface.m_mesh.m_vertex_to_edge_map[i].size(); ++edge) {
               int edgeID = surface.m_mesh.m_vertex_to_edge_map[i][edge];
               int otherVert = surface.m_mesh.m_edges[edgeID][0] == i? surface.m_mesh.m_edges[edgeID][1] : surface.m_mesh.m_edges[edgeID][0];
               sum += vertex_curvatures[otherVert];
            }
            temp_curvatures[i] = 0.5f * vertex_curvatures[i] + 0.5f * sum / (double)surface.m_mesh.m_vertex_to_edge_map[i].size();
         }
         vertex_curvatures = temp_curvatures;
      }

      // Create space to store the computed curvatures at interface crossings
      face_curvatures.resize(mesh.edges.size());
   }


   // Free surface boundary condition choice...
   // Clamping the surface position to 0.999 effectively puts zero at the neighbour pressure (1st order)
   // Clamping it slightly away from zero (0.001f) lets the surface lie at the true position (2nd order), 
   // or nudges it slightly to avoid division by zero.
   static const float theta_clamp = 0.01f;            

   // Size the matrix data
   SparseMatrixd matrix( mesh.vertices.size() );
   std::vector<double> rhs( mesh.vertices.size(), 0 );
   std::vector<double> pressure( mesh.vertices.size(), 0 );

   // Consider each Voronoi site, set up its row in the linear system

   for(unsigned int i = 0; i < mesh.vertices.size(); ++i) 
   {
      
      // Skip free surface regions
      int region_ID = regions[i];
      if(densities[region_ID] == 0) continue; //density == 0 indicates a free surface region, so skip it

      // Check if the Voronoi cell is entirely in the solid, and if so, skip processing it.
      const std::vector<unsigned int>& incident_faces = mesh.vert_to_edge_map[i];
      bool empty = true;
      for(unsigned int j = 0; j < incident_faces.size(); ++j) 
      {
         int face_index = incident_faces[j];
         if(solid_weights[face_index] > 1e-8)
            empty = false;
      }
      if ( empty ) 
         continue;     

      float diagonal_sum = 0;
      float self_phi = liquid_phi[i];
      for(unsigned int j = 0; j < incident_faces.size(); ++j) 
      {
         // Consider the Voronoi cell on the other side of the current face
         int face_index = incident_faces[j];
         int neighbour_index = mesh.edges[face_index][0]==i? mesh.edges[face_index][1] : mesh.edges[face_index][0];

         // There is no neighboring cell. (We're near the edge of the mesh, so skip it to be safe - probably shouldn't happen.)
         if ( neighbour_index == -1 ) 
            continue;

         // Get the distance between the Voronoi sites (needed for finite differences gradient estimates)
         float dist = mesh.edge_lengths[face_index];         
         dist = max( dist, 1e-7f );

         // Get the data on the neighbour region: distance to the surface and region ID
         float neighbour_phi = liquid_phi[neighbour_index];
         int nbr_region_ID = regions[neighbour_index];
         
         if( densities[nbr_region_ID] != 0)
         {
            //Not a free surface, so need entries for both cells

            float face_density = 1;
            if(region_ID == nbr_region_ID) { 
               //We are in the same region, so use the same density.
               face_density = densities[region_ID];
            }
            else {
               // Different region - multiphase case of ghost fluid method.
               // Compute the interpolated interface position from unsigned distances and
               // use it to interpolate density (per e.g. [Boyd/Bridson 2011] or Losasso [2006])
               float interface_theta = abs(self_phi) / (abs(self_phi) + abs(neighbour_phi)); 
               interface_theta = clamp(interface_theta, 0.01f, 0.99f);
               face_density = interface_theta  * densities[region_ID] + (1-interface_theta)*densities[nbr_region_ID];

               if(surface_tension_coeff > 0) {
                  // Estimate the mean curvature at the crossing point, and use it to incorporate surface tension on the RHS.
                  
                  //Always orient the ray from lower region ID to higher one for consistency
                  int lower_index = region_ID < nbr_region_ID? i : neighbour_index;
                  int higher_index = lower_index == i ? neighbour_index : i;
                  face_curvatures[face_index] = get_surface_curvature(mesh, surface, vertex_curvatures, lower_index, higher_index);
                  
                  float sign_flip = (lower_index == i? 1.0f : -1.0f);
                  rhs[i] += sign_flip * solid_weights[face_index] * mesh.voronoi_face_areas[face_index] * surface_tension_coeff * face_curvatures[face_index] / dist / face_density;
               }
            }

            // Add diagonal and off-diagonal entries.
            matrix.add_to_element(i, neighbour_index, -solid_weights[face_index] * mesh.voronoi_face_areas[face_index] / dist / face_density);
            diagonal_sum += abs(-solid_weights[face_index] * mesh.voronoi_face_areas[face_index] / dist / face_density);      
         }
         else  
         {
            // Free surface case, so use ghost fluid boundary condition
            float face_density = densities[region_ID]; //Use the interior density.
            float theta = max( theta_clamp, std::abs(self_phi) / (std::abs(self_phi) + std::abs(neighbour_phi))); //Determine the interface position
            
            // Increment the diagonal entry accordingly.
            diagonal_sum += solid_weights[face_index] * mesh.voronoi_face_areas[face_index] / dist / theta / face_density;

            // Add surface tension
            if(surface_tension_coeff > 0) {
               // Estimate the mean curvature at the crossing point, and use it to incorporate surface tension on the RHS.
               face_curvatures[face_index] = get_surface_curvature(mesh, surface, vertex_curvatures, i, neighbour_index);
               rhs[i] += solid_weights[face_index] * mesh.voronoi_face_areas[face_index] * surface_tension_coeff * face_curvatures[face_index] / dist / theta / face_density;
            }
         }

         // Now increment the divergence (i.e., the RHS)...

         // Account for face orientation of the divergence stencil by flipping the sign.
         float sign = mesh.edges[face_index][0]==i ? 1.0f : -1.0f;
         rhs[i] -= solid_weights[face_index] * face_velocities[face_index] * sign * mesh.voronoi_face_areas[face_index];

         // Debugging check for good measure.
         assert ( rhs[i] == rhs[i] );

      }

      // Set the diagonal entry with the accumulated value
      // This ensures (weak) diagonal dominance (which was previously not guaranteed due to some roundoff issues.)
      matrix.set_element(i, i, diagonal_sum);

   }

   /*
   // Debugging
   std::ofstream outfile("output.txt");
   matrix.write_matlab(outfile, "matrix");
   write_matlab(outfile, rhs, "rhs");
   outfile.close();
   */

   printf("Solving matrix\n");

   
   CGSolver<double> solver;
   solver.tolerance_factor = 1e-12;
   solver.max_iterations = 1000;

   double residual;
   int iterations;
   solver.solve(matrix, rhs, pressure, residual, iterations);

   printf("Iterations: %d\n", iterations);

   // Apply the gradient of pressure to the velocity field.
   // Consider each face in the mesh, and the two cells on either side.
   for(unsigned int i = 0; i < mesh.edges.size(); ++i) 
   {
      Vec2st verts = mesh.edges[i]; // Already ordered correctly for the next step

      if(solid_weights[i] > 0) 
      {
         // Get distances for each side.
         float phi0 = liquid_phi[verts[0]];
         float phi1 = liquid_phi[verts[1]];

         // Determine region ID's.
         int region0 = regions[verts[0]];
         int region1 = regions[verts[1]];

         // Check for free surface regions.
         bool isFS0 = (densities[region0] == 0);
         bool isFS1 = (densities[region1] == 0);

         // Skip free-surface faces (both sides are free surface regions).
         if(isFS0 && isFS1) 
         {
            face_velocities[i] = 0.0f;
            continue; 
         }

         double p0 = pressure[verts[0]];
         double p1 = pressure[verts[1]];
         float theta = 1;
         float face_density = 1;

         if(isFS0) 
         {
            p0 = 0;
            if(surface_tension_coeff > 0)
               p0 = surface_tension_coeff * face_curvatures[i]; //Use the stored curvature from before.
            theta = max(theta_clamp, std::abs(phi1) / (std::abs(phi1) + std::abs(phi0)));
            face_density = densities[region1];
         }
         else if(isFS1) 
         {
            p1 = 0;
            if(surface_tension_coeff > 0)
               p1 = surface_tension_coeff * face_curvatures[i]; //Use the stored curvature from before.
            theta = max(theta_clamp, std::abs(phi0) / (std::abs(phi0) + std::abs(phi1)));
            face_density = densities[region0];
         }
         else { // Neither side is a free surface.
            theta = 1;

            if(region0 != region1) { // There's an interface between the two regions
               //Multiphase ghost fluid again, compute the interpolated density.
               float interface_theta = abs(phi0) / (abs(phi0) + abs(phi1));
               interface_theta = clamp(interface_theta, 0.01f, 0.99f);
               face_density = interface_theta  * densities[region0] + (1-interface_theta )*densities[region1];

               //Increment one of the pressures to account for surface tension pressure jump.
               if(surface_tension_coeff > 0) {
                  //Make sure we incorporate surface tension in the proper direction as before.
                  int lower_index = region0 < region1? verts[0] : verts[1];
               
                  float sign_flip = (lower_index == verts[0] ? 1.0f : -1.0f);
                  p0 += sign_flip * surface_tension_coeff * face_curvatures[i];
               }
            }
            else { // It's all one fluid so do nothing special.
               face_density = densities[region0];
            }

         }
         
         // Get the distance between the pressure samples for the gradient finite difference.
         float dist = mesh.edge_lengths[i];
         dist = max( dist, 1e-7f );

         // Compute and subtract the discrete gradient.
         face_velocities[i] -= (float)(p1 - p0) / dist /  theta / face_density; 
      }
      else 
      {  
         // Just set this to zero for safety -> Will be overwritten by extrapolation.
         face_velocities[i] = 0.0f;
      }
   }

   return pressure;

}


