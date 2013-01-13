
#ifndef DUALFLUIDSIM3D_H
#define DUALFLUIDSIM3D_H

#include <ccd_wrapper.h>
#include <surftrack.h>
#include "tetmesh.h"
#include "collisionqueries.h"

class VelocityFunctor3D;

// ---------------------------------------------------------
///
/// Voronoi cell fluid simulation.
///
// ---------------------------------------------------------


class DualFluidSim3D
{
   
public:
   
   DualFluidSim3D( const std::vector<ElTopo::Vec3d>& surface_vertices, 
                   const std::vector<ElTopo::Vec3ui>& surface_triangles, 
                   const std::vector<double>& surface_vertex_masses,
                   const ElTopo::SurfTrackInitializationParameters& initial_parameters );
   
   ~DualFluidSim3D();
   
private:
   DualFluidSim3D();
   DualFluidSim3D( const DualFluidSim3D& );
   DualFluidSim3D& operator=( const DualFluidSim3D& );
   
public:
   
   //
   // init and re-init functions
   //
   
   void initialize();
   
   void compute_solids();
   
   
   //
   // helpers for the main fluid sim steps
   //

   float get_cfl_limit();

   float max_velocity();
   
   inline void trace_rk2( const ElTopo::Vec3f& start, ElTopo::Vec3f& end, float dt, const VelocityFunctor3D& get_velocity );

   //void surface_laplacian_smoothing( double coefficient );
   
   void correct_volume( );
   
   // barycentric interpolation
   ElTopo::Vec3f get_velocity_from_tet_vertices( const ElTopo::Vec3f& point );

   // Whitney-style interpolation
   ElTopo::Vec3f get_velocity_from_tet_edges( const ElTopo::Vec3f& point );
   
   // sharper barycentric interpolation
   ElTopo::Vec3f get_sharper_barycentric_velocity( const ElTopo::Vec3f& point );

   //generalized barycentric interpolation over Voronoi region
   ElTopo::Vec3f get_generalized_barycentric_velocity( const ElTopo::Vec3f& point );

   // helpers for barycentric interpolation
   void reset_vertex_velocities();
   void tet_edge_to_vertex_velocities();
   void extrapolate_tet_vertex_velocities();
   void clamp_solid_tet_vertex_velocities();
   
   //helpers for sharper barycentric interpolation
   void tet_edge_to_circumcentre_velocities();
   void construct_full_voronoi_face_velocities();
   void extrapolate_all_velocities();
   void clamp_solid_all_velocities();

   // helpers for Whitney interpolation
   void extrapolate_tet_edge_velocities();
   
   void add_thermal_buoyancy( float dt );
   void add_local_force( float dt );

   float get_distance_to_surface_triangle( unsigned int triangle_index, const ElTopo::Vec3f& point );

   
   //
   // the main fluid sim steps
   //
      
   void reconstruct_and_extrapolate_velocities();
   
   void advance_surface( float dt );
   
   void remesh_and_advect_semilagrangian( float dt );
   
   void semi_lagrangian_advection( float dt );
   
   void add_forces( float dt );
   
   void compute_liquid_phi();
   
   void extrapolate_liquid_phi_into_solid();
   
   void solve_pressure();
   
   void advance( float dt, unsigned int num_surface_substeps = 1 );

   void test_interpolation( float dt );
  
   bool is_liquid(int vert_index) {
     //return liquid_phi[vert_index] < 0;
     return densities[region_IDs[vert_index]] != 0;
   }

   //
   // data members
   //
   
   // simulation domain
   ElTopo::Vec3f domain_low;
   ElTopo::Vec3f domain_high;
   float characteristic_distance;
   
   float surface_tension_coefficient;

   ElTopo::Vec3f solid_low;
   ElTopo::Vec3f solid_high;
   
   TetMesh* mesh;
   
   ElTopo::SurfTrack *surface_tracker;
   
   bool free_surface;
   bool should_remesh;
   
   bool volume_correction;
   double initial_volume;
   
   bool allow_solid_overlap;     // allow the surface mesh to penetrate the solid (e.g. to extend a surface into a solid boundary)
   
   std::vector<ElTopo::Vec3f> markers;

   ElTopo::Vec3f gravity;
   
   // on tet vertices / voronoi centres
   std::vector<float> liquid_phi;         

   // multiphase support
   std::vector<int> region_IDs; //identify which fluid, located at vertices / voronoi centres
   std::vector<float> densities; //the density of the particular fluid - 0 implies free surface region. One per region ID.

   // Lagrange multipliers used during pressure projection (for debugging/visualization)
   std::vector<double> pressures;
   
   // tangential velocities on tet edges
   // (dual: normal velocities on voronoi faces)
   std::vector<float> tet_edge_velocities;
   
   // reconstructed velocity vectors on tet vertices
   // (dual: reconstructed velocity vectors on voronoi sites )
   std::vector<ElTopo::Vec3f> tet_vertex_velocities;

   // reconstructed velocity vectors on voronoi vertices
   // (dual: reconstructed velocity vectors on tet circumcentre )
   std::vector<ElTopo::Vec3f> voronoi_vertex_velocities;

   //full vectors on voronoi face (tet edges) 
   std::vector<ElTopo::Vec3f> full_voronoi_face_velocities;

   //validity flags for extrapolation
   std::vector<bool> tet_vertex_velocity_is_valid;
   std::vector<bool> voronoi_vertex_velocity_is_valid;
   std::vector<bool> full_voronoi_face_velocity_is_valid;

   std::vector<bool> edge_is_valid;
   
   std::vector<float> solid_phi;          // on tet circumcentres / voronoi cell vertices
   std::vector<float> solid_weights;      // on tet edges / voronoi faces
   
   double total_remesh_time;
   double total_semilagrangian_time;
   
   enum { BARYCENTRIC, GENERALIZED_BARYCENTRIC, IMPROVED_BARYCENTRIC, WHITNEY, NUM_INTERPOLATION_SCHEMES };
   unsigned int interpolation_scheme;
   
};


// ---------------------------------------------------------
///
/// Abstract velocity interpolation function object.
///
// ---------------------------------------------------------

class VelocityFunctor3D 
{
public:
   virtual ~VelocityFunctor3D() {}
   virtual ElTopo::Vec3f operator()(const ElTopo::Vec3f& point) const = 0;
};


// ---------------------------------------------------------
///
/// Assuming we have 3d velocity vectors at all tet vertices, use barycentric interpolation.
///
// ---------------------------------------------------------

class BarycentricTetVelocityFunctor : public VelocityFunctor3D 
{
   DualFluidSim3D& sim;
public:
   virtual ~BarycentricTetVelocityFunctor() {}
   
   BarycentricTetVelocityFunctor( DualFluidSim3D& sim_ ) : sim(sim_) {}
   
   ElTopo::Vec3f operator()(const ElTopo::Vec3f& pt) const 
   {
      return sim.get_velocity_from_tet_vertices(pt);
   }
};


// ---------------------------------------------------------
///
/// Use Whitney-element-style interpolation
///
// ---------------------------------------------------------

class WhitneyEdgeVelocityFunctor : public VelocityFunctor3D 
{
   DualFluidSim3D& sim;
public:
   virtual ~WhitneyEdgeVelocityFunctor() {}
   
   WhitneyEdgeVelocityFunctor( DualFluidSim3D& sim_ ) : sim(sim_) {}
   
   ElTopo::Vec3f operator()(const ElTopo::Vec3f& pt) const 
   {
      return sim.get_velocity_from_tet_edges(pt);
   }
};

// ---------------------------------------------------------
///
/// Assuming we have 3d velocity vectors at all tet vertices, voronoi vertices, and voronoi face centres,
/// use sharper barycentric interpolation
///
// ---------------------------------------------------------

class SharperBarycentricVelocityFunctor : public VelocityFunctor3D 
{
   DualFluidSim3D& sim;
public:
   virtual ~SharperBarycentricVelocityFunctor() {}
   
   SharperBarycentricVelocityFunctor( DualFluidSim3D& sim_ ) : sim(sim_) {}
   
   ElTopo::Vec3f operator()(const ElTopo::Vec3f& pt) const 
   {
      return sim.get_sharper_barycentric_velocity(pt);
   }
};

// ---------------------------------------------------------
///
/// Assuming we have 3d velocity vectors at all voronoi vertices
/// use generalized barycentric velocity interpolation
///
// ---------------------------------------------------------

class GeneralizedBarycentricVelocityFunctor : public VelocityFunctor3D 
{
   DualFluidSim3D& sim;
public:
   virtual ~GeneralizedBarycentricVelocityFunctor() {}
   
   GeneralizedBarycentricVelocityFunctor( DualFluidSim3D& sim_ ) : sim(sim_) {}
   
   ElTopo::Vec3f operator()(const ElTopo::Vec3f& pt) const 
   {
      return sim.get_generalized_barycentric_velocity(pt);
   }
};

// ---------------------------------------------------------
// Inline function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------

inline void DualFluidSim3D::trace_rk2( const ElTopo::Vec3f& start, ElTopo::Vec3f& end, float dt, const VelocityFunctor3D& get_velocity ) 
{
   
   // advance to midpoint
   ElTopo::Vec3f vel = get_velocity(start);
   ElTopo::Vec3f mid = start + 0.5f * dt * vel;
   
   // get velocity at midpoint
   vel = get_velocity(mid);
   
   // use it to compute the final position
   end = start + dt * vel;
   
}


// ---------------------------------------------------------
///
/// 
///
// ---------------------------------------------------------

inline float DualFluidSim3D::get_distance_to_surface_triangle( unsigned int triangle_index, const ElTopo::Vec3f& point )
{
   double dist;
   //unsigned int dummy_index = surface_tracker->get_num_vertices();
   const ElTopo::Vec3ui& tri = surface_tracker->m_mesh.m_tris[triangle_index];
   
   ElTopo::check_point_triangle_proximity(ElTopo::Vec3d(point), 
      ElTopo::Vec3d(surface_tracker->get_position(tri[0])),
      ElTopo::Vec3d(surface_tracker->get_position(tri[1])),
      ElTopo::Vec3d(surface_tracker->get_position(tri[2])),
      dist);
   
     /* ElTopo::point_triangle_distance( ElTopo::Vec3d(point), dummy_index,
                            ElTopo::Vec3d(surface_tracker->get_position(tri[0])), tri[0],
                            ElTopo::Vec3d(surface_tracker->get_position(tri[1])), tri[1],
                            ElTopo::Vec3d(surface_tracker->get_position(tri[2])), tri[2],
                            dist );*/
   
   return (float)dist;
   
}

#endif
