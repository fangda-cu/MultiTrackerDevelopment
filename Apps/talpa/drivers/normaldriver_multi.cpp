// ---------------------------------------------------------
//
//  normaldriver_multi.cpp
//  Christopher Batty 2013
//
//  Mesh driver for motion in the normal direction using vertex 
//  normals (area-weighted average of incident triangle normals).
//  Extended to multiphase case by Christopher, originally written
//  by Tyson Brochu 2008.
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include "normaldriver_multi.h"
#include "faceoff.h"
#include "../geometryinit.h"

#include <meshsmoother.h>
#include <nondestructivetrimesh.h>
#include <runstats.h>
#include <surftrack.h>

using namespace ElTopo;

// ---------------------------------------------------------
// Global externs
// ---------------------------------------------------------

// ---------------------------------------------------------
// Local constants, typedefs, macros
// ---------------------------------------------------------

// ---------------------------------------------------------
// Static function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Assign the normal at each vertex as the vertex's velocity
///
// ---------------------------------------------------------


void NormalDriverMulti::set_predicted_vertex_positions( const SurfTrack& surf, std::vector<Vec3d>& predicted_positions, double current_t, double& adaptive_dt )
{
    std::vector<Vec3d> displacements( surf.get_num_vertices(), Vec3d(0,0,0) );
    std::vector<Vec3d> velocities( surf.get_num_vertices(), Vec3d(0,0,0) );
    
    for ( unsigned int i = 0; i < surf.get_num_vertices(); ++i )
    {
        if ( surf.m_mesh.m_vertex_to_triangle_map[i].empty() ) 
        { 
            displacements[i] = Vec3d(0,0,0);
            continue;
        }

        Vec3d normal(0,0,0);
        double sum_areas = 0.0;
        double speed;
        bool any_moving = false;
        for ( unsigned int j = 0; j < surf.m_mesh.m_vertex_to_triangle_map[i].size(); ++j )
        {
           size_t tri = surf.m_mesh.m_vertex_to_triangle_map[i][j];
           Vec2i label = surf.m_mesh.get_triangle_label(tri);
           double area = surf.get_triangle_area( tri );
           if(label[0] == expanding_surface || label[1] == expanding_surface) {
              Vec3d tri_normal = surf.get_triangle_normal(tri);   
              speed = speed_matrix[label[0]][label[1]]; //let's assume all speeds are the same.
              normal += tri_normal * area * speed;
              sum_areas += area;
              any_moving = true;
           }
        }

        normal /= sum_areas;
        if(any_moving) {
           //normal /= mag(normal);

           velocities[i] = normal;
           displacements[i] = adaptive_dt * velocities[i];
        }
        else {
           velocities[i] = Vec3d(0,0,0);
           displacements[i] = adaptive_dt * velocities[i];
        }
    }
    
    double capped_dt = MeshSmoother::compute_max_timestep_quadratic_solve( surf.m_mesh.get_triangles(), surf.get_positions(), displacements, false );
    
    adaptive_dt = min( adaptive_dt, capped_dt );
    
    for ( unsigned int i = 0; i < surf.get_num_vertices(); ++i )
    {
        predicted_positions[i] = surf.get_position(i) + adaptive_dt * velocities[i];
    }
    
}


// ---------------------------------------------------------
///
/// Compute and output L1 and L_inf error measures
///
// ---------------------------------------------------------

namespace ElTopo{
extern RunStats g_stats; 
}

void NormalDriverMulti::compute_error( const SurfTrack& surf, double current_t )
{
  
}




