// ---------------------------------------------------------
//
//  meancurvature.cpp
//  Tyson Brochu 2008
//
//  Mesh driver for motion in the normal direction scaled by mean curvature
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include "meancurvature_multi.h"

#include <array3_utils.h>
#include <fstream>
#include <iomesh.h>
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
// Static and nonmember function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Constructor
///
// ---------------------------------------------------------

MeanCurvatureMultiDriver::MeanCurvatureMultiDriver(double in_curvature_multiplier) :
   curvature_multiplier( in_curvature_multiplier )
{
}

// ---------------------------------------------------------
///
/// Compute mean curvature times normal at a vertex and return the sum of weights used (for computing the time step restriction)
///
// ---------------------------------------------------------

void MeanCurvatureMultiDriver::add_triangle_contribution_of_mean_curvature_normal(size_t triangle_index, const SurfTrack & surf, Vec3d & vert0, Vec3d & vert1, Vec3d & vert2)
{
//    Vec3d mean_curvature_normal( 0, 0, 0 );
//    weight_sum = 0;
//    
//    double edge_length_sum = 0.0;
//    
//    for ( size_t i = 0; i < surf.m_mesh.m_vertex_to_edge_map[vertex_index].size(); ++i )
//    {
//        size_t e = surf.m_mesh.m_vertex_to_edge_map[vertex_index][i];
//        const Vec2st& curr_edge = surf.m_mesh.m_edges[e];
//        Vec3d edge_vector;
//        if ( curr_edge[0] == vertex_index )
//        {
//            edge_vector = surf.get_position( curr_edge[1] ) - surf.get_position( vertex_index );
//        }
//        else
//        {
//            assert( curr_edge[1] == vertex_index );
//            edge_vector = surf.get_position( curr_edge[0] ) - surf.get_position( vertex_index );
//        }
//        
//        edge_length_sum += mag( edge_vector );
//        
//        if ( surf.m_mesh.m_edge_to_triangle_map[e].size() != 2 )
//        {
//            // TODO: properly handle more than 2 incident triangles
//            continue;
//        }
//        
//        size_t tri0 = surf.m_mesh.m_edge_to_triangle_map[e][0];
//        size_t tri1 = surf.m_mesh.m_edge_to_triangle_map[e][1];
//        
//        size_t third_vertex_0 = surf.m_mesh.get_third_vertex( curr_edge[0], curr_edge[1], surf.m_mesh.get_triangle(tri0) );
//        size_t third_vertex_1 = surf.m_mesh.get_third_vertex( curr_edge[0], curr_edge[1], surf.m_mesh.get_triangle(tri1) );
//        
//        Vec3d v00 = surf.get_position( curr_edge[0] ) - surf.get_position( third_vertex_0 );
//        Vec3d v10 = surf.get_position( curr_edge[1] ) - surf.get_position( third_vertex_0 );
//        
//        double cross_0 = mag( cross( v00, v10 ) );
//        if ( cross_0 < 1e-10 )
//        {
//            continue;
//        }
//        double cot_0 = dot(v00, v10) / cross_0;
//        
//        Vec3d v01 = surf.get_position( curr_edge[0] ) - surf.get_position( third_vertex_1 );
//        Vec3d v11 = surf.get_position( curr_edge[1] ) - surf.get_position( third_vertex_1 );
//        
//        double cross_1 = mag( cross( v01, v11 ) );
//        if ( cross_1 < 1e-10 )
//        {
//            continue;
//        }
//        
//        double cot_1 = dot(v01, v11) / cross_1;
//        
//        double weight = cot_0 + cot_1;
//        weight_sum += weight;
//        
//        mean_curvature_normal += weight * edge_vector;
//        
//    }
//    
//    double vertex_area = 0.0;
//    for ( size_t i = 0; i < surf.m_mesh.m_vertex_to_triangle_map[vertex_index].size(); ++i )
//    {
//        vertex_area += mixed_area( vertex_index, surf.m_mesh.m_vertex_to_triangle_map[vertex_index][i], surf );
//    }
//    
//    double coeff = 1.0 / (2.0 * vertex_area);
//    
//    weight_sum *= coeff;
//    
//    out = coeff * mean_curvature_normal;
//

    const Vec3st & t = surf.m_mesh.get_triangle(triangle_index);
    
    Vec3d v01 = surf.get_position(t[1]) - surf.get_position(t[0]);
    Vec3d v20 = surf.get_position(t[0]) - surf.get_position(t[2]);
    Vec3d v12 = surf.get_position(t[2]) - surf.get_position(t[1]);
    
    Vec3d A = cross(v01, -v20);
    double Anorm = mag(A);
    Vec3d mul = A / Anorm;
    Vec3d out2 = curvature_multiplier * cross(v01, mul);
    Vec3d out1 = curvature_multiplier * cross(mul, -v20);
    Vec3d out0 = -(out1 + out2);
    
    vert0 += out0;
    vert1 += out1;
    vert2 += out2;
}


// ---------------------------------------------------------
///
/// Set velocities on each mesh vertex
///
// ---------------------------------------------------------

void MeanCurvatureMultiDriver::set_predicted_vertex_positions(const SurfTrack & surf, std::vector<Vec3d> & predicted_positions, double current_t, double & adaptive_dt)
{
    const NonDestructiveTriMesh & mesh = surf.m_mesh;
    std::vector<Vec3d> v(mesh.nv(), Vec3d(0, 0, 0));
    
    for (size_t i = 0; i < mesh.nt(); i++)
    {
        const Vec3st & t = mesh.get_triangle(i);
        add_triangle_contribution_of_mean_curvature_normal(i, surf, v[t[0]], v[t[1]], v[t[2]]);
    }
    
    double global_max_dt = 1.0;
    for (size_t i = 0; i < mesh.ne(); i++)
    {
        size_t v0 = mesh.m_edges[i][0];
        size_t v1 = mesh.m_edges[i][1];
        Vec3d edge = surf.get_position(v1) - surf.get_position(v0);
        Vec3d relv = v[v1] - v[v0];
        double approaching_velocity = -dot(relv, edge) / dot(edge, edge);
        
        double max_dt = 0;
        if (approaching_velocity <= 0)
            max_dt = std::numeric_limits<double>::infinity();
        else
            max_dt = 1 / approaching_velocity * 0.8;    // allow the edge to be shrinked by 80% in one time step at most
        
        if (max_dt < global_max_dt)
            global_max_dt = max_dt;
    }
    
    adaptive_dt = std::min(adaptive_dt, global_max_dt);
    
    predicted_positions.resize(mesh.nv());
    for (size_t i = 0; i < mesh.nv(); i++)
    {
        predicted_positions[i] = surf.get_position(i) + adaptive_dt * v[i];
        if (surf.get_position(i)[0] == 0) predicted_positions[i][0] = 0;
        if (surf.get_position(i)[1] == 0) predicted_positions[i][1] = 0;
        if (surf.get_position(i)[2] == 0) predicted_positions[i][2] = 0;
        if (surf.get_position(i)[0] == 1) predicted_positions[i][0] = 1;
        if (surf.get_position(i)[1] == 1) predicted_positions[i][1] = 1;
        if (surf.get_position(i)[2] == 1) predicted_positions[i][2] = 1;
    }
}
