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

void MeanCurvatureMultiDriver::initialize(ElTopo::SurfTrack & st)
{
    st.m_solid_vertices_callback = this;
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


int onBBWall(const Vec3d & pos)
{
    int walls = 0;
    if (pos[0] == 0) walls |= (1 << 0);
    if (pos[1] == 0) walls |= (1 << 1);
    if (pos[2] == 0) walls |= (1 << 2);
    if (pos[0] == 1) walls |= (1 << 3);
    if (pos[1] == 1) walls |= (1 << 4);
    if (pos[2] == 1) walls |= (1 << 5);
    
    return walls;
}

Vec3d enforceBBWallConstraint(const Vec3d & input, int constraints)
{
    Vec3d output = input;
    if (constraints & (1 << 0)) output[0] = 0;
    if (constraints & (1 << 1)) output[1] = 0;
    if (constraints & (1 << 2)) output[2] = 0;
    if (constraints & (1 << 3)) output[0] = 1;
    if (constraints & (1 << 4)) output[1] = 1;
    if (constraints & (1 << 5)) output[2] = 1;
    
    return output;
}

bool MeanCurvatureMultiDriver::generate_collapsed_position(ElTopo::SurfTrack & st, size_t v0, size_t v1, ElTopo::Vec3d & pos)
{
    ElTopo::Vec3d x0 = st.get_position(v0);
    ElTopo::Vec3d x1 = st.get_position(v1);
    
    int label0 = onBBWall(x0);
    int label1 = onBBWall(x1);
    
    if (label0 == label1)
    {
        // on the same wall(s), prefer the one with higher max edge valence
        size_t maxedgevalence0 = 0;
        size_t maxedgevalence1 = 0;
        for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v0].size(); i++)
            if (st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v0][i]].size() > maxedgevalence0)
                maxedgevalence0 = st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v0][i]].size();
        for (size_t i = 0; i < st.m_mesh.m_vertex_to_edge_map[v1].size(); i++)
            if (st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v1][i]].size() > maxedgevalence1)
                maxedgevalence1 = st.m_mesh.m_edge_to_triangle_map[st.m_mesh.m_vertex_to_edge_map[v1][i]].size();
        
        if (maxedgevalence0 == maxedgevalence1) // same max edge valence, use their midpoint
            pos = (x0 + x1) / 2;
        else if (maxedgevalence0 < maxedgevalence1)
            pos = x1;
        else
            pos = x0;
        
        return true;
        
    } else if ((label0 & ~label1) == 0)
    {
        // label0 is a proper subset of label1 (since label0 != label1)
        pos = x1;
        
        return true;
        
    } else if ((label1 & ~label0) == 0)
    {
        // label1 is a proper subset of label0
        pos = x0;
        
        return true;
        
    } else
    {
        // label0 and label1 are not subset of each other
        int newlabel = label0 | label1;
        assert(label0 != newlabel); // not subset of each other
        assert(label1 != newlabel);
        
        assert(!((label0 & (1 << 0)) != 0 && (label0 & (1 << 3)) != 0)); // can't have conflicting constraints in label0 and label1 already
        assert(!((label0 & (1 << 1)) != 0 && (label0 & (1 << 4)) != 0));
        assert(!((label0 & (1 << 2)) != 0 && (label0 & (1 << 5)) != 0));
        assert(!((label1 & (1 << 0)) != 0 && (label1 & (1 << 3)) != 0));
        assert(!((label1 & (1 << 1)) != 0 && (label1 & (1 << 4)) != 0));
        assert(!((label1 & (1 << 2)) != 0 && (label1 & (1 << 5)) != 0));
        
        bool conflict = false;
        if ((newlabel & (1 << 0)) != 0 && (newlabel & (1 << 3)) != 0) conflict = true;
        if ((newlabel & (1 << 1)) != 0 && (newlabel & (1 << 4)) != 0) conflict = true;
        if ((newlabel & (1 << 2)) != 0 && (newlabel & (1 << 5)) != 0) conflict = true;
        
        if (conflict)
        {
            // the two vertices are on opposite walls (conflicting constraints). Can't collapse this edge (which shouldn't have become a collapse candidate in the first place)
            return false;
        }
        
        pos = (x0 + x1) / 2;
        pos = enforceBBWallConstraint(pos, newlabel);
        
        return true;
    }
    
    return false;
}

bool MeanCurvatureMultiDriver::generate_split_position(ElTopo::SurfTrack & st, size_t v0, size_t v1, ElTopo::Vec3d & pos)
{
    pos = (st.get_position(v0) + st.get_position(v1)) / 2;
    
    return true;
}

ElTopo::Vec3c MeanCurvatureMultiDriver::generate_collapsed_solid_label(ElTopo::SurfTrack & st, size_t v0, size_t v1, const ElTopo::Vec3c & label0, const ElTopo::Vec3c & label1)
{
    ElTopo::Vec3d x0 = st.get_position(v0);
    ElTopo::Vec3d x1 = st.get_position(v1);
    
    int constraint0 = onBBWall(x0);
    int constraint1 = onBBWall(x1);
    
    assert(((constraint0 & (1 << 0)) || (constraint0 & (1 << 3))) == (bool)label0[0]);
    assert(((constraint0 & (1 << 1)) || (constraint0 & (1 << 4))) == (bool)label0[1]);
    assert(((constraint0 & (1 << 2)) || (constraint0 & (1 << 5))) == (bool)label0[2]);
    assert(((constraint1 & (1 << 0)) || (constraint1 & (1 << 3))) == (bool)label1[0]);
    assert(((constraint1 & (1 << 1)) || (constraint1 & (1 << 4))) == (bool)label1[1]);
    assert(((constraint1 & (1 << 2)) || (constraint1 & (1 << 5))) == (bool)label1[2]);
    
    ElTopo::Vec3c result;  // if either endpoint is constrained, the collapsed point shold be constrained. more specifically it should be on all the walls any of the two endpoints is on (implemented in generate_collapsed_position())
    int result_constraint = (constraint0 | constraint1);
    result[0] = ((result_constraint & (1 << 0)) || (result_constraint & (1 << 3)));
    result[1] = ((result_constraint & (1 << 1)) || (result_constraint & (1 << 4)));
    result[2] = ((result_constraint & (1 << 2)) || (result_constraint & (1 << 5)));
    
    return result;
}

ElTopo::Vec3c MeanCurvatureMultiDriver::generate_split_solid_label(ElTopo::SurfTrack & st, size_t v0, size_t v1, const ElTopo::Vec3c & label0, const ElTopo::Vec3c & label1)
{
    ElTopo::Vec3d x0 = st.get_position(v0);
    ElTopo::Vec3d x1 = st.get_position(v1);
    
    int constraint0 = onBBWall(x0);
    int constraint1 = onBBWall(x1);
    
    assert(((constraint0 & (1 << 0)) || (constraint0 & (1 << 3))) == (bool)label0[0]);
    assert(((constraint0 & (1 << 1)) || (constraint0 & (1 << 4))) == (bool)label0[1]);
    assert(((constraint0 & (1 << 2)) || (constraint0 & (1 << 5))) == (bool)label0[2]);
    assert(((constraint1 & (1 << 0)) || (constraint1 & (1 << 3))) == (bool)label1[0]);
    assert(((constraint1 & (1 << 1)) || (constraint1 & (1 << 4))) == (bool)label1[1]);
    assert(((constraint1 & (1 << 2)) || (constraint1 & (1 << 5))) == (bool)label1[2]);
    
    ElTopo::Vec3c result;  // the splitting midpoint has a positive constraint label only if the two endpoints are on a same wall (sharing a bit in their constraint bitfield representation)
    int result_constraint = (constraint0 & constraint1);
    result[0] = ((result_constraint & (1 << 0)) || (result_constraint & (1 << 3)));
    result[1] = ((result_constraint & (1 << 1)) || (result_constraint & (1 << 4)));
    result[2] = ((result_constraint & (1 << 2)) || (result_constraint & (1 << 5)));
    
    return result;
    
}

bool MeanCurvatureMultiDriver::generate_edge_popped_positions(ElTopo::SurfTrack & st, size_t oldv, const ElTopo::Vec2i & cut, ElTopo::Vec3d & pos_upper, ElTopo::Vec3d & pos_lower)
{
    int original_constraint = onBBWall(st.get_position(oldv));
    
    pos_upper = enforceBBWallConstraint(pos_upper, original_constraint);
    pos_lower = enforceBBWallConstraint(pos_lower, original_constraint);
    
    return true;
}

bool MeanCurvatureMultiDriver::generate_vertex_popped_positions(ElTopo::SurfTrack & st, size_t oldv, int A, int B, ElTopo::Vec3d & pos_a, ElTopo::Vec3d & pos_b)
{
    int original_constraint = onBBWall(st.get_position(oldv));
    
    pos_a = enforceBBWallConstraint(pos_a, original_constraint);
    pos_b = enforceBBWallConstraint(pos_b, original_constraint);
    
    return true;
}

bool MeanCurvatureMultiDriver::solid_edge_is_feature(const ElTopo::SurfTrack & st, size_t e)
{
    int constraint0 = onBBWall(st.get_position(st.m_mesh.m_edges[e][0]));
    int constraint1 = onBBWall(st.get_position(st.m_mesh.m_edges[e][1]));
    
    if (constraint0 & constraint1)  // edge is completely inside a wall
        return true;
    else
        return false;
}

