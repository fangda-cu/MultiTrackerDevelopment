// ---------------------------------------------------------
//
//  t1transition.cpp
//  Fang Da 2013
//  
//  Functions handling T1 transitions (edge popping and vertex popping).
//
// ---------------------------------------------------------

#include <queue>
#include <set>

#include <t1transition.h>
#include <broadphase.h>
#include <collisionqueries.h>
#include <runstats.h>
#include <subdivisionscheme.h>
#include <surftrack.h>
#include <trianglequality.h>


// ---------------------------------------------------------
//  Extern globals
// ---------------------------------------------------------

namespace ElTopo {
    
extern RunStats g_stats;

struct T1Transition::InteriorStencil
{
    Vec3st vertices;
    Vec3st vertex_indices;
    Vec2st edges;
    Vec2st edge_indices;
};

    

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Constructor.  Active SurfTrack object must be supplied.
///
// ---------------------------------------------------------

T1Transition::T1Transition(SurfTrack & surf, bool remesh_boundaries) :
    m_remesh_boundaries(remesh_boundaries),
    m_surf(surf)
{

}


// --------------------------------------------------------
///
/// Perform edge popping
///
// --------------------------------------------------------

bool T1Transition::pop_edges()
{
    if ( m_surf.m_verbose )
        std::cout << "---------------------- T1 Transition: Edge popping ----------------------" << std::endl;
    
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;
    bool pop_occurred = false;
    
    // find all the X-junction edges (higher valence junctions are not considered. hopefully they are rare.)
    // TODO: robustly handle higher valence junctions
    std::vector<size_t> xjunctions;
    for (size_t i = 0; i < mesh.m_edges.size(); i++)
        if (mesh.m_edge_to_triangle_map[i].size() == 4)
            xjunctions.push_back(i);
    
    // arrange the X-junction edges into groups of two consecutive X-junction edges with the same cut as much as possible. 
    // unpaired edge will be popped by splitting in the middle first. pairing them as much as possible helps reducing unnecessary splitting.
    std::vector<std::pair<Vec2st, Mat2i> > xjgroups;    // each element is a pair of x-junction edges, along with the id of the two regions that would end up adjacent after pull-apart, and the id of the other two regions
    for (size_t i = 0; i < xjunctions.size(); i++)
    {
        Mat2i cut = cut_x_junction_edge(xjunctions[i]);
        bool paired = false;
        for (size_t j = 0; j < xjgroups.size(); j++)
        {
            if (cut != xjgroups[j].second)  // different cut, or different surrounding environment
                continue;
            
            if (xjgroups[j].first[1] < mesh.ne())   // this pair is already completed
                continue;
            
            assert(xjgroups[j].first[0] < mesh.ne());
            if (mesh.get_common_vertex(xjunctions[i], xjgroups[j].first[0]) < mesh.nv())
            {
                // put this edge into this group
                paired = true;
                xjgroups[j].first[1] = xjunctions[i];
            }
        }
        
        if (!paired)    // new group
            xjgroups.push_back(std::pair<Vec2st, Mat2i>(Vec2st(xjunctions[i], static_cast<size_t>(~0)), cut));
    }
    
    // process the groups one after another
    for (size_t i = 0; i < xjgroups.size(); i++)
    {
        if (xjgroups[i].first[1] >= mesh.ne())
        {
            // if a group has only one edge, split the edge in half to create two consecutive edges
            // the use of midpoint should be automatically collision safe (same assumtion in EdgeSplitter)
            size_t edge = xjgroups[i].first[0];
            
            size_t v0 = mesh.m_edges[edge][0];
            size_t v1 = mesh.m_edges[edge][1];
            
            size_t midpoint = m_surf.add_vertex(0.5 * (m_surf.get_position(v0) + m_surf.get_position(v1)), 0.5 * (m_surf.m_masses[v0] + m_surf.m_masses[v1]));
            
            std::vector<size_t> faces_to_delete;
            std::vector<Vec3st> faces_to_create;
            std::vector<Vec2i> face_labels_to_create;
            std::vector<size_t> faces_created;
            
            for (size_t j = 0; j < mesh.m_edge_to_triangle_map[edge].size(); j++)
            {
                size_t triangle = mesh.m_edge_to_triangle_map[edge][j];
                size_t third_vertex = mesh.get_third_vertex(v0, v1, mesh.get_triangle(triangle));
                
                faces_to_delete.push_back(triangle);
                
                Vec2i label = mesh.get_triangle_label(triangle);
                if (mesh.oriented(v0, v1, mesh.get_triangle(triangle)))
                {
                    faces_to_create.push_back(Vec3st(v0, midpoint, third_vertex));
                    faces_to_create.push_back(Vec3st(midpoint, v1, third_vertex));
                } else
                {
                    faces_to_create.push_back(Vec3st(v1, midpoint, third_vertex));
                    faces_to_create.push_back(Vec3st(midpoint, v0, third_vertex));
                }
                face_labels_to_create.push_back(label);
                face_labels_to_create.push_back(label);
            }
            
            // apply the changes
            for (size_t j = 0; j < faces_to_delete.size(); j++)
                m_surf.remove_triangle(faces_to_delete[j]);
            
            assert(faces_to_create.size() == face_labels_to_create.size());
            for (size_t j = 0; j < faces_to_create.size(); j++)
            {
                size_t nf = m_surf.add_triangle(faces_to_create[j]);
                mesh.set_triangle_label(nf, face_labels_to_create[j]);
                faces_created.push_back(nf);
            }
            
            // update the group
            xjgroups[i].first[0] = mesh.get_edge_index(v0, midpoint);
            xjgroups[i].first[1] = mesh.get_edge_index(midpoint, v1);
            
            // add a history item
            MeshUpdateEvent split(MeshUpdateEvent::EDGE_SPLIT);
            split.m_v0 = v0;
            split.m_v1 = v1;
            split.m_vert_position = m_surf.get_newposition(midpoint);
            split.m_created_verts.push_back(midpoint);
            split.m_created_vert_data.push_back(split.m_vert_position);
            split.m_deleted_tris = faces_to_delete;
            split.m_created_tris = faces_created;
            split.m_created_tri_data = faces_to_create;
            split.m_created_tri_labels = face_labels_to_create;
            
            m_surf.m_mesh_change_history.push_back(split);
        }
        
        // now the group has exactly two edges. we try to pull apart the interior vertex, ensuring collision safety
        Mat2i cut_regions = xjgroups[i].second;
        Vec2i cut = cut_regions[0];        
        
        size_t edge0 = xjgroups[i].first[0];
        size_t edge1 = xjgroups[i].first[1];
        assert(edge0 < mesh.ne());
        assert(edge1 < mesh.ne());
        assert(edge0 != edge1);
        
        size_t v = mesh.get_common_vertex(edge0, edge1);
        assert(v < mesh.nv());
        
        size_t v0 = (v == mesh.m_edges[edge0][0] ? mesh.m_edges[edge0][1] : mesh.m_edges[edge0][0]);
        size_t v1 = (v == mesh.m_edges[edge1][0] ? mesh.m_edges[edge1][1] : mesh.m_edges[edge1][0]);
        assert(v0 != v1);
        
        bool edge0_oriented = (v == mesh.m_edges[edge0][1]);
        bool edge1_oriented = (v == mesh.m_edges[edge1][0]);

        bool original_constraint = m_surf.m_mesh.get_vertex_constraint_label(v);
        Vec3d original_position = m_surf.get_position(v);
        
        int upper_region = -1;  // upper means the region is on the top when looking from vertex v0 to vertex v (may or may not be the direction of edge0), with region cut.x() on the left and cut.y() on the right.
        int lower_region = -1;
        std::vector<size_t> region0faces;
        for (size_t k = 0; k < mesh.m_edge_to_triangle_map[edge0].size(); k++)
            if (mesh.get_triangle_label(mesh.m_edge_to_triangle_map[edge0][k])[0] == cut[0] ||
                mesh.get_triangle_label(mesh.m_edge_to_triangle_map[edge0][k])[1] == cut[0])
                region0faces.push_back(mesh.m_edge_to_triangle_map[edge0][k]);
        assert(region0faces.size() == 2);
        
        Vec2i label0 = mesh.get_triangle_label(region0faces[0]);
        upper_region = (label0[0] == cut[0] ? label0[1] : label0[0]);
        Vec2i label1 = mesh.get_triangle_label(region0faces[1]);
        lower_region = (label1[0] == cut[0] ? label1[1] : label1[0]);
        
        if (( mesh.oriented(v0, v, mesh.get_triangle(region0faces[0])) && label0[1] == cut[0]) ||
            (!mesh.oriented(v0, v, mesh.get_triangle(region0faces[0])) && label0[0] == cut[0]))
        {
            std::swap(upper_region, lower_region);
            assert(( mesh.oriented(v0, v, mesh.get_triangle(region0faces[1])) && label1[0] == cut[0]) ||
                   (!mesh.oriented(v0, v, mesh.get_triangle(region0faces[1])) && label1[1] == cut[0]));
        }
        assert(upper_region >= 0);
        assert(lower_region >= 0);
        
        std::vector<Vec3d> upper_vertices;
        std::vector<Vec3d> lower_vertices;
        
        for (size_t k = 0; k < mesh.m_edge_to_triangle_map[edge0].size(); k++)
        {
            size_t triangle = mesh.m_edge_to_triangle_map[edge0][k];
            size_t v2 = mesh.get_third_vertex(edge0, triangle);
            
            Vec2i label = mesh.get_triangle_label(triangle);
            assert(label[0] == upper_region || label[1] == upper_region || label[0] == lower_region || label[1] == lower_region);
            ((label[0] == upper_region || label[1] == upper_region) ? upper_vertices : lower_vertices).push_back(m_surf.get_position(v2));
        }
        
        for (size_t k = 0; k < mesh.m_edge_to_triangle_map[edge1].size(); k++)
        {
            size_t triangle = mesh.m_edge_to_triangle_map[edge1][k];
            size_t v2 = mesh.get_third_vertex(edge1, triangle);
            
            Vec2i label = mesh.get_triangle_label(triangle);
            assert(label[0] == upper_region || label[1] == upper_region || label[0] == lower_region || label[1] == lower_region);
            ((label[0] == upper_region || label[1] == upper_region) ? upper_vertices : lower_vertices).push_back(m_surf.get_position(v2));
        }
        
        Vec3d upper_vertices_mean(0, 0, 0);
        for (size_t i = 0; i < upper_vertices.size(); i++)
            upper_vertices_mean += upper_vertices[i];
        upper_vertices_mean /= upper_vertices.size();
        
        Vec3d lower_vertices_mean(0, 0, 0);
        for (size_t i = 0; i < lower_vertices.size(); i++)
            lower_vertices_mean += lower_vertices[i];
        lower_vertices_mean /= lower_vertices.size();
        
        Vec3d pull_apart_offset = (upper_vertices_mean - lower_vertices_mean);
        pull_apart_offset /= mag(pull_apart_offset);
        pull_apart_offset *= mag(m_surf.get_position(v1) - m_surf.get_position(v0));
        
        // compute the desired destination positions for the new vertices
        Vec3d upper_junction_desired_position = original_position + pull_apart_offset * 0.1;
        Vec3d lower_junction_desired_position = original_position - pull_apart_offset * 0.1;
        
        // enforce constraints
        if (original_constraint)
        {
            assert(m_surf.m_constrained_vertices_callback);
            m_surf.m_constrained_vertices_callback->generate_edge_popped_positions(m_surf, v, cut, upper_junction_desired_position, lower_junction_desired_position);
        }
        
        // collision test
        // sort the incident faces and edges into those that go with nv0, and those that go with nv1 (the two groups are not necessarily disjoint)
        std::vector<size_t> upper_faces;
        std::vector<size_t> upper_edges;
        
        for (size_t j = 0; j < mesh.m_vertex_to_triangle_map[v].size(); j++)
        {
            size_t triangle = mesh.m_vertex_to_triangle_map[v][j];
            
            Vec2i label = mesh.get_triangle_label(triangle);
            assert(label[0] == upper_region || label[1] == upper_region || label[0] == lower_region || label[1] == lower_region);
            if (label[0] == upper_region || label[1] == upper_region)
                upper_faces.push_back(triangle);
        }
        
        for (size_t j = 0; j < mesh.m_vertex_to_edge_map[v].size(); j++)
        {
            size_t edge = mesh.m_vertex_to_edge_map[v][j];
            
            bool upper = false;
            for (size_t k = 0; k < mesh.m_edge_to_triangle_map[edge].size(); k++)
            {
                Vec2i label = mesh.get_triangle_label(mesh.m_edge_to_triangle_map[edge][k]);
                if (label[0] == upper_region || label[1] == upper_region)
                    upper = true;
            }
            
            if (upper)
                upper_edges.push_back(edge);
        }
        
        // move v to lower_junction_desired_position
        if (vertex_pseudo_motion_introduces_collision(v, original_position, lower_junction_desired_position))
        {
            if (m_surf.m_verbose)
                std::cout << "Edge popping: collision introduced." << std::endl;
            continue;
        }
        
        // move the upper junction from v to upper_junction_desired_position
        m_surf.set_position(v, lower_junction_desired_position);
        if (vertex_pseudo_motion_introduces_collision(v, lower_junction_desired_position, upper_junction_desired_position, upper_faces, upper_edges))
        {
            if (m_surf.m_verbose)
                std::cout << "Edge popping: collision introduced." << std::endl;
            continue;
        }
        
        // check intersection in the final configuration
        const std::vector<Vec3d> & x = m_surf.get_positions();
        bool collision = false;
        
        // point-tet
        for (size_t j = 0; j < upper_faces.size(); j++)
        {
            Vec3st t = mesh.get_triangle(upper_faces[j]);
            
            Vec3d low, high;
            minmax(x[t[0]], x[t[1]], x[t[2]], upper_junction_desired_position, low, high);

            std::vector<size_t> overlapping_vertices;
            m_surf.m_broad_phase->get_potential_vertex_collisions(low, high, true, true, overlapping_vertices);
            
            for (size_t k = 0; k < overlapping_vertices.size(); k++) 
            { 
                size_t ov = overlapping_vertices[k];
                if (ov == t[0] || ov == t[2] || ov == t[1])
                    continue;
                
                if (point_tetrahedron_intersection(x[ov], ov, x[t[0]], t[0], x[t[1]], t[1], x[t[2]], t[2], upper_junction_desired_position, mesh.nv()))
                    collision = true;
            }
        }
        
        // edge-triangle
        for (size_t j = 0; j < upper_faces.size(); j++)
        {
            Vec3st t = mesh.get_triangle(upper_faces[j]);
            if (t[1] == v) std::swap(t[0], t[1]);
            if (t[2] == v) std::swap(t[0], t[2]);
            
            Vec3d low, high;
            minmax(x[t[1]], x[t[2]], upper_junction_desired_position, low, high);
            
            std::vector<size_t> overlapping_edges;
            m_surf.m_broad_phase->get_potential_edge_collisions(low, high, true, true, overlapping_edges);
            
            for (size_t k = 0; k < overlapping_edges.size(); k++) 
            { 
                const Vec2st & e = mesh.m_edges[overlapping_edges[k]];
                if (e[0] == t[1] || e[1] == t[1] || e[0] == t[2] || e[1] == t[2])
                    continue;
                
                bool incident = false;
                for (size_t l = 0; l < upper_edges.size(); l++)
                    if (upper_edges[l] == overlapping_edges[k])
                    {
                        incident = true;
                        break;
                    }
                if (incident)
                    continue;
                
                if (segment_triangle_intersection(x[e[0]], e[0], x[e[1]], e[1], x[t[1]], t[1], x[t[2]], t[2], upper_junction_desired_position, mesh.nv(), true))
                    collision = true;
            }
        }
        
        // triangle-edge
        for (size_t j = 0; j < upper_edges.size(); j++)
        {
            Vec2st e = mesh.m_edges[upper_edges[j]];
            if (e[1] == v) std::swap(e[0], e[1]);
            
            Vec3d low, high;
            minmax(x[e[1]], upper_junction_desired_position, low, high);
            
            std::vector<size_t> overlapping_triangles;
            m_surf.m_broad_phase->get_potential_triangle_collisions(low, high, true, true, overlapping_triangles);
            
            for (size_t k = 0; k < overlapping_triangles.size(); k++)
            {
                const Vec3st & t = mesh.get_triangle(overlapping_triangles[k]);
                if (e[1] == t[0] || e[1] == t[1] || e[1] == t[2])
                    continue;
                
                bool incident = false;
                for (size_t l = 0; l < upper_faces.size(); l++)
                    if (upper_faces[l] == overlapping_triangles[k])
                    {
                        incident = true;
                        break;
                    }
                if (incident)
                    continue;
                
                if (segment_triangle_intersection(x[e[1]], e[1], upper_junction_desired_position, mesh.nv(), x[t[0]], t[0], x[t[1]], t[1], x[t[2]], t[2], true))
                    collision = true;
            }
        }
                
        if (collision)
        {
            if (m_surf.m_verbose)
                std::cout << "Edge popping: collision introduced." << std::endl;
            continue;
        }        
        
        // collision test give green light. pull apart the interior vertex
        size_t nv0 = m_surf.add_vertex(upper_junction_desired_position, m_surf.m_masses[v]);
        size_t nv1 = m_surf.add_vertex(lower_junction_desired_position, m_surf.m_masses[v]);
        
        m_surf.set_remesh_velocity(nv0, m_surf.get_remesh_velocity(v));
        m_surf.set_remesh_velocity(nv1, m_surf.get_remesh_velocity(v));
        mesh.set_vertex_constraint_label(nv0, original_constraint);
        mesh.set_vertex_constraint_label(nv1, original_constraint);
        
        std::vector<size_t> verts_to_delete;
        std::vector<Vec3d> verts_to_create;
        std::vector<size_t> verts_created;
        
        verts_to_delete.push_back(v);
        verts_created.push_back(nv0);
        verts_created.push_back(nv1);
        verts_to_create.push_back(upper_junction_desired_position);
        verts_to_create.push_back(lower_junction_desired_position);
        
        // retriangulate around the pull-apart regions
        std::vector<size_t> faces_to_delete;
        std::vector<Vec3st> faces_to_create;
        std::vector<Vec2i> face_labels_to_create;
        std::vector<size_t> faces_created;
        
        // first, update the faces incident to the X-junction edges
        for (size_t j = 0; j < 2; j++)
        {
            size_t edge = (j == 0 ? edge0 : edge1);
            bool edge_oriented = (j == 0 ? edge0_oriented : edge1_oriented);
            
            for (size_t k = 0; k < mesh.m_edge_to_triangle_map[edge].size(); k++)
            {
                size_t triangle = mesh.m_edge_to_triangle_map[edge][k];
                faces_to_delete.push_back(triangle);
                
                size_t v2 = mesh.get_third_vertex(edge, triangle);
                
                Vec2i label = mesh.get_triangle_label(triangle);
                assert(label[0] == upper_region || label[1] == upper_region || label[0] == lower_region || label[1] == lower_region);
                
                size_t ev0 = (j == 0 ? v0 : ((label[0] == upper_region || label[1] == upper_region) ? nv0 : nv1));
                size_t ev1 = (j == 1 ? v1 : ((label[0] == upper_region || label[1] == upper_region) ? nv0 : nv1));
                
                if (mesh.oriented(mesh.m_edges[edge][0], mesh.m_edges[edge][1], mesh.get_triangle(triangle)) == edge_oriented)
                    faces_to_create.push_back(Vec3st(ev0, ev1, v2));
                else
                    faces_to_create.push_back(Vec3st(ev1, ev0, v2));
                face_labels_to_create.push_back(label);
            }
        }     
        
        // second, update the faces incident only to the interior vertex v but not to any X-junction edge
        for (size_t k = 0; k < mesh.m_vertex_to_triangle_map[v].size(); k++)
        {
            size_t triangle = mesh.m_vertex_to_triangle_map[v][k];
            
            // skip if the triangle contains edge0 or edge1
            if (mesh.triangle_contains_edge(mesh.get_triangle(triangle), mesh.m_edges[edge0]) ||
                mesh.triangle_contains_edge(mesh.get_triangle(triangle), mesh.m_edges[edge1]))
                continue;
            
            faces_to_delete.push_back(triangle);
            
            // find the edge in triangle triangle that's opposite to vertex v
            size_t l = 0;
            size_t edge2 = static_cast<size_t>(~0);
            for (l = 0; l < 3; l++)
            {
                if (mesh.m_edges[mesh.m_triangle_to_edge_map[triangle][l]][0] != v &&
                    mesh.m_edges[mesh.m_triangle_to_edge_map[triangle][l]][1] != v)
                    edge2 = mesh.m_triangle_to_edge_map[triangle][l];
            }
            assert(edge2 < mesh.ne());
            size_t v1 = mesh.m_edges[edge2][0];
            size_t v2 = mesh.m_edges[edge2][1];
            
            Vec2i label = mesh.get_triangle_label(triangle);
            assert(label[0] == upper_region || label[1] == upper_region || label[0] == lower_region || label[1] == lower_region);
            
            size_t nv = ((label[0] == upper_region || label[1] == upper_region) ? nv0 : nv1);
            
            if (mesh.oriented(v1, v2, mesh.get_triangle(triangle)))
                faces_to_create.push_back(Vec3st(nv, v1, v2));
            else
                faces_to_create.push_back(Vec3st(nv, v2, v1));
            face_labels_to_create.push_back(label);
        }
        
        // triangulate the new interface in the middle
        for (size_t j = 0; j < 2; j++)
        {
            Vec2st left_junction  = (j == 0 ? Vec2st(v0, v0) : Vec2st(nv0, nv1));
            Vec2st right_junction = (j == 1 ? Vec2st(v1, v1) : Vec2st(nv0, nv1));
            if (left_junction[0] != left_junction[1] && right_junction[0] != right_junction[1])
            {
                // quad
                faces_to_create.push_back(Vec3st(left_junction[0],  left_junction[1],  right_junction[0]));
                faces_to_create.push_back(Vec3st(right_junction[1], right_junction[0], left_junction[1] ));
                face_labels_to_create.push_back(cut);
                face_labels_to_create.push_back(cut);
            } else if (left_junction[0] != left_junction[1])
            {
                // triangle, pointing to the right
                faces_to_create.push_back(Vec3st(left_junction[0],  left_junction[1],  right_junction[0]));
                face_labels_to_create.push_back(cut);
            } else if (right_junction[0] != right_junction[1])
            {
                // triangle, pointing to the left
                faces_to_create.push_back(Vec3st(left_junction[0],  right_junction[1], right_junction[0]));
                face_labels_to_create.push_back(cut);
            } else
            {
                // no hole to begin with, nothing to do
            }
        }
        
        // apply the deletion/creation
        assert(faces_to_create.size() == face_labels_to_create.size());
        for (size_t j = 0; j < faces_to_create.size(); j++)
        {
            size_t nf = m_surf.add_triangle(faces_to_create[j]);
            mesh.set_triangle_label(nf, face_labels_to_create[j]);
            faces_created.push_back(nf);
        }
        
        for (size_t j = 0; j < faces_to_delete.size(); j++)
            m_surf.remove_triangle(faces_to_delete[j]);
        
        // Add to new history log
        MeshUpdateEvent edgepop(MeshUpdateEvent::EDGE_POP);
        edgepop.m_deleted_tris = faces_to_delete;
        edgepop.m_created_tris = faces_created;
        edgepop.m_created_tri_data = faces_to_create;
        edgepop.m_created_tri_labels = face_labels_to_create;
        edgepop.m_deleted_verts = verts_to_delete;
        edgepop.m_created_verts = verts_created;
        edgepop.m_created_vert_data = verts_to_create;
        m_surf.m_mesh_change_history.push_back(edgepop);
        
        pop_occurred = true;
    }
    
    return pop_occurred;
}


// --------------------------------------------------------
///
/// Perform vertex popping
///
// --------------------------------------------------------

bool T1Transition::pop_vertices()
{
    if ( m_surf.m_verbose )
        std::cout << "---------------------- T1 Transition: Vertex popping ----------------------" << std::endl;
    
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;
    bool pop_occurred = false;
    
    // Pull apart the X-junciton vertices
    // In fact this function attempts to make the region adjacency graph complete on ever vertex in the mesh (see explanation below), or in
    //  other words, so that every region incident to a vertex is adjacent to every other region on that vertex through a face. No two regions
    //  can be adjacent only through a vertex, as long as the surface tension force at the vertex is tensile instead of compressive.
    
    // find the region count
    int max_region = -1;
    for (size_t i = 0; i < mesh.nt(); i++)
    {
        Vec2i label = mesh.get_triangle_label(i);
        assert(label[0] >= 0);
        assert(label[1] >= 0);
        if (label[0] > max_region) max_region = label[0];
        if (label[1] > max_region) max_region = label[1];
    }
    int nregion = max_region + 1;
    
    // this is the incident matrix for a region graph. not every col/row need to be filled for a particular vertex becuase the vertex may not be incident to all regions.
    std::vector<std::vector<int> > region_graph(nregion, std::vector<int>(nregion, 0));
    
    // prepare an initial list of unprocessed vertices
    std::vector<size_t> vertices_to_process;
    for (size_t i = 0; i < mesh.nv(); i++)
        vertices_to_process.push_back(i);
    
    // process the vertices one after another
    while (vertices_to_process.size() > 0)
    {
        size_t xj = vertices_to_process.back();
        vertices_to_process.pop_back();
        
        std::set<int> vertex_regions_set;
        for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[xj].size(); i++)
        {
            Vec2i label = mesh.get_triangle_label(mesh.m_vertex_to_triangle_map[xj][i]);
            vertex_regions_set.insert(label[0]);
            vertex_regions_set.insert(label[1]);
        }
        
        std::vector<int> vertex_regions;
        vertex_regions.assign(vertex_regions_set.begin(), vertex_regions_set.end());
        
        // cull away the interior vertices of plateau borders (which are the majority)
        if (vertex_regions_set.size() < 4)
            continue;
        
        //
        // Pull apart strategy: in order to cope with various configurations (such as a region-valence-5 vertex being a triple junction, or
        //  junctions on the BB walls), the following strategy is adopted:
        //
        // The regions incident on the vertex form a graph, with each region being a node and an edge exists between two regions iff the two
        //  regions share a triangle face in the mesh. The objective at the end of this processing, is to make sure this graph is complete for
        //  every vertex. For example the original region-valence-4 X-junction vertex (the "hourglass" junction) forms a graph like this:
        //
        //  A o--o B
        //    | /|
        //    |/ |
        //  C o--o D
        //
        //  where the edge between node A and D (the two bulbs of the hourglass) is missing because the two regions only meet at one vertex
        //  (the hourglass neck).
        //
        // The processing here makes an incomplete graph complete by pulling apart vertices. If two nodes (e.g. A and D in the figure above)
        //  do not share an edge, it means the center vertex is the only interface between the two regions, and it needs to be pulled apart.
        //  Region A and D each keep one of the two duplicates of the vertex. If a region has edges with both regions A and D, such as region
        //  C and B, will contain both duplicates (the region's shape is turned from a cone into a flat screwdriver's tip). Now look at the
        //  resulting graphs for the two duplicate vertices. For the duplicate that goes with region A, region D is no longer incident and
        //  thus removed from the graph:
        //
        //  A o--o B
        //    | /
        //    |/
        //  C o
        //
        //  and this is now a complete graph. Similarly the graph for the duplicate that goes with region D will not have node A, and it is 
        //  complete too. Of course in more complex scenarios the graphs may not directly become complete immediately after we process one
        //  missing edge. We will visit the two resulting vertices again as if they are a regular vertex that may need to be pulled apart too.
        //
        // The algorithm pseudocode is as following:
        //
        //  Pop a vertex from the stack of vertices to be processed, construct a graph of region adjacency
        //  If the graph is already complete, skip this vertex
        //  Pick an arbitrary pair of unconnected nodes A and B from the graph (there may be more than one pair)
        //  Pull apart the vertex into two (a and b, corresponding to region A and B respectively), initialized to have the same coordinates and then pulled apart
        //  For each face incident to the center vertex in the original mesh
        //    If it is incident to region A, update it to use vertex a
        //    If it is incident to region B, update it to use vertex b
        //    Otherwise, update it to use vertex b
        //  For each edge incident to the center vertex in the original mesh, and incident to any face with A label
        //    Create a new triangle with both a and b
        //  Push vertex a and b on the stack to be visited next.
        //
        // The first loop essentially pulls region A from everything else around the vertex, leaving a ring of blank around the vertex. Then
        //  the next loop fills the blank. Note that this algorithm works even in presence of unresolved X-junction edges, if any (due to 
        //  collision etc). Another interpretation is: region A and region B remain intact; all the other triangles (incident to neither A
        //  nor B) form a number of triangle fans, which may share edges/triangles between them. There are a number of "head" edges, which 
        //  are located at the head of the fans from region A to region B. These edges are the edges incident to region A. The second loop 
        //  basically sweeps these head edges into a triangles (from a to b).
        //
        // Collision test is performed before the operation to ensure collision safety. First the center vertex is moved to the position of 
        //  new vertex b without changing connectivity. Collision is checked for this motion. Then vertex a is separated from vertex b, 
        //  bringing all region A faces along with it. Collision is checked for this motion, using only region A faces as incident faces and
        //  region A edges as incident edges. This CCD only detects the collision of this motion against the rest of the mesh, excluding all
        //  faces/edges connected to both a and b (b faces/edges are excluded too because initially a and b coincide). Then an instantaneous 
        //  intersection test is perform that finds intersections in the final configuration.
        //  
        //
        
        // construct the region graph
        region_graph.assign(nregion, std::vector<int>(nregion, 0));
        for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[xj].size(); i++)
        {
            Vec2i label = mesh.get_triangle_label(mesh.m_vertex_to_triangle_map[xj][i]);
            region_graph[label[0]][label[1]] = 1;
            region_graph[label[1]][label[0]] = 1;
        }
        
        // find a missing edge
        int A = -1;
        int B = -1;
        Vec3d pull_apart_direction;
        for (size_t i = 0; i < vertex_regions.size() && A < 0; i++)
        {
            for (size_t j = i + 1; j < vertex_regions.size() && A < 0; j++)
            {
                if (region_graph[vertex_regions[i]][vertex_regions[j]] == 0)
                {
                    if (should_pull_vertex_apart(xj, vertex_regions[i], vertex_regions[j], pull_apart_direction))
                    {
                        A = vertex_regions[i];
                        B = vertex_regions[j];
                        break;
                    }
                }
            }
        }
        
        // skip the vertex if the graph is complete already
        if (A < 0)
            continue;
        
        bool original_constraint = m_surf.m_mesh.get_vertex_constraint_label(xj);
        Vec3d original_position = m_surf.get_position(xj);
        
        double mean_edge_length = 0;
        int edge_count = 0;
        for (size_t i = 0; i < mesh.m_vertex_to_edge_map[xj].size(); i++)
        {
            size_t v0 = mesh.m_edges[mesh.m_vertex_to_edge_map[xj][i]][0];
            size_t v1 = mesh.m_edges[mesh.m_vertex_to_edge_map[xj][i]][1];
            mean_edge_length += mag(m_surf.get_position(v1) - m_surf.get_position(v0));
            edge_count++;
        }
        assert(edge_count > 0);
        mean_edge_length /= edge_count;
        
        Vec3d pull_apart_offset = pull_apart_direction * mean_edge_length;
        
        // compute the desired destination positions, enforcing constraints
        Vec3d a_desired_position = original_position + pull_apart_offset * 0.1;
        Vec3d b_desired_position = original_position - pull_apart_offset * 0.1;
        size_t a = static_cast<size_t>(~0);
        size_t b = static_cast<size_t>(~0);

        if (original_constraint)
        {
            assert(m_surf.m_constrained_vertices_callback);            
            m_surf.m_constrained_vertices_callback->generate_vertex_popped_positions(m_surf, xj, A, B, a_desired_position, b_desired_position);
        }
        
        // collision test
        // sort the incident faces and edges into those that go with nv0, and those that go with nv1 (the two groups are not necessarily disjoint)
        std::vector<size_t> A_faces;
        std::vector<size_t> A_edges;
        
        for (size_t j = 0; j < mesh.m_vertex_to_triangle_map[xj].size(); j++)
        {
            size_t triangle = mesh.m_vertex_to_triangle_map[xj][j];
            
            Vec2i label = mesh.get_triangle_label(triangle);
            if (label[0] == A || label[1] == A)
                A_faces.push_back(triangle);
        }
        
        for (size_t j = 0; j < mesh.m_vertex_to_edge_map[xj].size(); j++)
        {
            size_t edge = mesh.m_vertex_to_edge_map[xj][j];
            
            bool adjA = false;
            for (size_t k = 0; k < mesh.m_edge_to_triangle_map[edge].size(); k++)
            {
                Vec2i label = mesh.get_triangle_label(mesh.m_edge_to_triangle_map[edge][k]);
                if (label[0] == A || label[1] == A)
                    adjA = true;
            }
            
            if (adjA)
                A_edges.push_back(edge);
        }
        
        if (vertex_pseudo_motion_introduces_collision(xj, original_position, b_desired_position))
        {
            if (m_surf.m_verbose)
                std::cout << "Vertex popping: pulling vertex " << xj << " apart introduces collision." << std::endl;
            continue;
        }
        
        m_surf.set_position(xj, b_desired_position);
        if (vertex_pseudo_motion_introduces_collision(xj, b_desired_position, a_desired_position, A_faces, A_edges))
        {
            if (m_surf.m_verbose)
                std::cout << "Vertex popping: pulling vertex " << xj << " apart introduces collision." << std::endl;
            continue;
        }
        
        // check intersection in the final configuration
        const std::vector<Vec3d> & x = m_surf.get_positions();
        bool collision = false;
        
        // point-tet
        for (size_t j = 0; j < A_faces.size(); j++)
        {
            Vec3st t = mesh.get_triangle(A_faces[j]);
            
            Vec3d low, high;
            minmax(x[t[0]], x[t[1]], x[t[2]], a_desired_position, low, high);
            
            std::vector<size_t> overlapping_vertices;
            m_surf.m_broad_phase->get_potential_vertex_collisions(low, high, true, true, overlapping_vertices);
            
            for (size_t k = 0; k < overlapping_vertices.size(); k++) 
            { 
                size_t ov = overlapping_vertices[k];
                if (mesh.m_vertex_to_triangle_map[ov].size() == 0)
                    continue;
                if (ov == t[0] || ov == t[2] || ov == t[1])
                    continue;
                
                if (point_tetrahedron_intersection(x[ov], ov, x[t[0]], t[0], x[t[1]], t[1], x[t[2]], t[2], a_desired_position, mesh.nv()))
                    collision = true;
            }
        }
        
        // edge-triangle
        for (size_t j = 0; j < A_faces.size(); j++)
        {
            Vec3st t = mesh.get_triangle(A_faces[j]);
            if (t[1] == xj) std::swap(t[0], t[1]);
            if (t[2] == xj) std::swap(t[0], t[2]);
            
            Vec3d low, high;
            minmax(x[t[1]], x[t[2]], a_desired_position, low, high);
            
            std::vector<size_t> overlapping_edges;
            m_surf.m_broad_phase->get_potential_edge_collisions(low, high, true, true, overlapping_edges);
            
            for (size_t k = 0; k < overlapping_edges.size(); k++) 
            { 
                const Vec2st & e = mesh.m_edges[overlapping_edges[k]];
                if (e[0] == e[1])
                    continue;
                if (e[0] == t[1] || e[1] == t[1] || e[0] == t[2] || e[1] == t[2])
                    continue;
                
                bool incident = false;
                for (size_t l = 0; l < A_edges.size(); l++)
                    if (A_edges[l] == overlapping_edges[k])
                    {
                        incident = true;
                        break;
                    }
                if (incident)
                    continue;
                
                if (segment_triangle_intersection(x[e[0]], e[0], x[e[1]], e[1], x[t[1]], t[1], x[t[2]], t[2], a_desired_position, mesh.nv(), true))
                    collision = true;
            }
        }
        
        // triangle-edge
        for (size_t j = 0; j < A_edges.size(); j++)
        {
            Vec2st e = mesh.m_edges[A_edges[j]];
            if (e[1] == xj) std::swap(e[0], e[1]);
            
            Vec3d low, high;
            minmax(x[e[1]], a_desired_position, low, high);
            
            std::vector<size_t> overlapping_triangles;
            m_surf.m_broad_phase->get_potential_triangle_collisions(low, high, true, true, overlapping_triangles);
            
            for (size_t k = 0; k < overlapping_triangles.size(); k++)
            {
                const Vec3st & t = mesh.get_triangle(overlapping_triangles[k]);
                if (t[0] == t[1] || t[0] == t[2] || t[1] == t[2])
                    continue;
                if (e[1] == t[0] || e[1] == t[1] || e[1] == t[2])
                    continue;
                
                bool incident = false;
                for (size_t l = 0; l < A_faces.size(); l++)
                    if (A_faces[l] == overlapping_triangles[k])
                    {
                        incident = true;
                        break;
                    }
                if (incident)
                    continue;
                
                if (segment_triangle_intersection(x[e[1]], e[1], a_desired_position, mesh.nv(), x[t[0]], t[0], x[t[1]], t[1], x[t[2]], t[2], true))
                    collision = true;
            }
        }
        
        if (collision)
        {
            if (m_surf.m_verbose)
                std::cout << "Vertex popping: collision introduced." << std::endl;
            continue;
        }        
        
        // pull apart
        std::vector<size_t> verts_to_delete;
        std::vector<Vec3d> verts_to_create;
        std::vector<size_t> verts_created;
        
        a = m_surf.add_vertex(a_desired_position, m_surf.m_masses[xj]);
        b = m_surf.add_vertex(b_desired_position, m_surf.m_masses[xj]);
        
        m_surf.set_remesh_velocity(a, m_surf.get_remesh_velocity(xj));
        m_surf.set_remesh_velocity(b, m_surf.get_remesh_velocity(xj));
        mesh.set_vertex_constraint_label(a, original_constraint);
        mesh.set_vertex_constraint_label(b, original_constraint);
        
        verts_to_delete.push_back(xj);
        verts_to_create.push_back(a_desired_position);
        verts_to_create.push_back(b_desired_position);
        verts_created.push_back(a);
        verts_created.push_back(b);
        
        // update the face connectivities
        std::vector<size_t> faces_to_delete;
        std::vector<Vec3st> faces_to_create;
        std::vector<Vec2i> face_labels_to_create;
        std::vector<size_t> faces_created;
        
        triangulate_popped_vertex(xj, A, B, a, b, faces_to_delete, faces_to_create, face_labels_to_create);
        
        // apply the deletion/addition
        assert(faces_to_create.size() == face_labels_to_create.size());
        for (size_t i = 0; i < faces_to_create.size(); i++)
        {
            size_t nf = m_surf.add_triangle(faces_to_create[i]);
            mesh.set_triangle_label(nf, face_labels_to_create[i]);
            faces_created.push_back(nf);
        }
        
        for (size_t i = 0; i < faces_to_delete.size(); i++)
            m_surf.remove_triangle(faces_to_delete[i]);
        
        // mark the two new vertices a and b as dirty
        vertices_to_process.push_back(a);
        vertices_to_process.push_back(b);
        
        // Add to new history log
        MeshUpdateEvent vertpop(MeshUpdateEvent::VERTEX_POP);
        vertpop.m_deleted_tris = faces_to_delete;
        vertpop.m_created_tris = faces_created;
        vertpop.m_created_tri_data = faces_to_create;
        vertpop.m_created_tri_labels = face_labels_to_create;
        vertpop.m_deleted_verts = verts_to_delete;
        vertpop.m_created_verts = verts_created;
        vertpop.m_created_vert_data = verts_to_create;
        m_surf.m_mesh_change_history.push_back(vertpop);
        
        pop_occurred = true;
        
    }
    
    return pop_occurred;
}

    
    
// --------------------------------------------------------
///
/// Decide the cut direction on an X-junction edge
///
// --------------------------------------------------------
    
Mat2i T1Transition::cut_x_junction_edge(size_t e)
{
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;

    // for now use the angles to decide cut direction
    assert(mesh.m_edge_to_triangle_map[e].size() == 4);
    
    size_t v0 = mesh.m_edges[e][0];
    size_t v1 = mesh.m_edges[e][1];
    
    Vec3d x0 = m_surf.get_position(v0);
    Vec3d x1 = m_surf.get_position(v1);
    
    // find all the regions around this X junction edge e, in CCW order when looking down the edge
    std::vector<int> regions;
    size_t tmp = mesh.m_edge_to_triangle_map[e][0];
    if (mesh.oriented(v0, v1, mesh.get_triangle(tmp)))
    {
        regions.push_back(mesh.get_triangle_label(tmp)[1]);
        regions.push_back(mesh.get_triangle_label(tmp)[0]);
    } else
    {
        regions.push_back(mesh.get_triangle_label(tmp)[0]);
        regions.push_back(mesh.get_triangle_label(tmp)[1]);
    }
    
    while (true)
    {
        size_t s = regions.size();
        for (size_t i = 0; i < mesh.m_edge_to_triangle_map[e].size(); i++)
        {
            Vec2i label = mesh.get_triangle_label(mesh.m_edge_to_triangle_map[e][i]);
            if (label[0] == regions.back() && label[1] != *(regions.rbegin() + 1) && label[1] != regions.front())
            {
                regions.push_back(label[1]);
                break;
            }
            if (label[1] == regions.back() && label[0] != *(regions.rbegin() + 1) && label[0] != regions.front())
            {
                regions.push_back(label[0]);
                break;
            }
        }
        assert(s == regions.size() || s + 1 == regions.size());
        if (s == regions.size())
            break;
    }
    assert(regions.size() == 4);
    
    // create an arbitrary 2D frame in the cross section plane of the edge e
    Vec3d t = x1 - x0;
    t /= mag(t);
    Vec3d u = (::fabs(t[0]) <= ::fabs(t[1]) && ::fabs(t[0]) <= ::fabs(t[2]) ? Vec3d(0, t[2], -t[1]) : (::fabs(t[1]) <= ::fabs(t[2]) ? Vec3d(-t[2], 0, t[0]) : Vec3d(t[1], -t[0], 0)));
    u /= mag(u);
    Vec3d v = cross(u, t);
    v /= mag(v);
    
    // compute the angle of each region
    std::vector<double> regionangles;
    for (size_t i = 0; i < regions.size(); i++)
    {
        std::vector<size_t> regionfaces;
        for (size_t j = 0; j < mesh.m_edge_to_triangle_map[e].size(); j++)
        {
            Vec2i label = mesh.get_triangle_label(mesh.m_edge_to_triangle_map[e][j]);
            if (label[0] == regions[i] || label[1] == regions[i])
                regionfaces.push_back(mesh.m_edge_to_triangle_map[e][j]);
        }
        assert(regionfaces.size() == 2);
        
        Vec2i label = mesh.get_triangle_label(regionfaces[0]);
        if ((label[0] == regions[i] ? label[1] : label[0]) == regions[(i + 1) % 4])
            std::swap(regionfaces[0], regionfaces[1]);  // this ensures CCW ordering of the two faces
        
        size_t v20 = mesh.get_third_vertex(e, regionfaces[0]);
        size_t v21 = mesh.get_third_vertex(e, regionfaces[1]);
        
        Vec3d x20 = m_surf.get_position(v20);
        Vec3d x21 = m_surf.get_position(v21);
        
        Vec2d p20 = Vec2d(dot(x20 - x0, u), dot(x20 - x0, v));
        Vec2d p21 = Vec2d(dot(x21 - x0, u), dot(x21 - x0, v));
        
        double theta0 = atan2(p20[1], p20[0]);
        double theta1 = atan2(p21[1], p21[0]);
        
        while (theta0 > theta1)
            theta1 += M_PI * 2;
        
        regionangles.push_back(theta1 - theta0);
    }
    
    // pick the pair of opposing regions with larger summed angles
    Vec2i pair0 = regions[0] < regions[2] ? Vec2i(regions[0], regions[2]) : Vec2i(regions[2], regions[0]); // order the two regions with ascending region label
    Vec2i pair1 = regions[1] < regions[3] ? Vec2i(regions[1], regions[3]) : Vec2i(regions[3], regions[1]);  
    Mat2i ret;
    if (regionangles[0] + regionangles[2] > regionangles[1] + regionangles[3])
    {
        ret[0] = pair0; // cut regions
        ret[1] = pair1; // the other two regions
    } else
    {
        ret[0] = pair1;
        ret[1] = pair0;
    }
    
    // if the X junction should not be cut, return (-1, -1).
    return ret;
}

    
    
// --------------------------------------------------------
///
/// Decide whether to cut an X-junction vertex between two given regions
///
// --------------------------------------------------------

bool T1Transition::should_pull_vertex_apart(size_t xj, int A, int B, Vec3d & pull_apart_direction)
{
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;

    // compute surface tension pulling force to see if this vertex pair needs to be pulled open.
    std::vector<Vec3d> vertsA;
    std::vector<Vec3d> vertsB;
    for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[xj].size(); i++)
    {
        size_t triangle = mesh.m_vertex_to_triangle_map[xj][i];
        
        // find the edge in triangle triangle that's opposite to vertex xj
        size_t l = 0;
        size_t other_edge = static_cast<size_t>(~0);
        for (l = 0; l < 3; l++)
        {
            if (mesh.m_edges[mesh.m_triangle_to_edge_map[triangle][l]][0] != xj &&
                mesh.m_edges[mesh.m_triangle_to_edge_map[triangle][l]][1] != xj)
                other_edge = mesh.m_triangle_to_edge_map[triangle][l];
        }
        assert(other_edge < mesh.ne());
        size_t v0 = mesh.m_edges[other_edge][0];
        size_t v1 = mesh.m_edges[other_edge][1];
        
        Vec3d x0 = m_surf.get_position(v0);
        Vec3d x1 = m_surf.get_position(v1);
        
        Vec2i label = mesh.get_triangle_label(triangle);
        if (label[0] == A || label[1] == A)
        {
            vertsA.push_back(x0);
            vertsA.push_back(x1);
        } else if (label[0] == B || label[1] == B)
        {
            vertsB.push_back(x0);
            vertsB.push_back(x1);
        }
    }
    assert(vertsA.size() > 0);
    assert(vertsB.size() > 0);
    
    Vec3d centroidA(0, 0, 0);
    Vec3d centroidB(0, 0, 0);
    for (size_t i = 0; i < vertsA.size(); i++) 
        centroidA += vertsA[i];
    centroidA /= vertsA.size();
    for (size_t i = 0; i < vertsB.size(); i++) 
        centroidB += vertsB[i];
    centroidB /= vertsB.size();
    
    pull_apart_direction = (centroidA - centroidB);
    pull_apart_direction /= mag(pull_apart_direction);
    
    Vec3d xxj = m_surf.get_position(xj);
    Vec3d force_a(0);
    Vec3d force_b(0);
    
    std::vector<size_t> faces_to_delete;
    std::vector<Vec3st> faces_to_create;
    std::vector<Vec2i> face_labels_to_create;
    
    size_t a = mesh.nv() + 1;   // placeholder vertices
    size_t b = mesh.nv() + 2;
    triangulate_popped_vertex(xj, A, B, a, b, faces_to_delete, faces_to_create, face_labels_to_create);
    
    for (size_t i = 0; i < faces_to_create.size(); i++)
    {
        Vec3st t = sort_triangle(faces_to_create[i]);
        
        assert(t[0] < mesh.nv());
        assert(t[2] >= mesh.nv());
        if (t[1] < mesh.nv())
        {
            // a regular triangle containing only either a or b, not both -- this triangle is not degenerate, so we can use the area gradient
            Vec3d x0 = (t[0] > mesh.nv() ? xxj : m_surf.get_position(t[0]));
            Vec3d x1 = (t[1] > mesh.nv() ? xxj : m_surf.get_position(t[1]));
            
            Vec3d foot = dot(xxj - x0, x1 - x0) / dot(x1 - x0, x1 - x0) * (x1 - x0) + x0;
            Vec3d force = (foot - xxj); // the surface tension force, i.e. area gradient
            force /= mag(force);
            force *= mag(x1 - x0);
            
            if (t[2] == a)
                force_a += force;
            if (t[2] == b)
                force_b += force;
        } else
        {
            // a triangle that contains both a and b -- this triangle is degenerate. we must pretent a and b are moved apart infinitesimally along pull_apart_direction
            force_a += -pull_apart_direction * mag(cross(m_surf.get_position(t[0]) - xxj, pull_apart_direction));
            force_b += pull_apart_direction * mag(cross(m_surf.get_position(t[0]) - xxj, pull_apart_direction));
        }
    }
    
    double tensile_force = dot(force_a - force_b, pull_apart_direction);
    
    return tensile_force > 0;
}
    
void T1Transition::triangulate_popped_vertex(size_t xj, int A, int B, size_t a, size_t b, std::vector<size_t> & faces_to_delete, std::vector<Vec3st> & faces_to_create, std::vector<Vec2i> & face_labels_to_create)
{
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;
    
    for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[xj].size(); i++)
    {
        size_t triangle = mesh.m_vertex_to_triangle_map[xj][i];
        
        faces_to_delete.push_back(triangle);
        
        // find the edge in triangle triangle that's opposite to vertex xj
        size_t l = 0;
        size_t edge2 = static_cast<size_t>(~0);
        for (l = 0; l < 3; l++)
        {
            size_t e = mesh.m_triangle_to_edge_map[triangle][l];
            if (mesh.m_edges[e][0] != xj && mesh.m_edges[e][1] != xj)
                edge2 = e;
        }
        assert(edge2 < mesh.ne());
        size_t v0 = mesh.m_edges[edge2][0];
        size_t v1 = mesh.m_edges[edge2][1];
        
        if (!mesh.oriented(v0, v1, mesh.get_triangle(triangle)))
            std::swap(v0, v1);
        
        Vec2i label = mesh.get_triangle_label(triangle);
        if (label[0] == A || label[1] == A)
        {
            faces_to_create.push_back(Vec3st(a, v0, v1));
            face_labels_to_create.push_back(label);
        } else
        {
            faces_to_create.push_back(Vec3st(b, v0, v1));
            face_labels_to_create.push_back(label);                
        }
    }
    
    std::vector<size_t> A_edges;
    for (size_t j = 0; j < mesh.m_vertex_to_edge_map[xj].size(); j++)
    {
        size_t edge = mesh.m_vertex_to_edge_map[xj][j];
        for (size_t k = 0; k < mesh.m_edge_to_triangle_map[edge].size(); k++)
        {
            Vec2i label = mesh.get_triangle_label(mesh.m_edge_to_triangle_map[edge][k]);
            if (label[0] == A || label[1] == A)
            {
                A_edges.push_back(edge);
                break;
            }
        }
    }

    // sweep A region edges
    for (size_t i = 0; i < A_edges.size(); i++)
    {
        size_t edge = A_edges[i];
        size_t v2 = (mesh.m_edges[edge][0] == xj ? mesh.m_edges[edge][1] : mesh.m_edges[edge][0]);
        
        int upper_region = -1;  // the region on the top when looking down the edge from xj to v2, with region B on the right
        int lower_region = -1;
        for (size_t j = 0; j < mesh.m_edge_to_triangle_map[edge].size(); j++)
        {
            size_t triangle = mesh.m_edge_to_triangle_map[edge][j];
            bool oriented = mesh.oriented(xj, v2, mesh.get_triangle(triangle));
            
            Vec2i label = mesh.get_triangle_label(triangle);
            
            if ((label[0] == A &&  oriented) ||
                (label[1] == A && !oriented))
            {
                upper_region = (label[0] == A ? label[1] : label[0]);
            }
            if ((label[0] == A && !oriented) ||
                (label[1] == A &&  oriented))
            {
                lower_region = (label[0] == A ? label[1] : label[0]);
            }
        }
        
        if (upper_region >= 0 && lower_region >= 0) // if this is not true, then the neighborhood around this edge is not complete, which can oly happen on the boundary.
        {
            if (upper_region == lower_region)
            {
                // this means either this edge is just a manifold edge between region A faces (thus pulling apart xj doesn't affect this edge), or it is an X junction edge with the same region above and below (pulling this edge apart creates a tunnel connecting them)
                // in either case, no triangle should be created.
            } else
            {
                faces_to_create.push_back(Vec3st(a, b, v2));
                face_labels_to_create.push_back(Vec2i(lower_region, upper_region));
            }
        }
    }

}

bool T1Transition::pulling_vertex_apart_introduces_collision(size_t v, const Vec3d & oldpos, const Vec3d & newpos0, const Vec3d & newpos1)
{
    bool collision0 = vertex_pseudo_motion_introduces_collision(v, oldpos, newpos0);
    bool collision1 = vertex_pseudo_motion_introduces_collision(v, oldpos, newpos1);
    
    return collision0 || collision1;
}
    
bool T1Transition::vertex_pseudo_motion_introduces_collision(size_t v, const Vec3d & oldpos, const Vec3d & newpos)
{
    // code adapted from EdgeSplitter::split_edge_pseudo_motion_introduces_intersection()
    
    if (!m_surf.m_collision_safety)
        return false;
    
    NonDestructiveTriMesh & m_mesh = m_surf.m_mesh;
    
    const std::vector<Vec3d> & x = m_surf.get_positions();
    
    std::vector<size_t> & tris = m_surf.m_mesh.m_vertex_to_triangle_map[v];
    std::vector<size_t> & edges = m_surf.m_mesh.m_vertex_to_edge_map[v];
    std::vector<size_t> edge_other_endpoints(edges.size());
    
    for (size_t i = 0; i < edges.size(); i++)
        edge_other_endpoints[i] = (m_mesh.m_edges[edges[i]][0] == v ? m_mesh.m_edges[edges[i]][1] : m_mesh.m_edges[edges[i]][0]);
    
    // new point vs all triangles
    {
        
        Vec3d aabb_low, aabb_high;
        minmax(oldpos, newpos, aabb_low, aabb_high);
        
        aabb_low  -= m_surf.m_aabb_padding * Vec3d(1,1,1);
        aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);
        
        std::vector<size_t> overlapping_triangles;
        m_surf.m_broad_phase->get_potential_triangle_collisions(aabb_low, aabb_high, true, true, overlapping_triangles);
        
        for (size_t i = 0; i < overlapping_triangles.size(); i++)
        {
            // exclude incident triangles
            if (m_mesh.get_triangle(overlapping_triangles[i])[0] == v ||
                m_mesh.get_triangle(overlapping_triangles[i])[1] == v ||
                m_mesh.get_triangle(overlapping_triangles[i])[2] == v)
                continue;
            
            Vec3st sorted_triangle = sort_triangle(m_mesh.get_triangle(overlapping_triangles[i]));
            size_t a = sorted_triangle[0];
            size_t b = sorted_triangle[1];
            size_t c = sorted_triangle[2];
            
            double t_zero_distance;
            check_point_triangle_proximity(oldpos, x[a], x[b], x[c], t_zero_distance);
            if (t_zero_distance < m_surf.m_improve_collision_epsilon)
                return true;
            
            if (point_triangle_collision(oldpos, newpos, v, x[a], x[a], a, x[b], x[b], b, x[c], x[c], c))
            {
//                if (m_surf.m_verbose)
                    std::cout << "Popping collision: point triangle: with triangle " << overlapping_triangles[i] << std::endl;
                return true;
            }
        }
        
    }
    
    // new edges vs all edges
    {
        Vec3d edge_aabb_low, edge_aabb_high;
        
        // do one big query into the broad phase for all new edges
        minmax(oldpos, newpos, edge_aabb_low, edge_aabb_high);
        for (size_t i = 0; i < edge_other_endpoints.size(); ++i)
            update_minmax(m_surf.get_position(edge_other_endpoints[i]), edge_aabb_low, edge_aabb_high);
        
        edge_aabb_low  -= m_surf.m_aabb_padding * Vec3d(1,1,1);
        edge_aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);
        
        std::vector<size_t> overlapping_edges;
        m_surf.m_broad_phase->get_potential_edge_collisions(edge_aabb_low, edge_aabb_high, true, true, overlapping_edges);
        
        for (size_t i = 0; i < overlapping_edges.size(); i++)
        {
            if (m_mesh.m_edges[overlapping_edges[i]][0] == m_mesh.m_edges[overlapping_edges[i]][1])
                continue;
            
            for (size_t j = 0; j < edges.size(); j++)
            {
                // exclude adjacent edges
                if (m_mesh.get_common_vertex(edges[j], overlapping_edges[i]) < m_mesh.nv())
                    continue;
                
                size_t n = edge_other_endpoints[j];
                size_t e0 = m_mesh.m_edges[overlapping_edges[i]][0];
                size_t e1 = m_mesh.m_edges[overlapping_edges[i]][1];
                if (e0 > e1)
                    std::swap(e0, e1);
                
                double t_zero_distance;
                check_edge_edge_proximity(oldpos, x[n], x[e0], x[e1], t_zero_distance);
                if (t_zero_distance < m_surf.m_improve_collision_epsilon)
                    return true;

                bool collision = (n < v ?
                                  segment_segment_collision(x[n], x[n], n, oldpos, newpos, v, x[e0], x[e0], e0, x[e1], x[e1], e1) :
                                  segment_segment_collision(oldpos, newpos, v, x[n], x[n], n, x[e0], x[e0], e0, x[e1], x[e1], e1));
                
                if (collision)
                {
//                    if (m_surf.m_verbose)
                        std::cout << "Popping collision: edge edge: edge other vertex = " << edge_other_endpoints[j] << " edge = " << overlapping_edges[i] << std::endl;
                    return true;
                }
            }
        }      
    }
    
    // new triangles vs all points
    {
        Vec3d triangle_aabb_low, triangle_aabb_high;
        
        // do one big query into the broad phase for all new triangles
        minmax(oldpos, newpos, triangle_aabb_low, triangle_aabb_high);
        for (size_t i = 0; i < edge_other_endpoints.size(); ++i)
            update_minmax(m_surf.get_position(edge_other_endpoints[i]), triangle_aabb_low, triangle_aabb_high);
        
        triangle_aabb_low  -= m_surf.m_aabb_padding * Vec3d(1,1,1);
        triangle_aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);
        
        std::vector<size_t> overlapping_vertices;
        m_surf.m_broad_phase->get_potential_vertex_collisions(triangle_aabb_low, triangle_aabb_high, true, true, overlapping_vertices);
        
        for (size_t i = 0; i < overlapping_vertices.size(); i++)
        {
            if (m_mesh.m_vertex_to_triangle_map[overlapping_vertices[i]].empty()) 
                continue; 
            
            const Vec3d & vert = m_surf.get_position(overlapping_vertices[i]);
            
            for (size_t j = 0; j < tris.size(); j++)
            {
                // exclude incident triangles
                if (m_mesh.get_triangle(tris[j])[0] == overlapping_vertices[i] ||
                    m_mesh.get_triangle(tris[j])[1] == overlapping_vertices[i] ||
                    m_mesh.get_triangle(tris[j])[2] == overlapping_vertices[i])
                    continue;
                
                Vec3st sorted_triangle = sort_triangle(m_mesh.get_triangle(tris[j]));
                size_t a = sorted_triangle[0];
                size_t b = sorted_triangle[1];
                size_t c = sorted_triangle[2];                
                
                Vec3d oldxa = (a == v ? oldpos : x[a]);
                Vec3d newxa = (a == v ? newpos : x[a]);
                Vec3d oldxb = (b == v ? oldpos : x[b]);
                Vec3d newxb = (b == v ? newpos : x[b]);
                Vec3d oldxc = (c == v ? oldpos : x[c]);
                Vec3d newxc = (c == v ? newpos : x[c]);
                
                double t_zero_distance;
                check_point_triangle_proximity(vert, oldxa, oldxb, oldxc, t_zero_distance);
                if (t_zero_distance < m_surf.m_improve_collision_epsilon)
                    return true;
                
                if (point_triangle_collision(vert, vert, overlapping_vertices[i], oldxa, newxa, a, oldxb, newxb, b, oldxc, newxc, c))
                {
//                    if (m_surf.m_verbose)
                        std::cout << "Popping collision: triangle point: with triangle " << tris[j] << " with vertex " << overlapping_vertices[i] << std::endl;
                    return true;
                }
            }
        }
    }
    
    return false;
}
    
bool T1Transition::vertex_pseudo_motion_introduces_collision(size_t v, const Vec3d & oldpos, const Vec3d & newpos, const std::vector<size_t> & tris, const std::vector<size_t> & edges)
{
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;
    
    const std::vector<Vec3d> & x = m_surf.get_positions();
    
    // sanity check: all tris triangles and all edges edges are actually incident to v
    for (size_t i = 0; i < tris.size(); i++)
    {
        const Vec3st & t = mesh.get_triangle(tris[i]);
        assert(t[0] == v || t[1] == v || t[2] == v);
    }
    
    for (size_t i = 0; i < edges.size(); i++)
    {
        const Vec2st & e = mesh.m_edges[edges[i]];
        assert(e[0] == v || e[1] == v);
    }
    
    // new point vs all triangles
    {
        Vec3d aabb_low, aabb_high;
        minmax(oldpos, newpos, aabb_low, aabb_high);
        
        aabb_low  -= m_surf.m_aabb_padding * Vec3d(1,1,1);
        aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);
        
        std::vector<size_t> overlapping_triangles;
        m_surf.m_broad_phase->get_potential_triangle_collisions(aabb_low, aabb_high, true, true, overlapping_triangles);
        
        for (size_t i = 0; i < overlapping_triangles.size(); i++)
        {
            const Vec3st & t = mesh.get_triangle(overlapping_triangles[i]);
            
            // exclude incident triangles
            if (t[0] == v || t[1] == v || t[2] == v)
                continue;
            
            Vec3st sorted_triangle = sort_triangle(t);
            size_t a = sorted_triangle[0];
            size_t b = sorted_triangle[1];
            size_t c = sorted_triangle[2];
            
            double t_zero_distance;
            check_point_triangle_proximity(oldpos, x[a], x[b], x[c], t_zero_distance);
            if (t_zero_distance < m_surf.m_improve_collision_epsilon)
                return true;
            
            if (point_triangle_collision(oldpos, newpos, v, x[a], x[a], a, x[b], x[b], b, x[c], x[c], c))
            {
//                if (m_surf.m_verbose)
                    std::cout << "Popping collision: point triangle: with triangle " << overlapping_triangles[i] << std::endl;
                return true;
            }
        }
        
    }
    
    // new edges vs all edges
    {
        Vec3d edge_aabb_low, edge_aabb_high;
        
        // do one big query into the broad phase for all new edges
        minmax(oldpos, newpos, edge_aabb_low, edge_aabb_high);
        for (size_t i = 0; i < edges.size(); ++i)
        {
            size_t other_endpoint = (mesh.m_edges[edges[i]][0] == v ? mesh.m_edges[edges[i]][1] : mesh.m_edges[edges[i]][0]);
            update_minmax(m_surf.get_position(other_endpoint), edge_aabb_low, edge_aabb_high);
        }
        
        edge_aabb_low  -= m_surf.m_aabb_padding * Vec3d(1,1,1);
        edge_aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);
        
        std::vector<size_t> overlapping_edges;
        m_surf.m_broad_phase->get_potential_edge_collisions(edge_aabb_low, edge_aabb_high, true, true, overlapping_edges);
        
        for (size_t i = 0; i < overlapping_edges.size(); i++)
        {
            const Vec2st & e = mesh.m_edges[overlapping_edges[i]];
            
            if (e[0] == e[1])
                continue;
            
            for (size_t j = 0; j < edges.size(); j++)
            {
                // exclude adjacent edges
                if (mesh.get_common_vertex(edges[j], overlapping_edges[i]) < mesh.nv())
                    continue;
                
                size_t n = (mesh.m_edges[edges[j]][0] == v ? mesh.m_edges[edges[j]][1] : mesh.m_edges[edges[j]][0]);
                size_t e0 = e[0];
                size_t e1 = e[1];
                if (e0 > e1)
                    std::swap(e0, e1);
                
                double t_zero_distance;
                check_edge_edge_proximity(oldpos, x[n], x[e0], x[e1], t_zero_distance);
                if (t_zero_distance < m_surf.m_improve_collision_epsilon)
                    return true;
                
                bool collision = (n < v ?
                                  segment_segment_collision(x[n], x[n], n, oldpos, newpos, v, x[e0], x[e0], e0, x[e1], x[e1], e1) :
                                  segment_segment_collision(oldpos, newpos, v, x[n], x[n], n, x[e0], x[e0], e0, x[e1], x[e1], e1));
                
                if (collision)
                {
//                    if (m_surf.m_verbose)
                        std::cout << "Popping collision: edge edge: edge other vertex = " << n << " edge = " << overlapping_edges[i] << std::endl;
                    return true;
                }
            }
        }      
    }
    
    // new triangles vs all points
    {
        Vec3d triangle_aabb_low, triangle_aabb_high;
        
        // do one big query into the broad phase for all new triangles
        minmax(oldpos, newpos, triangle_aabb_low, triangle_aabb_high);
        for (size_t i = 0; i < tris.size(); ++i)
        {
            const Vec3st & t = mesh.get_triangle(tris[i]);
            Vec2st other_vertices;
            if (t[0] == v)
                other_vertices = Vec2st(t[1], t[2]);
            else if (t[1] == v)
                other_vertices = Vec2st(t[2], t[0]);
            else if (t[2] == v)
                other_vertices = Vec2st(t[0], t[1]);
            else
                assert(!"triangle in tris does not contain vertex v");
            update_minmax(m_surf.get_position(other_vertices[0]), triangle_aabb_low, triangle_aabb_high);
            update_minmax(m_surf.get_position(other_vertices[1]), triangle_aabb_low, triangle_aabb_high);
        }
        
        triangle_aabb_low  -= m_surf.m_aabb_padding * Vec3d(1,1,1);
        triangle_aabb_high += m_surf.m_aabb_padding * Vec3d(1,1,1);
        
        std::vector<size_t> overlapping_vertices;
        m_surf.m_broad_phase->get_potential_vertex_collisions(triangle_aabb_low, triangle_aabb_high, true, true, overlapping_vertices);
        
        for (size_t i = 0; i < overlapping_vertices.size(); i++)
        {
            size_t & ov = overlapping_vertices[i];
            
            if (mesh.m_vertex_to_triangle_map[ov].empty()) 
                continue; 
            
            const Vec3d & vert = m_surf.get_position(ov);
            
            for (size_t j = 0; j < tris.size(); j++)
            {
                const Vec3st & t = mesh.get_triangle(tris[j]);
                
                // exclude incident triangles
                if (t[0] == ov || t[1] == ov || t[2] == ov)
                    continue;
                
                Vec3st sorted_triangle = sort_triangle(t);
                size_t a = sorted_triangle[0];
                size_t b = sorted_triangle[1];
                size_t c = sorted_triangle[2];
                
                Vec3d oldxa = (a == v ? oldpos : x[a]);
                Vec3d newxa = (a == v ? newpos : x[a]);
                Vec3d oldxb = (b == v ? oldpos : x[b]);
                Vec3d newxb = (b == v ? newpos : x[b]);
                Vec3d oldxc = (c == v ? oldpos : x[c]);
                Vec3d newxc = (c == v ? newpos : x[c]);
                
                double t_zero_distance;
                check_point_triangle_proximity(vert, oldxa, oldxb, oldxc, t_zero_distance);
                if (t_zero_distance < m_surf.m_improve_collision_epsilon)
                    return true;
                
                if (point_triangle_collision(vert, vert, overlapping_vertices[i], oldxa, newxa, a, oldxb, newxb, b, oldxc, newxc, c))
                {
//                    if (m_surf.m_verbose)
                        std::cout << "Popping collision: triangle point: with triangle " << tris[j] << " with vertex " << ov << std::endl;
                    return true;
                }
            }
        }
    }
    
    return false;
    
    
    
}

    

}
