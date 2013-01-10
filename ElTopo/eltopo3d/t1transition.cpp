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
    
    // sort the X-junction edges into connected groups with the same cut
    std::vector<std::pair<std::vector<size_t>, Mat2i> > xjgroups;   // each element is a list of consecutive x junction edges, along with the id of the two regions that should end up adjacent after pull-apart, and the id of the other two regions

    for (size_t i = 0; i < xjunctions.size(); i++)
    {
        Mat2i cut = cut_x_junction_edge(xjunctions[i]);
        std::vector<size_t> found_groups;
        for (size_t j = 0; j < xjgroups.size(); j++)
        {
            if (cut != xjgroups[j].second)
                continue;
            
            for (size_t k = 0; k < xjgroups[j].first.size(); k++)
            {
                if (mesh.get_common_vertex(xjunctions[i], xjgroups[j].first[k]) < mesh.nv())
                {
                    found_groups.push_back(j);
                    break;
                }
            }
        }
        
        if (found_groups.size() == 0)
        {
            // new group
            xjgroups.push_back(std::pair<std::vector<size_t>, Mat2i>(std::vector<size_t>(1, xjunctions[i]), cut));
        } else
        {
            // joining all the groups in found_groups together
            std::vector<size_t> newgroup;
            for (size_t j = 0; j < found_groups.size(); j++)
                newgroup.insert(newgroup.end(), xjgroups[found_groups[j]].first.begin(), xjgroups[found_groups[j]].first.end());
            newgroup.push_back(xjunctions[i]);
            
            std::sort(found_groups.rbegin(), found_groups.rend());
            for (size_t j = 0; j < found_groups.size(); j++)
                xjgroups.erase(xjgroups.begin() + found_groups[j]);
            
            xjgroups.push_back(std::pair<std::vector<size_t>, Mat2i>(newgroup, cut));
        }
    }
    
    // process the groups one after another
    for (size_t i = 0; i < xjgroups.size(); i++)
    {
        // TODO: merge the following two cases into one codepath
        if (xjgroups[i].first.size() == 1)
        {
            // if a group has only one edge, split the edge in half and pull the midpoint vertex apart.
            size_t edge = xjgroups[i].first.front();
            Mat2i cut_regions = xjgroups[i].second;
            Vec2i cut = cut_regions[0];
            
            size_t v0 = mesh.m_edges[edge][0];
            size_t v1 = mesh.m_edges[edge][1];
            
            size_t nv0 = m_surf.add_vertex(m_surf.get_position(v0), m_surf.m_masses[v0]);
            size_t nv1 = m_surf.add_vertex(m_surf.get_position(v0), m_surf.m_masses[v0]);
            
            m_surf.set_newposition(nv0, (m_surf.get_position(v0) + m_surf.get_position(v1)) / 2);  // the two new vertices will be pulled apart later
            m_surf.set_newposition(nv1, (m_surf.get_position(v0) + m_surf.get_position(v1)) / 2);
            m_surf.set_remesh_velocity(nv0, (m_surf.get_remesh_velocity(v0) + m_surf.get_remesh_velocity(v1)) / 2);
            m_surf.set_remesh_velocity(nv1, (m_surf.get_remesh_velocity(v0) + m_surf.get_remesh_velocity(v1)) / 2);
            
            int upper_region = -1;  // upper means the region is on the top when looking down the edge "edge", with region cut.x() on the left and cut.y() on the right.
            int lower_region = -1;
            std::vector<size_t> region0faces;
            for (size_t j = 0; j < mesh.m_edge_to_triangle_map[edge].size(); j++)
                if (mesh.get_triangle_label(mesh.m_edge_to_triangle_map[edge][j])[0] == cut[0] ||
                    mesh.get_triangle_label(mesh.m_edge_to_triangle_map[edge][j])[1] == cut[0])
                    region0faces.push_back(mesh.m_edge_to_triangle_map[edge][j]);
            assert(region0faces.size() == 2);
            
            Vec2i label0 = mesh.get_triangle_label(region0faces[0]);
            upper_region = (label0[0] == cut[0] ? label0[1] : label0[0]);
            Vec2i label1 = mesh.get_triangle_label(region0faces[1]);
            lower_region = (label1[0] == cut[0] ? label1[1] : label1[0]);
            
            if (( mesh.oriented(v0, v1, mesh.get_triangle(region0faces[0])) && label0[1] == cut[0]) ||
                (!mesh.oriented(v0, v1, mesh.get_triangle(region0faces[0])) && label0[0] == cut[0]))
            {
                std::swap(upper_region, lower_region);
                assert(( mesh.oriented(v0, v1, mesh.get_triangle(region0faces[1])) && label1[0] == cut[0]) ||
                       (!mesh.oriented(v0, v1, mesh.get_triangle(region0faces[1])) && label1[1] == cut[0]));
            }                
            assert(upper_region >= 0);
            assert(lower_region >= 0);
            
            std::vector<size_t> faces_to_delete;
            std::vector<Vec3st> faces_to_create;
            std::vector<Vec2i> face_labels_to_create;
            std::vector<size_t> faces_created;
            
            std::vector<size_t> verts_to_delete;
            std::vector<Vec3d> verts_to_create;
            std::vector<size_t> verts_created;
            
            std::vector<Vec3d> upper_neighbors;
            std::vector<Vec3d> lower_neighbors;
            
            Vec3d upper_desired_position;
            Vec3d lower_desired_position;
            
            for (size_t j = 0; j < mesh.m_edge_to_triangle_map[edge].size(); j++)
            {
                size_t triangle = mesh.m_edge_to_triangle_map[edge][j];
                faces_to_delete.push_back(triangle);
                
                size_t v2 = mesh.get_third_vertex(edge, triangle);
                
                Vec2i label = mesh.get_triangle_label(triangle);
                assert(label[0] == upper_region || label[1] == upper_region || label[0] == lower_region || label[1] == lower_region);
                size_t nv = ((label[0] == upper_region || label[1] == upper_region) ? nv0 : nv1); // nv0 for upper region faces, nv1 for lower region faces
                
                if (mesh.oriented(v0, v1, mesh.get_triangle(triangle)))
                {
                    faces_to_create.push_back(Vec3st(v0, nv, v2));
                    faces_to_create.push_back(Vec3st(nv, v1, v2));
                } else {
                    faces_to_create.push_back(Vec3st(nv, v0, v2));
                    faces_to_create.push_back(Vec3st(v1, nv, v2));
                }
                face_labels_to_create.push_back(label);
                face_labels_to_create.push_back(label);
                
                ((label[0] == upper_region || label[1] == upper_region) ? upper_neighbors : lower_neighbors).push_back(m_surf.get_position(v2));
            }
            
            // walls
            faces_to_create.push_back(Vec3st(v0, nv1, nv0));
            faces_to_create.push_back(Vec3st(nv0, nv1, v1));
            
            face_labels_to_create.push_back(cut);
            face_labels_to_create.push_back(cut);
            
            // compute desired positions for new vertices
            Vec3d upper_neighbors_mean(0, 0, 0);
            for (size_t j = 0; j < upper_neighbors.size(); j++)
                upper_neighbors_mean += upper_neighbors[j];
            upper_neighbors_mean /= upper_neighbors.size();
            
            Vec3d lower_neighbors_mean(0, 0, 0);
            for (size_t j = 0; j < lower_neighbors.size(); j++)
                lower_neighbors_mean += lower_neighbors[j];
            lower_neighbors_mean /= lower_neighbors.size();
            
            Vec3d wall_breadth = (upper_neighbors_mean - lower_neighbors_mean);
            wall_breadth /= mag(wall_breadth);
            wall_breadth *= mag(m_surf.get_position(v1) - m_surf.get_position(v0));

            // compute the desired destination positions, enforcing constraints
            bool original_constraint = mesh.get_vertex_constraint_label(v0) && mesh.get_vertex_constraint_label(v1);
            upper_desired_position = m_surf.get_newposition(nv0) + wall_breadth * 0.1;
            lower_desired_position = m_surf.get_newposition(nv0) - wall_breadth * 0.1;
            
            if (original_constraint)
            {
                assert(m_surf.m_constrained_vertices_callback);
                m_surf.m_constrained_vertices_callback->generate_edge_popped_positions(m_surf, v0, cut, upper_desired_position, lower_desired_position);
            }
            
            mesh.set_vertex_constraint_label(nv0, original_constraint);
            mesh.set_vertex_constraint_label(nv1, original_constraint);
            
            // apply the deletion
            for (size_t j = 0; j < faces_to_delete.size(); j++)
                m_surf.remove_triangle(faces_to_delete[j]);    //&&&& recursive deletion
            
//            m_obj->deleteEdge(edge, false);   //&&&&
            
//            // triangulate with the new vertices  //&&&&
//            m_obj->addEdge(v0, nv0);
//            m_obj->addEdge(nv0, v1);
//            m_obj->addEdge(v0, nv1);
//            m_obj->addEdge(nv1, v1);
//            m_obj->addEdge(nv0, nv1);
            
            assert(faces_to_create.size() == face_labels_to_create.size());
            for (size_t j = 0; j < faces_to_create.size(); j++)
            {
                size_t nf = m_surf.add_triangle(faces_to_create[j]);
                mesh.set_triangle_label(nf, face_labels_to_create[j]);
                faces_created.push_back(nf);
            }
            
            //pull the two new vertices (nv0, nv1) apart (before this the two vertices have the same position)
            m_surf.set_newposition(nv0, upper_desired_position);
            m_surf.set_newposition(nv1, lower_desired_position);
            
            m_surf.set_position(nv0, m_surf.get_newposition(nv0));
            m_surf.set_position(nv1, m_surf.get_newposition(nv1));
            
            // vertex deletion/creation logging
            verts_to_create.push_back(m_surf.get_newposition(nv0));
            verts_to_create.push_back(m_surf.get_newposition(nv1));
            
            verts_created.push_back(nv0);
            verts_created.push_back(nv1);
            
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
            
        } else
        {
            // if a group has at least two edges forming a polyline (this line should have no branching but could be a closed loop), pull the interior vertices apart.
            assert(xjgroups[i].first.size() > 1);
            
            Mat2i cut_regions = xjgroups[i].second;
            Vec2i cut = cut_regions[0];
            
            // sort the edges into an ordered polyline
            std::deque<size_t> ordered_edges;
            ordered_edges.push_back(xjgroups[i].first.back());
            xjgroups[i].first.pop_back();
            while (xjgroups[i].first.size() > 0)
            {
                size_t s = xjgroups[i].first.size();
                for (size_t j = 0; j < xjgroups[i].first.size(); j++)
                {
                    if (mesh.get_common_vertex(ordered_edges.back(), xjgroups[i].first[j]) < mesh.nv())
                    {
                        ordered_edges.push_back(xjgroups[i].first[j]);
                        xjgroups[i].first.erase(xjgroups[i].first.begin() + j);
                        break;
                    } else if (mesh.get_common_vertex(ordered_edges.front(), xjgroups[i].first[j]) < mesh.nv())
                    {
                        ordered_edges.push_front(xjgroups[i].first[j]);
                        xjgroups[i].first.erase(xjgroups[i].first.begin() + j);
                        break;
                    }
                }
                assert(xjgroups[i].first.size() < s);
            }
            size_t ne = ordered_edges.size();
            
            // TODO: handle the closed loop case. for now, assert the edges do not form a closed loop.
            if (ne == 2)
                assert(ordered_edges.front() != ordered_edges.back());
            else
                assert(mesh.get_common_vertex(ordered_edges.front(), ordered_edges.back()) >= mesh.nv());
            
            // duplicate the interior vertices, and compute their desired pull-apart positions
            std::vector<size_t> upper_junctions(ne + 1);
            std::vector<size_t> lower_junctions(ne + 1);
            std::vector<Vec3d> upper_junction_desired_positions(ne + 1);
            std::vector<Vec3d> lower_junction_desired_positions(ne + 1);
            
            int upper_region = -1;
            int lower_region = -1;
            
            std::vector<Vec3d> pull_apart_offsets(ne - 1);
            std::vector<bool> edge_oriented(ne);
            
            for (size_t j = 0; j < ne - 1; j++)
            {
                size_t edge0 = ordered_edges[j];
                size_t edge1 = ordered_edges[(j + 1) % ne];
                
                size_t v = mesh.get_common_vertex(edge0, edge1);
                assert(v < mesh.nv());
                
                size_t v0 = (v == mesh.m_edges[edge0][0] ? mesh.m_edges[edge0][1] : mesh.m_edges[edge0][0]);
                size_t v1 = (v == mesh.m_edges[edge1][0] ? mesh.m_edges[edge1][1] : mesh.m_edges[edge1][0]);
                
                edge_oriented[j]     = (v == mesh.m_edges[edge0][1]);
                edge_oriented[j + 1] = (v == mesh.m_edges[edge1][0]);
                
                size_t nv0 = m_surf.add_vertex(m_surf.get_position(v), m_surf.m_masses[v]);
                size_t nv1 = m_surf.add_vertex(m_surf.get_position(v), m_surf.m_masses[v]);
                
                m_surf.set_newposition(nv0, m_surf.get_position(v));  // the two new vertices will be pulled apart later
                m_surf.set_newposition(nv1, m_surf.get_position(v));
                m_surf.set_remesh_velocity(nv0, m_surf.get_remesh_velocity(v));
                m_surf.set_remesh_velocity(nv1, m_surf.get_remesh_velocity(v));
                
                upper_junctions[j + 1] = nv0;
                lower_junctions[j + 1] = nv1;
                if (j == 0)
                {
                    upper_junctions[0] = v0;
                    lower_junctions[0] = v0;
                }
                if (j == ne - 2)
                {
                    upper_junctions[ne] = v1;
                    lower_junctions[ne] = v1;
                }
                
                int edge_upper_region = -1;  // upper means the region is on the top when looking from vertex v0 to vertex v (may or may not be the direction of edge0), with region cut.x() on the left and cut.y() on the right.
                int edge_lower_region = -1;
                std::vector<size_t> region0faces;
                for (size_t k = 0; k < mesh.m_edge_to_triangle_map[edge0].size(); k++)
                    if (mesh.get_triangle_label(mesh.m_edge_to_triangle_map[edge0][k])[0] == cut[0] ||
                        mesh.get_triangle_label(mesh.m_edge_to_triangle_map[edge0][k])[1] == cut[0])
                        region0faces.push_back(mesh.m_edge_to_triangle_map[edge0][k]);
                assert(region0faces.size() == 2);
                
                Vec2i label0 = mesh.get_triangle_label(region0faces[0]);
                edge_upper_region = (label0[0] == cut[0] ? label0[1] : label0[0]);
                Vec2i label1 = mesh.get_triangle_label(region0faces[1]);
                edge_lower_region = (label1[0] == cut[0] ? label1[1] : label1[0]);
                
                if (( mesh.oriented(v0, v, mesh.get_triangle(region0faces[0])) && label0[1] == cut[0]) ||
                    (!mesh.oriented(v0, v, mesh.get_triangle(region0faces[0])) && label0[0] == cut[0]))
                {
                    std::swap(edge_upper_region, edge_lower_region);
                    assert(( mesh.oriented(v0, v, mesh.get_triangle(region0faces[1])) && label1[0] == cut[0]) ||
                           (!mesh.oriented(v0, v, mesh.get_triangle(region0faces[1])) && label1[1] == cut[0]));
                }
                assert(edge_upper_region >= 0);
                assert(edge_lower_region >= 0);
                
                assert(upper_region < 0 || upper_region == edge_upper_region);
                upper_region = edge_upper_region;
                assert(lower_region < 0 || lower_region == edge_lower_region);
                lower_region = edge_lower_region;
                
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
                
                pull_apart_offsets[j] = (upper_vertices_mean - lower_vertices_mean);
                pull_apart_offsets[j] /= mag(pull_apart_offsets[j]);
                pull_apart_offsets[j] *= mag(m_surf.get_position(v1) - m_surf.get_position(v0));
                
                // compute the desired destination positions for nv0 and nv1, enforcing constraints
                bool original_constraint = m_surf.m_mesh.get_vertex_constraint_label(v);
                upper_junction_desired_positions[j + 1] = m_surf.get_newposition(nv0) + pull_apart_offsets[j] * 0.1;
                lower_junction_desired_positions[j + 1] = m_surf.get_newposition(nv1) - pull_apart_offsets[j] * 0.1;
                
                if (original_constraint)
                {
                    assert(m_surf.m_constrained_vertices_callback);
                    m_surf.m_constrained_vertices_callback->generate_edge_popped_positions(m_surf, v, cut, upper_junction_desired_positions[j + 1], lower_junction_desired_positions[j + 1]);
                }
                
                mesh.set_vertex_constraint_label(nv0, original_constraint);
                mesh.set_vertex_constraint_label(nv1, original_constraint);

            }
            
            std::vector<size_t> faces_to_delete;
            std::vector<Vec3st> faces_to_create;
            std::vector<Vec2i> face_labels_to_create;
            std::vector<size_t> faces_created;
            
            std::vector<size_t> verts_to_delete;
            std::vector<Vec3d> verts_to_create;
            std::vector<size_t> verts_created;
            
            // update the faces incident to the X-junction edges
            for (size_t j = 0; j < ne; j++)
            {
                size_t edge = ordered_edges[j];
                //&&&&
//                edges_to_delete.push_back(edge);
//                edges_to_create.push_back(Vec2st(upper_junctions[j + 0], upper_junctions[j + 1]));
//                edges_to_create.push_back(Vec2st(lower_junctions[j + 0], lower_junctions[j + 1]));
                
                for (size_t k = 0; k < mesh.m_edge_to_triangle_map[edge].size(); k++)
                {
                    size_t triangle = mesh.m_edge_to_triangle_map[edge][k];
                    faces_to_delete.push_back(triangle);
                    
                    size_t v2 = mesh.get_third_vertex(edge, triangle);
                    
                    Vec2i label = mesh.get_triangle_label(triangle);
                    assert(label[0] == upper_region || label[1] == upper_region || label[0] == lower_region || label[1] == lower_region);
                    
                    size_t v0 = ((label[0] == upper_region || label[1] == upper_region) ? upper_junctions[j + 0] : lower_junctions[j + 0]);
                    size_t v1 = ((label[0] == upper_region || label[1] == upper_region) ? upper_junctions[j + 1] : lower_junctions[j + 1]);
                    
                    if (mesh.oriented(mesh.m_edges[edge][0], mesh.m_edges[edge][1], mesh.get_triangle(triangle)) == edge_oriented[j])
                        faces_to_create.push_back(Vec3st(v0, v1, v2));
                    else
                        faces_to_create.push_back(Vec3st(v1, v0, v2));
                    face_labels_to_create.push_back(label);
                    
                }
            }     
            
            // update the faces incident only to the interior vertices but not to any X-junction edge
            for (size_t j = 0; j < ne - 1; j++)
            {
                size_t edge0 = ordered_edges[j];
                size_t edge1 = ordered_edges[(j + 1) % ne];
                
                size_t v = mesh.get_common_vertex(edge0, edge1);
                assert(v < mesh.nv());
                
                verts_to_delete.push_back(v);
                
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
                    
                    size_t nv = ((label[0] == upper_region || label[1] == upper_region) ? upper_junctions[j + 1] : lower_junctions[j + 1]);
                    
                    if (mesh.oriented(v1, v2, mesh.get_triangle(triangle)))
                        faces_to_create.push_back(Vec3st(nv, v1, v2));
                    else
                        faces_to_create.push_back(Vec3st(nv, v2, v1));
                    face_labels_to_create.push_back(label);
                    
                }

                //&&&& recursive deletion
//                for (VertexEdgeIterator veit = m_obj->ve_iter(v); veit; ++veit)
//                {
//                    if (*veit == edge0 || *veit == edge1)
//                        continue;
//                    
//                    edges_to_delete.push_back(*veit);
//                }
                
            }
            
            // triangulate the new interface between cut.x() and cut.y()
            for (size_t j = 0; j < ne; j++)
            {
                if (j == 0)
                {
                    faces_to_create.push_back(Vec3st(upper_junctions[j], lower_junctions[j + 1], upper_junctions[j + 1]));
                } else if (j == ne - 1)
                {
                    faces_to_create.push_back(Vec3st(upper_junctions[j], lower_junctions[j], upper_junctions[j + 1]));
                } else
                {
                    faces_to_create.push_back(Vec3st(upper_junctions[j], lower_junctions[j], upper_junctions[j + 1]));
                    faces_to_create.push_back(Vec3st(lower_junctions[j + 1], upper_junctions[j + 1], lower_junctions[j]));
                }
                face_labels_to_create.push_back(cut);
                if (j != 0 && j != ne - 1)
                    face_labels_to_create.push_back(cut);
            }

            // apply the deletion
            for (size_t j = 0; j < faces_to_delete.size(); j++)
                m_surf.remove_triangle(faces_to_delete[j]);    //&&&& recursive deletion
  
            //&&&&
//            for (size_t j = 0; j < edges_to_delete.size(); j++)
//                m_obj->deleteEdge(edges_to_delete[j], false);
//            
//            for (size_t j = 0; j < vertices_to_delete.size(); j++)
//                m_obj->deleteVertex(vertices_to_delete[j]);
            
            // triangulate with the new vertices
            //&&&&
//            for (size_t j = 0; j < edges_to_create.size(); j++)
//                m_obj->addEdge(edges_to_create[j].x(), edges_to_create[j].y());
            
            assert(faces_to_create.size() == face_labels_to_create.size());
            for (size_t j = 0; j < faces_to_create.size(); j++)
            {
                size_t nf = m_surf.add_triangle(faces_to_create[j]);
                mesh.set_triangle_label(nf, face_labels_to_create[j]);
                faces_created.push_back(nf);
            }
            
            // pull the two new vertices (nv0, nv1) apart (before this the two vertices have the same position)
            for (size_t j = 0; j < ne - 1; j++)
            {
                size_t edge0 = ordered_edges[j];
                size_t edge1 = ordered_edges[(j + 1) % ne];
                
                size_t v = mesh.get_common_vertex(edge0, edge1);
                assert(v < mesh.nv());

                size_t nv0 = upper_junctions[j + 1];
                size_t nv1 = lower_junctions[j + 1];
                
                m_surf.set_newposition(nv0, upper_junction_desired_positions[j + 1]);
                m_surf.set_newposition(nv1, lower_junction_desired_positions[j + 1]);
                
                m_surf.set_position(nv0, m_surf.get_newposition(nv0));
                m_surf.set_position(nv1, m_surf.get_newposition(nv1));                
            }
            
            // vertex deletion/creation logging
            for (size_t j = 0; j < ne - 1; j++)
            {
                size_t nv0 = upper_junctions[j + 1];
                size_t nv1 = lower_junctions[j + 1];
                
                verts_to_create.push_back(m_surf.get_newposition(nv0));
                verts_to_create.push_back(m_surf.get_newposition(nv1));
            
                verts_created.push_back(nv0);
                verts_created.push_back(nv1);
            }
            
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
        //    Otherwise, if it is touching region A through an edge, it will be updated into a quad spanning both vertex a and b, otherwise, update it to use vertex b
        //  Push vertex a and b on the stack to be visited next.
        //
        // Notes: 
        //  The algorithm above assumes that region A and B must always be connected by triangle fans, with one end of the fan touching region
        //  A and the other touching B. Pulling apart the vertex makes this a stripe, thus at least one of the triangles in the fan need to be
        //  turned into a quad. The algorithm picks the head of the fan triangle sequence (the one touching A) for this task. In the general
        //  setting, there is one more possibility: region A and B are adjacent through only an edge. However this is not possible in our
        //  workflow because performT1Transition() should have been called resolving all X-junction edges. (TODO: what if performT1Transition()
        //  fails due to collision?)
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
        for (size_t i = 0; i < vertex_regions.size(); i++)
        {
            for (size_t j = i + 1; j < vertex_regions.size(); j++)
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
        
        // pull apart
        size_t a = m_surf.add_vertex(m_surf.get_position(xj), m_surf.m_masses[xj]);
        size_t b = m_surf.add_vertex(m_surf.get_position(xj), m_surf.m_masses[xj]);
        
        // set the position/velocity of new vertices
        m_surf.set_newposition(a, m_surf.get_position(xj));
        m_surf.set_newposition(b, m_surf.get_position(xj));
        m_surf.set_remesh_velocity(a, m_surf.get_remesh_velocity(xj));
        m_surf.set_remesh_velocity(b, m_surf.get_remesh_velocity(xj));
        
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
        bool xj_constraint = m_surf.m_mesh.get_vertex_constraint_label(xj);
        Vec3d a_desired_position = m_surf.get_newposition(a) + pull_apart_offset * 0.1;
        Vec3d b_desired_position = m_surf.get_newposition(b) - pull_apart_offset * 0.1;

        if (xj_constraint)
        {
            assert(m_surf.m_constrained_vertices_callback);            
            m_surf.m_constrained_vertices_callback->generate_vertex_popped_positions(m_surf, xj, A, B, a_desired_position, b_desired_position);
        }
        
        mesh.set_vertex_constraint_label(a, xj_constraint);
        mesh.set_vertex_constraint_label(b, xj_constraint);
        
        // update the face connectivities
        std::vector<size_t> faces_to_delete;
        std::vector<Vec3st> faces_to_create;
        std::vector<Vec2i> face_labels_to_create;
        std::vector<size_t> faces_created;
        
        std::vector<size_t> verts_to_delete;
        std::vector<Vec3d> verts_to_create;
        std::vector<size_t> verts_created;
        
        //&&&&
//        std::vector<size_t> edges_to_delete;
        
        for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[xj].size(); i++)
        {
            size_t triangle = mesh.m_vertex_to_triangle_map[xj][i];
            
            faces_to_delete.push_back(triangle);
            
            // find the edge in triangle triangle that's opposite to vertex xj
            size_t l = 0;
            size_t edge0 = static_cast<size_t>(~0);
            size_t edge1 = static_cast<size_t>(~0);
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
            
            for (l = 0; l < 3; l++)
            {
                size_t e = mesh.m_triangle_to_edge_map[triangle][l];
                if ((mesh.m_edges[e][0] == xj && mesh.m_edges[e][1] == v0) ||
                    (mesh.m_edges[e][1] == xj && mesh.m_edges[e][0] == v0))
                    edge0 = e;
                else if ((mesh.m_edges[e][0] == xj && mesh.m_edges[e][1] == v1) ||
                         (mesh.m_edges[e][1] == xj && mesh.m_edges[e][0] == v1))
                    edge1 = e;
            }
            assert(edge0 < mesh.ne());
            assert(edge1 < mesh.ne());
            
            assert(v0 == mesh.m_edges[edge0][0] || v0 == mesh.m_edges[edge0][1]);
            assert(v1 == mesh.m_edges[edge1][0] || v1 == mesh.m_edges[edge1][1]);
            
            if (!mesh.oriented(v0, v1, mesh.get_triangle(triangle)))
            {
                std::swap(v0, v1);
                std::swap(edge0, edge1);
            }
            
            Vec2i label = mesh.get_triangle_label(triangle);
            if (label[0] == A || label[1] == A)
            {
                faces_to_create.push_back(Vec3st(a, v0, v1));
                face_labels_to_create.push_back(label);
                
            } else if (label[0] == B || label[1] == B)
            {
                faces_to_create.push_back(Vec3st(b, v0, v1));
                face_labels_to_create.push_back(label);
                
            } else
            {
                // find out which one out of xj-v0 and xj-v1 is next to region A and to region B (or neither), in order to determine triangulation
                size_t ev0 = edge0;
                size_t ev1 = edge1;
                
                bool v0adjA = false;
                bool v1adjA = false;
                for (size_t j = 0; j < mesh.m_edge_to_triangle_map[ev0].size(); j++)
                {
                    Vec2i label = mesh.get_triangle_label(mesh.m_edge_to_triangle_map[ev0][j]);
                    if (label[0] == A || label[1] == A)
                        v0adjA = true;
                }
                for (size_t j = 0; j < mesh.m_edge_to_triangle_map[ev1].size(); j++)
                {
                    Vec2i label = mesh.get_triangle_label(mesh.m_edge_to_triangle_map[ev1][j]);
                    if (label[0] == A || label[1] == A)
                        v1adjA = true;
                }
                
                assert(!v0adjA || !v1adjA);
                
                if (v0adjA)
                {
                    faces_to_create.push_back(Vec3st(a, v0, v1));
                    faces_to_create.push_back(Vec3st(a, v1, b));
                    face_labels_to_create.push_back(label);
                    face_labels_to_create.push_back(label);
                } else if (v1adjA)
                {
                    faces_to_create.push_back(Vec3st(a, v0, v1));
                    faces_to_create.push_back(Vec3st(a, b, v0));
                    face_labels_to_create.push_back(label);
                    face_labels_to_create.push_back(label);
                } else
                {
                    faces_to_create.push_back(Vec3st(b, v0, v1));
                    face_labels_to_create.push_back(label);
                }
                
            }
            
        }

        //&&&& recursive deletion
//        for (size_t i = 0; i < mesh.m_vertex_to_edge_map[xj].size(); i++)
//        {
//            edges_to_delete.push_back(mesh.m_vertex_to_edge_map[xj][i]);
//        }
        
        // prune flap triangles
        for (size_t i = 0; i < faces_to_create.size(); i++)
        {
            for (size_t j = i + 1; j < faces_to_create.size(); j++)
            {
                Vec3st & f0 = faces_to_create[i];
                Vec3st & f1 = faces_to_create[j];
                
                if (mesh.triangle_has_these_verts(f0, f1))
                {
                    // f0 and f1 have the same vertices
                    
                    Vec2i l0 = face_labels_to_create[i];
                    Vec2i l1 = face_labels_to_create[j];
                    
                    assert(l0[0] == l1[0] || l0[0] == l1[1] || l0[1] == l1[0] || l0[1] == l1[1]);
                    
                    Vec2i newlabel; // newlabel has the same orientation with l0
                    if (l0[0] == l1[0] || l0[0] == l1[1])
                        newlabel = Vec2i(l0[0] == l1[0] ? l1[1] : l1[0], l0[1]);
                    else
                        newlabel = Vec2i(l0[0], l0[1] == l1[0] ? l1[1] : l1[0]);
                    
                    face_labels_to_create[i] = newlabel;
                    
                    faces_to_create.erase(faces_to_create.begin() + j);
                    face_labels_to_create.erase(face_labels_to_create.begin() + j);
                    break;
                }
            }
        }
        
        // apply the deleteion/addition
        for (size_t i = 0; i < faces_to_delete.size(); i++)
            m_surf.remove_triangle(faces_to_delete[i]);

        //&&&& recursive deletion
//        for (size_t i = 0; i < edges_to_delete.size(); i++)
//            m_obj->deleteEdge(edges_to_delete[i], false);
//        m_obj->deleteVertex(xj);
        
        assert(faces_to_create.size() == face_labels_to_create.size());
        for (size_t i = 0; i < faces_to_create.size(); i++)
        {
            size_t nf = m_surf.add_triangle(faces_to_create[i]);
            mesh.set_triangle_label(nf, face_labels_to_create[i]);
            faces_created.push_back(nf);
        }
        
        // mark the two new vertices a and b as dirty
        vertices_to_process.push_back(a);
        vertices_to_process.push_back(b);
        
        // set the vertices new positions
        m_surf.set_newposition(a, a_desired_position);
        m_surf.set_newposition(b, b_desired_position);
        
        m_surf.set_position(a, m_surf.get_newposition(a));
        m_surf.set_position(b, m_surf.get_newposition(b));
        
        // vertex deletion/creation logging
        verts_to_delete.push_back(xj);
        
        verts_to_create.push_back(m_surf.get_newposition(a));
        verts_to_create.push_back(m_surf.get_newposition(b));
        
        verts_created.push_back(a);
        verts_created.push_back(b);

        // Add to new history log
        MeshUpdateEvent edgepop(MeshUpdateEvent::VERTEX_POP);
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
    if (mesh.oriented(v0, v1, mesh.get_triangle(tmp)) > 0)
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
    for (size_t i = 0; i < mesh.m_vertex_to_triangle_map[xj].size(); i++)
    {
        size_t triangle = mesh.m_vertex_to_triangle_map[xj][i];
        
        // find the edge in triangle triangle that's opposite to vertex xj
        size_t l = 0;
        size_t edge0 = static_cast<size_t>(~0);
        size_t edge1 = static_cast<size_t>(~0);
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
        
        for (l = 0; l < 3; l++)
        {
            size_t e = mesh.m_triangle_to_edge_map[triangle][l];
            if ((mesh.m_edges[e][0] == xj && mesh.m_edges[e][1] == v0) ||
                (mesh.m_edges[e][1] == xj && mesh.m_edges[e][0] == v0))
                edge0 = e;
            else if ((mesh.m_edges[e][0] == xj && mesh.m_edges[e][1] == v1) ||
                     (mesh.m_edges[e][1] == xj && mesh.m_edges[e][0] == v1))
                edge1 = e;
        }
        assert(edge0 < mesh.ne());
        assert(edge1 < mesh.ne());
        
        assert(v0 == mesh.m_edges[edge0][0] || v0 == mesh.m_edges[edge0][1]);
        assert(v1 == mesh.m_edges[edge1][0] || v1 == mesh.m_edges[edge1][1]);
        
        Vec3d x0 = m_surf.get_position(v0);
        Vec3d x1 = m_surf.get_position(v1);
        
        Vec3d foot = dot(xxj - x0, x1 - x0) / dot(x1 - x0, x1 - x0) * (x1 - x0) + x0;
        Vec3d force = (foot - xxj);
        force /= mag(force);
        force *= mag(x1 - x0);
        
        Vec2i label = mesh.get_triangle_label(triangle);
        if (label[0] == A || label[1] == A)
        {
            force_a += force;
        } else if (label[0] == B || label[1] == B)
        {
            force_b += force;
        } else
        {
            size_t ev0 = edge0;
            size_t ev1 = edge1;
            
            bool v0adjA = false;
            bool v1adjA = false;
            for (size_t j = 0; j < mesh.m_edge_to_triangle_map[ev0].size(); j++)
            {
                Vec2i label = mesh.get_triangle_label(mesh.m_edge_to_triangle_map[ev0][j]);
                if (label[0] == A || label[1] == A)
                    v0adjA = true;
            }
            for (size_t j = 0; j < mesh.m_edge_to_triangle_map[ev1].size(); j++)
            {
                Vec2i label = mesh.get_triangle_label(mesh.m_edge_to_triangle_map[ev1][j]);
                if (label[0] == A || label[1] == A)
                    v1adjA = true;
            }
            
            assert(!v0adjA || !v1adjA);
            
            if (v0adjA)
            {
                force_a += force;
                force_a += -pull_apart_direction * mag(x1 - xxj);
                force_b += pull_apart_direction * mag(x1 - xxj);
            } else if (v1adjA)
            {
                force_a += force;
                force_a += -pull_apart_direction * mag(x0 - xxj);
                force_b += pull_apart_direction * mag(x0 - xxj);
            } else
            {
                force_b += force;
            }
        }
    }
    
    double tensile_force = dot(force_a - force_b, pull_apart_direction);
    
    return tensile_force > 0;
}

bool T1Transition::pulling_vertex_apart_introduces_collision(size_t v, const Vec3d & oldpos, const Vec3d & newpos0, const Vec3d & newpos1)
{
    
    return true;
}

    

}
