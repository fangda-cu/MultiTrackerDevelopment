// ---------------------------------------------------------
//
//  t1transition.cpp
//  Fang Da 2013
//  
//  Functions handling T1 transitions (edge popping and vertex popping).
//
// ---------------------------------------------------------

#include <queue>

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
                if (mesh.get_common_vertex(xjunctions[i], xjgroups[i].first[k]) < mesh.nv())
                {
                    xjgroups[j].first.push_back(xjunctions[i]);
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
            
            size_t nv0 = mesh.nondestructive_add_vertex();
            size_t nv1 = mesh.nondestructive_add_vertex();
            
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
                assert((mesh.oriented(v0, v1, mesh.get_triangle(region0faces[1])) && label1[0] == cut[0]) ||
                       (mesh.oriented(v0, v1, mesh.get_triangle(region0faces[1])) && label1[1] == cut[0]));
            }                
            assert(upper_region >= 0);
            assert(lower_region >= 0);
            
            std::vector<size_t> faces_to_delete;
            std::vector<Vec3st> faces_to_create;
            std::vector<Vec2i> face_labels_to_create;
            
            std::vector<Vec3d> upper_neighbors;
            std::vector<Vec3d> lower_neighbors;
            
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
            
            Vec3d upper_neighbors_mean(0, 0, 0);
            for (size_t j = 0; j < upper_neighbors.size(); j++)
                upper_neighbors_mean += upper_neighbors[j];
            upper_neighbors_mean /= upper_neighbors.size();
            
            Vec3d lower_neighbors_mean(0, 0, 0);
            for (size_t j = 0; j < lower_neighbors.size(); j++)
                lower_neighbors_mean += lower_neighbors[j];
            lower_neighbors_mean /= lower_neighbors.size();
            
            // apply the deletion
            for (size_t j = 0; j < faces_to_delete.size(); j++)
                mesh.nondestructive_remove_triangle(faces_to_delete[j]);    //&&&& recursive deletion
            
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
                size_t nf = mesh.nondestructive_add_triangle(faces_to_create[j]);
                mesh.set_triangle_label(nf, face_labels_to_create[j]);
            }
            
            size_t wall0 = mesh.nondestructive_add_triangle(Vec3st(v0, nv1, nv0));
            size_t wall1 = mesh.nondestructive_add_triangle(Vec3st(nv0, nv1, v1));

            mesh.set_triangle_label(wall0, cut);
            mesh.set_triangle_label(wall1, cut);
            
            // pull the two new vertices (nv0, nv1) apart (before this the two vertices have the same position)
            Vec3d wall_breadth = (upper_neighbors_mean - lower_neighbors_mean);
            wall_breadth /= mag(wall_breadth);
            wall_breadth *= mag(m_surf.get_position(v1) - m_surf.get_position(v0));
            m_surf.set_newposition(nv0, m_surf.get_newposition(nv0) + wall_breadth * 0.1);
            m_surf.set_newposition(nv1, m_surf.get_newposition(nv1) - wall_breadth * 0.1);
            
            m_surf.set_position(nv0, m_surf.get_newposition(nv0));
            m_surf.set_position(nv1, m_surf.get_newposition(nv1));
            
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
                
                size_t nv0 = mesh.nondestructive_add_vertex();
                size_t nv1 = mesh.nondestructive_add_vertex();
                
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
            }
            
            std::vector<size_t> faces_to_delete;
            std::vector<Vec3st> faces_to_create;
            std::vector<Vec2i> face_labels_to_create;
            
            std::vector<size_t> edges_to_delete;
            std::vector<Vec2st> edges_to_create;
            
            std::vector<size_t> vertices_to_delete;
            
            // update the faces incident to the X-junction edges
            for (size_t j = 0; j < ne; j++)
            {
                size_t edge = ordered_edges[j];
                edges_to_delete.push_back(edge);
                edges_to_create.push_back(Vec2st(upper_junctions[j + 0], upper_junctions[j + 1]));
                edges_to_create.push_back(Vec2st(lower_junctions[j + 0], lower_junctions[j + 1]));
                
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
                
                vertices_to_delete.push_back(v);
                
                for (size_t k = 0; k < mesh.m_vertex_to_triangle_map[v].size(); k++)
                {
                    size_t triangle = mesh.m_vertex_to_triangle_map[v][k];
                    
                    // skip if the triangle contains edge0 or edge1
                    Vec2ui dummy;
                    if (mesh.index_in_triangle(mesh.get_triangle(triangle), mesh.m_edges[edge0][0], dummy) < 3 &&
                        mesh.index_in_triangle(mesh.get_triangle(triangle), mesh.m_edges[edge0][1], dummy) < 3)
                        continue;
                    if (mesh.index_in_triangle(mesh.get_triangle(triangle), mesh.m_edges[edge1][0], dummy) < 3 &&
                        mesh.index_in_triangle(mesh.get_triangle(triangle), mesh.m_edges[edge1][1], dummy) < 3)
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
            
            // apply the deletion
            for (size_t j = 0; j < faces_to_delete.size(); j++)
                mesh.nondestructive_remove_triangle(faces_to_delete[j]);    //&&&& recursive deletion
  
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
                size_t nf = mesh.nondestructive_add_triangle(faces_to_create[j]);
                mesh.set_triangle_label(nf, face_labels_to_create[j]);
            }
            
            // triangulate the new interface between cut.x() and cut.y()
            for (size_t j = 0; j < ne; j++)
            {
                size_t wall0, wall1;
                if (j == 0)
                {
                    wall0 = mesh.nondestructive_add_triangle(Vec3st(upper_junctions[j], lower_junctions[j + 1], upper_junctions[j + 1]));
                } else if (j == ne - 1)
                {
                    wall0 = mesh.nondestructive_add_triangle(Vec3st(upper_junctions[j], lower_junctions[j], upper_junctions[j + 1]));
                } else
                {
                    wall0 = mesh.nondestructive_add_triangle(Vec3st(upper_junctions[j], lower_junctions[j], upper_junctions[j + 1]));
                    wall1 = mesh.nondestructive_add_triangle(Vec3st(lower_junctions[j + 1], upper_junctions[j + 1], lower_junctions[j]));
                }

                mesh.set_triangle_label(wall0, cut);
                if (j != 0 && j != ne - 1)
                    mesh.set_triangle_label(wall1, cut);
            }
            
            // pull the two new vertices (nv0, nv1) apart (before this the two vertices have the same position)
            for (size_t j = 0; j < ne - 1; j++)
            {
                size_t nv0 = upper_junctions[j + 1];
                size_t nv1 = lower_junctions[j + 1];
                
                m_surf.set_newposition(nv0, m_surf.get_newposition(nv0) + pull_apart_offsets[j] * 0.1);
                m_surf.set_newposition(nv1, m_surf.get_newposition(nv1) - pull_apart_offsets[j] * 0.1);
                
                m_surf.set_position(nv0, m_surf.get_newposition(nv0));
                m_surf.set_position(nv1, m_surf.get_newposition(nv1));
            }   
            
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
    
    
    
    return pop_occurred;
}

// --------------------------------------------------------
///
/// Decide the cut direction on an X-junction edge
///
// --------------------------------------------------------
    
Mat2i T1Transition::cut_x_junction_edge(size_t e)
{
    Mat2i cut;
    return cut;
}
    


}
