#ifndef CGAL_WRAPPER_H
#define CGAL_WRAPPER_H

#include "vec.h"
#include <vector>
#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Delaunay_triangulation_3.h"
#include "CGAL/Triangulation_vertex_base_with_info_3.h"
#include "CGAL/Triangulation_cell_base_with_info_3.h"
#include "CGAL/Location_policy.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Triangulation_cell_base_with_info_3<int, K> Cb;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, K> Vb;

typedef CGAL::Triangulation_data_structure_3<Vb,Cb> Tds;

typedef CGAL::Delaunay_triangulation_3<K, Tds, CGAL::Fast_location>      Triangulation;

typedef Triangulation::Vertex_handle  Vertex_handle;
typedef Triangulation::Point          Point;
typedef Triangulation::Locate_type    Locate_type;
typedef Triangulation::Cell_handle    Cell_handle;

void compute_delaunay_CGAL(const std::vector<ElTopo::Vec2f>& points, std::vector<ElTopo::Vec3st>& tris);

void compute_delaunay_CGAL(const std::vector<ElTopo::Vec3f>& points, std::vector<ElTopo::Vec4st>& tets);
void compute_delaunay_CGAL(const std::vector<ElTopo::Vec3f>& points, std::vector<ElTopo::Vec4st>& tets, Triangulation& T);


void compute_regular_CGAL(const std::vector<ElTopo::Vec2f>& points, const std::vector<float>& weights, std::vector<ElTopo::Vec3st>& tris, std::vector<ElTopo::Vec2f>& circums);

#endif