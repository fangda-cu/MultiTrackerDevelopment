#include "CGAL_wrapper.h"

#include "CGAL/Exact_predicates_inexact_constructions_kernel.h"
#include "CGAL/Delaunay_triangulation_2.h"
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>

#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <map>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K2D;

typedef CGAL::Delaunay_triangulation_2<K2D>      Triangulation2D;

typedef Triangulation2D::Vertex_handle  Vertex_handle2D;
typedef Triangulation2D::Point          Point2D;

using namespace ElTopo;

void compute_delaunay_CGAL(const std::vector<Vec2f>& points, std::vector<Vec3ui>& tris) {
   //construct point list in CGAL format
   Triangulation2D T;
   std::map<Vertex_handle2D, int> vert_handles;
   for(unsigned int i = 0; i < points.size(); ++i) {
      Point2D p(points[i][0], points[i][1]);
      Vertex_handle2D vh = T.insert(p);
      vert_handles[vh] = i;
   }
   
   //iterate over all tris now
   Triangulation2D::Face_iterator fit = T.faces_begin();
   for(;fit != T.faces_end(); ++fit) {
      //get the face
      Triangulation2D::Face curFace = *fit;
      Vertex_handle2D va = curFace.vertex(0);
      Vertex_handle2D vb = curFace.vertex(1);
      Vertex_handle2D vc = curFace.vertex(2);
      if(va == T.infinite_vertex()) continue;
      if(vb == T.infinite_vertex()) continue;
      if(vc == T.infinite_vertex()) continue;

      Vec3ui myFace(vert_handles[va], vert_handles[vb], vert_handles[vc]);
      
      tris.push_back(myFace);
   }

}

typedef CGAL::Regular_triangulation_filtered_traits_2<K>  Traits;
typedef CGAL::Regular_triangulation_2<Traits> Regular_triangulation;

void compute_regular_CGAL(const std::vector<Vec2f>& points, const std::vector<float>& weights, std::vector<Vec3ui>& tris, std::vector<Vec2f>& weighted_circumcentres) {
   
   //construct point list in CGAL format
   Regular_triangulation T;
   std::map<Regular_triangulation::Vertex_handle, int> vert_handles;
   for(unsigned int i = 0; i < points.size(); ++i) {
      Regular_triangulation::Point p(points[i][0], points[i][1]);
      Regular_triangulation::Weighted_point wp(p, Regular_triangulation::Weight(weights[i]));
      Regular_triangulation::Vertex_handle vh = T.insert(wp);
      vert_handles[vh] = i;
   }

   //iterate over all tris now
   Regular_triangulation::Face_iterator fit = T.faces_begin();
   for(;fit != T.faces_end(); ++fit) {
      //get the face
      Regular_triangulation::Face curFace = *fit;
      Regular_triangulation::Vertex_handle va = curFace.vertex(0);
      Regular_triangulation::Vertex_handle vb = curFace.vertex(1);
      Regular_triangulation::Vertex_handle vc = curFace.vertex(2);
      if(va == T.infinite_vertex()) continue;
      if(vb == T.infinite_vertex()) continue;
      if(vc == T.infinite_vertex()) continue;
      
      Vec3ui myFace(vert_handles[va], vert_handles[vb], vert_handles[vc]);
      
      tris.push_back(myFace);
      Regular_triangulation::Point p = T.dual(fit);
      Vec2f circumcentre((float)p[0], (float)p[1]);
      weighted_circumcentres.push_back(circumcentre);
   }

   if(T.number_of_hidden_vertices() > 0)
      std::cout << "Non-zero hidden vertex count\n";

   
}
