#include "CGAL_wrapper.h"


#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include <vector>
#include <map>

using namespace ElTopo;

void compute_delaunay_CGAL(const std::vector<Vec3f>& points, std::vector<Vec4ui>& tets) {
  Triangulation T;
  compute_delaunay_CGAL(points, tets, T);
}

void compute_delaunay_CGAL(const std::vector<Vec3f>& points, std::vector<Vec4ui>& tets, Triangulation& T) {
   //construct point list in CGAL format
  T.clear();
  std::map<Vertex_handle, int> vert_handles;
   for(unsigned int i = 0; i < points.size(); ++i) {
      Point p(points[i][0], points[i][1], points[i][2]);
      Vertex_handle vh = T.insert(p);
      vert_handles[vh] = i;
   }
   
   //iterate over all tris now
   Triangulation::Cell_iterator cit = T.cells_begin();
   for(;cit != T.cells_end(); ++cit) {
      //get the face
      Triangulation::Cell curCell = *cit;
      cit->info() = -1;
      Vertex_handle va = curCell.vertex(0);
      Vertex_handle vb = curCell.vertex(1);
      Vertex_handle vc = curCell.vertex(2);
      Vertex_handle vd = curCell.vertex(3);
      if(va == T.infinite_vertex()) continue;
      if(vb == T.infinite_vertex()) continue;
      if(vc == T.infinite_vertex()) continue;
      if(vd == T.infinite_vertex()) continue;

      Vec4ui myTet(vert_handles[va], vert_handles[vb], vert_handles[vc], vert_handles[vd]);
      
      cit->info() = tets.size(); //store the tet's index
      tets.push_back(myTet);
   }
   

}