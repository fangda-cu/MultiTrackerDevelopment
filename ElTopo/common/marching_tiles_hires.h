#ifndef MARCHING_TILES_HIRES_H
#define MARCHING_TILES_HIRES_H

#include "array3.h"
#include "hashtable.h"
#include "vec.h"

namespace ElTopo {
struct MarchingTilesHiRes
{
   ElTopo::Array1<ElTopo::Vec3st> tri;
   ElTopo::Array1<ElTopo::Vec3d> x;
   ElTopo::Array1<ElTopo::Vec3d> normal;
   ElTopo::Vec3d origin;
   double dx;
   const ElTopo::Array3d& phi;

   MarchingTilesHiRes(const ElTopo::Vec3d &origin_, double dx_, const ElTopo::Array3d& phi_) : 
      tri(0), x(0), normal(0),
      origin(origin_), dx(dx_), phi(phi_),
      edge_cross()
   {}

   void contour(void);
   void improve_mesh(void);
   void estimate_normals(void);

   private:
   ElTopo::HashTable<ElTopo::Vec6i,unsigned int> edge_cross; // stores vertices that have been created already at given edge crossings

   double eval(double i, double j, double k); // interpolate if non-integer coordinates given
   void eval_gradient(double i, double j, double k, ElTopo::Vec3d& grad);
   void contour_tile(int i, int j, int k); // add triangles for contour in the given tile (starting at grid point (4*i,4*j,4*k))
   void contour_tet(const ElTopo::Vec3i& x0, const ElTopo::Vec3i& x1, const ElTopo::Vec3i& x2, const ElTopo::Vec3i& x3, double p0, double p1, double p2, double p3);
   int find_edge_cross(const ElTopo::Vec3i& x0, const ElTopo::Vec3i& x1, double p0, double p1);
};
}

#endif
