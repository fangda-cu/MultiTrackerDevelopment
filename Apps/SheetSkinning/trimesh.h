#ifndef TRIMESH_H
#define TRIMESH_H

#include <vector>
#include "vec.h"

// for general triangle meshes
struct TriMesh{
   // list of triangles: the fundamental data
   std::vector<Vec3ui> tri;
   unsigned int nv; // the highest vertex number tri contains (or at least an upper bound)

   // auxiliary structures computed from triangles
   std::vector<std::vector<unsigned int> > vtinc; // sorted triangle incidence lists (for each vertex, which triangles it's in)
   std::vector<std::vector<unsigned int> > vadj; // sorted vertex adjacency lists
   std::vector<std::vector<unsigned int> > tadj; // sorted triangle adjacency lists
   std::vector<Vec2ui> edge; // list of edges

   TriMesh()
   {}

   void clear(void);
   void clear_auxiliaries(void);
   void update_nv(void);
   void update_vtinc(void); // requires nv
   void update_vadj(void); // requires nv
   void update_tadj(void); // require vtinc
   void update_edge(void); // requires vadj
};

#endif
