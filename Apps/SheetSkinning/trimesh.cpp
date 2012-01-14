#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include "trimesh.h"

void TriMesh::
clear(void)
{
   tri.clear();
   clear_auxiliaries();
}

void TriMesh::
clear_auxiliaries(void)
{
   nv=0;
   vtinc.clear();
   vadj.clear();
   tadj.clear();
   edge.clear();
}

void TriMesh::
update_nv(void)
{
   nv=0;
   for(unsigned int t=0; t<tri.size(); ++t)
      nv=max(nv, tri[t][0], tri[t][1], tri[t][2]);
   ++nv;
}

void TriMesh::
update_vtinc(void)
{
   vtinc.clear();
   vtinc.resize(nv);
   for(unsigned int t=0; t<tri.size(); ++t){
      unsigned int i, j, k; assign(tri[t], i, j, k);
      vtinc[i].push_back(t);
      vtinc[j].push_back(t);
      vtinc[k].push_back(t);
   }
   /*
   for(unsigned int i=0; i<vtinc.size(); ++i)
      sort(vtinc[i].begin(), vtinc[i].end());
   */
}

// requires nv
void TriMesh::
update_vadj(void)
{
   vadj.clear();
   vadj.resize(nv);
   for(unsigned int t=0; t<tri.size(); ++t){
      unsigned int i, j, k; assign(tri[t], i, j, k);
      add_unique(vadj[i], j);
      add_unique(vadj[i], k);
      add_unique(vadj[j], i);
      add_unique(vadj[j], k);
      add_unique(vadj[k], i);
      add_unique(vadj[k], j);
   }
   for(unsigned int i=0; i<vadj.size(); ++i)
      sort(vadj[i].begin(), vadj[i].end());
}

static void add_unique_intersection(std::vector<unsigned int> &a, const std::vector<unsigned int> &b,
                                    const std::vector<unsigned int> &c, unsigned int excluded)
{
   unsigned int r=0, s=0;
   while(r<b.size() && s<c.size()){
      if(b[r]==c[s]){
         if(b[r]!=excluded) add_unique(a, b[r]);
         ++r;
         ++s;
      }else if(b[r]<c[s])
         ++r;
      else
         ++s;
   }
}

// requires vtinc
void TriMesh::
update_tadj(void)
{
   tadj.clear();
   tadj.resize(tri.size());
   for(unsigned int t=0; t<tri.size(); ++t){
      unsigned int i, j, k; assign(tri[t], i, j, k);
      add_unique_intersection(tadj[t], vtinc[i], vtinc[j], t);
      add_unique_intersection(tadj[t], vtinc[i], vtinc[k], t);
      add_unique_intersection(tadj[t], vtinc[j], vtinc[k], t);
   }
   for(unsigned int t=0; t<tadj.size(); ++t)
      sort(tadj[t].begin(), tadj[t].end());
}

// requires vadj
void TriMesh::
update_edge(void)
{
   edge.clear();
   unsigned int degreesum=0;
   for(unsigned int i=0; i<vadj.size(); ++i)
      degreesum+=(unsigned int)vadj[i].size();
   edge.reserve(degreesum/2);
   for(unsigned int i=0; i<vadj.size(); ++i)
      for(unsigned int j=0; j<vadj[i].size(); ++j)
         if(i<vadj[i][j])
            edge.push_back(Vec2ui(i,vadj[i][j]));
}

