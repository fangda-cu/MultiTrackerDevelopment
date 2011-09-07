// ---------------------------------------------------------
//
//  broadphasegrid.cpp
//  Tyson Brochu 2008
//  
//  Broad phase collision detection culling using three regular, volumetric grids.
//
// ---------------------------------------------------------

// ---------------------------------------------------------
// Includes
// ---------------------------------------------------------

#include "broadphasegrid.hh"

// ---------------------------------------------------------
// Global externs
// ---------------------------------------------------------

// ---------------------------------------------------------
// Local constants, typedefs, macros
// ---------------------------------------------------------

// ---------------------------------------------------------
// Static function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

// --------------------------------------------------------
///
/// Construct one grid from the given set of AABBs, using the given length scale as the cell size, with the given padding
///
// --------------------------------------------------------
namespace ElTopo {

Vec3d toElTopo(BASim::Vec3d vec) {
   return Vec3d(vec[0], vec[1], vec[2]);
}

void BroadPhaseGrid::build_acceleration_grid( AccelerationGrid& grid, 
                                              std::vector<Vec3d>& xmins, 
                                              std::vector<Vec3d>& xmaxs, 
                                              double length_scale, 
                                              double grid_padding )
{

   Vec3d xmax = xmaxs[0];
   Vec3d xmin = xmins[0];
   double maxdistance = 0;
   
   unsigned int n = xmins.size();
   
   for(unsigned int i = 0; i < n; i++)
   {
      update_minmax(xmins[i], xmin, xmax);
      update_minmax(xmaxs[i], xmin, xmax);
      maxdistance = std::max(maxdistance, mag(xmaxs[i] - xmins[i]));
   }
   
   for(unsigned int i = 0; i < 3; i++)
   {
      xmin[i] -= 2*maxdistance + grid_padding;
      xmax[i] += 2*maxdistance + grid_padding;
   }
   
   Vec3ui dims(1,1,1);
          
   if(mag(xmax-xmin) > grid_padding)
   {
      for(unsigned int i = 0; i < 3; i++)
      {
         unsigned int d = (unsigned int)std::ceil((xmax[i] - xmin[i])/length_scale);
         
         if(d < 1) d = 1;
         if(d > n) d = n;
         dims[i] = d;
      }
   }
      
   grid.set(dims, xmin, xmax);
   
   for(unsigned int i = n; i > 0; i--)
   {
      unsigned int index = i - 1;
      
      // don't add inside-out AABBs
      if ( xmins[index][0] > xmaxs[index][0] )  { continue; }
      
      grid.add_element(index, xmins[index], xmaxs[index]);
   }
}


// --------------------------------------------------------
///
/// Rebuild acceleration grids according to the given triangle mesh
///
// --------------------------------------------------------

void BroadPhaseGrid::update_broad_phase_static( const BASim::TopologicalObject& m_obj, const BASim::VertexProperty<BASim::Vec3d>& vertices, double proximity_epsilon  )
{
   
   double sum = 0;
   int count = 0;
   for(BASim::EdgeIterator iter = m_obj.edges_begin(); iter != m_obj.edges_end(); ++iter) {
      BASim::EdgeHandle edge = *iter;
      BASim::VertexHandle from = m_obj.fromVertex(edge);
      BASim::VertexHandle to = m_obj.toVertex(edge);
      BASim::Vec3d fromPos = vertices.operator[](from);
      BASim::Vec3d toPos = vertices.operator[](to);
      double len = (fromPos - toPos).norm();
      sum += len;
      ++count;
   }
   double grid_scale = sum / (double) count;
   
   {
      unsigned int num_vertices = m_obj.nv();
      std::vector<Vec3d> vertex_xmins(num_vertices), vertex_xmaxs(num_vertices);
      
      for(BASim::VertexIterator iter = m_obj.vertices_begin(); iter != m_obj.vertices_end(); ++iter) {
         BASim::VertexHandle vert = *iter;
         BASim::Vec3d pos = vertices[vert];
         int i = vert.idx(); //using internal index is a bit dangerous, if internal relabelling is ever allowed...
         if(vertex_xmins.size() <= i) {
           vertex_xmins.resize(i+1);
           vertex_xmaxs.resize(i+1);
         }
         vertex_xmins[i] = Vec3d(pos[0] - proximity_epsilon, pos[1] - proximity_epsilon, pos[2] - proximity_epsilon);
         vertex_xmaxs[i] = Vec3d(pos[0] + proximity_epsilon, pos[1] + proximity_epsilon, pos[2] + proximity_epsilon);;
      }      
      build_acceleration_grid( m_vertex_grid, vertex_xmins, vertex_xmaxs, grid_scale, proximity_epsilon );
   }
   
   
   {
      unsigned int num_edges = m_obj.ne();
      Vec3d offset(proximity_epsilon, proximity_epsilon, proximity_epsilon);
      std::vector<Vec3d> edge_xmins(num_edges), edge_xmaxs(num_edges);
      
      for(BASim::EdgeIterator iter = m_obj.edges_begin(); iter != m_obj.edges_end(); ++iter) {
         BASim::EdgeHandle edge = *iter;
         Vec3d pos0 = toElTopo(vertices[m_obj.fromVertex(edge)]), 
            pos1 = toElTopo(vertices[m_obj.toVertex(edge)]);
         Vec3d min_v, max_v;
         minmax(pos0, pos1, min_v, max_v);
         int i = edge.idx(); //using internal index is a bit dangerous, if internal relabelling is ever allowed...
         if(edge_xmins.size() <= i) {
           edge_xmins.resize(i+1);
           edge_xmaxs.resize(i+1);
         }
         edge_xmins[i] = min_v - offset;
         edge_xmaxs[i] = max_v + offset;
      }      
      build_acceleration_grid( m_edge_grid, edge_xmins, edge_xmaxs, grid_scale, proximity_epsilon );
   }
   
   {
      unsigned int num_triangles = m_obj.nf(); 
      Vec3d offset(proximity_epsilon, proximity_epsilon, proximity_epsilon);
      std::vector<Vec3d> tri_xmins(num_triangles), tri_xmaxs(num_triangles);
      
      for(BASim::FaceIterator iter = m_obj.faces_begin(); iter != m_obj.faces_end(); ++iter) {
         BASim::FaceHandle face = *iter;
         Vec3d min_v, max_v;
         int c = 0;
         for(BASim::FaceVertexIterator fvit = m_obj.fv_iter(face); fvit; ++fvit) {
           BASim::VertexHandle v = *fvit;
            Vec3d pos0 = toElTopo(vertices[v]);
            if(c == 0) {
              min_v = pos0; 
              max_v = pos0;
            }
            else {
              update_minmax(pos0, min_v, max_v);
            }
         }
         int i = face.idx(); //using internal index is a bit dangerous, if internal relabelling is ever allowed...
         if(tri_xmins.size() <= i) {
           tri_xmins.resize(i+1);
           tri_xmaxs.resize(i+1);
         }
         tri_xmins[i] = min_v - offset;
         tri_xmaxs[i] = max_v + offset;
      }  
      build_acceleration_grid( m_triangle_grid, tri_xmins, tri_xmaxs, grid_scale, proximity_epsilon );  
   }
   
}



// --------------------------------------------------------
///
/// Rebuild acceleration grids according to the given triangle mesh
///
// --------------------------------------------------------
/*
void BroadPhaseGrid::update_broad_phase_continuous( const BASim::TopologicalObject& m_obj, const BASim::VertexProperty<BASim::Vec3d>& vertices  )
{
   double grid_scale = surface.get_average_edge_length();
   
   {
      unsigned int num_vertices = surface.m_positions.size();
      std::vector<Vec3d> vertex_xmins(num_vertices), vertex_xmaxs(num_vertices);
      for(unsigned int i = 0; i < num_vertices; i++)
      {           
         surface.vertex_continuous_bounds(i, vertex_xmins[i], vertex_xmaxs[i]);
      }
      build_acceleration_grid( m_vertex_grid, vertex_xmins, vertex_xmaxs, grid_scale, surface.m_proximity_epsilon );
   }
   
   {
      unsigned int num_edges = surface.m_mesh.m_edges.size();
      std::vector<Vec3d> edge_xmins(num_edges), edge_xmaxs(num_edges);
      for(unsigned int i = 0; i < num_edges; i++)
      {
         surface.edge_continuous_bounds(i, edge_xmins[i], edge_xmaxs[i]);
      }
      build_acceleration_grid( m_edge_grid, edge_xmins, edge_xmaxs, grid_scale, surface.m_proximity_epsilon );
   }
   
   {
      unsigned int num_triangles = surface.m_mesh.m_tris.size();
      std::vector<Vec3d> tri_xmins(num_triangles), tri_xmaxs(num_triangles);
      for(unsigned int i = 0; i < num_triangles; i++)
      {            
         surface.triangle_continuous_bounds(i, tri_xmins[i], tri_xmaxs[i]);
      }
      build_acceleration_grid( m_triangle_grid, tri_xmins, tri_xmaxs, grid_scale, surface.m_proximity_epsilon );  
   }
   
}
*/
}