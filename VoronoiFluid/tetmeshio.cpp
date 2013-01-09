
#include "tetmeshio.h"

#include <fstream>

using namespace ElTopo;
// ---------------------------------------------------------

bool load_ele_file( const char* ele_file_name, std::vector<Vec4ui>& tets )
{
   std::ifstream ele_file;
   ele_file.open( ele_file_name );
   if ( !ele_file.good() )
   {
      std::cout << ".ele file not good" << std::endl;
      return false;
   }
   
   int num_tets, nodes_per_tet, boundary_markers;
   ele_file >> num_tets;
   ele_file >> nodes_per_tet;
   ele_file >> boundary_markers;
   
   if ( nodes_per_tet != 4 )
   {
      std::cout << "error: expected 4 nodes per tet. nodes_per_tet: " << nodes_per_tet << std::endl;
      return false;
   }
   
   if ( boundary_markers != 0 )
   {
      std::cout << "error: expected no boundary_markers. boundary_markers: " << boundary_markers << std::endl;
      return false;
   }
   
   for ( int i = 0; i < num_tets; ++i )
   {
      int index;
      ele_file >> index;
      assert( index == i + 1 );
      int a, b, c, d;
      ele_file >> a >> b >> c >> d;
      tets.push_back( Vec4ui(a,b,c,d) - Vec4ui(1) );
   }
   
   return true;
}

// ---------------------------------------------------------

bool write_ele_file( const char* ele_file_name, std::vector<Vec4ui>& tets )
{
   std::ofstream ele_file;
   ele_file.open( ele_file_name );
   if ( !ele_file.good() )
   {
      std::cout << ".ele file not good" << std::endl;
      return false;
   }
   
   int num_tets = tets.size();
   int nodes_per_tet = 4;
   int boundary_markers = 0;
   
   ele_file << num_tets;
   ele_file << nodes_per_tet;
   ele_file << boundary_markers;
      
   for ( int i = 0; i < num_tets; ++i )
   {
      ele_file << i+1;
      int a = tets[i][0] + 1, 
          b = tets[i][1] + 1, 
          c = tets[i][2] + 1, 
          d = tets[i][3] + 1;
      
      ele_file << a << b << c << d;
   }
   
   return true;
}

// ---------------------------------------------------------

bool load_node_file( const char* node_file_name, std::vector<Vec3f>& xs )
{
   std::ifstream node_file;
   node_file.open( node_file_name );
   if ( !node_file.good() )
   {
      return false;
   }
   
   int num_nodes, dim, attributes, boundary_markers;
   node_file >> num_nodes;
   node_file >> dim;
   if ( dim != 3 )
   {
      std::cout << "error: expected 3d tet mesh vertices.  dim: " << dim << std::endl;
      return false;
   }
   
   node_file >> attributes;
   if ( attributes != 0 )
   {
      std::cout << "error: expected no attributes.  attributes: " << attributes << std::endl;
      return false;
   }
   
   node_file >> boundary_markers;
   if ( boundary_markers != 0 )
   {
      std::cout << "error: expected no boundary_markers.  boundary_markers: " << attributes << std::endl;
      return false;
   }
   
   for ( int i = 0; i < num_nodes; ++i )
   {
      int index;
      node_file >> index;
      assert( index == i + 1 );
      float x, y, z;
      node_file >> x >> y >> z;
      xs.push_back( Vec3f(x,y,z) );
   }
   
   return true;
}


// ---------------------------------------------------------

bool write_node_file( const char* node_file_name, std::vector<Vec3f>& xs )
{
   std::ofstream node_file;
   node_file.open( node_file_name );
   if ( !node_file.good() )
   {
      return false;
   }
   
   int num_nodes = xs.size();
   int dim = 3;
   int attributes = 0;
   int boundary_markers = 0;
   
   node_file << num_nodes << dim << attributes << boundary_markers;
   
   for ( int i = 0; i < num_nodes; ++i )
   {
      node_file << i+1;
      float x = xs[i][0], 
            y = xs[i][1],
            z = xs[i][2];      
      node_file << x << y << z;
   }
   
   return true;
}
