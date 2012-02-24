// ---------------------------------------------------------
//
//  iomesh.h
//
//  Non-member functions for reading and writing various mesh file formats.
//
// ---------------------------------------------------------

#ifndef EL_TOPO_IOMESH_H
#define EL_TOPO_IOMESH_H

// ---------------------------------------------------------
// Nested includes
// ---------------------------------------------------------

#include <fstream>
#include <vec.h>
#include <vector>

// ---------------------------------------------------------
//  Forwards and typedefs
// ---------------------------------------------------------

namespace ElTopo {

class NonDestructiveTriMesh;

namespace Gluvi
{
    struct Target3D;
}

// ---------------------------------------------------------
//  Function declarations
// ---------------------------------------------------------

/// write mesh in our own binary format
///
bool write_binary_file(const NonDestructiveTriMesh &mesh,  
                       const std::vector<Vec3d> &x, 
                       const std::vector<double> &masses, 
                       double curr_t, const char *filename_format, ...);

/// For each vertex in a surface, write which surface it belongs to (useful for rendering multiple surfaces)
///
bool write_surface_ids(const std::vector<size_t> &ids,
                       const char *filename_format, ... );

/// Read mesh in binary format
///
bool read_binary_file( NonDestructiveTriMesh &mesh, std::vector<Vec3d> &x, std::vector<double> &masses, double& curr_t, const char *filename_format, ...);

/// For each vertex in a surface, read in which surface it belongs to (useful for rendering multiple surfaces)
///
bool read_surface_ids( std::vector<unsigned int> &ids,
                      const char *filename_format, ... );

/// Write an STL vector to an ASCII file.  Not really mesh-related, but useful.
///
template<class T> inline void dump_vector_to_file( const char* filename, const std::vector<T, std::allocator<T> >& vec )
{
    std::ofstream outfile( filename, std::ios::out|std::ios::trunc );
    for ( unsigned int i = 0; i < vec.size(); ++i )
    {
        outfile << vec[i] << std::endl;
    }         
    outfile.close();
}

}

#endif
