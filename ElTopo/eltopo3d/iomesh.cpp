// ---------------------------------------------------------
//
//  iomesh.cpp
//
//  Non-member functions for reading and writing various mesh file formats.
//
// ---------------------------------------------------------


#include <iomesh.h>

#include <bfstream.h>
#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <nondestructivetrimesh.h>

#ifndef NO_GUI
#include <gluvi.h>
#endif

#define LINESIZE 1024 // maximum line size when reading .OBJ files

// ---------------------------------------------------------
///
/// Write mesh in binary format
///
// ---------------------------------------------------------

bool write_binary_file( const NonDestructiveTriMesh &mesh, 
                       const std::vector<Vec3d> &x,
                       const std::vector<double> &masses, 
                       double curr_t, 
                       const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);   
    bofstream outfile( filename_format, ap );
    va_end(ap);
    
    outfile.write_endianity();
    
    outfile << curr_t;
    
    unsigned int nverts = static_cast<unsigned int>( x.size() );
    outfile << nverts;
    for ( unsigned int i = 0; i < x.size(); ++i )
    {
        outfile << x[i][0];
        outfile << x[i][1];
        outfile << x[i][2];
    }
    
    assert( x.size() == masses.size() );
    
    for ( unsigned int i = 0; i < masses.size(); ++i )
    {
        outfile << masses[i];
    }
    
    unsigned int ntris = static_cast<unsigned int>( mesh.num_triangles() );
    outfile << ntris;
    
    for ( unsigned int t = 0; t < mesh.num_triangles(); ++t )
    {
        const Vec3ui& tri = Vec3ui( mesh.get_triangle(t) );
        outfile << tri[0];
        outfile << tri[1];
        outfile << tri[2];      
    }
    
    outfile.close();
    
    return outfile.good();
}

// ---------------------------------------------------------
///
/// For each vertex in a surface, write which surface it belongs to (useful for rendering multiple surfaces)
///
// ---------------------------------------------------------

bool write_surface_ids( const std::vector<size_t> &ids,
                       const char *filename_format, ... )
{
    va_list ap;
    va_start(ap, filename_format);   
    bofstream outfile( filename_format, ap );
    va_end(ap);
    
    outfile.write_endianity();
    
    size_t nids = ids.size(); 
    outfile << nids;
    
    for ( unsigned int i = 0; i < ids.size(); ++i )
    {
        unsigned int curr_id = static_cast<unsigned int>(ids[i]);
        outfile << curr_id;
    }
    
    outfile.close();
    
    return outfile.good();
}


// ---------------------------------------------------------
///
/// Read mesh in binary format
///
// ---------------------------------------------------------

bool read_binary_file( NonDestructiveTriMesh &mesh, 
                      std::vector<Vec3d> &x, 
                      std::vector<double> &masses, 
                      double& curr_t, 
                      const char *filename_format, ...)
{
    va_list ap;
    va_start(ap, filename_format);   
    bifstream infile( filename_format, ap );
    va_end(ap);
    
    assert( infile.good() );
    
    infile.read_endianity();
    
    infile >> curr_t;
    
    unsigned int nverts;
    infile >> nverts;
    x.resize( nverts );
    for ( unsigned int i = 0; i < nverts; ++i )
    {
        infile >> x[i][0];
        infile >> x[i][1];
        infile >> x[i][2];  
        mesh.nondestructive_add_vertex();
    }
    
    masses.resize( nverts );
    for ( unsigned int i = 0; i < nverts; ++i )
    {
        infile >> masses[i];
    }
    
    unsigned int ntris;
    infile >> ntris;
    
    for ( unsigned int t = 0; t < ntris; ++t )
    {
        Vec3st tri;
        infile >> tri[0];
        infile >> tri[1];
        infile >> tri[2];
        mesh.nondestructive_add_triangle( tri );
    }
    
    infile.close();
    
    return infile.good();
}


// ---------------------------------------------------------
///
/// For each vertex in a surface, read in which surface it belongs to (useful for rendering multiple surfaces)
///
// ---------------------------------------------------------

bool read_surface_ids( std::vector<unsigned int>& ids,
                      const char *filename_format, ... )
{
    
    va_list ap;
    va_start(ap, filename_format);   
    bifstream infile( filename_format, ap );
    va_end(ap);
    
    infile.read_endianity();
    
    unsigned int n;
    infile >> n;
    
    ids.resize(n);
    for ( unsigned int i = 0; i < ids.size(); ++i )
    {
        infile >> ids[i];
    }
    
    infile.close();
    
    return infile.good();
    
}





