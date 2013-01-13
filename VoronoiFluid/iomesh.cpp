
#include "iomesh.h"

#include <cstdarg>
#include <cstdlib>
#include <cmath>
#include <fstream>

#include <nondestructivetrimesh.h>
#include <gluvi.h>
#include <bfstream.h>

#define LINESIZE 1024 // maximum line size when reading .OBJ files

using namespace ElTopo;

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
   
   outfile << x.size();
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
   
   outfile << mesh.num_triangles();
   
   for ( unsigned int t = 0; t < mesh.num_triangles(); ++t )
   {
      Vec3st tri = mesh.get_triangle(t);
      outfile << tri[0];
      outfile << tri[1];
      outfile << tri[2];      
   }
   
   outfile.close();
   
   return outfile.good();
}

// ---------------------------------------------------------
///
/// Write mesh in binary format, with per-vertex velocities
///
// ---------------------------------------------------------

bool write_binary_file_with_velocities( const NonDestructiveTriMesh &mesh, 
                                       const std::vector<Vec3d> &x,
                                       const std::vector<double> &masses,                                       
                                       const std::vector<Vec3d> &v,
                                       double curr_t, 
                                       const char *filename_format, ...)
{
   
   va_list ap;
   va_start(ap, filename_format);   
   bofstream outfile( filename_format, ap );
   va_end(ap);
   
   outfile.write_endianity();
   
   outfile << curr_t;
   
   outfile << x.size();
   
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
   
   for ( unsigned int i = 0; i < v.size(); ++i )
   {
      outfile << v[i][0];
      outfile << v[i][1];
      outfile << v[i][2];
   }
   
   outfile << mesh.num_triangles();
   
   for ( unsigned int t = 0; t < mesh.num_triangles(); ++t )
   {
      Vec3st tri = mesh.get_triangle(t);
      outfile << tri[0];
      outfile << tri[1];
      outfile << tri[2];      
   }
   
   outfile.close();
   
   return outfile.good();
}


// ---------------------------------------------------------
///
///
///
// ---------------------------------------------------------

bool write_binary_file_with_newpositions( const NonDestructiveTriMesh &mesh, 
                                          const std::vector<Vec3d> &x, 
                                          const std::vector<double> &masses, 
                                          const std::vector<Vec3d> &new_positions, 
                                          double curr_t, 
                                          const char *filename_format, ...)
{
   
   va_list ap;
   va_start(ap, filename_format);   
   bofstream outfile( filename_format, ap );
   va_end(ap);
   
   outfile.write_endianity();
   
   outfile << curr_t;
   
   outfile << x.size();
   
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
   
   for ( unsigned int i = 0; i < new_positions.size(); ++i )
   {
      outfile << new_positions[i][0];
      outfile << new_positions[i][1];
      outfile << new_positions[i][2];
   }
   
   outfile << mesh.num_triangles();
   
   for ( unsigned int t = 0; t < mesh.num_triangles(); ++t )
   {
      Vec3st tri = mesh.get_triangle(t);
      outfile << tri[0];
      outfile << tri[1];
      outfile << tri[2];      
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
   }
   
   masses.resize( nverts );
   for ( unsigned int i = 0; i < masses.size(); ++i )
   {
      infile >> masses[i];
   }
   
   unsigned int ntris;
   infile >> ntris;
   mesh.m_tris.resize( ntris );
   for ( unsigned int t = 0; t < ntris; ++t )
   {
      infile >> mesh.m_tris[t][0];
      infile >> mesh.m_tris[t][1];
      infile >> mesh.m_tris[t][2];
   }
   
   infile.close();
   
   return infile.good();
}


// ---------------------------------------------------------
///
/// Read mesh in binary format, with per-vertex velocities
///
// ---------------------------------------------------------

bool read_binary_file_with_velocities( NonDestructiveTriMesh &mesh, 
                                      std::vector<Vec3d> &x, 
                                      std::vector<double> &masses,
                                      std::vector<Vec3d> &v, 
                                      double& curr_t, 
                                      const char *filename_format, ...)
{
   va_list ap;
   va_start(ap, filename_format);   
   bifstream infile( filename_format, ap );
   va_end(ap);
   
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
   }
   
   masses.resize( nverts );
   for ( unsigned int i = 0; i < masses.size(); ++i )
   {
      infile >> masses[i];
   }
   
   v.resize( nverts );
   for ( unsigned int i = 0; i < nverts; ++i )
   {
      infile >> v[i][0];
      infile >> v[i][1];
      infile >> v[i][2];      
   }
   
   unsigned int ntris;
   infile >> ntris;
   mesh.m_tris.resize( ntris );
   for ( unsigned int t = 0; t < ntris; ++t )
   {
      infile >> mesh.m_tris[t][0];
      infile >> mesh.m_tris[t][1];
      infile >> mesh.m_tris[t][2];
   }
   
   infile.close();
   
   return infile.good();
}


// ---------------------------------------------------------
///
/// 
///
// ---------------------------------------------------------

bool read_binary_file_with_newpositions( NonDestructiveTriMesh &mesh, 
                                      std::vector<Vec3d> &x, 
                                      std::vector<double> &masses,
                                      std::vector<Vec3d> &new_positions, 
                                      double& curr_t, 
                                      const char *filename_format, ...)
{
   va_list ap;
   va_start(ap, filename_format);   
   bifstream infile( filename_format, ap );
   va_end(ap);
   
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
   }
   
   masses.resize( nverts );
   for ( unsigned int i = 0; i < masses.size(); ++i )
   {
      infile >> masses[i];
   }
   
   new_positions.resize( nverts );
   for ( unsigned int i = 0; i < nverts; ++i )
   {
      infile >> new_positions[i][0];
      infile >> new_positions[i][1];
      infile >> new_positions[i][2];      
   }
   
   unsigned int ntris;
   infile >> ntris;
   mesh.m_tris.resize( ntris );
   for ( unsigned int t = 0; t < ntris; ++t )
   {
      infile >> mesh.m_tris[t][0];
      infile >> mesh.m_tris[t][1];
      infile >> mesh.m_tris[t][2];
   }
   
   infile.close();
   
   return infile.good();
}

// ---------------------------------------------------------
///
/// Write mesh in Wavefront OBJ format
///
// ---------------------------------------------------------

bool write_objfile(const NonDestructiveTriMesh &mesh, const std::vector<Vec3d> &x, const char *filename_format, ...)
{
   va_list ap;
   va_start(ap, filename_format);
#ifdef WIN32
   int len=_vscprintf(filename_format, ap) +1;// _vscprintf doesn't count terminating '\0'
   char *filename=new char[len];
   vsprintf(filename, filename_format, ap);
#else
   char *filename;
   vasprintf(&filename, filename_format, ap);
#endif
   std::cout << "Writing " << filename << std::endl;
   
   std::ofstream output(filename, std::ofstream::binary);
#ifdef WIN32
   delete [] filename;
#else
   std::free(filename);
#endif
   va_end(ap);

   if(!output.good()) return false;

   output<<"# generated by editmesh"<<std::endl;
   for(unsigned int i=0; i<x.size(); ++i)
      output<<"v "<<x[i]<<std::endl;
   for(unsigned int t=0; t<mesh.m_tris.size(); ++t)
      output<<"f "<<mesh.m_tris[t][0]+1<<' '<<mesh.m_tris[t][1]+1<<' '<<mesh.m_tris[t][2]+1<<std::endl; // correct for 1-based indexing in OBJ files
   return output.good();
}

// ---------------------------------------------------------
///
/// Helper for reading OBJ file
///
// ---------------------------------------------------------

static bool read_int(const char *s, int &value, bool &leading_slash, int &position)
{
   leading_slash=false;
   for(position=0; s[position]!=0; ++position){
      switch(s[position]){
         case '/':
            leading_slash=true;
            break;
	 case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
            goto found_int;
      }
   }
   return false;
   
   found_int:
   value=0;
   for(;; ++position){
      switch(s[position]){
	 case '0': case '1': case '2': case '3': case '4': case '5': case '6': case '7': case '8': case '9':
            value=10*value+s[position]-'0';
            break;
         default:
            return true;
      }
   }

#ifndef _MSC_VER //this yields an annoying warning on VC++, so strip it.
   return true; // should never get here, but keeps compiler happy
#endif
}


// ---------------------------------------------------------
///
/// Helper for reading OBJ file
///
// ---------------------------------------------------------

static void read_face_list(const char *s, std::vector<int> &vertex_list)
{
   vertex_list.clear();
   int v, skip;
   bool leading_slash;
   for(int i=0;;){
      if(read_int(s+i, v, leading_slash, skip)){
         if(!leading_slash)
            vertex_list.push_back(v-1); // correct for 1-based index
         i+=skip;
      }else
         break;
   }
}

// ---------------------------------------------------------
///
/// Read mesh in Wavefront OBJ format
///
// ---------------------------------------------------------

bool read_objfile(NonDestructiveTriMesh &mesh, std::vector<Vec3d> &x, const char *filename_format, ...)
{
   va_list ap;
   va_start(ap, filename_format);

#ifdef WIN32
   int len=_vscprintf(filename_format, ap) +1;// _vscprintf doesn't count // terminating '\0'
   char *filename=new char[len];
   vsprintf(filename, filename_format, ap);
#else
   char *filename;
   vasprintf(&filename, filename_format, ap);
#endif

   std::ifstream input(filename, std::ifstream::binary);

#ifdef WIN32
   delete [] filename;
#else
   std::free(filename);
#endif

   va_end(ap);

   if(!input.good()) return false;

   x.clear();
   mesh.clear();

   char line[LINESIZE];
   std::vector<int> vertex_list;
   while(input.good()){
      input.getline(line, LINESIZE);
      switch(line[0]){
         case 'v': // vertex data
            if(line[1]==' '){
               Vec3d new_vertex;
               std::sscanf(line+2, "%lf %lf %lf", &new_vertex[0], &new_vertex[1], &new_vertex[2]);
               x.push_back(new_vertex);
            }
            break;
         case 'f': // face data
            if(line[1]==' '){
               read_face_list(line+2, vertex_list);
               for(int j=0; j<(int)vertex_list.size()-2; ++j)
                  mesh.m_tris.push_back(Vec3ui(vertex_list[0], vertex_list[j+1], vertex_list[j+2]));
            }
            break;
      }
   }
   return true;
}

// ---------------------------------------------------------
///
/// Write mesh in Renderman RIB format.
///
// ---------------------------------------------------------

bool write_ribfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, const char *filename_format, ...)
{
   va_list ap;
   va_start(ap, filename_format);
   
#ifdef WIN32
   int len=_vscprintf(filename_format, ap) +1;// _vscprintf doesn't count // terminating '\0'
   char *filename=new char[len];
   vsprintf(filename, filename_format, ap);
#else
   char *filename;
   vasprintf(&filename, filename_format, ap);
#endif
   
   std::ofstream output(filename, std::ofstream::binary);
#ifdef WIN32
   delete [] filename;
#else
   std::free(filename);
#endif
   va_end(ap);

   if(!output.good()) return false;
    return write_ribfile(mesh, x, output);
}

// ---------------------------------------------------------
///
/// Write mesh in Renderman RIB format.
///
// ---------------------------------------------------------

bool write_ribfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, std::ostream &output)
{
   output<<"# generated by editmesh"<<std::endl;
   output<<"PointsPolygons"<<std::endl;
   output<<" [ ";
   for(unsigned int i=0; i<mesh.m_tris.size(); ++i){
      output<<"3 ";
      if(i%38==37 && i!=mesh.m_tris.size()-1) output<<std::endl;
   }
   output<<"]"<<std::endl;
   output<<" [ ";
   for(unsigned int i=0; i<mesh.m_tris.size(); ++i){
      output<<mesh.m_tris[i]<<"  ";
      if(i%6==5 && i!=mesh.m_tris.size()-1) output<<std::endl;
   }
   output<<"]"<<std::endl;
   output<<" \"P\" [";
   for(unsigned int i=0; i<x.size(); ++i){
      output<<x[i]<<"  ";
      if(i%4==3 && i!=x.size()-1) output<<std::endl;
   }
   output<<"]"<<std::endl;
   
   return output.good();
}

// ---------------------------------------------------------
///
/// Write mesh in Renderman RIB format.
///
// ---------------------------------------------------------

bool write_ribfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, FILE *output)
{
   fprintf( output, "# generated by editmesh\n" );
   fprintf( output, "PointsPolygons\n" );
   fprintf( output, " [ " );
   for(unsigned int i=0; i<mesh.m_tris.size(); ++i){
      fprintf( output, "3 " );
      if(i%38==37 && i!=mesh.m_tris.size()-1) fprintf( output, "\n" );
   }
   fprintf( output, "]\n" );
   fprintf( output, " [ " );
   for(unsigned int i=0; i<mesh.m_tris.size(); ++i){
      fprintf( output, " %d %d %d ", mesh.m_tris[i][0], mesh.m_tris[i][1], mesh.m_tris[i][2] );
      if(i%6==5 && i!=mesh.m_tris.size()-1) fprintf( output, "\n" ); 
   }
   fprintf( output, "]\n" );
   fprintf( output, " \"P\" [" );
   for(unsigned int i=0; i<x.size(); ++i){
      fprintf( output, " %f ", x[i] );
      if(i%4==3 && i!=x.size()-1) fprintf( output, "\n" ); 
   }
   fprintf( output, "]\n" );
   
   return true; 
}


// ---------------------------------------------------------
///
/// Write an RIB file for the shadow map for the given light
///
// ---------------------------------------------------------

bool output_shadow_rib( Gluvi::Target3D& light, const std::vector<Vec3d>& positions, const NonDestructiveTriMesh& mesh, const char *filename_format, ...)
{
   va_list ap;
   va_start(ap, filename_format);
   
#ifdef WIN32
   int len=_vscprintf(filename_format, ap) +1;// _vscprintf doesn't count // terminating '\0'
   char *filename=new char[len];
   vsprintf(filename, filename_format, ap);
#else
   char *filename;
   vasprintf(&filename, filename_format, ap);
#endif
   
   std::ofstream out;
   out.open( filename );
   

   
   if( out == NULL )
   {
      return false;
   }
   
#ifdef WIN32
   len=_vscprintf("track%04d_shadow.tiff", ap) +1;// _vscprintf doesn't count // terminating '\0'
   delete[] filename;
   filename=new char[len];
   vsprintf(filename, "track%04d_shadow.tiff", ap);
#else
   vasprintf(&filename, "track%04d_shadow.tiff", ap);
#endif

   
#ifdef WIN32
   delete [] filename;
#else
   std::free(filename);
#endif

   va_end(ap);
   
   // flatten
   std::vector<float> xs;
   for ( unsigned int i = 0; i < positions.size(); ++i )
   {
      xs.push_back( (float) positions[i][0] );
      xs.push_back( (float) positions[i][1] );
      xs.push_back( (float) positions[i][2] );
   }
   
   out << "Display \"" << filename << "\" \"file\" \"z\"" << std::endl;
   delete[] filename;
   
   // next line: image format (width and height in pixels, pixel aspect ratio)
   out << "Format " << "1024 1024 1" << std::endl;
   out << "PixelFilter \"box\" 1 1 " << std::endl;
   
   // then write out the camera specification
   light.export_rib(out);
   
   // start the scene
   out << "WorldBegin\n";
   
   out << "AttributeBegin\n";
   out << "  Color [0.6 0.6 0.2]\n";
   out << "  Opacity [1 1 1]\n";
   out << "  Surface \"matte\" \"Kd\" 1\n";
   
   const float plane_limit = 50.0f;
   char buf[256];
   sprintf( buf, "  Polygon \"P\" [-%f -%f -%f %f -%f -%f  %f %f -%f  -%f %f -%f]\n", 
           plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit, plane_limit );
   out << buf;
   
   out << "AttributeEnd\n";
   
   out << "Color [0.7 0.7 0.9]\n";
   out << "Surface \"matte\" \n";
   
   write_ribfile( mesh, xs, out);
   
   // finish the scene
   out << "WorldEnd\n";
   
   out.flush();
   out.close();
   
   return true;
}   

// ---------------------------------------------------------
///
/// Write a render-ready RIB file.
///
// ---------------------------------------------------------

bool output_rib( const std::vector<Vec3d>& positions, const NonDestructiveTriMesh& mesh, const char *filename_format, ...)
{
   va_list ap;
   va_start(ap, filename_format);
   
#ifdef WIN32
   int len=_vscprintf(filename_format, ap) +1;// _vscprintf doesn't count // terminating '\0'
   char *filename=new char[len];
   vsprintf(filename, filename_format, ap);
#else
   char *filename;
   vasprintf(&filename, filename_format, ap);
#endif
    
   
   std::ofstream out;
   out.open( filename );
   
   delete[] filename;
   
   if( out == NULL )
   {
      return false;
   }
   
   // first line: what image file this RIB file should produce
#ifdef WIN32
   len=_vscprintf("track%04d.tiff", ap) +1;// _vscprintf doesn't count // terminating '\0'
   delete[] filename;
   filename=new char[len];
   vsprintf(filename, "track%04d.tiff", ap);
#else
   vasprintf(&filename, "track%04d.tiff", ap);
#endif
   

   char *shadow_filename;
#ifdef WIN32
   len=_vscprintf("track%04d_shadow.tiff", ap) +1;// _vscprintf doesn't count // terminating '\0'
   shadow_filename=new char[len];
   vsprintf(shadow_filename, "track%04d_shadow.tiff", ap);
#else
   vasprintf(&shadow_filename, "track%04d_shadow.tiff", ap);
#endif
   
   va_end(ap);
   
   // flatten
   std::vector<float> xs;
   for ( unsigned int i = 0; i < positions.size(); ++i )
   {
      xs.push_back( (float) positions[i][0] );
      xs.push_back( (float) positions[i][1] );
      xs.push_back( (float) positions[i][2] );
   }
   
   std::vector<Vec3f> normals;
   for ( unsigned int i = 0; i < positions.size(); ++i )
   {
      Vec3f n(0,0,0);
      for ( unsigned int j = 0; j < mesh.m_vertex_to_triangle_map[i].size(); ++j )
      {
         Vec3d u = positions[ mesh.m_tris[ mesh.m_vertex_to_triangle_map[i][j] ][1] ] - positions[ mesh.m_tris[ mesh.m_vertex_to_triangle_map[i][j] ][0] ];
         Vec3d v = positions[ mesh.m_tris[ mesh.m_vertex_to_triangle_map[i][j] ][2] ] - positions[ mesh.m_tris[ mesh.m_vertex_to_triangle_map[i][j] ][0] ];
         Vec3d tn = normalized(cross(u, v));          
         n += Vec3f( (float)tn[0], (float)tn[1], (float)tn[2] ); 
      }
      normals.push_back( n / (float) mesh.m_vertex_to_triangle_map[i].size() );
   }
   
   
   out << "Display \"" << filename << "\" \"file\" \"rgb\"" << std::endl;
   delete[] filename;
   
   // next line: image format (width and height in pixels, pixel aspect ratio)
   out << "Format " << Gluvi::winwidth << " " << Gluvi::winheight << " 1" << std::endl;
   out << "PixelSamples 2 2" << std::endl;
   out << "Exposure 1.0 2.2" << std::endl;
   
   // then write out the camera specification
   Gluvi::camera->export_rib(out);
   
   // start the scene
   out << "WorldBegin\n";
   
   out << "LightSource \"ambientlight\" 1 \"intensity\" .3 \"lightcolor\" [1 1 1]\n";
   
   /*
    for ( unsigned int i = 0; i < lights.size(); ++i )
    {
    // compute location of light
    Vec3f light_pos( 0, 0, lights[i].dist );
    rotate( light_pos, lights[i].pitch, Vec3f(1,0,0) );
    rotate( light_pos, lights[i].heading, Vec3f(0,1,0) );
    light_pos += Vec3f( lights[i].target[0], lights[i].target[1], lights[i].target[2] );
    
    std::cout << "light position: " << light_pos << std::endl;
    
    out << "LightSource \"singleshadowpoint\" 2 \"intensity\" 30 \"from\" [" << light_pos << "] \"shadowmap\" \"" << shadow_filename << "\"\n";
    }
    */
   
   //out << "LightSource \"distantlight\" 3 \"intensity\" 0.3 \"from\" [-5 -10 20] \"to\" [0 0 0]\n";
   
   out << "AttributeBegin\n";
   out << "  Color [0.6 0.6 0.2]\n";
   out << "  Opacity [1 1 1]\n";
   out << "  Surface \"matte\" \"Kd\" 1\n";
   
   const float plane_limit = 50.0f;
   const float plane_distance = 10.0f;
   char buf[256];
   sprintf( buf, "  Polygon \"P\" [-%f -%f -%f %f -%f -%f  %f %f -%f  -%f %f -%f]\n", 
           plane_limit, plane_limit, plane_distance, 
           plane_limit, plane_limit, plane_distance, 
           plane_limit, plane_limit, plane_distance, 
           plane_limit, plane_limit, plane_distance );
   out << buf;
   
   out << "AttributeEnd\n";
   
   out << "Color [0.3 0.3 0.9]\n";
   out << "Surface \"matte\" \n";
   
   write_ribfile( mesh, xs, out);
   
   out << " \"N\" [";
   for(unsigned int i=0; i<normals.size(); ++i)
   {
      out << normals[i] << "  ";
   }
   out << "]" << std::endl;
   
   // finish the scene
   out << "WorldEnd\n";
   
   out.flush();
   out.close();
   
   return true;
}

// ---------------------------------------------------------
//
// Write mesh in PBRT format
//
// ---------------------------------------------------------

bool write_pbrtfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, const char *filename_format, ...)
{
   va_list ap;
   va_start(ap, filename_format);
#ifdef WIN32
   int len=_vscprintf(filename_format, ap) +1;// _vscprintf doesn't count // terminating '\0'
   char *filename=new char[len];
   vsprintf(filename, filename_format, ap);
#else
   char *filename;
   vasprintf(&filename, filename_format, ap);
#endif
   std::ofstream output(filename, std::ofstream::binary);
   
#ifdef WIN32
   delete [] filename;
#else
   std::free(filename);
#endif
   
   va_end(ap);

   if(!output.good()) return false;
    return write_ribfile(mesh, x, output);
}

// ---------------------------------------------------------
//
// Write mesh in PBRT format
//
// ---------------------------------------------------------

bool write_pbrtfile(const NonDestructiveTriMesh &mesh, const std::vector<float> &x, std::ostream &output)
{
   output<<"# generated by editmesh"<<std::endl;

   //output<<"\"integer nlevels\" [3]"<<std::endl;
   output<<"\"point P\" ["<<std::endl;
   for(unsigned int i=0; i<x.size(); ++i){
      output<<x[i]<<"  ";
      if(i%4==3 && i!=x.size()-1) output<<std::endl;
   }
   output<<"]"<<std::endl;
   output<<" \"integer indices\" ["<<std::endl;
   for(unsigned int i=0; i<mesh.m_tris.size(); ++i){
      output<<mesh.m_tris[i]<<"  ";
      if(i%6==5 && i!=mesh.m_tris.size()-1) output<<std::endl;
   }
   output<<"]"<<std::endl;

   return output.good();
}





