#include "makelevelset3.h"
#include "marchingtets.h"
#include "array3_utils.h"

#include <fstream>
#include <string>
#include <iostream>
#include <limits>

class ShellMarchingTets : public MarchingTets {
public:
  Array3f& m_data;
  Vec3f m_data_origin;
  float m_data_dx;

  ShellMarchingTets(Array3f& data, const Vec3f& data_origin, float data_dx, 
                                   const Vec3f& mt_origin, float mt_dx):
    m_data(data),m_data_origin(data_origin), m_data_dx(data_dx), MarchingTets(mt_origin, mt_dx) 
    {}

  float eval(int i, int j, int k) {
    //interpolate in the distance field

    //compute query point (using marching tets origin/dx)
    Vec3f pos = origin + Vec3f(i*dx, j*dx, k*dx);
    
    //determine corresponding coordinates in data field (using data origin/dx)
    Vec3f grid_coords = (pos-m_data_origin)/m_data_dx;
    
    //do the lerp
    return interpolate_value(grid_coords, m_data);
  }

  float eval(Vec3f query_pos) {
    //interpolate in the distance field

    //determine corresponding coordinates in data field (using data origin/dx)
    Vec3f grid_coords = (query_pos-m_data_origin)/m_data_dx;

    //do the lerp
    return interpolate_value(grid_coords, m_data);
  }

};

void main(int argc, char* argv[]) {

  if(argc < 2) std::cout << "ERROR: Must provide PLY file as first command line arg.\n";
  std::string input_file(argv[1]);

  float scale_factor = 2.4f;

  std::vector<Vec3f> positions;
  std::vector<Vec3ui> faces;
  std::vector<double> thicknesses;

  float min_thickness =std::numeric_limits<float>::max();
  float max_thickness = 0;
  Vec3f min_bound(std::numeric_limits<float>::max(),std::numeric_limits<float>::max(),std::numeric_limits<float>::max()) , 
    max_bound(-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max(),-std::numeric_limits<float>::max());

  //load ply data
  std::string data;
  
  std::ifstream reader(input_file.c_str());
  if(!reader) {
    std::cout << "File open failed:" << input_file << "\n";
    exit(-1);
  }

  std::cout << "Reading data from input PLY file: " << input_file << std::endl;

  getline(reader, data); //ply
  getline(reader, data); //format
  getline(reader, data); //comment
  
  int num_verts;
  reader >> data >> data >> num_verts;

  std::cout << "Num verts: " << num_verts << std::endl;

  getline(reader, data); //eoln
  getline(reader, data); //x
  getline(reader, data); //y
  getline(reader, data); //z
  getline(reader, data); //thickness

  int num_faces;
  reader >> data >> data >> num_faces;

  std::cout << "Faces: " << num_faces << std::endl;

  getline(reader, data); //eoln
  getline(reader, data); //property list
  getline(reader, data); //end header
   
  positions.resize(num_verts);
  thicknesses.resize(num_verts);
  faces.resize(num_faces);

  std::cout << "Reading vertices.\n";
  for(int i = 0; i < num_verts; ++i) {
    Vec3f cur_pos; 
    float thickness;
    reader >> cur_pos[0] >> cur_pos[1] >> cur_pos[2];
    reader >> thickness;
    positions[i] = cur_pos;
    thicknesses[i] = thickness;
    min_thickness = min(thickness, min_thickness);
    max_thickness = max(max_thickness, thickness);
    update_minmax(cur_pos, min_bound, max_bound);
  }

  std::cout << "Reading faces.\n";
  for(int i = 0; i < num_faces; ++i) {
    Vec3ui cur_pos; 
    int dummy;
    reader >> dummy >> cur_pos[0] >> cur_pos[1] >> cur_pos[2];
    faces[i] = cur_pos;
  }
    
  //make signed distance field
  float dx = min_thickness / scale_factor; 
  //pad a bit
  float increase = 2*max_thickness; 
  min_bound -= increase*Vec3f(1,1,1);
  max_bound += increase*Vec3f(1,1,1);
  Vec3ui resolutions((unsigned int)((max_bound[0]-min_bound[0])/dx), 
                    (unsigned int)((max_bound[1]-min_bound[1])/dx),
                    (unsigned int)((max_bound[2]-min_bound[2])/dx));
  Array3f distance_data;
  int exact_band = 2;
  std::cout << "Lower bound: " << min_bound << std::endl;
  std::cout << "Upper bound: " << max_bound << std::endl;
  std::cout << "Resolution: " << resolutions << std::endl;
  std::cout << "Cell size:" << dx << std::endl;
  std::cout << "Building signed distance field.\n";
  
  //float cull_bound = dx*max(max(resolutions[0],resolutions[1]), resolutions[2]);
  float cull_bound = max_thickness + 3*dx; //far enough away that interpolation at the surface will not use invalid data
  
  make_distance_field3(faces, positions, thicknesses, min_bound, dx, 
                        resolutions[0], resolutions[1], resolutions[2], distance_data, cull_bound, exact_band);

  //mesh with marching tets
  //TODO Use different (higher) marching tets resolution as compared to distance field resolution,
  //for better results.  Requires editing ShellMarchingTets to distinguish interpolation from meshing.
  std::cout << "Building mesh surface.\n";
  ShellMarchingTets shellMarch(distance_data, min_bound, dx, min_bound, dx);

  for(unsigned int k = 0; k < resolutions[2]; ++k) {
    for(unsigned int j = 0; j < resolutions[1]; ++j) {
      for(unsigned int i = 0; i < resolutions[0]; ++i) {
        if(shellMarch.eval(i,j,k) < 2*min_thickness) //don't bother if we're far from the surface
          shellMarch.contour_cube(i,j,k);
      }
    }
  }

  
  std::ofstream out_file("output.OBJ");
  if(!out_file) {
    std::cout << "Failed to open output file for writing.\n";
    exit(-1);
  }
  std::cout << "Writing out corresponding OBJ file." << std::endl;
  std::cout << "Writing " << shellMarch.x.size() << " vertices.\n";
  for(unsigned int i = 0; i < shellMarch.x.size(); ++i) {
    out_file << "v " << shellMarch.x[i][0] << " " << shellMarch.x[i][1] << " " << shellMarch.x[i][2] << std::endl;
  }
  std::cout << "Writing " << shellMarch.mesh.tri.size() << " faces.\n";
  for(unsigned int i = 0; i < shellMarch.mesh.tri.size(); ++i) {
    out_file << "f " << shellMarch.mesh.tri[i][0]+1 << " " << shellMarch.mesh.tri[i][1]+1 << " " << shellMarch.mesh.tri[i][2]+1 << std::endl;
  }

  std::cout << "Thanks for playing!\n";

}