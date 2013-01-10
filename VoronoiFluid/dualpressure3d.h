#ifndef DUALPRESSURE3D_H
#define DUALPRESSURE3D_H

#include <vector>

class TetMesh;
namespace ElTopo{
class DynamicSurface;
}

//single phase flow
std::vector<double> pressure_solve_voronoi( TetMesh& mesh, 
                                            ElTopo::DynamicSurface& surface,
                                            float surfaceTensionCoeff,
                                            std::vector<float>& face_velocities,      
                                            const std::vector<float>& solid_weights,  
                                            const std::vector<float>& liquid_phi,     
                                            const std::vector<float>& wall_velocities);

//multiphase version of the above
std::vector<double> pressure_solve_multi( TetMesh& mesh, 
  ElTopo::DynamicSurface& surface,
  float surfaceTensionCoeff,
  std::vector<float>& face_velocities,      
  const std::vector<float>& solid_weights,  
  const std::vector<float>& liquid_phi,
  const std::vector<int>& regions,
  const std::vector<float>& densities);


#endif

