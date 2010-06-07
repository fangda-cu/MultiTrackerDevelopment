/**
 * \file HairyBall.hh
 *
 * \author smith@cs.columbia.edu
 * \date 05/11/2010
 */

#ifndef HAIRYBALL_HH
#define HAIRYBALL_HH

#include "ProblemBase.hh"
#include "BASim/src/Physics/ElasticRods/BridsonStepper.hh"
#include "BASim/src/Physics/ElasticRods/ParallelStepper.hh"
#include "BASim/src/Physics/ElasticRods/RodForce.hh"
#include <vector>
#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"
#include "BASim/src/Core/TriangleMesh.hh"
#include "BASim/src/IO/ObjParser.hh"

#include <fstream>

class HairyBall : public Problem
{
public:
  
  HairyBall();
  virtual ~HairyBall();
  
protected:
  
  void Setup();
  void AtEachTimestep();
  void AtEachTimestepOLD();
  
private:
  
  void loadHandCodedVectors( std::vector<Vec3d>& normals );
  void loadNormalsFile( const std::string& filename, std::vector<Vec3d>& normals );
  
  // These assume the objects are cenetered at the origin. Easy to relax this, I just haven't.
  void transformTriangleObject( TriangleMesh* tri_mesh, const Mat3d& transformation );
  void transformRodRoot( ElasticRod* rod, const Mat3d& transformation );

  void generateStraightHair( const Vec3d& initnorm, const Vec3d& startpoint, const double& dL, const int& nv, std::vector<Vec3d>& vertices );
  void generateCurlyHair( const Vec3d& initnorm, const Vec3d& startpoint, const double& dL, const int& nv, std::vector<Vec3d>& vertices );

  std::vector<TriangleMesh*> m_tri_meshes;
  std::vector<ElasticRod*> m_rods;
  std::vector<RodTimeStepper*> m_steppers;
  ParallelStepper* m_pr_stepper;
  
};

#endif // HAIRYBALL_HH
