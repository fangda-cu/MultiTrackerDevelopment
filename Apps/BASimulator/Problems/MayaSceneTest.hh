/**
 * \file HandTest.hh
 *
 * \author 
 * \date 
 */

#ifndef MayaSceneTest_HH
#define MayaSceneTest_HH

#include "ProblemBase.hh"
#include "BASim/src/Core/TriangleMesh.hh"
#include "BASim/src/IO/ObjParser.hh"

#include "BASim/src/Physics/ElasticRods/BridsonStepper.hh"
#include "BASim/src/Physics/ElasticRods/AdaptiveBinaryStepper.hh"
#include "BASim/src/Physics/ElasticRods/RodCollisionTimeStepper.hh"

#include <BASim/src/Collisions/CollisionMeshData.hh>

/**
 * 
 */


class MayaSceneTest : public Problem
{
  enum TimeStepperType { WETA, COLUMBIA };
      
public:

  MayaSceneTest();
  virtual ~MayaSceneTest();

  void reverse_test();
  void reverse_timestep();
  void reverse_solve();
  
  bool impulse_enabled;
  
protected:

  void Setup();
  void AtEachTimestep();

  FILE* readNumberOfRodsFromFile( const std::string i_cacheFilename, size_t& o_numRodsInFile);

private:
  
  TimeStepperType collision_method;
  
  void LoadNextFrameFromCache();
  void InitializeMeshes();
  void ReadMeshFile(int frame, size_t &numberOfVerts, size_t &numberOfFaces, std::vector<Vec3d> &verts, std::vector<uint> &faces);
  void UpdateCollisionMesh( std::vector<Vec3d> &points, std::vector<uint> &faces );
  void initialiseCollisionMesh( BASim::CollisionMeshData *collisionMeshData, size_t id );
  void MakeTriangleMesh( std::vector<Vec3d> &points, std::vector<uint> &faces );
  void UpdateTriangleMesh( std::vector<Vec3d> &points );
  
  int nSubsteps;
  int current_frame;
  int current_step;
  int fps;
  
  int mesh_frame;
  
  bool i_collisionsEnabled;
  
  bool is_selected_running;
  std::vector <int> selected_rods;
  
  std::vector<VecXd> prevBoundaryCondition;
  std::vector<VecXd> nextBoundaryCondition;

  std::vector<VecXd> prevMayaPosition;
  std::vector<VecXd> nextMayaPosition;
  
  CollisionMeshDataHashMap m_collisionMeshMap;
  CollisionMeshData* m_collisionMeshData;
      
  std::vector<ElasticRod*> m_rods;
  std::vector<RodTimeStepper*> m_steppers;
  std::vector<RodTimeStepper*> m_controllers;
//  BridsonStepper* m_br_stepper;
  
  std::vector<RodCollisionTimeStepper*> m_collision_controllers;

  std::vector<TriangleMesh*> m_tri_objs;
  TriangleMesh* br_mesh;
  
  //HandController* hand_controller;
 
};

#endif // HandTest_HH
