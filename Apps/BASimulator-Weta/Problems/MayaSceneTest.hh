/**
 * \file HandTest.hh
 *
 * \author 
 * \date 
 */

#ifndef MayaSceneTest_HH
#define MayaSceneTest_HH

//#define COLUMBIA_COLLISION

#include "ProblemBase.hh"
#include "BASim/src/Core/TriangleMesh.hh"
#include "BASim/src/IO/ObjParser.hh"

#include "BASim/src/Physics/ElasticRods/BridsonStepper.hh"
#include "BASim/src/Physics/ElasticRods/AdaptiveBinaryStepper.hh"

#ifndef COLUMBIA_COLLISION
#include "BASim/src/Physics/ElasticRods/RodCollisionTimeStepper.hh"
#include <BASim/src/Collisions/CollisionMeshData.hh>
#endif

#include "BASim/src/Physics/ElasticRods/RodPenaltyForce.hh"

// 


/**
 * 
 */
#ifdef COLUMBIA_COLLISION
class MeshController : public ObjectControllerBase
{
  public:
    MeshController( TriangleMesh* mesh );

    bool execute();
    
    void updateMesh(  std::vector<Vec3d> & prevPosition, std::vector<Vec3d> & nextMeshPosition, double s );
    
    TriangleMesh* m_mesh;
    std::vector<Vec3d> m_verts;

};
#endif


class MayaSceneTest : public Problem
{
public:

  MayaSceneTest();
  virtual ~MayaSceneTest();

  void reverse_test();
  void reverse_timestep();
  void reverse_solve();
  void reverse_findudf(uint id);
  
  void reverse_findudf_list();
  void reverse_setup_list();
  
  bool impulse_enabled;
  
//  RodGroupManager *rodGroupManager;
  
protected:

  void Setup();
  void AtEachTimestep();

  FILE* readNumberOfRodsFromFile( const std::string i_cacheFilename, size_t& o_numRodsInFile);

private:
  void LoadNextFrameFromCache();
  void InitializeMeshes();
  void ReadMeshFile(int frame, size_t &numberOfVerts, size_t &numberOfFaces, std::vector<Vec3d> &verts, std::vector<uint> &faces);
  void UpdateCollisionMesh( std::vector<Vec3d> &points, std::vector<uint> &faces );
  void MakeTriangleMesh( std::vector<Vec3d> &points, std::vector<uint> &faces );
  void UpdateTriangleMesh( std::vector<Vec3d> &points );
  
  void makeMultipleRodStepper();
  
#ifdef COLUMBIA_COLLISION
#else COLUMBIA_COLLISION
  void initialiseCollisionMesh( BASim::CollisionMeshData *collisionMeshData, size_t id );
#endif COLUMBIA_COLLISION

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

#ifdef COLUMBIA_COLLISION
  BridsonStepper* m_br_stepper;
  std::vector<MeshController*> m_mesh_controllers;
  std::vector<ObjectControllerBase*> m_obj_controllers;
  
  std::vector<Vec3d> prevMeshPosition;
  std::vector<Vec3d> nextMeshPosition;
  
#else
  CollisionMeshDataHashMap m_collisionMeshMap;
  CollisionMeshData* m_collisionMeshData;
  std::vector<RodCollisionTimeStepper*> m_collision_controllers;
#endif
      
  std::vector<ElasticRod*> m_rods;
//  std::vector<RodTimeStepper*> m_steppers;
  std::vector<RodTimeStepper*> m_controllers;

  std::vector<TriangleMesh*> m_tri_objs;
  TriangleMesh* br_mesh;
  
  MultipleRodTimeStepper* multi_stepper;
  bool is_multi_stepper;
  
//  RodGroupSpringForce* group_spring;
  
  std::vector<int> nvs;
  
  //HandController* hand_controller;
 
};

#endif // HandTest_HH
