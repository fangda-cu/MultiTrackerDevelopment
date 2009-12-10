#ifndef RODTREE_H
#define RODTREE_H

#include <BASim/Physics>
//#include <BASim/src/Physics/ElasticRods/ElasticRod.hh>

//#include <vector>

using namespace std;

namespace BASim {

class RodTree 
{
public:
  RodTree(ElasticRod* parentRod, ElasticRod* rod, Vec3d rootPoint=Vec3d(0,0,0));
  ~RodTree();
  void addChildNode(RodTree* child);
  vector<RodTree*>& getChildNodes();
  void enforceParentConnections();
  void draw();
  void getRodMaterialFrame( ElasticRod* rod, const int edge, Vec3d& m1, Vec3d& m2, Vec3d&t );
  
private:
  ElasticRod* m_rod;
  ElasticRod* m_parentRod;
  Vec3d m_rootPoint;
  vector<RodTree*> m_childNodes;

  Vec3d m_material1;
  Vec3d m_material2;
  Vec3d m_tangent;

  Quaternion m_material1Rot;
  Quaternion m_material2Rot;
  Quaternion m_rotationFromParentFrame;

  Vec3d testAxis;
  Quaternion testRot;
};

}

#endif
