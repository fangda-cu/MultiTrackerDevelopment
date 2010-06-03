// RodRodPenaltyForce.hh
//

#ifndef RODRodPenaltyForce_HH
#define RODRodPenaltyForce_HH

#include "ElasticRod.hh"
#include "RodExternalForce.hh"
#include "../../Collisions/CollisionMeshData.hh"

//#include <ext/hash_map>
#include <map>

namespace BASim {

//typedef __gnu_cxx::hash_map<int, std::pair<CollisionMeshData *, int> > VertexObjectMap;
typedef std::map<int, std::pair<CollisionMeshData *, int> > VertexObjectMap;
typedef VertexObjectMap::iterator VertexObjectMapIterator;

//typedef __gnu_cxx::hash_map<int, std::pair<ElasticRod *, int> > EdgeRodMap;
typedef std::map<int, std::pair<ElasticRod *, int> > EdgeRodMap;
typedef EdgeRodMap::iterator EdgeRodMapIterator;

typedef std::map<int, std::pair<ElasticRod *, int> > VertexRodMap;
typedef EdgeRodMap::iterator VertexRodMapIterator;

typedef std::vector<Scalar> RealArray;

class RodPenaltyForce : public RodExternalForce
{
public:
  RodPenaltyForce();
  ~RodPenaltyForce();

	void clearCollisions();
	
  //void computeEnergy(Real& e)

  virtual void computeForce(const ElasticRod& rod, VecXd& F);

  virtual void computeForceDX(const ElasticRod& rod, MatrixBase& J);
  virtual void computeForceDV(const ElasticRod& rod, MatrixBase& J);
  
  void addRodPenaltyForce(int vertex, CollisionMeshData *cmData, int triangle);
  void addRodPenaltyForce(int edge, ElasticRod *rod, int otherEdge);

  void addRodClumpingForce(int vertex, ElasticRod *rod, int otherVertex);

  void clearPenaltyForces();
	
  void setClumping(bool flag, Real coeff = 0.0) 
  {
    clumping_enbld = flag;
    clumping_coeff = coeff;
  }
	
	
protected:
  VertexObjectMap _vertexObjects;
  EdgeRodMap _edgeRods;
  
  VertexRodMap _clumpingVerts;

  Real getClosestPointsVertexTriangle(const Vec3d& v0, const Vec3d& v1,
                                      const Vec3d& v2, const Vec3d& v3,
                                      Real &t1, Real &t2, Real &t3) const;

  Real getClosestPointsEdgeEdge(const Vec3d& e11, const Vec3d& e12,
                                const Vec3d& e21, const Vec3d& e22,
                                Real &s, Real &t) const;
                                
  void localJacobian(MatXd& J, const Scalar stiffness, const Vec3d& normal);
	
  bool clumping_enbld;
  Real clumping_coeff;
                                
};

}

#endif
