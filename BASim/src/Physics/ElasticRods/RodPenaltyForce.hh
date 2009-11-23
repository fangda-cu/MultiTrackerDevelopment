// RodRodPenaltyForce.hh
//

#ifndef RODRodPenaltyForce_HH
#define RODRodPenaltyForce_HH

#include <BASim/BASim>
#include "ElasticRod.hh"
#include "RodExternalForce.hh"

#include <ext/hash_map>

namespace BASim {

typedef __gnu_cxx::hash_map<int, std::pair<CollisionMeshData *, int> > VertexObjectMap;
typedef __gnu_cxx::hash_map<int, std::pair<CollisionMeshData *, int> >::iterator VertexObjectMapIterator;

typedef __gnu_cxx::hash_map<int, std::pair<ElasticRod *, int> > EdgeRodMap;
typedef __gnu_cxx::hash_map<int, std::pair<ElasticRod *, int> >::iterator EdgeRodMapIterator;

typedef std::vector<Scalar> RealArray;

class RodPenaltyForce : public RodExternalForce
{
public:
    RodPenaltyForce();
    ~RodPenaltyForce();

    //void computeEnergy(Real& e)
    
    virtual void computeForce(const ElasticRod& rod, VecXd& F);
    virtual void computeForceDX(const ElasticRod& rod, MatrixBase& J);
    
    void addRodPenaltyForce(int vertex, CollisionMeshData *cmData, int triangle);
    void addRodPenaltyForce(int edge, ElasticRod *rod, int otherEdge);

protected:
    VertexObjectMap _vertexObjects;
    EdgeRodMap _edgeRods;

    Real getClosestPointsVertexTriangle(const Vec3d& v0, const Vec3d& v1,
                                        const Vec3d& v2, const Vec3d& v3,
                                        Real &t1, Real &t2, Real &t3) const;

    Real getClosestPointsEdgeEdge(const Vec3d& e11, const Vec3d& e12,
                                  const Vec3d& e21, const Vec3d& e22,
                                  Real &s, Real &t) const;
};

}

#endif

