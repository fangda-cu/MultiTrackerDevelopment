// RodRodPenaltyForce.hh
//

#ifndef RODRodPenaltyForce_HH
#define RODRodPenaltyForce_HH

#include "ElasticRod.hh"
#include "RodExternalForce.hh"
#include "../../Collisions/BVHAABB.hh"

//#include <ext/hash_map>
#include <map>

namespace BASim {

class RodPenaltyForce : public RodExternalForce
{
public:
  RodPenaltyForce();
  ~RodPenaltyForce();

	void clearCollisions();
	
  //void computeEnergy(Real& e)

  virtual void computeForce(const ElasticRod& rod, VecXd& F);

  void addRodPenaltyForce(int vertex, VertexFaceImplicitPenaltyCollision cllsn);

  void clearPenaltyForces();

  virtual void computeForceDX(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J);
  virtual void computeForceDV(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) {}

protected:

//  double getClosestPointsVertexTriangle(const Vec3d& v0, const Vec3d& v1,
//                                      const Vec3d& v2, const Vec3d& v3,
//                                      double &t1, double & &t2, double & &t3) const;
                                
  void localJacobian(MatXd& J, const Scalar stiffness, const Vec3d& normal);

	std::vector <int> vidx;
	std::vector <VertexFaceImplicitPenaltyCollision> vertex_face_collisions;

};

}

#endif
