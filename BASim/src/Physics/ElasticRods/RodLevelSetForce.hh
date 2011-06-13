// RodRodPenaltyForce.hh
//

#ifndef RodPenaltyForce_HH
#define RodPenaltyForce_HH

#include "ElasticRod.hh"
#include "RodExternalForce.hh"
//#include "../../Collisions/BVHAABB.hh"
#include "../../Collisions/Collision.hh"

//#include <ext/hash_map>
#include <map>

namespace BASim
{

class RodPenaltyForce: public RodExternalForce
{
public:
    RodPenaltyForce(double penaltyThicknessFraction = 0.1);
    ~RodPenaltyForce();

    virtual void computeForce(const ElasticRod& rod, VecXd& F) const;

    void registerProximityCollision(int vertex, VertexFaceProximityCollision* vfpcol);

    bool cleared() const;

    void clearProximityCollisions();

    virtual void computeForceDX(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) const;

    virtual void computeForceDV(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) const
    {
    }

protected:
    void localJacobian(MatXd& J, const Scalar stiffness, const Vec3d& normal) const;

    const double m_penaltyThicknessFraction; // fraction of proximity threshold to use for penalty thickness

    std::vector<int> vidx;

    std::vector<const VertexFaceProximityCollision*> vertex_face_collisions;

};

}

#endif
