// RodLevelSetForce.hh
//

#ifndef RodLevelSetForce_HH
#define RodLevelSetForce_HH

#include "ElasticRod.hh"
#include "RodExternalForce.hh"

namespace BASim
{

class RodLevelSetForce: public RodExternalForce
{
public:
    RodLevelSetForce(double LevelSetThicknessFraction = 0.1);
    ~RodLevelSetForce();

    virtual void computeForce(const ElasticRod& rod, VecXd& F) const;

    virtual void computeForceDX(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) const;

    virtual void computeForceDV(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) const
    {
    }

protected:
    void localJacobian(MatXd& J, const Scalar stiffness, const Vec3d& normal) const;

};

}

#endif
