// RodLevelSetForce.hh
//


#ifndef RodLevelSetForce_HH
#define RodLevelSetForce_HH

#include "ElasticRod.hh"
#include "RodExternalForce.hh"

namespace BASim
{

class LevelSet;


class RodLevelSetForce: public RodExternalForce
{
public:
    RodLevelSetForce( LevelSet* levelSet );

    ~RodLevelSetForce();

    virtual void computeForce(const ElasticRod& rod, VecXd& F) const;

    virtual void computeForceDX(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) const;

    virtual void computeForceDV(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) const
    {
    }

    virtual Scalar computeEnergy(const ElasticRod& rod) const;
 
    virtual void computeForceEnergy(const ElasticRod& rod, VecXd& force, Scalar& energy) const;

protected:
    void localJacobian(MatXd& J, const Scalar stiffness, const Vec3d& normal) const;

    LevelSet* m_levelSet;

};

}

#endif
