// RodLevelSetForce.hh
//


#ifndef RodLevelSetForce_HH
#define RodLevelSetForce_HH

#include "ElasticRod.hh"
#include "RodExternalConservativeForce.hh"

namespace BASim
{

class LevelSet;


class RodLevelSetForce: public RodExternalConservativeForce
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

    void setStiffness( Scalar stiffness )
    {
	m_stiffness = stiffness;
    }

protected:
    void localJacobian(MatXd& J, const Scalar stiffness, const Vec3d& normal) const;

    LevelSet* m_levelSet;

    Scalar m_stiffness;
};

}

#endif
