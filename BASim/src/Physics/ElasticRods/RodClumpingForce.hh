/*
 * RodClumpingForce.hh
 *
 *  Created on: 30/06/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef RODCLUMPINGFORCE_HH_
#define RODCLUMPINGFORCE_HH_

#include "RodExternalConservativeForce.hh"
#include "ElasticRod.hh"

namespace BASim
{

class RodClumpingForce: public BASim::RodExternalConservativeForce
{
public:
    RodClumpingForce();
    virtual ~RodClumpingForce();

    virtual Scalar computeEnergy(const ElasticRod& rod) const;

    virtual void computeForce(const ElasticRod& rod, VecXd& force) const;

    virtual void computeForceEnergy(const ElasticRod& rod, VecXd& force, Scalar& energy) const;

    virtual void computeForceDX(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) const;

    virtual void computeForceDV(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) const;

    void setCharge(Scalar charge)
    {
        q = charge;
    }

    void setPower(Scalar power)
    {
        r = power;
    }

    /*void setVertexPowerMap(Scalar power)
    {
        r = power;
    }*/

    void setDistance(Scalar dist)
    {
        rho2 = dist*dist;
    }

    const double getCharge() const
    {
        return q;
    }

    const double getPower() const
    {
        return r;
    }

    const double getDistance() const
    {
        return sqrt(rho2);
    }

private:
    Scalar computeEnergy(const ElasticRod& rod, const ElasticRod& other) const;

    void computeForce(const ElasticRod& rod, const ElasticRod& other, VecXd& force) const;

    void computeForceDX(int baseindex, const ElasticRod& rod, const ElasticRod& other, Scalar scale, MatrixBase& J) const;

    void computeLocalForceDX(const ElasticRod& rod, int vidx, const ElasticRod& other, Mat3d& localJ) const;

    double q;
    double r;
    double rho2;
    double vertexPowerMap [];
};

}

#endif /* RODCLUMPINGFORCE_HH_ */
