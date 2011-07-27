/*
 * ClumpingForce.hh
 *
 *  Created on: 27/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef CLUMPINGFORCE_HH_
#define CLUMPINGFORCE_HH_

#include "ForceBase.hh"

namespace strandsim
{

class ClumpingForce: public strandsim::ForceBase
{
public:
    ClumpingForce();

    virtual ~ClumpingForce();

    virtual std::string getName()
    {
        return "clumping";
    }

    virtual void accumulateEFJ( StrandGeometry& geometry, const ElasticStrand& strand ) const;

    virtual void accumulateEF( StrandGeometry& geometry, const ElasticStrand& strand ) const;

    virtual void accumulateJ( StrandGeometry& geometry, const ElasticStrand& strand ) const;

private:
    void accumulateMutualEF( Scalar& energy, VecXd& force, const ElasticStrand& strand,
            const ElasticStrand& other ) const;

    void accumulateMutualEFJ( Scalar& energy, VecXd& force, JacobianMatrixType& Jacobian,
            const ElasticStrand& strand, const ElasticStrand& other ) const;

    void accumulateMutualJ( JacobianMatrixType& Jacobian, const ElasticStrand& strand,
            const ElasticStrand& other ) const;

    void computeLocalJacobian( Mat3d& localJ, const Vec3d& x, const Vec3d& y, const Scalar normep2,
            const Scalar charge ) const;

    Scalar m_charge;
    Scalar m_power;
    Scalar m_rho2;
    std::vector<Scalar> m_vertexQMap;

};

}

#endif /* CLUMPINGFORCE_HH_ */
