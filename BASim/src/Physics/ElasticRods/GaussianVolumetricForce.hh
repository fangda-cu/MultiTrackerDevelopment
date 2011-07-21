/*
 * GaussianVolumetricForce.hh
 *
 *  Created on: 21/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef GAUSSIANVOLUMETRICFORCE_HH_
#define GAUSSIANVOLUMETRICFORCE_HH_

#include "RodExternalConservativeForce.hh"

namespace BASim
{

/**
 * A conservative force with potential E(x) = charge * exp(-1/2 ( x - center)^T Sigma^{-1} (x - center) /(scale^2))
 */
class GaussianVolumetricForce: public BASim::RodExternalConservativeForce
{
public:
    GaussianVolumetricForce( const Scalar charge = 1.0, const Scalar scale = 1.0,
            const Vec3d& center = Vec3d(), const Mat3d& invSigma = Mat3d::Identity() );

    virtual ~GaussianVolumetricForce();

    void setCharge( const Scalar charge );

    void setScale( const Scalar scale );

    // Compute mean and variance of the given points and assign its result to m_center and m_invSigma (and m_scaledInvSigma)
    void setupSigma( const Eigen::Matrix<Scalar, 3, Eigen::Dynamic>& points );

    virtual Scalar computeEnergy( const ElasticRod& rod ) const;

    virtual void computeForce(const ElasticRod& rod, VecXd& force) const;

    virtual void computeForceEnergy( const ElasticRod& rod, VecXd& force, Scalar& energy ) const;

    virtual void
            computeForceDX( int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J ) const;

    virtual void
            computeForceDV( int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J ) const;

private:
    Scalar m_charge;
    Scalar m_scale;
    Vec3d m_center;
    Mat3d m_invSigma;
    Mat3d m_scaledInvSigma;
};

}

#endif /* GAUSSIANVOLUMETRICFORCE_HH_ */
