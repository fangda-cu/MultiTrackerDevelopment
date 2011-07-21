/*
 * GaussianVolumetricForce.cc
 *
 *  Created on: 21/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#include "GaussianVolumetricForce.hh"
#include "ElasticRod.hh"

namespace BASim
{

GaussianVolumetricForce::GaussianVolumetricForce( const Scalar charge, const Scalar scale,
        const Vec3d& center, const Mat3d& invSigma ) :
    m_charge( charge ), m_scale( scale ), m_center( center ), m_invSigma( invSigma )
{
    m_scaledInvSigma = invSigma / ( scale * scale );
}

GaussianVolumetricForce::~GaussianVolumetricForce()
{
}

void GaussianVolumetricForce::setScale( const Scalar scale )
{
    m_scale = scale;
    m_scaledInvSigma = m_invSigma / ( m_scale * m_scale );
}

void GaussianVolumetricForce::setCharge( const Scalar charge )
{
    m_charge = charge;
}

void GaussianVolumetricForce::setupSigma( const Eigen::Matrix<Scalar, 3, Eigen::Dynamic>& points )
{
    int n = points.cols();
    m_center = points.rowwise().sum() / n;

    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> centeredPoints( 3, n );
    for ( int i = 0; i < n; ++i )
        centeredPoints.block<3, 1> ( 0, i ) = points.block<3, 1> ( 0, i ) - m_center;

    m_invSigma = ( centeredPoints * centeredPoints.transpose() ).inverse();
    m_scaledInvSigma = m_invSigma / ( m_scale * m_scale );
}

Scalar GaussianVolumetricForce::computeEnergy( const ElasticRod& rod ) const
{
    Scalar energy = 0.0;

    for ( int vidx = 0; vidx < rod.nv(); ++vidx )
    {
        const Vec3d& x = rod.getVertex( vidx );

        energy += exp(
                -0.5 * ( ( x - m_center ).transpose() * m_scaledInvSigma * ( x - m_center ) )[0] );
    }

    return m_charge * energy;

}
void GaussianVolumetricForce::computeForce( const ElasticRod& rod, VecXd& force ) const
{
    for ( int vidx = 0; vidx < rod.nv(); ++vidx )
    {
        const Vec3d& x = rod.getVertex( vidx );
        const Vec3d& y = m_scaledInvSigma * ( x - m_center );
        Scalar localE = m_charge * exp( -0.5 * ( ( x - m_center ).transpose() * y )[0] );

        force.segment<3> ( rod.vertIdx( vidx, 0 ) ) += y * localE;
    }
}

void GaussianVolumetricForce::computeForceEnergy( const ElasticRod& rod, VecXd& force,
        Scalar& energy ) const
{
    for ( int vidx = 0; vidx < rod.nv(); ++vidx )
    {
        const Vec3d& x = rod.getVertex( vidx );
        const Vec3d& y = m_scaledInvSigma * ( x - m_center );
        Scalar localE = m_charge * exp( -0.5 * ( ( x - m_center ).transpose() * y )[0] );

        energy += localE;
        force.segment<3> ( rod.vertIdx( vidx, 0 ) ) += y * localE;
    }
}

void GaussianVolumetricForce::computeForceDX( int baseindex, const ElasticRod& rod, Scalar scale,
        MatrixBase& J ) const
{
}

void GaussianVolumetricForce::computeForceDV( int baseindex, const ElasticRod& rod, Scalar scale,
        MatrixBase& J ) const
{
}

}
