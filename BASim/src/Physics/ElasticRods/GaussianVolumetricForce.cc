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

void GaussianVolumetricForce::setCenter( const Vec3d& center )
{
    m_center = center;
}

void GaussianVolumetricForce::setCovariance( const Mat3d& sigma )
{
    m_invSigma = sigma.inverse();
    m_scaledInvSigma = m_invSigma / ( m_scale * m_scale );
}

void GaussianVolumetricForce::setupSigma( const Eigen::Matrix<Scalar, 3, Eigen::Dynamic>& points )
{
    int n = points.cols();
    m_center = points.rowwise().sum() / n;

    Eigen::Matrix<Scalar, 3, Eigen::Dynamic> centeredPoints( 3, n );
    for ( int i = 0; i < n; ++i )
        centeredPoints.block<3, 1> ( 0, i ) = points.block<3, 1> ( 0, i ) - m_center;

    const Mat3d& Sigma = centeredPoints * centeredPoints.transpose();
    std::cout << "Sigma = " << Sigma << '\n';

    m_invSigma = Sigma.inverse();
    m_scaledInvSigma = m_invSigma / ( m_scale * m_scale );

    m_charge = 100000.0*pow( m_scaledInvSigma.determinant(), 0.5 ); // L^1 scaling

    std::cout << "charge = " << m_charge << '\n';
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
        const Scalar localE = m_charge * exp( -0.5 * ( ( x - m_center ).transpose() * y )[0] );

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
        const Scalar localE = m_charge * exp( -0.5 * ( ( x - m_center ).transpose() * y )[0] );

        energy += localE;
        force.segment<3> ( rod.vertIdx( vidx, 0 ) ) += y * localE;
    }
}

void GaussianVolumetricForce::computeForceDX( int baseindex, const ElasticRod& rod, Scalar scale,
        MatrixBase& J ) const
{
    Mat3d localJ;
    for ( int vidx = 0; vidx < rod.nv(); ++vidx )
    {
        computeLocalForceDX( rod, vidx, localJ );
        localJ *= scale;
        J.pointStencilAdd( rod.vertIdx( vidx, 0 ) + baseindex, localJ );
    }
}

void GaussianVolumetricForce::computeLocalForceDX( const ElasticRod& rod, int vidx, Mat3d& localJ ) const
{
    const Vec3d& x = rod.getVertex( vidx );
    const Vec3d& y = m_scaledInvSigma * ( x - m_center );
    const Scalar localE = m_charge * exp( -0.5 * ( ( x - m_center ).transpose() * y )[0] );
    localJ = localE * ( m_scaledInvSigma - y * y.transpose() );
}

void GaussianVolumetricForce::computeForceDV( int baseindex, const ElasticRod& rod, Scalar scale,
        MatrixBase& J ) const
{
}

}
