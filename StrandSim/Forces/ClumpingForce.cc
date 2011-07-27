/*
 * ClumpingForce.cc
 *
 *  Created on: 27/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#include "ClumpingForce.hh"
#include "../ElasticStrand.hh"

namespace strandsim
{

ClumpingForce::ClumpingForce() :
    m_charge( 0.0 ), m_power( 1.0 ), m_rho2( 0.000025 )
{
    // TODO Auto-generated constructor stub

}

ClumpingForce::~ClumpingForce()
{
    // TODO Auto-generated destructor stub
}

void ClumpingForce::accumulateEF( StrandGeometry& geometry, const ElasticStrand& strand ) const
{
    const std::list<ElasticStrand*>& attractors = strand.getClumpingAttractors();

    for ( std::list<ElasticStrand*>::const_iterator attractor = attractors.begin(); attractor
            != attractors.end(); ++attractor )
        accumulateMutualEF( geometry.m_totalEnergy, geometry.m_totalForce, strand, *( *attractor ) );
}

void ClumpingForce::accumulateEFJ( StrandGeometry& geometry, const ElasticStrand& strand ) const
{
    const std::list<ElasticStrand*>& attractors = strand.getClumpingAttractors();

    for ( std::list<ElasticStrand*>::const_iterator attractor = attractors.begin(); attractor
            != attractors.end(); ++attractor )
        accumulateMutualEFJ( geometry.m_totalEnergy, geometry.m_totalForce,
                *( geometry.m_totalJacobian ), strand, *( *attractor ) );
}

void ClumpingForce::accumulateJ( StrandGeometry& geometry, const ElasticStrand& strand ) const
{
    const std::list<ElasticStrand*>& attractors = strand.getClumpingAttractors();

    for ( std::list<ElasticStrand*>::const_iterator attractor = attractors.begin(); attractor
            != attractors.end(); ++attractor )
        accumulateMutualJ( *( geometry.m_totalJacobian ), strand, *( *attractor ) );
}

void ClumpingForce::accumulateMutualEF( Scalar& energy, VecXd& force, const ElasticStrand& strand,
        const ElasticStrand& other ) const
{
    Scalar chargeMultiplier = 1.0;

    for ( int vtx = 0; vtx < strand.m_numVertices; ++vtx )
    {
        if ( vtx < m_vertexQMap.size() )
            chargeMultiplier = m_vertexQMap[vtx];

        const Vec3d& x = strand.getVertex( vtx );
        const Vec3d& y = other.closestPoint( x );
        const Scalar normep2 = ( x[0] - y[0] ) * ( x[0] - y[0] ) + ( x[1] - y[1] ) * ( x[1] - y[1] )
                + ( x[2] - y[2] ) * ( x[2] - y[2] );

        energy -= m_charge * chargeMultiplier * normep2 / pow( normep2 + m_rho2, 1 + 0.5 * m_power );

        force.segment<3> ( 4 * vtx ) += m_charge * chargeMultiplier * ( x - y ) * pow(
                normep2 + m_rho2, -1 - 0.5 * m_power ) * ( 2 - ( 2 + m_power ) * normep2
                / ( normep2 + m_rho2 ) );
    }
}

void ClumpingForce::accumulateMutualEFJ( Scalar& energy, VecXd& force,
        JacobianMatrixType& Jacobian, const ElasticStrand& strand, const ElasticStrand& other ) const
{
    Scalar chargeMultiplier = 1.0;

    for ( int vtx = 0; vtx < strand.m_numVertices; ++vtx )
    {
        if ( vtx < m_vertexQMap.size() )
            chargeMultiplier = m_vertexQMap[vtx];

        const Vec3d& x = strand.getVertex( vtx );
        const Vec3d& y = other.closestPoint( x );
        const Scalar normep2 = ( x[0] - y[0] ) * ( x[0] - y[0] ) + ( x[1] - y[1] ) * ( x[1] - y[1] )
                + ( x[2] - y[2] ) * ( x[2] - y[2] );

        energy -= m_charge * chargeMultiplier * normep2 / pow( normep2 + m_rho2, 1 + 0.5 * m_power );

        force.segment<3> ( 4 * vtx ) += m_charge * chargeMultiplier * ( x - y ) * pow(
                normep2 + m_rho2, -1 - 0.5 * m_power ) * ( 2 - ( 2 + m_power ) * normep2
                / ( normep2 + m_rho2 ) );

        Mat3d localJ;
        computeLocalJacobian( localJ, x, y, normep2, chargeMultiplier );
        Jacobian.localStencilAdd<3> ( 4 * vtx, localJ );
    }

}

void ClumpingForce::accumulateMutualJ( JacobianMatrixType& Jacobian, const ElasticStrand& strand,
        const ElasticStrand& other ) const
{
    Scalar chargeMultiplier = 1.0;

    for ( int vtx = 0; vtx < strand.m_numVertices; ++vtx )
    {
        if ( vtx < m_vertexQMap.size() )
            chargeMultiplier = m_vertexQMap[vtx];

        const Vec3d& x = strand.getVertex( vtx );
        const Vec3d& y = other.closestPoint( x );
        const Scalar normep2 = ( x[0] - y[0] ) * ( x[0] - y[0] ) + ( x[1] - y[1] ) * ( x[1] - y[1] )
                + ( x[2] - y[2] ) * ( x[2] - y[2] );

        Mat3d localJ;
        computeLocalJacobian( localJ, x, y, normep2, chargeMultiplier );
        Jacobian.localStencilAdd<3> ( 4 * vtx, localJ );
    }

}

void ClumpingForce::computeLocalJacobian( Mat3d& localJ, const Vec3d& x, const Vec3d& y,
        const Scalar normep2, const Scalar chargeMultiplier ) const
{
    const Scalar normep4 = normep2 * normep2;

    for ( int i = 0; i < 3; i++ )
    {
        const double xx2 = ( x[i] - y[i] ) * ( x[i] - y[i] );
        localJ( i, i ) = m_charge * chargeMultiplier * pow( normep2 + m_rho2, -3 - m_power / 2. )
                * ( 2 * m_rho2 * ( -2 * ( 2 + m_power ) * xx2 + m_rho2 ) + normep2 * ( 2 * m_power
                        * xx2 + m_power * m_power * xx2 + 2 * m_rho2 - m_power * m_rho2 ) - normep4
                        * m_power );

        for ( int j = 0; j < i; j++ )
            localJ( i, j ) = localJ( j, i ) = m_charge * chargeMultiplier * ( 2 + m_power )
                    * ( x[i] - y[i] ) * ( x[j] - y[j] ) * ( normep2 * m_power - 4 * m_rho2 ) * pow(
                    normep2 + m_rho2, -3 - m_power / 2. );
    }
}

}
