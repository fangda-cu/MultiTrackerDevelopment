/*
 * BendingForce.cc
 *
 *  Created on: 12/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#include "BendingForce.hh"

namespace strandsim
{

BendingForce::BendingForce()
{
    // TODO Auto-generated constructor stub

}

BendingForce::~BendingForce()
{
    // TODO Auto-generated destructor stub
}

Scalar BendingForce::localEnergy( const ElasticStrand& strand, const StrandGeometry& geometry,
        const IndexType vtx )
{
    const Mat2d& B = strand.m_bendingMatrix;// strand.m_bendingMatrices[vtx];// std::cout << "B = " << B <<'\n';
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];// std::cout << "ilen = " << ilen <<'\n';
    const Vec2d& kappa = geometry.m_kappa[vtx];// std::cout << "kappa = " << kappa <<'\n';
    const Vec2d& kappaBar = strand.m_restBends[vtx];// std::cout << "kappaBar = " << kappaBar <<'\n';

    return 0.5 * ilen * ( kappa - kappaBar ).dot( B * ( kappa - kappaBar ) );
}

void BendingForce::computeLocalForce( BendingForce::LocalForceType& localF,
        const ElasticStrand& strand, const StrandGeometry& geometry, const IndexType vtx )
{
    // std::cout << "Local bending force (vertex " << vtx << ")\n";
    const Mat2d& B = strand.m_bendingMatrix;// strand.m_bendingMatrices[vtx];// std::cout << "B = " << B <<'\n';
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];// std::cout << "ilen = " << ilen <<'\n';
    const Vec2d& kappa = geometry.m_kappa[vtx];// std::cout << "kappa = " << kappa <<'\n';
    const Vec2d& kappaBar = strand.m_restBends[vtx];// std::cout << "kappaBar = " << kappaBar <<'\n';
    const Eigen::Matrix<Scalar, 11, 2>& gradKappa = geometry.m_gradKappa[vtx];// std::cout << "gradKappa = " << gradKappa <<'\n';

    //    localF = -ilen * gradKappa * B * ( kappa - kappaBar );
    localF = -ilen * B( 0, 0 ) * ( gradKappa * ( kappa - kappaBar ) ); // Okay if B is a constant x identity i.e. if the rods are isotropic.
}

//BendingForce::LocalJacobianType BendingForce::localJacobian( const ElasticStrand& strand,
//        const StrandGeometry& geometry, const IndexType vtx )
void BendingForce::computeLocalJacobian( BendingForce::LocalJacobianType& localJ,
        const ElasticStrand& strand, const StrandGeometry& geometry, const IndexType vtx )
{
    const Mat2d& B = strand.m_bendingMatrix;// strand.m_bendingMatrices[vtx];
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Vec2d& kappa = geometry.m_kappa[vtx];
    const Vec2d& kappaBar = strand.m_restBends[vtx];
    const Eigen::Matrix<Scalar, 11, 2>& gradKappa = geometry.m_gradKappa[vtx];

    symBProduct( localJ, B, gradKappa );
    localJ *= -ilen;

    std::pair<LocalJacobianType, LocalJacobianType> hessKappa;
    geometry.computeHessKappa( hessKappa, vtx );
    const Vec2d temp = -ilen * ( kappa - kappaBar ).transpose() * B;
    localJ += temp( 0 ) * hessKappa.first + temp( 1 ) * hessKappa.second;
}

void BendingForce::addInPosition( ForceVectorType& globalForce, const IndexType vtx,
        const LocalForceType& localForce )
{
    globalForce.segment<11> ( 4 * ( vtx - 1 ) ) += localForce;
}

void BendingForce::addInPosition( JacobianMatrixType& globalJacobian, const IndexType vtx,
        const LocalJacobianType& localJacobian )
{
    globalJacobian.localStencilAdd<11> ( 4 * ( vtx - 1 ), localJacobian );
}

}
