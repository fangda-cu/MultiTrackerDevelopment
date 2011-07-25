/*
 * TwistingForce.cc
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#include "TwistingForce.hh"

namespace strandsim
{

TwistingForce::TwistingForce()
{
    // TODO Auto-generated constructor stub

}

TwistingForce::~TwistingForce()
{
    // TODO Auto-generated destructor stub
}

Scalar TwistingForce::localEnergy(const ElasticStrand& strand, const StrandGeometry& geometry, const IndexType vtx)
{
    const Scalar kt = strand.m_parameters.m_kt;
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Scalar undefTwist = strand.m_restTwists[vtx];
    const Scalar twist = geometry.m_twists[vtx];

    return 0.5 * kt * square(twist - undefTwist) * ilen;
}

TwistingForce::LocalForceType TwistingForce::localForce(const ElasticStrand& strand, const StrandGeometry& geometry,
        const IndexType vtx)
{
    const Scalar kt = strand.m_parameters.m_kt;
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Scalar undefTwist = strand.m_restTwists[vtx];
    const Scalar twist = geometry.m_twists[vtx];

    return -kt * ilen * (twist - undefTwist) * geometry.m_gradTwists[vtx];
}

TwistingForce::LocalJacobianType TwistingForce::localJacobian(const ElasticStrand& strand, const StrandGeometry& geometry,
        const IndexType vtx)
{
    const Scalar kt = strand.m_parameters.m_kt;
    const Scalar ilen = strand.m_invVoronoiLengths[vtx];
    const Scalar twist = geometry.m_twists[vtx];
    const Scalar undeformedTwist = strand.m_restTwists[vtx];

    const LocalForceType& gradTwist = geometry.m_gradTwists[vtx];
    LocalJacobianType hessTwist;
    geometry.computeHessTwist(hessTwist, vtx );

    return -kt * ilen * ((twist - undeformedTwist) * hessTwist + gradTwist * gradTwist.transpose());
}

void TwistingForce::addInPosition(ForceVectorType& globalForce, const IndexType vtx, const LocalForceType& localForce)
{
    globalForce.segment<11> (4 * (vtx - 1)) += localForce;
}

void TwistingForce::addInPosition(JacobianMatrixType& globalJacobian, const IndexType vtx,
        const LocalJacobianType& localJacobian)
{
    globalJacobian.localStencilAdd<11> (4 * (vtx - 1), localJacobian);
}

}
