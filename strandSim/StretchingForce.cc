/*
 * StretchingForce.cc
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#include "StretchingForce.hh"
#include "ElasticStrand.hh"

namespace strandsim
{

StretchingForce::StretchingForce()
{
    // TODO Auto-generated constructor stub

}

StretchingForce::~StretchingForce()
{
    // TODO Auto-generated destructor stub
}

Scalar StretchingForce::localEnergy(const ElasticStrand& strand, const StrandGeometry& geometry, const IndexType vtx)
{
    const Scalar ks = strand.m_parameters.m_ks;
    const Scalar restLength = strand.m_restLengths[vtx];
    const Scalar length = geometry.m_lengths[vtx];

    return 0.5 * ks * square(length / restLength - 1.0) * restLength;
}

StretchingForce::LocalForceType StretchingForce::localForce(const ElasticStrand& strand, const StrandGeometry& geometry,
        const IndexType vtx)
{
    LocalForceType force;
    const Scalar ks = strand.m_parameters.m_ks;
    const Scalar restLength = strand.m_restLengths[vtx];
    const Scalar length = geometry.m_lengths[vtx];

    Vec3d f = ks * (length / restLength - 1.0) * geometry.getEdgeVector(vtx).normalized();
    force.segment<3> (0) = f;
    force.segment<3> (3) = -f;

    return force;
}

StretchingForce::LocalJacobianType StretchingForce::localJacobian(const ElasticStrand& strand, const StrandGeometry& geometry,
        const IndexType vtx)
{
    LocalJacobianType Jacobian;

    const Scalar ks = strand.m_parameters.m_ks;
    const Scalar restLength = strand.m_restLengths[vtx];
    const Scalar length = geometry.m_lengths[vtx];
    const Vec3d& edge = geometry.getEdgeVector(vtx);
    const Mat3d M = ks * ((1.0 / restLength - 1.0 / length) * Mat3d::Identity() + 1.0 / length * edge * edge.transpose()
            / square(length));

    Jacobian.block<3, 3> (0, 0) = Jacobian.block<3, 3> (3, 3) = -M;
    Jacobian.block<3, 3> (0, 3) = Jacobian.block<3, 3> (3, 0) = M;
    assert(isSymmetric(Jacobian));

    return Jacobian;
}

void StretchingForce::addInPosition(ForceVectorType& globalForce, const IndexType vtx, const LocalForceType& localForce)
{
    globalForce.segment<3> (4 * vtx) += localForce.segment<3> (0);
    globalForce.segment<3> (4 * (vtx + 1)) += localForce.segment<3> (3);
}

void StretchingForce::addInPosition(JacobianMatrixType& globalJacobian, const IndexType vtx,
        const LocalJacobianType& localJacobian)
{
    globalJacobian.edgeStencilAdd<6> (4 * vtx, localJacobian);
}

}
