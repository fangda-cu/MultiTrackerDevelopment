/*
 * GravitationForce.cc
 *
 *  Created on: 14/07/2011
 *      Author: jaubry
 */

#include "GravitationForce.hh"

namespace strandsim
{

GravitationForce::GravitationForce()
{
    // TODO Auto-generated constructor stub
}

GravitationForce::~GravitationForce()
{
    // TODO Auto-generated destructor stub
}

Scalar GravitationForce::localEnergy( const ElasticStrand& strand, const StrandGeometry& geometry,
        const IndexType vtx )
{
    return -strand.m_vertexMasses[vtx] * geometry.getVertex( vtx ).dot( s_gravity );
}

void GravitationForce::computeLocalForce(GravitationForce::LocalForceType& localF, const ElasticStrand& strand,
        const StrandGeometry& geometry, const IndexType vtx )
{
    localF = strand.m_vertexMasses[vtx] * s_gravity;
}

void GravitationForce::computeLocalJacobian( GravitationForce::LocalJacobianType& localJ,
        const ElasticStrand& strand, const StrandGeometry& geometry, const IndexType vtx )
{
    localJ.setZero(); // Jacobian is zero
}

void GravitationForce::addInPosition( ForceVectorType& globalForce, const IndexType vtx,
        const LocalForceType& localForce )
{
    globalForce.segment<3> ( 4 * vtx ) += localForce;
}

void GravitationForce::addInPosition( JacobianMatrixType& globalJacobian, const IndexType vtx,
        const LocalJacobianType& localJacobian )
{
    globalJacobian.localStencilAdd<3> ( 4 * vtx, localJacobian );
}

Vec3d GravitationForce::s_gravity( 0.0, 0.0, -981.0 ); // Acceleration vector, in cm/s^2 to match Maya's units

}
