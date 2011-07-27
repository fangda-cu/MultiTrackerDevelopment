/*
 * GravitationForce.hh
 *
 *  Created on: 14/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef GRAVITATIONFORCE_HH_
#define GRAVITATIONFORCE_HH_

#include "../ElasticStrand.hh"

namespace strandsim
{

class GravitationForce
{
public:
    static const IndexType s_first = 0; // The first index on which this force can apply
    static const IndexType s_last = 0; // The last index (counting from the end)

    typedef Vec3d LocalForceType;
    typedef Mat3d LocalJacobianType;
    typedef VecXd ForceVectorType;

    GravitationForce();
    virtual ~GravitationForce();

    static std::string getName()
    {
        return "gravitation";
    }

    static Scalar localEnergy( const ElasticStrand& strand, const StrandGeometry& geometry,
            const IndexType vtx );

    static void computeLocalForce( LocalForceType& localF, const ElasticStrand& strand,
            const StrandGeometry& geometry, const IndexType vtx );

    static void computeLocalJacobian( LocalJacobianType& localJ, const ElasticStrand& strand,
            const StrandGeometry& geometry, const IndexType vtx );

    static void addInPosition( ForceVectorType& globalForce, const IndexType vtx,
            const LocalForceType& localForce );

    static void addInPosition( JacobianMatrixType& globalJacobian, const IndexType vtx,
            const LocalJacobianType& localJacobian );

    static void setGravity( const Vec3d& gravity )
    {
        s_gravity = gravity;
    }

private:
    static Vec3d s_gravity;
};

}

#endif /* GRAVITATIONFORCE_HH_ */
