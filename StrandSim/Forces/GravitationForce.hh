/*
 * GravitationForce.hh
 *
 *  Created on: 14/07/2011
 *      Author: jaubry
 */

#ifndef GRAVITATIONFORCE_HH_
#define GRAVITATIONFORCE_HH_

#include "ForceBase.hh"
#include "../ElasticStrand.hh"

namespace strandsim
{

class GravitationForce: public ForceBase<ElasticStrand>
{
public:
    static const IndexType s_first = 0; // The first index on which this force can apply
    static const IndexType s_last = 0; // The last index (counting from the end)

    typedef Vec3d LocalForceType;
    typedef Mat3d LocalJacobianType;
    typedef ElasticStrand::ForceVectorType ForceVectorType;
    typedef ElasticStrand::JacobianMatrixType JacobianMatrixType;

    GravitationForce();
    virtual ~GravitationForce();

    static Scalar localEnergy(const ElasticStrand& strand, const StrandGeometry& geometry, const IndexType vtx);
    static LocalForceType localForce(const ElasticStrand& strand, const StrandGeometry& geometry, const IndexType vtx);
    static LocalJacobianType localJacobian(const ElasticStrand& strand, const StrandGeometry& geometry, const IndexType vtx);

    static void addInPosition(ForceVectorType& globalForce, const IndexType vtx, const LocalForceType& localForce);
    static void addInPosition(JacobianMatrixType& globalJacobian, const IndexType vtx, const LocalJacobianType& localJacobian);

    static void setGravity(const Vec3d& gravity)
    {
        s_gravity = gravity;
    }

private:
    static Vec3d s_gravity;
};

}

#endif /* GRAVITATIONFORCE_HH_ */
