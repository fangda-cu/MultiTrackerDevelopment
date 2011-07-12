/*
 * TwistingForce.hh
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#ifndef TWISTINGFORCE_HH_
#define TWISTINGFORCE_HH_

#include "Force.hh"
#include "ElasticStrand.hh"

namespace strandsim
{

class TwistingForce: public Force<ElasticStrand>
{
public:
    static const IndexType s_first = 1; // The first index on which this force can apply

    typedef Eigen::Matrix<Scalar, 11, 1> LocalForceType;
    typedef Eigen::Matrix<Scalar, 11, 11> LocalJacobianType;
    typedef ElasticStrand::ForceVectorType ForceVectorType;
    typedef ElasticStrand::JacobianMatrixType JacobianMatrixType;

    TwistingForce();
    virtual ~TwistingForce();

    static Scalar localEnergy(const ElasticStrand& strand, const IndexType vtx);
    static LocalForceType localForce(const ElasticStrand& strand, const IndexType vtx);
    static LocalJacobianType localJacobian(const ElasticStrand& strand, const IndexType vtx);

    static void addInPosition(ForceVectorType& totalForce, const IndexType vtx, const LocalForceType& localForce);
    static void addInPosition(JacobianMatrixType& totalForce, const IndexType vtx, const LocalJacobianType& localJacobian);
};

}

#endif /* TWISTINGFORCE_HH_ */
