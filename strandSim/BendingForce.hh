/*
 * BendingForce.hh
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#ifndef BENDINGFORCE_HH_
#define BENDINGFORCE_HH_

#include "Force.hh"
#include "ElasticStrand.hh"

namespace strandsim
{

class BendingForce: public Force<ElasticStrand>
{
public:
    typedef Eigen::Matrix<Scalar, 11, 1> LocalForceType;
    typedef Eigen::Matrix<Scalar, 11, 11> LocalJacobianType;
    typedef ElasticStrand::ForceVectorType ForceVectorType;
    typedef ElasticStrand::JacobianMatrixType JacobianMatrixType;

    BendingForce();
    virtual ~BendingForce();

    static Scalar localEnergy(const ElasticStrand& strand, const IndexType vtx);
    static LocalForceType localForce(const ElasticStrand& strand, const IndexType vtx);
    static LocalJacobianType localJacobian(const ElasticStrand& strand, const IndexType vtx);

    static void addInPosition(ForceVectorType& totalForce, const IndexType vtx, const LocalForceType& localForce);
    static void addInPosition(JacobianMatrixType& totalForce, const IndexType vtx, const LocalJacobianType& localJacobian);
};

}

#endif /* BENDINGFORCE_HH_ */
