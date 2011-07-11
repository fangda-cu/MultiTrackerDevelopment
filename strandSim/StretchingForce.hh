/*
 * StretchingForce.hh
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#ifndef STRETCHINGFORCE_HH_
#define STRETCHINGFORCE_HH_

#include "Force.hh"
#include "ElasticStrand.hh"

namespace strandsim
{

class StretchingForce: public Force<ElasticStrand>
{
public:
    typedef Eigen::Matrix<Scalar, 6, 1> LocalForceType;
    typedef Eigen::Matrix<Scalar, 6, 6> LocalJacobianType;
    typedef ElasticStrand::ForceVectorType ForceVectorType;
    typedef ElasticStrand::JacobianMatrixType JacobianMatrixType;

    StretchingForce();
    virtual ~StretchingForce();

    static Scalar localEnergy(const ElasticStrand& strand, const IndexType vtx);
    static LocalForceType localForce(const ElasticStrand& strand, const IndexType vtx);
    static LocalJacobianType localJacobian(const ElasticStrand& strand, const IndexType vtx);

    static void addInPosition(ForceVectorType& totalForce, const IndexType vtx, const LocalForceType& localForce);
    static void addInPosition(JacobianMatrixType& totalForce, const IndexType vtx, const LocalJacobianType& localJacobian);
};

}

#endif /* STRETCHINGFORCE_HH_ */
