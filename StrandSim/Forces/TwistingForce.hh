/*
 * TwistingForce.hh
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#ifndef TWISTINGFORCE_HH_
#define TWISTINGFORCE_HH_

#include "ForceBase.hh"
#include "../ElasticStrand.hh"

namespace strandsim
{

class TwistingForce: public ForceBase<ElasticStrand>
{
public:
    static const IndexType s_first = 1; // The first index on which this force can apply
    static const IndexType s_last = 1; // The last index (counting from the end)

    typedef Eigen::Matrix<Scalar, 11, 1> LocalForceType;
    typedef Eigen::Matrix<Scalar, 11, 11> LocalJacobianType;
    typedef VecXd ForceVectorType;
    typedef ElasticStrand::JacobianMatrixType JacobianMatrixType;

    TwistingForce();
    virtual ~TwistingForce();

    static std::string getName() {return "twisting";}
    static Scalar localEnergy(const ElasticStrand& strand, const StrandGeometry& geometry, const IndexType vtx);
    static LocalForceType localForce(const ElasticStrand& strand, const StrandGeometry& geometry, const IndexType vtx);
    static LocalJacobianType localJacobian(const ElasticStrand& strand, const StrandGeometry& geometry, const IndexType vtx);

    static void addInPosition(ForceVectorType& globalForce, const IndexType vtx, const LocalForceType& localForce);
    static void addInPosition(JacobianMatrixType& globalJacobian, const IndexType vtx, const LocalJacobianType& localJacobian);
};

}

#endif /* TWISTINGFORCE_HH_ */
