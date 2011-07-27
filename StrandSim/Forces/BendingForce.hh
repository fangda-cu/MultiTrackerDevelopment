/*
 * BendingForce.hh
 *
 *  Created on: 12/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef BENDINGFORCE_HH_
#define BENDINGFORCE_HH_

#include "../ElasticStrand.hh"

namespace strandsim
{

class BendingForce
{
public:
    static const IndexType s_first = 1; // The first index on which this force can apply
    static const IndexType s_last = 1; // The last index (counting from the end)

    typedef Eigen::Matrix<Scalar, 11, 1> LocalForceType;
    typedef Eigen::Matrix<Scalar, 11, 11> LocalJacobianType;

    BendingForce();
    virtual ~BendingForce();

    static std::string getName()
    {
        return "bending";
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

private:

};

}

#endif /* BENDINGFORCE_HH_ */
