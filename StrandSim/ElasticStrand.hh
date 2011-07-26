/*
 * ElasticStrand.hh
 *
 *  Created on: 7/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef ELASTICSTRAND_HH_
#define ELASTICSTRAND_HH_

#include "StrandBase.hh"
#include "StrandGeometry.hh"
#include "ElasticStrandUtils.hh"
#include "BandMatrix.hh"
#include "ElasticStrandParameters.hh"
#include "Typelist.hh"

namespace strandsim
{

class StretchingForce;
class BendingForce;
class TwistingForce;
class GravitationForce;

typedef Typelist<StretchingForce, Typelist<TwistingForce, Typelist<BendingForce, Typelist<
        GravitationForce, NullType> > > > BuiltInForcesList;

class ElasticStrand: public StrandBase
{
public:
    typedef ElasticStrandParameters ParametersType;
    typedef strandsim::BandMatrix<Scalar, 10, 10> JacobianMatrixType; // TODO: replace this with a type that has built-in symmetry. Or at least squareness.

    ElasticStrand( const VecXd& dofs, const ParametersType& parameters );

    virtual ~ElasticStrand();

    const VecXd& getDegreesOfFreedom() const
    {
        return m_currentGeometry->m_degreesOfFreedom;
    }

    // Expose the future position storage so it can be changed by the stepper.
    VecXd& getNewDegreesOfFreedom()
    {
        m_futureGeometry->m_framesUpToDate = false;
        m_readyForExamining = false;

        return m_futureGeometry->m_degreesOfFreedom;
    }

    void filterNewGeometryLength();

    bool readyForSolving() const
    {
        return m_readyForSolving;
    }

    Scalar getTotalEnergy() const
    {
        return m_currentGeometry->m_totalEnergy;
    }

    const VecXd& getTotalForces() const
    {
        return m_currentGeometry->m_totalForces;
    }

    VecXd& getTotalForces()
    {
        return m_currentGeometry->m_totalForces;
    }

    const JacobianMatrixType& getTotalJacobian() const
    {
        return m_totalJacobian;
    }

    JacobianMatrixType& getTotalJacobian()
    {
        return m_totalJacobian;
    }

    Scalar getNewTotalEnergy() const
    {
        return m_futureGeometry->m_totalEnergy;
    }

    const VecXd& getNewTotalForces() const
    {
        return m_futureGeometry->m_totalForces;
    }

    VecXd& getNewTotalForces()
    {
        return m_futureGeometry->m_totalForces;
    }

    Scalar getMass( IndexType i ) const
    {
        return m_vertexMasses[i];
    }

    void prepareForSolving();
    void prepareForExamining();
    void acceptNewPositions();

public:
    // For testing only, otherwise private:
    void resizeInternals();
    void freezeRestShape();

    Mat2d computeBendingMatrix( const IndexType vtx ) const;

    template<typename ForceT>
    void accumulateEF( StrandGeometry* geometry ) const;

    template<typename ForceT>
    void accumulateJ( StrandGeometry* geometry ); // TODO: const once Jacobian in geometry

    template<typename ForceT>
    void accumulateEFJ( StrandGeometry* geometry ); // TODO: const once Jacobian in geometry

    /**
     * Member variables
     */

    // Size of the strand
    IndexType m_numVertices;

    // Current and future geometry
    StrandGeometry* m_currentGeometry;
    StrandGeometry* m_futureGeometry;

    // Rest shape and other physical parameters
    ParametersType m_parameters;
    std::vector<Scalar> m_vertexMasses;
    std::vector<Scalar> m_VoronoiLengths; // rest length around each vertex
    std::vector<Scalar> m_invVoronoiLengths; // their inverses
    std::vector<Scalar> m_restLengths;
    std::vector<Vec2d, Eigen::aligned_allocator<Vec2d> > m_restBends;
    std::vector<Scalar> m_restTwists;

    // Flags
    bool m_readyForExamining;
    bool m_readyForSolving;

    // Stuff that should go in the geometry
    JacobianMatrixType m_totalJacobian;
    // std::vector<Mat2d, Eigen::aligned_allocator<Mat2d> > m_bendingMatrices;
    Mat2d m_bendingMatrix;

public:
    friend class StretchingForce;
    //    friend class ForceAccumulator<StretchingForce> ;
    friend class BendingForce;
    //    friend class ForceAccumulator<BendingForce> ;
    friend class TwistingForce;
    //    friend class ForceAccumulator<TwistingForce> ;
    friend class GravitationForce;
    //    friend class ForceAccumulator<GravitationForce> ;
    friend std::ostream& operator<<( std::ostream& os, const ElasticStrand& strand );
};

}

#endif /* ELASTICSTRAND_HH_ */
