/*
 * ElasticStrand.hh
 *
 *  Created on: 7/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef ELASTICSTRAND_HH_
#define ELASTICSTRAND_HH_

#include <list>

#include "StrandBase.hh"
#include "StrandGeometry.hh"
#include "ElasticStrandUtils.hh"
#include "BandMatrix.hh"
#include "ElasticStrandParameters.hh"
#include "Typelist.hh"
#include "Forces/ForceBase.hh"

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

    Vec3d getVertex( const IndexType vtx ) const
    {
        return m_currentGeometry->getVertex( vtx );
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
        return m_currentGeometry->m_totalForce;
    }

    VecXd& getTotalForces()
    {
        return m_currentGeometry->m_totalForce;
    }

    const JacobianMatrixType& getTotalJacobian() const
    {
        return *( m_currentGeometry->m_totalJacobian );
    }

    JacobianMatrixType& getTotalJacobian()
    {
        return *( m_currentGeometry->m_totalJacobian );
    }

    Scalar getNewTotalEnergy() const
    {
        return m_futureGeometry->m_totalEnergy;
    }

    const VecXd& getNewTotalForces() const
    {
        return m_futureGeometry->m_totalForce;
    }

    VecXd& getNewTotalForces()
    {
        return m_futureGeometry->m_totalForce;
    }

    Scalar getMass( IndexType i ) const
    {
        return m_vertexMasses[i];
    }

    void prepareForSolving();
    void prepareForExamining();
    void acceptNewPositions();

    void addExternalForce( ForceBase* force )
    {
        m_externalForces.push_back( force );
    }

private:
    const std::list<ElasticStrand*>& getClumpingAttractors() const
    {
        return m_clumpingAttractors;
    }

    // For testing only, otherwise private:
    void resizeInternals();
    void freezeRestShape();

    Mat2d computeBendingMatrix( const IndexType vtx ) const;

    template<typename ForceT>
    void accumulateEF( StrandGeometry* geometry ) const;

    template<typename ForceT>
    void accumulateJ( StrandGeometry* geometry ) const;

    template<typename ForceT>
    void accumulateEFJ( StrandGeometry* geometry ) const;

    Vec3d closestPoint( const Vec3d& x ) const;

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

    // If we do anisotropic rods, the following will go in the geometry
    // std::vector<Mat2d, Eigen::aligned_allocator<Mat2d> > m_bendingMatrices;
    Mat2d m_bendingMatrix;

    // Flags
    bool m_readyForExamining;
    bool m_readyForSolving;

    // Forces that are not built-in
    std::list<ForceBase*> m_externalForces;

    std::list<ElasticStrand*> m_clumpingAttractors;

public:
    friend class StretchingForce;
    friend class BendingForce;
    friend class TwistingForce;
    friend class GravitationForce;
    friend class ClumpingForce;
    friend std::ostream& operator<<( std::ostream& os, const ElasticStrand& strand );
};

}

#endif /* ELASTICSTRAND_HH_ */
