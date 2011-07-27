/*
 * ElasticStrand.cc
 *
 *  Created on: 7/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#include "ElasticStrand.hh"
#include "Forces/StretchingForce.hh"
#include "Forces/TwistingForce.hh"
#include "Forces/BendingForce.hh"
#include "Forces/GravitationForce.hh"

namespace strandsim
{

using namespace BASim;

ElasticStrand::ElasticStrand( const VecXd& dofs, const ParametersType& parameters ) :
    m_parameters( parameters ), m_currentGeometry( new StrandGeometry( dofs ) ),
            m_futureGeometry( new StrandGeometry( dofs ) ), m_readyForSolving( false )
{
    // Allocate the Jacobian matrix and store it in a shared pointer.
    m_currentGeometry->m_totalJacobian = m_futureGeometry->m_totalJacobian = boost::shared_ptr<
            JacobianMatrixType>( new JacobianMatrixType );

    resizeInternals();

    m_currentGeometry->storeInitialFrames();
    freezeRestShape(); // for now the rest shape is the shape in which the strand is created, unless modified later on.
    m_bendingMatrix = computeBendingMatrix( 0 ); // to be replaced by m_bendingMatrices (evolving with geometry) if we do anisotropic rods

    m_futureGeometry->storeInitialFrames();

    prepareForSolving();
}

ElasticStrand::~ElasticStrand()
{
    delete m_currentGeometry;
    delete m_futureGeometry;
}

// To be called on creation
void ElasticStrand::resizeInternals()
{
    m_currentGeometry->resizeSelf();
    m_futureGeometry->resizeSelf();

    const IndexType ndofs = static_cast<IndexType> ( m_currentGeometry->m_degreesOfFreedom.size() );
    m_numVertices = ( ndofs + 1 ) / 4;

    m_restLengths.resize( m_numVertices - 1 );
    m_restBends.resize( m_numVertices - 1 );
    m_restTwists.resize( m_numVertices - 1 );

    // m_bendingMatrices.resize( m_numVertices - 1 ); // NB m_bendingMatrices[0] not used
    m_vertexMasses.resize( m_numVertices );
    m_VoronoiLengths.resize( m_numVertices );
    m_invVoronoiLengths.resize( m_numVertices );
}

// Take the current geometry as rest shape
void ElasticStrand::freezeRestShape()
{
    // Fix rest lengths
    for ( IndexType vtx = 0; vtx < m_numVertices - 1; ++vtx )
        m_restLengths[vtx] = m_currentGeometry->getEdgeVector( vtx ).norm();

    // Compute Voronoi lengths
    m_VoronoiLengths[0] = 0.5 * m_restLengths[0];
    for ( IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx )
        m_VoronoiLengths[vtx] = 0.5 * ( m_restLengths[vtx - 1] + m_restLengths[vtx] );
    m_VoronoiLengths[m_numVertices - 1] = 0.5 * m_restLengths[m_numVertices - 2];

    // Compute masses and inverse of Voronoi lengths
    for ( IndexType vtx = 0; vtx < m_numVertices; ++vtx )
    {
        m_vertexMasses[vtx] = m_parameters.m_density * m_VoronoiLengths[vtx] * M_PI
                * m_parameters.m_radius * m_parameters.m_radius;
        m_invVoronoiLengths[vtx] = 1.0 / m_VoronoiLengths[vtx];
    }

    for ( IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx )
    {
        m_currentGeometry->computeKappa( m_restBends[vtx], vtx );
        m_restTwists[vtx] = 0; // TODO: check that this is correct
    }
}

// For anisotropic strands this will become meaningful
Mat2d ElasticStrand::computeBendingMatrix( const IndexType vtx ) const
{
    return m_parameters.m_YoungsModulus * 0.25 * M_PI * square( square( m_parameters.m_radius ) )
            * Mat2d::Identity();
}

template<typename ForceT>
void ElasticStrand::accumulateEF( StrandGeometry* geometry ) const
{
    assert( geometry->m_framesUpToDate );

    typename ForceT::LocalForceType localF;

    for ( IndexType vtx = ForceT::s_first; vtx < m_numVertices - ForceT::s_last; ++vtx )
    {
        geometry->m_totalEnergy += ForceT::localEnergy( *this, *geometry, vtx );

        ForceT::computeLocalForce( localF, *this, *geometry, vtx );
        ForceT::addInPosition( geometry->m_totalForce, vtx, localF );
    }
    // std::cout << "With " << ForceT::getName() << ": " << geometry->m_totalForce << '\n';
}

template<typename ForceT>
void ElasticStrand::accumulateJ( StrandGeometry* geometry ) const
{
    assert( geometry->m_framesUpToDate );

    typename ForceT::LocalJacobianType localJ;

    for ( IndexType vtx = ForceT::s_first; vtx < m_numVertices - ForceT::s_last; ++vtx )
    {
        ForceT::computeLocalJacobian( localJ, *this, *geometry, vtx );
        ForceT::addInPosition( *( geometry->m_totalJacobian ), vtx, localJ );
    }
}

template<typename ForceT>
void ElasticStrand::accumulateEFJ( StrandGeometry* geometry ) const
{
    assert( geometry->m_framesUpToDate );

    typename ForceT::LocalForceType localF;
    typename ForceT::LocalJacobianType localJ;

    for ( IndexType vtx = ForceT::s_first; vtx < m_numVertices - ForceT::s_last; ++vtx )
    {
        geometry->m_totalEnergy += ForceT::localEnergy( *this, *geometry, vtx );

        ForceT::computeLocalForce( localF, *this, *geometry, vtx );
        ForceT::addInPosition( geometry->m_totalForce, vtx, localF );

        ForceT::computeLocalJacobian( localJ, *this, *geometry, vtx );
        ForceT::addInPosition( *( geometry->m_totalJacobian ), vtx, localJ );
    }
}

// Compute energy, force and Jacobian, based on current geometry
void ElasticStrand::prepareForSolving()
{
    if ( m_readyForSolving )
        return;

    m_currentGeometry->m_totalEnergy = 0.0;
    m_currentGeometry->m_totalForce.setZero();
    m_currentGeometry->m_totalJacobian->setZero();

    m_currentGeometry->updateFrames();
    // for ( IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx )
    //     m_bendingMatrices[vtx] = computeBendingMatrix( vtx );

    accumulateEFJ<StretchingForce> ( m_currentGeometry );
    accumulateEFJ<TwistingForce> ( m_currentGeometry );
    accumulateEFJ<BendingForce> ( m_currentGeometry );
    accumulateEFJ<GravitationForce> ( m_currentGeometry );

    for ( std::list<ForceBase*>::const_iterator force = m_externalForces.begin(); force
            != m_externalForces.end(); ++force )
        if ( *force )
            ( *force )->accumulateEFJ( *m_currentGeometry, *this );

    *( m_currentGeometry->m_totalJacobian ) *= -1.0; // To match BASim's sign conventions

    m_readyForSolving = true;
}

// Compute energy and force, based on tentative geometry
void ElasticStrand::prepareForExamining()
{
    if ( m_readyForExamining )
        return;

    m_futureGeometry->m_referenceFrames1 = m_currentGeometry->m_referenceFrames1; // We need the old ones to compute the new ones
    m_futureGeometry->m_previousTangents = m_currentGeometry->m_previousTangents; // Can we avoid the copis?
    m_futureGeometry->m_referenceTwists = m_currentGeometry->m_referenceTwists;

    m_futureGeometry->updateFrames();
    //  for ( IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx )
    //      m_bendingMatrices[vtx] = computeBendingMatrix( vtx );

    m_futureGeometry->m_totalEnergy = 0.0;
    m_futureGeometry->m_totalForce.setZero();

    accumulateEF<StretchingForce> ( m_futureGeometry );
    accumulateEF<TwistingForce> ( m_futureGeometry );
    accumulateEF<BendingForce> ( m_futureGeometry );
    accumulateEF<GravitationForce> ( m_futureGeometry );

    for ( std::list<ForceBase*>::const_iterator force = m_externalForces.begin(); force
            != m_externalForces.end(); ++force )
        if ( *force )
            ( *force )->accumulateEF( *m_futureGeometry, *this );

    m_readyForExamining = true;
}

void ElasticStrand::acceptNewPositions()
{
    assert( m_readyForExamining );

    std::swap( m_currentGeometry, m_futureGeometry );

    // Compute the Jacobian
    m_currentGeometry->m_totalJacobian->setZero();
    // for ( IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx )
    //     m_bendingMatrices[vtx] = computeBendingMatrix( vtx );

    accumulateJ<StretchingForce> ( m_currentGeometry );
    accumulateJ<TwistingForce> ( m_currentGeometry );
    accumulateJ<BendingForce> ( m_currentGeometry );
    accumulateJ<GravitationForce> ( m_currentGeometry );

    for ( std::list<ForceBase*>::const_iterator force = m_externalForces.begin(); force
            != m_externalForces.end(); ++force )
        if ( *force )
            ( *force )->accumulateJ( *m_currentGeometry, *this );

    *( m_currentGeometry->m_totalJacobian ) *= -1.0; // To match BASim's sign conventions

    m_readyForSolving = true;
}

void ElasticStrand::filterNewGeometryLength()
{
    m_futureGeometry->m_framesUpToDate = false; // Because we are changing stuff below
    m_readyForExamining = false;

    Vec3d xaP = m_currentGeometry->getVertex( 1 );
    Vec3d xaN = m_futureGeometry->getVertex( 1 );
    Vec3d disp;

    for ( int i = 2; i < m_numVertices; ++i )
    {
        const Vec3d& xbP = m_currentGeometry->getVertex( i );
        const Vec3d& xbN = m_futureGeometry->getVertex( i ) + disp;

        const Scalar lP = ( xbP - xaP ).norm(); // This was computed already, fetch it from where it is cached
        const Scalar lN = ( xbN - xaN ).norm();

        // compute and store revised delta
        const Vec3d xbNrev = xaN + ( xbN - xaN ) * lP / lN;
        m_futureGeometry->setVertex( i, xbNrev );
        disp += xbNrev - xbN;

        xaP = xbP;
        xaN = xbNrev;
    }
}

Vec3d ElasticStrand::closestPoint( const Vec3d& x ) const
{
    Scalar mindist = std::numeric_limits<Scalar>::max();
    Vec3d winner;

    for ( int vtx = 0; vtx < m_numVertices - 1; ++vtx )
    {
        Vec3d y = ClosestPtPointSegment( x, getVertex( vtx ), getVertex( vtx + 1 ) );
        Scalar dist = ( y - x ).squaredNorm();
        if ( dist < mindist )
        {
            mindist = dist;
            winner = y;
        }
    }

    return winner;
}

std::ostream& operator<<( std::ostream& os, const ElasticStrand& strand )
{
    os << '{';
    for ( int i = 0; i < strand.m_numVertices - 1; i++ )
    {
        os << '{' << strand.m_currentGeometry->m_degreesOfFreedom[4 * i] << ", "
                << strand.m_currentGeometry->m_degreesOfFreedom[4 * i + 1] << ", "
                << strand.m_currentGeometry->m_degreesOfFreedom[4 * i + 2] << "}, ";
        // os << strand.m_currentGeometry->m_degreesOfFreedom[4 * i + 3] << ', ';
    }
    os << '{' << strand.m_currentGeometry->m_degreesOfFreedom[4 * ( strand.m_numVertices - 1 )]
            << ", " << strand.m_currentGeometry->m_degreesOfFreedom[4 * ( strand.m_numVertices - 1 )
            + 1] << ", " << strand.m_currentGeometry->m_degreesOfFreedom[4 * ( strand.m_numVertices
            - 1 ) + 2] << '}';
    os << '}';

    return os;
}

}
