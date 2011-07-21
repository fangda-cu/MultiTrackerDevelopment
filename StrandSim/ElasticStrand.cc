/*
 * ElasticStrand.cc
 *
 *  Created on: 7/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#include "ElasticStrand.hh"

namespace strandsim
{

using namespace BASim;

ElasticStrand::ElasticStrand( VecXd& dofs, const ParametersType& parameters ) :
    m_parameters( parameters ), m_geometry( dofs ), m_dofsStorage(),
            m_newGeometry( m_dofsStorage ), m_readyForSolving( false )
{
    resizeInternals();

    m_geometry.storeInitialFrames();
    freezeRestShape(); // for now the rest shape is the shape in which the strand is created, unless modified later on.

    m_newGeometry.storeInitialFrames();

    prepareForSolving();
}

ElasticStrand::~ElasticStrand()
{
    // TODO Auto-generated destructor stub
}

// To be called on creation
void ElasticStrand::resizeInternals()
{
    m_geometry.resizeSelf();

    m_dofsStorage = m_geometry.m_degreesOfFreedom; // Copy of the initial position so m_newGeometry can compute its initial frames
    m_newGeometry.resizeSelf();

    const IndexType ndofs = static_cast<IndexType> ( m_geometry.m_degreesOfFreedom.size() );
    m_numVertices = ( ndofs + 1 ) / 4;

    m_restLengths.resize( m_numVertices - 1 );
    m_kappaBar.resize( m_numVertices - 1 );
    m_restTwists.resize( m_numVertices - 1 );

    m_bendingMatrices.resize( m_numVertices - 1 ); // NB m_bendingMatrices[0] not used
    m_vertexMasses.resize( m_numVertices );
    m_VoronoiLengths.resize( m_numVertices );
    m_invVoronoiLengths.resize( m_numVertices );

    m_totalForces.resize( ndofs );
    m_totalJacobian.resize( ndofs, ndofs );
    m_newTotalForces.resize( ndofs );
}

// Take the current geometry as rest shape
void ElasticStrand::freezeRestShape()
{
    // Fix rest lengths
    for ( IndexType vtx = 0; vtx < m_numVertices - 1; ++vtx )
        m_restLengths[vtx] = m_geometry.getEdgeVector( vtx ).norm();

    // Compute Voronoi lenghts
    m_VoronoiLengths[0] = 0.5 * m_restLengths[0];
    for ( IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx )
        m_VoronoiLengths[vtx] = 0.5 * ( m_restLengths[vtx - 1] + m_restLengths[vtx] );
    m_VoronoiLengths[m_numVertices - 1] = 0.5 * m_restLengths[m_numVertices - 2];

    // Compute masses and inverse of Voronoi lenghts
    for ( IndexType vtx = 0; vtx < m_numVertices; ++vtx )
    {
        m_vertexMasses[vtx] = m_parameters.m_density * m_VoronoiLengths[vtx] * M_PI
                * m_parameters.m_radius * m_parameters.m_radius;
        m_invVoronoiLengths[vtx] = 1.0 / m_VoronoiLengths[vtx];
    }

    for ( IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx )
    {
        m_kappaBar[vtx] = m_geometry.computeKappa( vtx );
        m_restTwists[vtx] = 0; // TODO: check that this is correct
    }
}

// For anisotropic strands this will become meaningful
Mat2d ElasticStrand::computeBendingMatrix( const IndexType vtx ) const
{
    return m_parameters.m_YoungsModulus * 0.25 * M_PI * square( square( m_parameters.m_radius ) )
            * Mat2d::Identity();
}

// Compute energy, force and Jacobian, based on current geometry
void ElasticStrand::prepareForSolving()
{
    if ( m_readyForSolving )
        return;

    m_totalEnergy = 0.0;
    m_totalForces.setZero();
    m_totalJacobian.setZero();

    m_geometry.updateFrames();
    for ( IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx )
        m_bendingMatrices[vtx] = computeBendingMatrix( vtx );

    ForceAccumulator<StretchingForce>::accumulate( m_totalEnergy, m_totalForces, m_totalJacobian,
            *this, m_geometry );
    ForceAccumulator<BendingForce>::accumulate( m_totalEnergy, m_totalForces, m_totalJacobian,
            *this, m_geometry );
    ForceAccumulator<TwistingForce>::accumulate( m_totalEnergy, m_totalForces, m_totalJacobian,
            *this, m_geometry );
    ForceAccumulator<GravitationForce>::accumulate( m_totalEnergy, m_totalForces, m_totalJacobian,
            *this, m_geometry );

    m_readyForSolving = true;
}

// Compute energy and force, based on tentative geometry
void ElasticStrand::prepareForExamining()
{
    if ( m_readyForExamining )
        return;

    m_newTotalEnergy = 0.0;
    m_newTotalForces.setZero();

    m_newGeometry.updateFrames(); // This is why we don't want to cache the 2nd derivatives, they won't be used here.
    for ( IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx )
        m_bendingMatrices[vtx] = computeBendingMatrix( vtx );

    ForceAccumulator<StretchingForce>::accumulate( m_newTotalEnergy, m_newTotalForces, *this,
            m_newGeometry );
    std::cout << "Stretching energy = " << m_newTotalEnergy << " forces norm = " << m_newTotalForces.norm() << '\n';
    ForceAccumulator<BendingForce>::accumulate( m_newTotalEnergy, m_newTotalForces, *this,
            m_newGeometry );
    std::cout << "Bending energy = " << m_newTotalEnergy << " forces norm = " << m_newTotalForces.norm() << '\n';
    ForceAccumulator<TwistingForce>::accumulate( m_newTotalEnergy, m_newTotalForces, *this,
            m_newGeometry );
    std::cout << "Twisting energy = " << m_newTotalEnergy << " forces norm = " << m_newTotalForces.norm() << '\n';
    ForceAccumulator<GravitationForce>::accumulate( m_newTotalEnergy, m_newTotalForces, *this,
            m_newGeometry );
    std::cout << "Grawity energy = " << m_newTotalEnergy << " forces norm = " << m_newTotalForces.norm() << '\n';

    m_readyForExamining = true;
}

void ElasticStrand::acceptNewPositions()
{
    assert( m_readyForExamining );

    m_geometry = m_newGeometry; // copy the new positions, frames, everything. TODO: copy only what's needed

    // Copy energy and forces from new
    m_totalEnergy = m_newTotalEnergy;
    m_totalForces = m_newTotalForces;

    // Compute the Jacobian
    m_totalJacobian.setZero();
    for ( IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx )
        m_bendingMatrices[vtx] = computeBendingMatrix( vtx );

    ForceAccumulator<StretchingForce>::accumulate( m_totalJacobian, *this, m_geometry );
    ForceAccumulator<BendingForce>::accumulate( m_totalJacobian, *this, m_geometry );
    ForceAccumulator<TwistingForce>::accumulate( m_totalJacobian, *this, m_geometry );
    ForceAccumulator<TwistingForce>::accumulate( m_totalJacobian, *this, m_geometry );

    m_readyForSolving = true;
}

void ElasticStrand::filterNewGeometryLength()
{
    m_newGeometry.m_framesUpToDate = false; // Because we are changing stuff below
    m_readyForExamining = false;

    Vec3d xaP = m_geometry.getVertex( 1 );// std::cout << "xaP = " << xaP << '\n';
    Vec3d xaN = m_newGeometry.getVertex( 1 );// std::cout << "xaN = " << xaN << '\n';
    Vec3d disp;

    for ( int i = 2; i < m_numVertices; ++i )
    {
        Vec3d xbP = m_geometry.getVertex( i );
        Vec3d xbN = m_newGeometry.getVertex( i ) + disp;

        Scalar lP = ( xbP - xaP ).norm(); // This was computed already, fetch it from where it is cached
        Scalar lN = ( xbN - xaN ).norm();

        // compute and store revised delta
        Vec3d xbNrev = xaN + ( xbN - xaN ) * lP / lN;
        m_newGeometry.setVertex( i, xbNrev );
        disp += xbNrev - xbN;

        xaP = xbP;// std::cout << "xaP = " << xaP << '\n';
        xaN = xbNrev;// std::cout << "xaN = " << xaN << '\n';
    }

}

std::ostream& operator<<( std::ostream& os, const ElasticStrand& strand )
{
    os << '{';
    for ( int i = 0; i < strand.m_numVertices - 1; i++ )
    {
        os << '{' << strand.m_geometry.m_degreesOfFreedom[4 * i] << ", "
                << strand.m_geometry.m_degreesOfFreedom[4 * i + 1] << ", "
                << strand.m_geometry.m_degreesOfFreedom[4 * i + 2] << "}, ";
        // os << strand.m_geometry.m_degreesOfFreedom[4 * i + 3] << ', ';
    }
    os << '{' << strand.m_geometry.m_degreesOfFreedom[4 * ( strand.m_numVertices - 1 )] << ", "
            << strand.m_geometry.m_degreesOfFreedom[4 * ( strand.m_numVertices - 1 ) + 1] << ", "
            << strand.m_geometry.m_degreesOfFreedom[4 * ( strand.m_numVertices - 1 ) + 2] << '}';
    os << '}';

    return os;
}

}
