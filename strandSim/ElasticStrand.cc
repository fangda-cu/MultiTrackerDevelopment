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

ElasticStrand::ElasticStrand(VecXd& dofs, const ParametersType& parameters) :
    m_parameters(parameters), m_geometry(dofs), m_dofsStorage(), m_newGeometry(m_dofsStorage), m_readyForSolving(false)
{
    resizeInternals();

    m_geometry.storeInitialFrames(); // Probably not needed because the frames in m_geometry will be copied from the ones in m_newGeometry
    m_newGeometry.storeInitialFrames();

    freezeRestShape(); // for now the rest shape is the shape in which the strand is created, unless modified later on.

    prepareForSolving(); // does nothing if the rest shape is the same as creation shape (see above).
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

    const IndexType ndofs = static_cast<IndexType> (m_geometry.m_degreesOfFreedom.size());
    m_numVertices = (ndofs + 1) / 4;

    m_restLengths.resize(m_numVertices - 1);
    m_kappaBar.resize(m_numVertices - 1);
    m_restTwists.resize(m_numVertices - 1);

    m_bendingMatrices.resize(m_numVertices - 1); // NB m_bendingMatrices[0] not used

    m_totalForces.resize(ndofs);
    m_totalJacobian.resize(ndofs, ndofs);
    m_newTotalForces.resize(ndofs);
}

// Take the current geometry as rest shape
void ElasticStrand::freezeRestShape()
{
    for (IndexType vtx = 0; vtx < m_numVertices - 1; ++vtx)
        m_restLengths[vtx] = m_geometry.getEdgeVector(vtx).norm();

    for (IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx)
    {
        m_kappaBar[vtx] = m_geometry.computeKappa(vtx);
        m_restTwists[vtx] = 0; // TODO: check that this is correct
    }
}

// For anisotropic strands this will become meaningful
Mat2d ElasticStrand::computeBendingMatrix(const IndexType vtx) const
{
    return m_parameters.m_YoungsModulus * 0.25 * M_PI * square(square(m_parameters.m_radius)) * Mat2d::Identity();
}

// Compute energy, force and Jacobian, based on current geometry
void ElasticStrand::prepareForSolving()
{
    if (m_readyForSolving)
        return;

    m_totalEnergy = 0.0;
    m_totalForces.setZero();
    m_totalJacobian.setZero();

    m_geometry.updateFrames();
    for (IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx)
        m_bendingMatrices[vtx] = computeBendingMatrix(vtx);

    ForceAccumulator<StretchingForce>::accumulate(m_totalEnergy, m_totalForces, m_totalJacobian, *this, m_geometry);
    ForceAccumulator<BendingForce>::accumulate(m_totalEnergy, m_totalForces, m_totalJacobian, *this, m_geometry);
    ForceAccumulator<TwistingForce>::accumulate(m_totalEnergy, m_totalForces, m_totalJacobian, *this, m_geometry);

    m_readyForSolving = true;
}

// Compute energy and force, based on tentative geometry
void ElasticStrand::prepareForExamining()
{
    if (m_readyForExamining)
        return;

    m_newTotalEnergy = 0.0;
    m_newTotalForces.setZero();

    m_newGeometry.updateFrames(); // This is why we don't want to cache the 2nd derivatives, they won't be used here.
    for (IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx)
        m_bendingMatrices[vtx] = computeBendingMatrix(vtx);

    ForceAccumulator<StretchingForce>::accumulate(m_newTotalEnergy, m_newTotalForces, *this, m_newGeometry);
    ForceAccumulator<BendingForce>::accumulate(m_newTotalEnergy, m_newTotalForces, *this, m_newGeometry);
    ForceAccumulator<TwistingForce>::accumulate(m_newTotalEnergy, m_newTotalForces, *this, m_newGeometry);

    m_readyForExamining = true;
}

void ElasticStrand::acceptNewPositions()
{
    assert(m_readyForExamining);

    m_geometry = m_newGeometry; // copy the new positions, frames, everything. TODO: copy only what's needed

    // Copy energy and forces from new
    m_totalEnergy = m_newTotalEnergy;
    m_totalForces = m_newTotalForces;

    // Compute the Jacobian
    m_totalJacobian.setZero();
    for (IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx)
        m_bendingMatrices[vtx] = computeBendingMatrix(vtx);

    ForceAccumulator<StretchingForce>::accumulate(m_totalJacobian, *this, m_geometry);
    ForceAccumulator<BendingForce>::accumulate(m_totalJacobian, *this, m_geometry);
    ForceAccumulator<TwistingForce>::accumulate(m_totalJacobian, *this, m_geometry);

    m_readyForSolving = true;
}

}
