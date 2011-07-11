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

ElasticStrand::ElasticStrand(VecXd& dofs, ParametersType parameters) :
    m_forceCachesUpToDate(false), m_degreesOfFreedom(dofs)
{
    // TODO: something with parameters

    resizeInternals();
    storePreviousTangents();
    updateForceCaches();
}

ElasticStrand::~ElasticStrand()
{
    // TODO Auto-generated destructor stub
}

void ElasticStrand::growVertex(const Vec3d& vertex, const Scalar torsion)
{
    assert(m_degreesOfFreedom.size() % 4 == 3);

    const IndexType oldNumDofs = static_cast<IndexType> (m_degreesOfFreedom.size());
    m_degreesOfFreedom.resize(oldNumDofs + 4); // TODO: check that resize preserves existing data.
    m_degreesOfFreedom[oldNumDofs] = torsion;
    m_degreesOfFreedom.segment<3> (oldNumDofs + 1) = vertex;

    resizeInternals();
}

void ElasticStrand::popVertex()
{
    assert(m_degreesOfFreedom.size() % 4 == 3);
    assert(m_degreesOfFreedom.size() > 7);

    const IndexType oldNumDofs = static_cast<IndexType> (m_degreesOfFreedom.size());
    m_degreesOfFreedom.resize(oldNumDofs - 4); // TODO: check that resize preserves existing data.

    resizeInternals();
}

void ElasticStrand::resizeInternals()
{
    assert(m_degreesOfFreedom.size() % 4 == 3); // dofs are 3 per vertex, one per edge
    assert(m_degreesOfFreedom.size() > 3); // minimum two vertices per rod

    const IndexType numDofs = static_cast<IndexType> (m_degreesOfFreedom.size());

    m_numVertices = (numDofs + 1) / 4;

    m_previousTangents.resize(3 * (m_numVertices - 1));

    m_referenceFrames1.resize(3 * (m_numVertices - 1));
    m_referenceFrames2.resize(3 * (m_numVertices - 1));
    m_materialFrames1.resize(3 * (m_numVertices - 1));
    m_materialFrames2.resize(3 * (m_numVertices - 1));

    m_totalForces.resize(numDofs);
    m_totalJacobian.resize(numDofs, numDofs);

    m_newDegreesOfFreedom.resize(numDofs);
}

void ElasticStrand::updateFrames()
{
    if (m_framesUpToDate)
        return;

    // Update reference frames by time-parallel transport
    for (IndexType vtx = 0; vtx < m_numVertices - 1; ++vtx)
    {
        const Vec3d& currentTangent = (getVertex(vtx + 1) - getVertex(vtx)).normalized();
        Vec3d u = parallelTransport(getReferenceFrame1(vtx), getPreviousTangent(vtx), currentTangent);
        u -= u.dot(currentTangent)*currentTangent; // superfluous?
        u.normalize();
        setReferenceFrame1(vtx, u);
        setReferenceFrame2(vtx, currentTangent.cross(u));
    }

    m_framesUpToDate = true;
}

void ElasticStrand::updateForceCaches()
{
    if (m_forceCachesUpToDate)
        return;

    m_totalEnergy = 0.0;
    m_totalForces.setZero();
    m_totalJacobian.setZero();

    accumulateInternalForces(); // Updates reference frames and internal energy, force, Jacobian.
    accumulateExternalForces(); // Updates external energy, force, Jacobian.

    m_forceCachesUpToDate = true;
}

void ElasticStrand::accumulateInternalForces()
{
    updateFrames(); // Necessary to compute internal forces.

    accumulateStretchingForces();
    accumulateBendingForces();
    accumulateTorsionForces();
}

void ElasticStrand::acceptNewPositions()
{
    m_framesUpToDate = false;
    m_forceCachesUpToDate = false;

    storePreviousTangents(); // for time-parallel transport
    m_degreesOfFreedom = m_newDegreesOfFreedom; // copy the new positions

    updateForceCaches();
}

}
