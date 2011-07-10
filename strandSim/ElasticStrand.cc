/*
 * ElasticStrand.cc
 *
 *  Created on: 7/07/2011
 *      Author: jaubry
 */

#include "ElasticStrand.hh"

namespace strandsim
{

using namespace BASim;

ElasticStrand::ElasticStrand(VecXd& dofs, ParametersType parameters) : m_forceCachesUpToDate(false),
    m_degreesOfFreedom(dofs)
{
    // TODO: do something with parameters

    resizeInternals();
    updateForceCaches();
}

ElasticStrand::~ElasticStrand()
{
    // TODO Auto-generated destructor stub
}

void ElasticStrand::growLastVertex(const Vec3d& vertex, const Scalar torsion)
{
    assert(m_degreesOfFreedom.size() % 4 == 3);

    const IndexType oldNumDofs = static_cast<IndexType> (m_degreesOfFreedom.size());
    m_degreesOfFreedom.resize(oldNumDofs + 4); // TODO: check that resize preserves existing data.
    m_degreesOfFreedom[oldNumDofs] = torsion;
    m_degreesOfFreedom.segment<3> (oldNumDofs + 1) = vertex;

    resizeInternals();
}

void ElasticStrand::popLastVertex()
{
    assert(m_degreesOfFreedom.size() % 4 == 3);
    assert(m_degreesOfFreedom.size() > 7);

    const IndexType oldNumDofs = static_cast<IndexType> (m_degreesOfFreedom.size());
    m_degreesOfFreedom.resize(oldNumDofs - 4); // TODO: check that resize preserves existing data.

    resizeInternals();
}

void ElasticStrand::resizeInternals()
{
    assert(m_degreesOfFreedom.size() % 4 == 3);

    const IndexType numDofs = static_cast<IndexType> (m_degreesOfFreedom.size());

    m_numVertices = (numDofs + 1) / 4;
    m_totalForces.resize(numDofs);
    m_totalJacobian.resize(numDofs, numDofs);
    m_newDegreesOfFreedom.resize(numDofs);
}

void ElasticStrand::updateReferenceFrames()
{
    if (m_referenceFramesUpToDate)
        return;

    // Do some stuff here.

    m_referenceFramesUpToDate = true;
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
    updateReferenceFrames(); // Necessary to compute internal forces.

    accumulateStretchingForces();
    accumulateBendingForces();
    accumulateTorsionForces();
}

void ElasticStrand::acceptNewPositions()
{
    m_referenceFramesUpToDate = false;
    m_forceCachesUpToDate = false;

    m_degreesOfFreedom = m_newDegreesOfFreedom;

    updateForceCaches();
}

}
