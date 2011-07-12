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
    m_parameters(parameters), m_forceCachesUpToDate(false), m_degreesOfFreedom(dofs)
{
    resizeInternals();
    freezeRestShape(); // for now the rest shape is the shape in which the strand is created, unless modified later on.
    storeInitialFrames();
    updateForceCaches();
}

ElasticStrand::~ElasticStrand()
{
    // TODO Auto-generated destructor stub
}

void ElasticStrand::resizeInternals()
{
    assert(m_degreesOfFreedom.size() % 4 == 3); // dofs are 3 per vertex, one per edge
    assert(m_degreesOfFreedom.size() > 3); // minimum two vertices per rod

    const IndexType numDofs = static_cast<IndexType> (m_degreesOfFreedom.size());

    m_numVertices = (numDofs + 1) / 4;

    m_restLengths.resize(m_numVertices - 1);
    m_kappaBar.resize(m_numVertices - 1);

    m_currentLengths.resize(m_numVertices - 1);
    m_kappa.resize(m_numVertices - 1);
    m_gradKappa.resize(m_numVertices - 1);
    m_bendingMatrices.resize(m_numVertices - 1); // NB m_bendingMatrices[0] not used

    m_previousTangents.resize(3 * (m_numVertices - 1));

    m_referenceFrames1.resize(3 * (m_numVertices - 1));
    m_referenceFrames2.resize(3 * (m_numVertices - 1));
    m_materialFrames1.resize(3 * (m_numVertices - 1));
    m_materialFrames2.resize(3 * (m_numVertices - 1));

    m_totalForces.resize(numDofs);
    m_totalJacobian.resize(numDofs, numDofs);

    m_newDegreesOfFreedom.resize(numDofs);
}

void ElasticStrand::freezeRestShape()
{
    for (IndexType vtx = 0; vtx < m_numVertices - 1; ++vtx)
        m_restLengths[vtx] = getEdgeVector(vtx).norm();

    // TODO: compute m_kappaBar here?
}

void ElasticStrand::storeInitialFrames() // Called only in constructor
{
    // Previous tangents vector initially contains normalised edges
    for (IndexType vtx = 0; vtx < m_numVertices - 1; ++vtx)
        setPreviousTangent(vtx, getEdgeVector(vtx).normalized());

    // Initial first reference frame is arbitrary
    setReferenceFrame1(0, findNormal(getPreviousTangent(0)));
    setReferenceFrame2(0, getPreviousTangent(0).cross(getReferenceFrame1(0)));
    assert(fabs(getReferenceFrame2(0).norm() - 1) < std::numeric_limits<Scalar>::epsilon());

    // Next initial reference frames are obtained by space-parallel transportation
    for (IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx)
    {
        setReferenceFrame1(vtx,
                parallelTransport(getReferenceFrame1(vtx - 1), getPreviousTangent(vtx - 1), getPreviousTangent(vtx)));
        setReferenceFrame2(vtx, getPreviousTangent(vtx).cross(getReferenceFrame1(vtx)));
    }
}

void ElasticStrand::updateFrames() // and related stuff
{
    if (m_framesUpToDate)
        return;

    // Update reference frames by time-parallel transportation
    for (IndexType vtx = 0; vtx < m_numVertices - 1; ++vtx)
    {
        m_currentLengths[vtx] = getEdgeVector(vtx).norm();
        assert(!isSmall(m_currentLengths[vtx]));

        const Vec3d& currentTangent = getEdgeVector(vtx).normalized();
        Vec3d u = parallelTransport(getReferenceFrame1(vtx), getPreviousTangent(vtx), currentTangent);
        u -= u.dot(currentTangent) * currentTangent; // superfluous?
        u.normalize();
        setReferenceFrame1(vtx, u);
        setReferenceFrame2(vtx, currentTangent.cross(u));

        setPreviousTangent(vtx, currentTangent); // TODO: make sure that previous tangents are not used elsewhere
    }

    // Update other cached frame-dependent quantities
    for (IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx)
    {
        m_bendingMatrices[vtx] = computeBendingMatrix(vtx);
        m_kappa[vtx] = computeKappa(vtx);
        m_gradKappa[vtx] = computeGradKappa(vtx);
    }

    m_framesUpToDate = true;
}

Mat2d ElasticStrand::computeBendingMatrix(const IndexType vtx) const
{
    // In the case of isotropic strands the bending matrix is a multiple of identity. When we allow non isotropic strands, elliptic section + rotation have to be taken into account.
    return m_parameters.m_YoungsModulus * 0.25 * M_PI * square(square(m_parameters.m_radius)) * Mat2d::Identity();
}

Vec2d ElasticStrand::computeKappa(const IndexType vtx) const
{
    // TODO: compute here.
}

Eigen::Matrix<Scalar, 11, 2> ElasticStrand::computeGradKappa(const IndexType vtx) const
{
    // TODO: compute here.
}

std::pair<Eigen::Matrix<Scalar, 11, 11>, Eigen::Matrix<Scalar, 11, 11> > ElasticStrand::computeHessKappa(const IndexType vtx) const
{
    // TODO: compute here.
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
    assert(m_framesUpToDate);

    ForceAccumulator<StretchingForce>::accumulate(*this);
    ForceAccumulator<BendingForce>::accumulate(*this);
    ForceAccumulator<TwistingForce>::accumulate(*this);
}

void ElasticStrand::acceptNewPositions()
{
    m_framesUpToDate = false;
    m_forceCachesUpToDate = false;

    m_degreesOfFreedom = m_newDegreesOfFreedom; // copy the new positions

    updateFrames(); // This also updates m_previousTangents
    updateForceCaches();
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

}
