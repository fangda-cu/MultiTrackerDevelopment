/*
 * ElasticStrand.hh
 *
 *  Created on: 7/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef ELASTICSTRAND_HH_
#define ELASTICSTRAND_HH_

#include "StrandBase.hh"
#include "ElasticStrandUtils.hh"
#include "Edge.hh"
#include "BandMatrix.hh"

namespace strandsim
{

class ElasticStrandParameters
{

};

class ElasticStrand: public StrandBase
{
public:
    typedef ElasticStrandParameters ParametersType;
    typedef VecXd ForceVectorType;
    typedef strandsim::BandMatrix<Scalar, 10, 10> JacobianMatrixType; // TODO: replace this with a type that has built-in symmetry

    ElasticStrand(VecXd& dofs, ParametersType parameters);

    virtual ~ElasticStrand();

    const Vec3d getVertex(const IndexType vtx) const
    {
        assert(vtx < m_numVertices);

        return m_degreesOfFreedom.segment<3> (4 * vtx);
    }

    Scalar getTorsion(const IndexType vtx) const
    {
        assert(vtx < m_numVertices - 1);

        return m_degreesOfFreedom[4 * vtx + 3];
    }

    const Edge<ElasticStrand> getEdge(const IndexType vtx) const
    {
        assert(vtx < m_numVertices - 1);

        return Edge<ElasticStrand> (*this, vtx);
    }

    const Vec3d getReferenceFrame1(const IndexType vtx) const
    {
        assert(vtx < m_numVertices - 1);

        return m_referenceFrames1.segment<3> (3 * vtx);
    }

    void setReferenceFrame1(const IndexType vtx, const Vec3d& vec)
    {
        assert(vtx < m_numVertices - 1);

        m_referenceFrames1.segment<3> (3 * vtx) = vec;
    }

    const Vec3d getReferenceFrame2(const IndexType vtx) const
    {
        assert(vtx < m_numVertices - 1);

        return m_referenceFrames2.segment<3> (3 * vtx);
    }

    void setReferenceFrame2(const IndexType vtx, const Vec3d& vec)
    {
        assert(vtx < m_numVertices - 1);

        m_referenceFrames2.segment<3> (3 * vtx) = vec;
    }

    const Vec3d getMaterialFrame1(const IndexType vtx) const
    {
        assert(vtx < m_numVertices - 1);

        return m_materialFrames1.segment<3> (3 * vtx);
    }

    void setMaterialFrame1(const IndexType vtx, const Vec3d& vec)
    {
        assert(vtx < m_numVertices - 1);

        m_materialFrames1.segment<3> (3 * vtx) = vec;
    }

    const Vec3d getMaterialFrame2(const IndexType vtx) const
    {
        assert(vtx < m_numVertices - 1);

        return m_materialFrames2.segment<3> (3 * vtx);
    }

    void setMaterialFrame2(const IndexType vtx, const Vec3d& vec)
    {
        assert(vtx < m_numVertices - 1);

        m_materialFrames2.segment<3> (3 * vtx) = vec;
    }

    void storePreviousTangents()
    {
        for (IndexType vtx = 0; vtx < m_numVertices - 1; ++vtx)
            m_previousTangents.segment<3> (vtx) = (getVertex(vtx + 1) - getVertex(vtx)).normalized();
    }

    const Vec3d getPreviousTangent(const IndexType vtx)
    {
        assert(vtx < m_numVertices - 1);

        return m_previousTangents.segment<3> (vtx);
    }

    void growVertex(const Vec3d& vertex, const Scalar torsion);
    void popVertex();

    const VecXd& getDegreesOfFreedom() const
    {
        return m_degreesOfFreedom;
    }

    // Expose the future position storage so it can be changed by the stepper.
    VecXd& getNewDegreesOfFreedom()
    {
        return m_newDegreesOfFreedom;
    }

    bool forceCachesUpToDate() const
    {
        return m_forceCachesUpToDate;
    }

    Scalar getTotalEnergy() const
    {
        return m_totalEnergy;
    }

    const ForceVectorType& getTotalForces() const
    {
        return m_totalForces;
    }

    const JacobianMatrixType& getTotalJacobian() const
    {
        return m_totalJacobian;
    }

    void acceptNewPositions();

private:
    void resizeInternals();

    void updateFrames();

    void updateForceCaches();
    void accumulateInternalForces();
    void accumulateExternalForces();

    void accumulateStretchingForces();
    void accumulateBendingForces();
    void accumulateTorsionForces();

    /**
     * Member variables
     */

    // Size of the strand
    IndexType m_numVertices;

    // References passed by Maya, who owns the original
    VecXd& m_degreesOfFreedom; // size = m_numVertices * 4 -1

    // Previous time storage, for time-parallel stepping
    VecXd m_previousTangents;

    // Reference and material frames
    bool m_framesUpToDate;
    VecXd m_referenceFrames1; // first vectors of reference frame
    VecXd m_referenceFrames2; // second vectors of reference frame
    VecXd m_materialFrames1;// first vectors of material frame
    VecXd m_materialFrames2; // second vectors of material frame

    // Force caching
    bool m_forceCachesUpToDate;
    Scalar m_totalEnergy;
    ForceVectorType m_totalForces;
    JacobianMatrixType m_totalJacobian;

    // Future position storage. Position solve write its result here.
    VecXd m_newDegreesOfFreedom;
};

}

#endif /* ELASTICSTRAND_HH_ */
