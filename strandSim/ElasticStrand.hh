/*
 * ElasticStrand.hh
 *
 *  Created on: 7/07/2011
 *      Author: jaubry
 */

#ifndef ELASTICSTRAND_HH_
#define ELASTICSTRAND_HH_

#include "StrandBase.hh"
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
        assert(0 <= vtx && vtx < m_numVertices);

        return m_degreesOfFreedom.segment<3> (4 * vtx);
    }

    Scalar getTorsion(const IndexType vtx) const
    {
        assert(0 <= vtx && vtx < m_numVertices - 1);

        return m_degreesOfFreedom[4 * vtx + 3];
    }

    const Edge<ElasticStrand> getEdge(const IndexType firstVtx) const
    {
        assert(0 <= vtx && vtx < m_numVertices - 1);

        return Edge<ElasticStrand> (*this, firstVtx);
    }

    void growLastVertex(const Vec3d& vertex, const Scalar torsion);
    void popLastVertex();

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

    void updateReferenceFrames();

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

    // Reference frames
    bool m_referenceFramesUpToDate;
    // other stuff here

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
