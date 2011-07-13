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
#include "BandMatrix.hh"
#include "ForceAccumulator.hh"
#include "ElasticStrandParameters.hh"

namespace strandsim
{

class StretchingForce;
class BendingForce;
class TwistingForce;

class ElasticStrand: public StrandBase
{
public:
    typedef ElasticStrandParameters ParametersType;
    typedef VecXd ForceVectorType;
    typedef strandsim::BandMatrix<Scalar, 10, 10> JacobianMatrixType; // TODO: replace this with a type that has built-in symmetry

    ElasticStrand(VecXd& dofs, ParametersType parameters);

    virtual ~ElasticStrand();

    void storeInitialFrames();

    const Vec3d getVertex(const IndexType vtx) const
    {
        assert(vtx < m_numVertices);

        return m_degreesOfFreedom.segment<3> (4 * vtx);
    }

    Scalar getTheta(const IndexType vtx) const
    {
        assert(vtx < m_numVertices - 1);

        return m_degreesOfFreedom[4 * vtx + 3];
    }

    const Vec3d getEdgeVector(const IndexType vtx) const
    {
        assert(vtx < m_numVertices - 1);

        return getVertex(vtx + 1) - getVertex(vtx);
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

    const Vec3d getPreviousTangent(const IndexType vtx)
    {
        assert(vtx < m_numVertices - 1);

        return m_previousTangents.segment<3> (3 * vtx);
    }

    void setPreviousTangent(const IndexType vtx, const Vec3d& vec)
    {
        assert(vtx < m_numVertices - 1);

        m_previousTangents.segment<3> (3 * vtx) = vec;
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
    typedef Eigen::Matrix<Scalar, 11, 1> Vec11d;
    typedef Eigen::Matrix<Scalar, 11, 11> Mat11d;
    typedef std::pair<Mat11d, Mat11d> Mat11dPair;

    void resizeInternals();
    void freezeRestShape();
    void updateFrames();

    Mat2d computeBendingMatrix(const IndexType vtx) const;
    Vec2d computeKappa(const IndexType vtx) const;
    Eigen::Matrix<Scalar, 11, 2> computeGradKappa(const IndexType vtx) const;
    Mat11dPair computeHessKappa(const IndexType vtx) const;

    Scalar computeReferenceTwist(const IndexType vtx) const;
    Scalar computeTwist(const IndexType vtx) const;
    Vec11d computeGradTwist(const IndexType vtx) const;
    Mat11d computeHessTwist(const IndexType vtx) const;

    void updateForceCaches();
    void accumulateInternalForces();
    void accumulateExternalForces();

    /**
     * Member variables
     */

    // Size of the strand
    IndexType m_numVertices;

    // References passed by Maya, who owns the original
    VecXd& m_degreesOfFreedom; // size = m_numVertices * 4 -1

    // Other physical parameters
    ParametersType m_parameters;

    // Rest shape
    std::vector<Scalar> m_restLengths;
    std::vector<Vec2d, Eigen::aligned_allocator<Vec2d> > m_kappaBar;
    std::vector<Scalar> m_restTwists;

    // Previous time storage, for time-parallel stepping
    VecXd m_previousTangents;

    // Reference frames, material frames and other geometric properties
    bool m_framesUpToDate;
    std::vector<Scalar> m_lengths; // lengths of edges, cached
    std::vector<Vec3d, Eigen::aligned_allocator<Vec3d> > m_tangents;
    VecXd m_referenceFrames1; // first vectors of reference frame
    VecXd m_referenceFrames2; // second vectors of reference frame
    VecXd m_materialFrames1;// first vectors of material frame
    VecXd m_materialFrames2; // second vectors of material frame

    // Caches related to bending
    std::vector<Mat2d, Eigen::aligned_allocator<Mat2d> > m_bendingMatrices;
    std::vector<Vec2d, Eigen::aligned_allocator<Vec2d> > m_kappa;
    std::vector<Eigen::Matrix<Scalar, 11, 2>, Eigen::aligned_allocator<Eigen::Matrix<Scalar, 11, 2> > > m_gradKappa;
    std::vector<Mat11dPair, Eigen::aligned_allocator<Mat11dPair> > m_HessKappa; // Maybe not

    // Caches related to twisting
    std::vector<Scalar> m_referenceTwists;
    std::vector<Scalar> m_twists;
    std::vector<Vec11d> m_gradTwists;
    std::vector<Mat11d> m_HessTwists;

    // Force caching
    bool m_forceCachesUpToDate;
    Scalar m_totalEnergy;
    ForceVectorType m_totalForces;
    JacobianMatrixType m_totalJacobian;

    // Future position storage. Position solve writes its result here.
    VecXd m_newDegreesOfFreedom;

public:
    friend class StretchingForce;
    friend class ForceAccumulator<StretchingForce> ;
    friend class BendingForce;
    friend class ForceAccumulator<BendingForce> ;
    friend class TwistingForce;
    friend class ForceAccumulator<TwistingForce> ;
};

}

#endif /* ELASTICSTRAND_HH_ */
