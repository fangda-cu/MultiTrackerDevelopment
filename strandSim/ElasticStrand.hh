/*
 * ElasticStrand.hh
 *
 *  Created on: 7/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef ELASTICSTRAND_HH_
#define ELASTICSTRAND_HH_

#include "StrandBase.hh"
#include "StrandGeometry.hh"
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

    ElasticStrand(VecXd& dofs, const ParametersType& parameters);

    virtual ~ElasticStrand();

    const VecXd& getDegreesOfFreedom() const
    {
        return m_geometry.m_degreesOfFreedom;
    }

    // Expose the future position storage so it can be changed by the stepper.
    VecXd& getNewDegreesOfFreedom()
    {
        return m_newGeometry.m_degreesOfFreedom;
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
    void freezeRestShape();

    Mat2d computeBendingMatrix(const IndexType vtx) const;
    void updateForceCaches();
    void accumulateInternalForces();
    void accumulateExternalForces();

    /**
     * Member variables
     */

    // Size of the strand
    IndexType m_numVertices;

    // Storage structure for a copy of the degrees of freedom vector, because m_newGeometry wants a reference.
    VecXd m_dofsStorage;

    // Current and tentative geometry
    StrandGeometry m_geometry;
    StrandGeometry m_newGeometry;

    // Other physical parameters
    ParametersType m_parameters;
    std::vector<Mat2d, Eigen::aligned_allocator<Mat2d> > m_bendingMatrices;

    // Rest shape
    std::vector<Scalar> m_restLengths;
    std::vector<Vec2d, Eigen::aligned_allocator<Vec2d> > m_kappaBar;
    std::vector<Scalar> m_restTwists;

    // Force caching
    bool m_forceCachesUpToDate;
    Scalar m_totalEnergy;
    ForceVectorType m_totalForces;
    JacobianMatrixType m_totalJacobian;

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
