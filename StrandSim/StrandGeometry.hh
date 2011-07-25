/*
 * StrandGeometry.hh
 *
 *  Created on: 14/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef STRANDGEOMETRY_HH_
#define STRANDGEOMETRY_HH_

#include "Definitions.hh"
#include "ElasticStrandUtils.hh"

namespace strandsim
{

// This class contains all the rod's variable data (i.e. changed by simulation)
class StrandGeometry
{
public:
    explicit StrandGeometry( const VecXd& dofs );

    explicit StrandGeometry( const size_t dofSize );

    ~StrandGeometry();

    void storeInitialFrames();

    const Vec3d getVertex( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices );

        return m_degreesOfFreedom.segment<3> ( 4 * vtx );
    }

    void setVertex( const IndexType vtx, const Vec3d coordinates )
    {
        assert( vtx < m_numVertices );

        m_degreesOfFreedom.segment<3> ( 4 * vtx ) = coordinates;
    }

    Scalar getTheta( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices - 1 );

        return m_degreesOfFreedom[4 * vtx + 3];
    }

    const Vec3d getEdgeVector( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices - 1 );

        return getVertex( vtx + 1 ) - getVertex( vtx );
    }

    const Vec3d getReferenceFrame1( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices - 1 );

        return m_referenceFrames1.segment<3> ( 3 * vtx );
    }

    void setReferenceFrame1( const IndexType vtx, const Vec3d& vec )
    {
        assert( vtx < m_numVertices - 1 );

        m_referenceFrames1.segment<3> ( 3 * vtx ) = vec;
    }

    const Vec3d getReferenceFrame2( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices - 1 );

        return m_referenceFrames2.segment<3> ( 3 * vtx );
    }

    void setReferenceFrame2( const IndexType vtx, const Vec3d& vec )
    {
        assert( vtx < m_numVertices - 1 );

        m_referenceFrames2.segment<3> ( 3 * vtx ) = vec;
    }

    const Vec3d getMaterialFrame1( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices - 1 );

        return m_materialFrames1.segment<3> ( 3 * vtx );
    }

    void setMaterialFrame1( const IndexType vtx, const Vec3d& vec )
    {
        assert( vtx < m_numVertices - 1 );

        m_materialFrames1.segment<3> ( 3 * vtx ) = vec;
    }

    const Vec3d getMaterialFrame2( const IndexType vtx ) const
    {
        assert( vtx < m_numVertices - 1 );

        return m_materialFrames2.segment<3> ( 3 * vtx );
    }

    void setMaterialFrame2( const IndexType vtx, const Vec3d& vec )
    {
        assert( vtx < m_numVertices - 1 );

        m_materialFrames2.segment<3> ( 3 * vtx ) = vec;
    }

    const Vec3d getPreviousTangent( const IndexType vtx )
    {
        assert( vtx < m_numVertices - 1 );

        return m_previousTangents.segment<3> ( 3 * vtx );
    }

    void setPreviousTangent( const IndexType vtx, const Vec3d& vec )
    {
        assert( vtx < m_numVertices - 1 );

        m_previousTangents.segment<3> ( 3 * vtx ) = vec;
    }

    void resizeSelf();
    void updateFrames();
    void updateKappaAndTwist();

    void computeKappa( Vec2d& kappa, const IndexType vtx ) const;
    void computeGradKappa( Eigen::Matrix<Scalar, 11, 2> & gradKappa, const IndexType vtx ) const;
    void computeHessKappa( Mat11dPair &, const IndexType vtx ) const;

    Scalar computeReferenceTwist( const IndexType vtx ) const;
    Scalar computeTwist( const IndexType vtx ) const;
    void computeGradTwist( Vec11d& gradTwist, const IndexType vtx ) const;
    void computeHessTwist( Mat11d &, const IndexType vtx ) const;

    IndexType m_numVertices;

    // The actual positions (3 per vertex + 1 per edge)
    VecXd m_degreesOfFreedom; // size = m_numVertices * 4 -1

    // Previous time storage, for time-parallel stepping
    VecXd m_previousTangents;

    // Reference frames, material frames and other geometric properties
    bool m_framesUpToDate;
    std::vector<Scalar> m_lengths; // lengths of edges, cached
    std::vector<Vec3d, Eigen::aligned_allocator<Vec3d> > m_tangents;
    std::vector<Vec3d, Eigen::aligned_allocator<Vec3d> > m_curvatureBinormals;
    VecXd m_referenceFrames1; // first vectors of reference frame
    VecXd m_referenceFrames2; // second vectors of reference frame
    VecXd m_materialFrames1;// first vectors of material frame
    VecXd m_materialFrames2; // second vectors of material frame

    // Caches related to bending
    std::vector<Vec2d, Eigen::aligned_allocator<Vec2d> > m_kappa;
    std::vector<Eigen::Matrix<Scalar, 11, 2>,
            Eigen::aligned_allocator<Eigen::Matrix<Scalar, 11, 2> > > m_gradKappa;

    // Caches related to twisting
    std::vector<Scalar> m_referenceTwists;
    std::vector<Scalar> m_twists;
    std::vector<Vec11d> m_gradTwists;

    // Energy, force, Jacobian
    Scalar m_totalEnergy;
    VecXd m_totalForces;

};

}

#endif /* STRANDGEOMETRY_HH_ */
