/*
 * StrandGeometry.cc
 *
 *  Created on: 14/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#include "StrandGeometry.hh"

namespace strandsim
{

StrandGeometry::StrandGeometry( const VecXd& dofs ) :
    m_degreesOfFreedom( dofs ), m_framesUpToDate( false )
{
}

StrandGeometry::StrandGeometry( const size_t dofSize ) :
    m_degreesOfFreedom( VecXd( dofSize ) ), m_framesUpToDate( false )
{
}

StrandGeometry::~StrandGeometry()
{
    // TODO Auto-generated destructor stub
}

const StrandGeometry& StrandGeometry::operator=( const StrandGeometry& newGeo )
{
    assert( m_numVertices == newGeo.m_numVertices );

    m_degreesOfFreedom = VecXd( newGeo.m_degreesOfFreedom ); // Copy of the values

    m_previousTangents = newGeo.m_previousTangents; // Probably useless, temporary.

    m_framesUpToDate = newGeo.m_framesUpToDate;

    m_lengths = newGeo.m_lengths; // Probably useless, temporary.
    m_tangents = newGeo.m_tangents; // Probably useless, temporary.
    m_referenceFrames1 = newGeo.m_referenceFrames1; // Probably useless, temporary.
    m_referenceFrames2 = newGeo.m_referenceFrames2; // Probably useless, temporary.
    m_materialFrames1 = newGeo.m_materialFrames1; // Probably useless, temporary.
    m_materialFrames2 = newGeo.m_materialFrames2; // Probably useless, temporary.

    // Caches related to bending
    for ( int i = 0; i < newGeo.m_kappa.size(); i++ )
        std::cout << newGeo.m_kappa[i] << ' ';
    std::cout << '\n';
    m_kappa = newGeo.m_kappa;
    m_gradKappa = newGeo.m_gradKappa;
    m_HessKappa = newGeo.m_HessKappa; // Maybe not

    // Caches related to twisting
    m_referenceTwists = newGeo.m_referenceTwists; // Probably useless, temporary.
    m_twists = newGeo.m_twists;
    m_gradTwists = newGeo.m_gradTwists;
    m_HessTwists = newGeo.m_HessTwists;

    return *this;
}

void StrandGeometry::resizeSelf()
{
    assert( m_degreesOfFreedom.size() % 4 == 3 ); // dofs are 3 per vertex, one per edge
    assert( m_degreesOfFreedom.size() > 3 ); // minimum two vertices per rod

    m_numVertices = static_cast<IndexType> ( m_degreesOfFreedom.size() + 1 ) / 4;

    m_lengths.resize( m_numVertices - 1 );
    m_tangents.resize( m_numVertices - 1 );

    m_kappa.resize( m_numVertices - 1 );
    m_gradKappa.resize( m_numVertices - 1 );
    m_HessKappa.resize( m_numVertices - 1 );

    m_referenceTwists.resize( m_numVertices - 1 );
    m_twists.resize( m_numVertices - 1 );
    m_gradTwists.resize( m_numVertices - 1 );
    m_HessTwists.resize( m_numVertices - 1 );

    m_previousTangents.resize( 3 * ( m_numVertices - 1 ) );
    m_referenceFrames1.resize( 3 * ( m_numVertices - 1 ) );
    m_referenceFrames2.resize( 3 * ( m_numVertices - 1 ) );
    m_materialFrames1.resize( 3 * ( m_numVertices - 1 ) );
    m_materialFrames2.resize( 3 * ( m_numVertices - 1 ) );

    m_totalForces.resize(m_degreesOfFreedom.size());
}

void StrandGeometry::storeInitialFrames()
{
    // Previous tangents vector initially contains normalised edges
    for ( IndexType vtx = 0; vtx < m_numVertices - 1; ++vtx )
        setPreviousTangent( vtx, getEdgeVector( vtx ).normalized() );

    // Initial first reference frame is arbitrary
    setReferenceFrame1( 0, findNormal( getPreviousTangent( 0 ) ) );
    setReferenceFrame2( 0, getPreviousTangent( 0 ).cross( getReferenceFrame1( 0 ) ) );
    assert( fabs( getReferenceFrame2( 0 ).norm() - 1 ) < std::numeric_limits<Scalar>::epsilon() );

    // Next initial reference frames are obtained by space-parallel transportation
    for ( IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx )
    {
        setReferenceFrame1(
                vtx,
                parallelTransport( getReferenceFrame1( vtx - 1 ), getPreviousTangent( vtx - 1 ),
                        getPreviousTangent( vtx ) ) );
        setReferenceFrame2( vtx, getPreviousTangent( vtx ).cross( getReferenceFrame1( vtx ) ) );
    }

    updateFrames();
}

void StrandGeometry::updateFrames() // and related stuff
{
    if ( m_framesUpToDate )
        return;

    // Update reference frames by time-parallel transportation
    for ( IndexType vtx = 0; vtx < m_numVertices - 1; ++vtx )
    {
        m_lengths[vtx] = getEdgeVector( vtx ).norm(); // std::cerr << "m_lengths[vtx] = " <<  m_lengths[vtx] <<'\n';
        assert( !isSmall( m_lengths[vtx] ) );

        m_tangents[vtx] = getEdgeVector( vtx ) / m_lengths[vtx];

        Vec3d u = parallelTransport( getReferenceFrame1( vtx ), getPreviousTangent( vtx ),
                m_tangents[vtx] );
        u -= u.dot( m_tangents[vtx] ) * m_tangents[vtx];
        u.normalize();
        const Vec3d& v = m_tangents[vtx].cross( u );
        // std::cout << "t = " << m_tangents[vtx] << '\n';
        // std::cout << "u = " << u << '\n';
        // std::cout << "v = " << v << '\n';
        setReferenceFrame1( vtx, u );
        setReferenceFrame2( vtx, v );

        const Scalar theta = m_degreesOfFreedom[4 * vtx + 3];
        const Scalar c = cos( theta );
        const Scalar s = sin( theta );
        setMaterialFrame1( vtx, c * u + s * v );
        setMaterialFrame2( vtx, -s * u + c * v );

        setPreviousTangent( vtx, m_tangents[vtx] ); // TODO: make sure that previous tangents are not used elsewhere
    }

    updateKappaAndTwist();

    m_framesUpToDate = true;
}

void StrandGeometry::updateKappaAndTwist()
{
    // Update other cached frame-dependent quantities
    for ( IndexType vtx = 1; vtx < m_numVertices - 1; ++vtx )
    {
        m_kappa[vtx] = computeKappa( vtx );
        m_gradKappa[vtx] = computeGradKappa( vtx );
        m_HessKappa[vtx] = computeHessKappa( vtx ); // This caching may not be necessary as HessKappa will be used only once.

        m_referenceTwists[vtx] = computeReferenceTwist( vtx ); // Idem: does this really need to be cached?
        m_twists[vtx] = computeTwist( vtx );
        m_gradTwists[vtx] = computeGradTwist( vtx );
        m_HessTwists[vtx] = computeHessTwist( vtx ); // Idem.
    }
}

Vec2d StrandGeometry::computeKappa( const IndexType vtx ) const
{
    const Vec3d& kb = discreteCurvatureBinormal( getEdgeVector( vtx - 1 ), getEdgeVector( vtx ) ); // TODO: cache
    // std::cout << "kb = " << kb << '\n';
    const Vec3d& m1e = getMaterialFrame1( vtx - 1 );
    // std::cout << "m1e = " << m1e << '\n';
    const Vec3d& m2e = getMaterialFrame2( vtx - 1 );
    // std::cout << "m2e = " << m2e << '\n';
    const Vec3d& m1f = getMaterialFrame1( vtx );
    // std::cout << "m1f = " << m1f << '\n';
    const Vec3d& m2f = getMaterialFrame2( vtx );
    // std::cout << "m2f = " << m2f << '\n';

    return Vec2d( 0.5 * kb.dot( m2e + m2f ), -0.5 * kb.dot( m1e + m1f ) );
}

Eigen::Matrix<Scalar, 11, 2> StrandGeometry::computeGradKappa( const IndexType vtx ) const
{
    Eigen::Matrix<Scalar, 11, 2> gradKappa;

    const Scalar norm_e = m_lengths[vtx - 1];
    const Scalar norm_f = m_lengths[vtx];

    const Vec3d& te = m_tangents[vtx - 1];
    const Vec3d& tf = m_tangents[vtx];

    const Vec3d& m1e = getMaterialFrame1( vtx - 1 );
    const Vec3d& m2e = getMaterialFrame2( vtx - 1 );
    const Vec3d& m1f = getMaterialFrame1( vtx );
    const Vec3d& m2f = getMaterialFrame2( vtx );

    const Scalar chi = 1.0 + te.dot( tf );
    const Vec3d& tilde_t = ( te + tf ) / chi;
    const Vec3d& tilde_d1 = ( m1e + m1f ) / chi;
    const Vec3d& tilde_d2 = ( m2e + m2f ) / chi;

    const Vec2d& kappa = m_kappa[vtx];

    const Vec3d& Dkappa0De = 1.0 / norm_e * ( -kappa[0] * tilde_t + tf.cross( tilde_d2 ) );
    const Vec3d& Dkappa0Df = 1.0 / norm_f * ( -kappa[0] * tilde_t - te.cross( tilde_d2 ) );
    const Vec3d& Dkappa1De = 1.0 / norm_e * ( -kappa[1] * tilde_t - tf.cross( tilde_d1 ) );
    const Vec3d& Dkappa1Df = 1.0 / norm_f * ( -kappa[1] * tilde_t + te.cross( tilde_d1 ) );

    gradKappa.block<3, 1> ( 0, 0 ) = -Dkappa0De;
    gradKappa.block<3, 1> ( 4, 0 ) = Dkappa0De - Dkappa0Df;
    gradKappa.block<3, 1> ( 8, 0 ) = Dkappa0Df;
    gradKappa.block<3, 1> ( 0, 1 ) = -Dkappa1De;
    gradKappa.block<3, 1> ( 4, 1 ) = Dkappa1De - Dkappa1Df;
    gradKappa.block<3, 1> ( 8, 1 ) = Dkappa1Df;

    const Vec3d& kb = discreteCurvatureBinormal( getEdgeVector( vtx - 1 ), getEdgeVector( vtx ) ); // TODO: cache

    gradKappa( 3, 0 ) = -0.5 * kb.dot( m1e );
    gradKappa( 7, 0 ) = -0.5 * kb.dot( m1f );
    gradKappa( 3, 1 ) = -0.5 * kb.dot( m2e );
    gradKappa( 7, 1 ) = -0.5 * kb.dot( m2f );

    return gradKappa;
}

Mat11dPair StrandGeometry::computeHessKappa( const IndexType vtx ) const
{
    Mat11dPair HessKappa;

    Mat11d& DDkappa1 = HessKappa.first;
    Mat11d& DDkappa2 = HessKappa.second;

    const Scalar norm_e = m_lengths[vtx - 1];
    const Scalar norm_f = m_lengths[vtx];
    const Scalar norm2_e = square( norm_e ); // That's bloody stupid, taking the square of a square root.
    const Scalar norm2_f = square( norm_f );

    const Vec3d& te = m_tangents[vtx - 1];
    const Vec3d& tf = m_tangents[vtx];

    const Vec3d& m1e = getMaterialFrame1( vtx - 1 );
    const Vec3d& m2e = getMaterialFrame2( vtx - 1 );
    const Vec3d& m1f = getMaterialFrame1( vtx );
    const Vec3d& m2f = getMaterialFrame2( vtx );

    const Scalar chi = 1.0 + te.dot( tf );
    const Vec3d& tilde_t = ( te + tf ) / chi;
    const Vec3d& tilde_d1 = ( m1e + m1f ) / chi;
    const Vec3d& tilde_d2 = ( m2e + m2f ) / chi;

    const Vec2d& kappa = m_kappa[vtx];

    const Vec3d& kb = discreteCurvatureBinormal( getEdgeVector( vtx - 1 ), getEdgeVector( vtx ) ); // TODO: cache

    const Mat3d& tt_o_tt = outerProd( tilde_t, tilde_t );
    const Mat3d& tf_c_d2t_o_tt = outerProd( tf.cross( tilde_d2 ), tilde_t );
    const Mat3d& tt_o_tf_c_d2t = tf_c_d2t_o_tt.transpose();
    const Mat3d& kb_o_d2e = outerProd( kb, m2e );
    const Mat3d& d2e_o_kb = kb_o_d2e.transpose();

    const Mat3d& Id = Mat3d::Identity();

    const Mat3d& D2kappa1De2 = 1.0 / norm2_e * ( 2 * kappa[0] * tt_o_tt - ( tf_c_d2t_o_tt
            + tt_o_tf_c_d2t ) ) - kappa[0] / ( chi * norm2_e ) * ( Id - outerProd( te, te ) ) + 1.0
            / ( 4.0 * norm2_e ) * ( kb_o_d2e + d2e_o_kb );

    const Mat3d& te_c_d2t_o_tt = outerProd( te.cross( tilde_d2 ), tilde_t );
    const Mat3d& tt_o_te_c_d2t = te_c_d2t_o_tt.transpose();
    const Mat3d& kb_o_d2f = outerProd( kb, m2f );
    const Mat3d& d2f_o_kb = kb_o_d2f.transpose();

    const Mat3d& D2kappa1Df2 = 1.0 / norm2_f * ( 2 * kappa[0] * tt_o_tt + ( te_c_d2t_o_tt
            + tt_o_te_c_d2t ) ) - kappa[0] / ( chi * norm2_f ) * ( Id - outerProd( tf, tf ) ) + 1.0
            / ( 4.0 * norm2_f ) * ( kb_o_d2f + d2f_o_kb );

    const Mat3d& D2kappa1DeDf = -kappa[0] / ( chi * norm_e * norm_f ) * ( Id + outerProd( te, tf ) )
            + 1.0 / ( norm_e * norm_f ) * ( 2 * kappa[0] * tt_o_tt - tf_c_d2t_o_tt + tt_o_te_c_d2t
                    - crossMat( tilde_d2 ) );
    const Mat3d& D2kappa1DfDe = D2kappa1DeDf.transpose();

    const Mat3d& tf_c_d1t_o_tt = outerProd( tf.cross( tilde_d1 ), tilde_t );
    const Mat3d& tt_o_tf_c_d1t = tf_c_d1t_o_tt.transpose();
    const Mat3d& kb_o_d1e = outerProd( kb, m1e );
    const Mat3d& d1e_o_kb = kb_o_d1e.transpose();

    const Mat3d& D2kappa2De2 = 1.0 / norm2_e * ( 2 * kappa[1] * tt_o_tt + ( tf_c_d1t_o_tt
            + tt_o_tf_c_d1t ) ) - kappa[1] / ( chi * norm2_e ) * ( Id - outerProd( te, te ) ) - 1.0
            / ( 4.0 * norm2_e ) * ( kb_o_d1e + d1e_o_kb );

    const Mat3d& te_c_d1t_o_tt = outerProd( te.cross( tilde_d1 ), tilde_t );
    const Mat3d& tt_o_te_c_d1t = te_c_d1t_o_tt.transpose();
    const Mat3d& kb_o_d1f = outerProd( kb, m1f );
    const Mat3d& d1f_o_kb = kb_o_d1f.transpose();

    const Mat3d& D2kappa2Df2 = 1.0 / norm2_f * ( 2 * kappa[1] * tt_o_tt - ( te_c_d1t_o_tt
            + tt_o_te_c_d1t ) ) - kappa[1] / ( chi * norm2_f ) * ( Id - outerProd( tf, tf ) ) - 1.0
            / ( 4.0 * norm2_f ) * ( kb_o_d1f + d1f_o_kb );

    const Mat3d& D2kappa2DeDf = -kappa[1] / ( chi * norm_e * norm_f ) * ( Id + outerProd( te, tf ) )
            + 1.0 / ( norm_e * norm_f ) * ( 2 * kappa[1] * tt_o_tt + tf_c_d1t_o_tt - tt_o_te_c_d1t
                    + crossMat( tilde_d1 ) );
    const Mat3d& D2kappa2DfDe = D2kappa2DeDf.transpose();

    const Scalar D2kappa1Dthetae2 = -0.5 * kb.dot( m2e );
    const Scalar D2kappa1Dthetaf2 = -0.5 * kb.dot( m2f );
    const Scalar D2kappa2Dthetae2 = 0.5 * kb.dot( m1e );
    const Scalar D2kappa2Dthetaf2 = 0.5 * kb.dot( m1f );

    const Vec3d& D2kappa1DeDthetae = 1.0 / norm_e * ( 0.5 * kb.dot( m1e ) * tilde_t - 1.0 / chi
            * tf.cross( m1e ) );
    const Vec3d& D2kappa1DeDthetaf = 1.0 / norm_e * ( 0.5 * kb.dot( m1f ) * tilde_t - 1.0 / chi
            * tf.cross( m1f ) );
    const Vec3d& D2kappa1DfDthetae = 1.0 / norm_f * ( 0.5 * kb.dot( m1e ) * tilde_t + 1.0 / chi
            * te.cross( m1e ) );
    const Vec3d& D2kappa1DfDthetaf = 1.0 / norm_f * ( 0.5 * kb.dot( m1f ) * tilde_t + 1.0 / chi
            * te.cross( m1f ) );
    const Vec3d& D2kappa2DeDthetae = 1.0 / norm_e * ( 0.5 * kb.dot( m2e ) * tilde_t - 1.0 / chi
            * tf.cross( m2e ) );
    const Vec3d& D2kappa2DeDthetaf = 1.0 / norm_e * ( 0.5 * kb.dot( m2f ) * tilde_t - 1.0 / chi
            * tf.cross( m2f ) );
    const Vec3d& D2kappa2DfDthetae = 1.0 / norm_f * ( 0.5 * kb.dot( m2e ) * tilde_t + 1.0 / chi
            * te.cross( m2e ) );
    const Vec3d& D2kappa2DfDthetaf = 1.0 / norm_f * ( 0.5 * kb.dot( m2f ) * tilde_t + 1.0 / chi
            * te.cross( m2f ) );

    DDkappa1.block<3, 3> ( 0, 0 ) = D2kappa1De2;
    DDkappa1.block<3, 3> ( 0, 4 ) = -D2kappa1De2 + D2kappa1DeDf;
    DDkappa1.block<3, 3> ( 4, 0 ) = -D2kappa1De2 + D2kappa1DfDe;
    DDkappa1.block<3, 3> ( 4, 4 ) = D2kappa1De2 - ( D2kappa1DeDf + D2kappa1DfDe ) + D2kappa1Df2;
    DDkappa1.block<3, 3> ( 0, 8 ) = -D2kappa1DeDf;
    DDkappa1.block<3, 3> ( 8, 0 ) = -D2kappa1DfDe;
    DDkappa1.block<3, 3> ( 4, 8 ) = D2kappa1DeDf - D2kappa1Df2;
    DDkappa1.block<3, 3> ( 8, 4 ) = D2kappa1DfDe - D2kappa1Df2;
    DDkappa1.block<3, 3> ( 8, 8 ) = D2kappa1Df2;
    DDkappa1( 3, 3 ) = D2kappa1Dthetae2;
    DDkappa1( 7, 7 ) = D2kappa1Dthetaf2;
    DDkappa1.block<3, 1> ( 0, 3 ) = -D2kappa1DeDthetae;
    DDkappa1.block<1, 3> ( 3, 0 ) = DDkappa1.block<3, 1> ( 0, 3 ).transpose();
    DDkappa1.block<3, 1> ( 4, 3 ) = D2kappa1DeDthetae - D2kappa1DfDthetae;
    DDkappa1.block<1, 3> ( 3, 4 ) = DDkappa1.block<3, 1> ( 4, 3 ).transpose();
    DDkappa1.block<3, 1> ( 8, 3 ) = D2kappa1DfDthetae;
    DDkappa1.block<1, 3> ( 3, 8 ) = DDkappa1.block<3, 1> ( 8, 3 ).transpose();
    DDkappa1.block<3, 1> ( 0, 7 ) = -D2kappa1DeDthetaf;
    DDkappa1.block<1, 3> ( 7, 0 ) = DDkappa1.block<3, 1> ( 0, 7 ).transpose();
    DDkappa1.block<3, 1> ( 4, 7 ) = D2kappa1DeDthetaf - D2kappa1DfDthetaf;
    DDkappa1.block<1, 3> ( 7, 4 ) = DDkappa1.block<3, 1> ( 4, 7 ).transpose();
    DDkappa1.block<3, 1> ( 8, 7 ) = D2kappa1DfDthetaf;
    DDkappa1.block<1, 3> ( 7, 8 ) = DDkappa1.block<3, 1> ( 8, 7 ).transpose();

    assert( isSymmetric( DDkappa1 ) );

    DDkappa2.block<3, 3> ( 0, 0 ) = D2kappa2De2;
    DDkappa2.block<3, 3> ( 0, 4 ) = -D2kappa2De2 + D2kappa2DeDf;
    DDkappa2.block<3, 3> ( 4, 0 ) = -D2kappa2De2 + D2kappa2DfDe;
    DDkappa2.block<3, 3> ( 4, 4 ) = D2kappa2De2 - ( D2kappa2DeDf + D2kappa2DfDe ) + D2kappa2Df2;
    DDkappa2.block<3, 3> ( 0, 8 ) = -D2kappa2DeDf;
    DDkappa2.block<3, 3> ( 8, 0 ) = -D2kappa2DfDe;
    DDkappa2.block<3, 3> ( 4, 8 ) = D2kappa2DeDf - D2kappa2Df2;
    DDkappa2.block<3, 3> ( 8, 4 ) = D2kappa2DfDe - D2kappa2Df2;
    DDkappa2.block<3, 3> ( 8, 8 ) = D2kappa2Df2;
    DDkappa2( 3, 3 ) = D2kappa2Dthetae2;
    DDkappa2( 7, 7 ) = D2kappa2Dthetaf2;
    DDkappa2.block<3, 1> ( 0, 3 ) = -D2kappa2DeDthetae;
    DDkappa2.block<1, 3> ( 3, 0 ) = DDkappa2.block<3, 1> ( 0, 3 ).transpose();
    DDkappa2.block<3, 1> ( 4, 3 ) = D2kappa2DeDthetae - D2kappa2DfDthetae;
    DDkappa2.block<1, 3> ( 3, 4 ) = DDkappa2.block<3, 1> ( 4, 3 ).transpose();
    DDkappa2.block<3, 1> ( 8, 3 ) = D2kappa2DfDthetae;
    DDkappa2.block<1, 3> ( 3, 8 ) = DDkappa2.block<3, 1> ( 8, 3 ).transpose();
    DDkappa2.block<3, 1> ( 0, 7 ) = -D2kappa2DeDthetaf;
    DDkappa2.block<1, 3> ( 7, 0 ) = DDkappa2.block<3, 1> ( 0, 7 ).transpose();
    DDkappa2.block<3, 1> ( 4, 7 ) = D2kappa2DeDthetaf - D2kappa2DfDthetaf;
    DDkappa2.block<1, 3> ( 7, 4 ) = DDkappa2.block<3, 1> ( 4, 7 ).transpose();
    DDkappa2.block<3, 1> ( 8, 7 ) = D2kappa2DfDthetaf;
    DDkappa2.block<1, 3> ( 7, 8 ) = DDkappa2.block<3, 1> ( 8, 7 ).transpose();

    assert( isSymmetric( DDkappa2 ) );

    return HessKappa;
}

Scalar StrandGeometry::computeReferenceTwist( const IndexType vtx ) const
{
    Scalar referenceTwist = m_referenceTwists[vtx];

    const Vec3d& u0 = getReferenceFrame1( vtx - 1 );
    const Vec3d& u1 = getReferenceFrame1( vtx );
    const Vec3d& tangent = m_tangents[vtx];

    // transport reference frame to next edge
    Vec3d ut = parallelTransport( u0, m_tangents[vtx - 1], tangent );

    // rotate by current value of reference twist
    rotateAxisAngle( ut, tangent, referenceTwist );

    // compute increment to reference twist to align reference frames
    referenceTwist += signedAngle( ut, u1, tangent );

    return referenceTwist;
}

Scalar StrandGeometry::computeTwist( const IndexType vtx ) const
{
    return m_referenceTwists[vtx] + getTheta( vtx ) - getTheta( vtx - 1 );
}

Vec11d StrandGeometry::computeGradTwist( const IndexType vtx ) const
{
    Vec11d Dtwist;

    const Vec3d& kb = discreteCurvatureBinormal( getEdgeVector( vtx - 1 ), getEdgeVector( vtx ) ); // TODO: cache

    Dtwist.segment<3> ( 0 ) = -0.5 / m_lengths[vtx - 1] * kb;
    Dtwist.segment<3> ( 8 ) = 0.5 / m_lengths[vtx] * kb;
    Dtwist.segment<3> ( 4 ) = -( Dtwist.segment<3> ( 0 ) + Dtwist.segment<3> ( 8 ) );
    Dtwist( 3 ) = -1;
    Dtwist( 7 ) = 1;

    return Dtwist;
}

Mat11d StrandGeometry::computeHessTwist( const IndexType vtx ) const
{
    Mat11d DDtwist;

    const Vec3d& te = m_tangents[vtx - 1];
    const Vec3d& tf = m_tangents[vtx];
    const Scalar norm_e = m_lengths[vtx - 1];
    const Scalar norm_f = m_lengths[vtx];
    const Vec3d& kb = discreteCurvatureBinormal( getEdgeVector( vtx - 1 ), getEdgeVector( vtx ) ); // TODO: cache

    const Scalar chi = 1 + te.dot( tf );
    const Vec3d& tilde_t = 1.0 / chi * ( te + tf );

    const Mat3d& D2mDe2 = -0.25 / square( norm_e ) * ( outerProd<3> ( kb, te + tilde_t )
            + outerProd<3> ( te + tilde_t, kb ) );
    const Mat3d& D2mDf2 = -0.25 / square( norm_f ) * ( outerProd<3> ( kb, tf + tilde_t )
            + outerProd<3> ( tf + tilde_t, kb ) );
    const Mat3d& D2mDeDf = 0.5 / ( norm_e * norm_f ) * ( 2.0 / chi * crossMat( te ) - outerProd(
            kb, tilde_t ) );
    const Mat3d& D2mDfDe = D2mDeDf.transpose();

    DDtwist.block<3, 3> ( 0, 0 ) = D2mDe2;
    DDtwist.block<3, 3> ( 0, 4 ) = -D2mDe2 + D2mDeDf;
    DDtwist.block<3, 3> ( 4, 0 ) = -D2mDe2 + D2mDfDe;
    DDtwist.block<3, 3> ( 4, 4 ) = D2mDe2 - ( D2mDeDf + D2mDfDe ) + D2mDf2;
    DDtwist.block<3, 3> ( 0, 8 ) = -D2mDeDf;
    DDtwist.block<3, 3> ( 8, 0 ) = -D2mDfDe;
    DDtwist.block<3, 3> ( 8, 4 ) = D2mDfDe - D2mDf2;
    DDtwist.block<3, 3> ( 4, 8 ) = D2mDeDf - D2mDf2;
    DDtwist.block<3, 3> ( 8, 8 ) = D2mDf2;
    assert( isSymmetric( DDtwist ) );

    return DDtwist;
}

}
