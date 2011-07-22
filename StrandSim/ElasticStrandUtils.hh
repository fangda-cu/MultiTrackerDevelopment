/*
 * ElasticStrandUtils.hh
 *
 *  Created on: 11/07/2011
 *      Author: jaubry
 */

#ifndef ELASTICSTRANDUTILS_HH_
#define ELASTICSTRANDUTILS_HH_

#include <limits>
#include "Definitions.hh"

namespace strandsim
{

static const Scalar SMALL_NUMBER = 1.0e-12;

inline Scalar square( const Scalar x )
{
    return x * x;
}

template<typename NormableT>
inline bool isClose( const NormableT& x1, const NormableT& x2 )
{
    return ( x1 - x2 ).norm() < SMALL_NUMBER;
}

inline bool isSmall( Scalar x )
{
    return fabs( x ) < SMALL_NUMBER;
}

template<typename MatrixT>
inline bool isSymmetric( const MatrixT& A )
{
    for ( int i = 0; i < A.rows(); ++i )
        for ( int j = i + 1; j < A.cols(); ++j )
            if ( A( i, j ) != A( j, i ) )
            {
                std::cerr << "isSymmetric failing by " << fabs( A( i, j ) - A( j, i ) ) << '\n';
                return false;
            }
    return true;
}

inline Vec3d parallelTransport( const Vec3d& u, const Vec3d& t1, const Vec3d& t2 )
{
    assert( isSmall( t1.norm() - 1 ) );
    assert( isSmall( t2.norm() - 1 ) );

    Vec3d b = t1.cross( t2 );
    if ( isSmall( b.norm() ) ) // vectors are colinear
        return u;

    b.normalize();
    const Vec3d n1 = t1.cross( b );
    const Vec3d n2 = t2.cross( b );
    return u.dot( t1 ) * t2 + u.dot( n1 ) * n2 + u.dot( b ) * b;
}

// Find an arbitrary normal (unit length) vector
template<int n>
inline Eigen::Matrix<Scalar, n, 1> findNormal( const Eigen::Matrix<Scalar, n, 1>& u )
{
    assert( u.norm() != 0 );
    Eigen::Matrix<Scalar, n, 1> v;

    int maxCoordinate = 0;
    for ( int i = 0; i < n; ++i )
    {
        if ( isSmall( u[i] ) )
        {
            v[i] = 1;
            goto finished;
        }
        if ( fabs( u[i] ) > fabs( u[maxCoordinate] ) )
            maxCoordinate = i;
    }
    {
        const int otherCoordinate = ( maxCoordinate + 1 ) % n;
        v[otherCoordinate] = u[maxCoordinate];
        v[maxCoordinate] = -u[otherCoordinate];
    }
    v.normalize();

    finished: assert( isSmall( u.dot( v ) ) );
    return v;
}

// Discrete curvature binormal between two consecutive tangent vectors (they don't need to be normalized). Note that the returned vector is not normalised but has magnitude 2 \tan(\phi/2), \phi being the angle between t1 and t2
inline Vec3d discreteCurvatureBinormal( const Vec3d& t1, const Vec3d& t2 )
{
    return 2.0 * t1.cross( t2 ) / ( t1.norm() * t2.norm() + t1.dot( t2 ) );
}

// Computes the signed angle from one vector to another given an orientation vector.
inline Scalar signedAngle( const Vec3d& u, const Vec3d& v, const Vec3d& n )
{
    Vec3d w = u.cross( v );
    Scalar angle = atan2( w.norm(), u.dot( v ) );
    if ( n.dot( w ) < 0 )
        return -angle;
    return angle;
}

// Rotates a vector about the given axis (unit vector) by the given angle.
inline void rotateAxisAngle( Vec3d& v, const Vec3d& z, const Scalar theta )
{
    assert( isSmall( z.norm() - 1.0 ) );

    if ( theta == 0 )
        return;

    const Scalar c = cos( theta );
    const Scalar s = sin( theta );

    v = c * v + s * z.cross( v ) + z.dot( v ) * ( 1.0 - c ) * z;
}

// Outer product of two vectors
template<int n>
inline Eigen::Matrix<Scalar, n, n> outerProd( const Eigen::Matrix<Scalar, n, 1>& a,
        const Eigen::Matrix<Scalar, n, 1>& b )
{
    return a * b.transpose();
}

// Matrix representation of the cross product operator
inline Mat3d crossMat( const Vec3d& a )
{
    Mat3d M;
    M << 0, -a( 2 ), a( 1 ), a( 2 ), 0, -a( 0 ), -a( 1 ), a( 0 ), 0;

    return M;
}

// Computes u B v^T, assuming B is symmetric 2x2 and u, v are 2x1 vectors.
inline Scalar BProduct( const Mat2d& B, const Vec2d& u, const Vec2d& v )
{
    assert( isSymmetric( B ) );

    // return u[0] * (B(0, 0) * v[0] + B(0, 1) * v[1]) + u[1] * (B(1, 0) * v[0] + B(1, 1) * v[1]); // Bad
    return B( 0, 0 ) * u[0] * v[0] + B( 0, 1 ) * ( u[0] * v[1] + u[1] * v[0] ) + B( 1, 1 ) * u[1]
            * v[1]; // Good
}

// Computes Q B Q^T, assuming B is symmetric 2x2 and Q is nx2. The result is then (exactly) symmetric nxn.
template<int n>
inline void symBProduct( Eigen::Matrix<Scalar, n, n>& result, const Mat2d& B,
        const Eigen::Matrix<Scalar, n, 2>& Q )
{
    assert( isSymmetric( B ) );
    assert( Q.rows() == n );
    assert( Q.cols() == 2 );

    for ( int i = 0; i < n; ++i )
    {
        result( i, i ) = BProduct( B, Q.row( i ), Q.row( i ) );
        for ( int j = 0; j < i; ++j )
            result( i, j ) = result( j, i ) = BProduct( B, Q.row( i ), Q.row( j ) );
    }
}

}
#endif /* ELASTICSTRANDUTILS_HH_ */
