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

template<typename ScalarT>
inline ScalarT square(const ScalarT x)
{
    return x * x;
}

template<typename NormableT>
inline bool isClose(const NormableT& x1, const NormableT& x2)
{
    return (x1 - x2).norm() < std::numeric_limits<ScalarT>::epsilon();
}

template<typename ScalarT>
inline bool isSmall(ScalarT x)
{
    return fabs(x) < std::numeric_limits<ScalarT>::epsilon();
}

template<typename MatrixT>
inline bool isSymmetric(const MatrixT& A)
{
    for (int i = 0; i < A.rows(); ++i)
        for (int j = i + 1; j < A.cols(); ++j)
            if (A(i, j) != A(j, i))
            {
                std::cerr << "isSymmetric failing by " << fabs(A(i, j) - A(j, i)) << '\n';
                return false;
            }
    return true;
}

inline Vec3d parallelTransport(const Vec3d& u, const Vec3d& t1, const Vec3d& t2)
{
    assert(isSmall(t1.norm() - 1));
    assert(isSmall(t2.norm() - 1));

    Vec3d b = t1.cross(t2);
    if (isSmall(b.norm())) // vectors are colinear
        return u;

    b.normalize();
    const Vec3d n1 = t1.cross(b);
    const Vec3d n2 = t2.cross(b);
    return u.dot(t1) * t2 + u.dot(n1) * n2 + u.dot(b) * b;
}

inline Vec3d findNormal(const Vec3d& u)
{
    assert(u.norm() != 0);
    Vec3d v;

    int maxCoordinate = 0;
    for (int i = 0; i < 3; ++i)
    {
        if (isSmall(u[i]))
        {
            v[i] = 1;
            goto finished;
        }
        if (fabs(u[i]) > fabs(u[maxCoordinate]))
            maxCoordinate = i;
    }
    {
        const int idx = (maxCoordinate + 1) % 3;
        v[idx] = u[maxCoordinate];
        v[maxCoordinate] = -u[idx];
    }

    finished: assert(isSmall(u.dot(v)));
    return v;
}

// Computes u B v^T, assuming B is symmetric 2x2 and u, v are 2x1 vectors.
inline Scalar BProduct(const Mat2d& B, const Vec2d& u, const Vec2d& v)
{
    assert(isSymmetric(B));

    // return u[0] * (B(0, 0) * v[0] + B(0, 1) * v[1]) + u[1] * (B(1, 0) * v[0] + B(1, 1) * v[1]); // Bad
    return B(0, 0) * u[0] * v[0] + B(0, 1) * (u[0] * v[1] + u[1] * v[0]) + B(1, 1) * u[1] * v[1]; // Good
}

// Computes Q B Q^T, assuming B is symmetric 2x2 and Q is nx2. The result is then (exactly) symmetric nxn.
template<int n>
inline void symBProduct(Eigen::Matrix<Scalar, n, n>& result, const Mat2d& B, const Eigen::Matrix<Scalar, n, 2>& Q)
{
    assert(isSymmetric(B));
    assert(Q.rows() == n);
    assert(Q.cols() == 2);

    for (int i = 0; i < n; ++i)
    {
        result(i, i) = BProduct(B, Q.row(i), Q.row(i));
        for (int j = 0; j < i; ++j)
            result(i, j) = result(j, i) = BProduct(B, Q.row(i), Q.row(j));
    }
}


}
#endif /* ELASTICSTRANDUTILS_HH_ */
