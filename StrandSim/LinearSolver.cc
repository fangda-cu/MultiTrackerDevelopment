/*
 * BandMatrixLinearSolver.cc
 *
 *  Created on: 18/07/2011
 *      Author: jaubry
 */

#include "LinearSolver.hh"
#include "BandMatrix.hh"
#include <mkl.h>

namespace strandsim
{

template<int kl, int ku>
BandMatrixLinearSolver<kl, ku>::BandMatrixLinearSolver()
{
}

template<int kl, int ku>
BandMatrixLinearSolver<kl, ku>::~BandMatrixLinearSolver()
{
}

template<int kl, int ku>
int BandMatrixLinearSolver<kl, ku>::solve( VecXd& x, const BandMatrix<double, kl, ku>& A,
        const VecXd& b )
{
    assert( A.rows() == A.cols() );
    assert( b.rows() == A.cols() );
    assert( x.rows() == A.cols() );

    static const int s_kl = kl; // Just because dgbsv wants a pointer
    static const int s_ku = ku;
    static const int ldab = 2 * kl + ku + 1;
    static const int nrhs = 1;
    const int n = A.cols();
    const int ldb = n;

    // Prepare data for LAPACK-style solver
    int info;
    m_ipiv.resize( n );
    m_ab.resize( ldab * n ); // Space must be made for an  additional kl super-diagonals for LU factorization

    // Copy the data from A to m_ab
    const std::vector<Scalar>& data = A.getData();
    int mj = kl;
    for ( int j = 0; j < n; ++j )
    {
        for ( int i = 0; i < kl + ku + 1; ++i )
            m_ab[i + mj] = data[i * n + j];
        mj += ldab;
    }

    // Copy the right-hand side
    x = b;

    // Solve it!
    dgbsv_( &n, &s_kl, &s_ku, &nrhs, &m_ab[0], &ldab, &m_ipiv[0], x.data(), &ldb, &info );

    return info;
}

template class BandMatrixLinearSolver<10, 10> ;

}
