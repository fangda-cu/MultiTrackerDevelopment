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
    const int n = A.rows();
    const int nrhs = 1;
    const int ldab = 2 * kl + ku + 1;
    const int ldb = n;

    // Prepare data for LAPACK-style solver
    int info;
    m_ipiv.resize( n );
    m_ab.resize( ldab * n ); // Space must be made for an  additional kl super-diagonals for LU factorization

    for ( int j = 0; j < n; ++j )
        for ( int i = std::max( 0, j - ku ); i < std::min( n, j + kl + 1 ); ++i )
        {
            int row = kl + ku + i - j;
            int col = j;
            m_ab[row + col * ldab] = A( i, j );
        }
    x = b;
    static const int sckl = kl; // Just because dgbsv wants a pointer
    static const int scku = ku;

    dgbsv_( &n, &sckl, &scku, &nrhs, &m_ab[0], &ldab, &m_ipiv[0], x.data(), &ldb, &info );

    return info;
}

template class BandMatrixLinearSolver<10, 10> ;

}
