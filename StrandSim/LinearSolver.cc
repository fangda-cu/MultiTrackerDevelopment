/*
 * LinearSolver.cc
 *
 *  Created on: 18/07/2011
 *      Author: jaubry
 */

#include "LinearSolver.hh"
#include "BandMatrix.hh"
#include <mkl.h>

namespace strandsim
{

template<typename MatrixT>
LinearSolver<MatrixT>::LinearSolver()
{
    // TODO Auto-generated constructor stub

}

template<typename MatrixT>
LinearSolver<MatrixT>::~LinearSolver()
{
    // TODO Auto-generated destructor stub
}

template<int kl, int ku>
void convert(double* ab, const BandMatrix<Scalar, kl, ku>& A, int n)
{
    int NUMROWS = 2 * kl + ku + 1;
    for (int j = 0; j < n; ++j)
    {
        for (int i = std::max(0, j - ku); i < std::min(n, j + kl + 1); ++i)
        {
            int row = kl + ku + i - j;
            int col = j;
            int offset = row + col * NUMROWS;
            ab[offset] = A(i, j);
        }
    }
}

template<>
int LinearSolver<BandMatrix<Scalar, 10, 10> >::solve(VecXd& x, const BandMatrix<Scalar, 10, 10>& A, const VecXd& b) const
{
    static  int kl = 10; // replace with templates
    static  int ku = 10; // replace with templates

    int n = A.rows();
    int nrhs = 1;
    int ldab = 2 * kl + ku + 1;
    int ldb = n;
    int info;
    int* ipiv = new int[n];

    // ab holds the entries of the matrix. Space must be made for an
    // additional kl super-diagonals for LU factorization
    Scalar* ab = new Scalar[(2 * kl + ku + 1) * n];

    convert<10, 10> (ab, A, n);
    x = b;

    dgbsv_(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, x.data(), &ldb, &info);

    return info;
}

template class LinearSolver<BandMatrix<Scalar, 10, 10> > ;

}
