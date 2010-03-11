/**
 * \file MKLLinearSolver.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 11/29/2009
 */

#include "MKLLinearSolver.hh"
#include "../BandMatrix.hh"
#ifdef HAVE_MKL
#include <mkl_lapack.h>
#else // HAVE_MKL
extern "C" void dgbsv_( int* n, int* kl, int* ku, int* nrhs, double* ab, int* ldab, int* ipiv, double* b, int* ldb, int* info );
#define dgbsv dgbsv_
#endif // HAVE_MKL

namespace BASim {

MKLLinearSolver::MKLLinearSolver(BandMatrix& A)
  : LinearSolverBase(A)
{
  assert(A.rows() == A.cols());
  int n = A.rows();
  int kl = smart_cast<BandMatrix&>(A).kl();
  int ku = smart_cast<BandMatrix&>(A).ku();
  ipiv = new int[n];

  // ab holds the entries of the matrix. Space must be made for an
  // additional kl super-diagonals for LU factorization
  ab = new double[(2 * kl + ku + 1) * n];
}

MKLLinearSolver::~MKLLinearSolver()
{
  delete [] ipiv;
  delete [] ab;
}

inline void convert(double* ab, BandMatrix& A, int kl, int ku, int n)
{
  int NUMROWS = 2 * kl + ku + 1;
  for (int j = 0; j < n; ++j) {
    for (int i = std::max(0, j - ku); i < std::min(n, j + kl + 1); ++i) {
      int row = kl + ku + i - j;
      int col = j;
      int offset = row + col * NUMROWS;
      ab[offset] = A(i, j);
    }
  }
}

int MKLLinearSolver::solve(VecXd& x, const VecXd& b)
{
  int n, kl, ku, nrhs, ldab, ldb, info;

  BandMatrix& A = smart_cast<BandMatrix&>(m_A);
  n = A.rows();
  kl = A.kl();
  ku = A.ku();
  nrhs = 1;
  ldab = 2 * kl + ku + 1;
  ldb = n;
  convert(ab, A, kl, ku, n);
  x = b;

  dgbsv(&n, &kl, &ku, &nrhs, ab, &ldab, ipiv, x.data(), &ldb, &info);

  // check return value for errors
  if (info < 0) {
    std::cerr << "Error in parameter " << -info << " in call to dgbsv"
              << std::endl;
    return -1;
  } else if (info > 0) {
    std::cerr << "Factor U is singular (" << info << std::endl;
    return -1;
  }

  return 0;
}

} // namespace BASim
