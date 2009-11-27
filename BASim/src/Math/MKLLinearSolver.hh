/**
 * \file MLLLinearSolver.hh
 *
 * \author acoull@wetafx.co.nz
 * \date 11/11/2009
 */
#ifndef MKLLINEARSOLVER_HH
#define MKLLINEARSOLVER_HH

#ifdef USING_INTEL_COMPILER

namespace BASim {

#include <mkl_solver.h>
#include <mkl_lapack.h>

/** Linear solver that uses MKL. */
class MKLLinearSolver : public LinearSolverBase
{
public:

  MKLLinearSolver(MatrixBase& A)
    : LinearSolverBase(A)
  {
    assert(m_A.rows() == m_A.cols());
  }

  ~MKLLinearSolver()
  {
  }

  int solve(VecXd& x, const VecXd& b)
  {

    vector<Scalar> c;
    size_t numElements = b.size();
    c.resize(numElements);
    for (size_t e=0; e<numElements; e++)
        c[e] = b[e];

    MKL_INT n    = c.size();
    MKL_INT nrhs = 1;
    MKL_INT lda  = n;
    MKL_INT ldb  = n;
    MKL_INT info = 0;
    vector<MKL_INT> ipiv(n);
    
    MKLMatrix& M = dynamic_cast<MKLMatrix&>(m_A);
    M.transpose();
    vector<Scalar>& a = M.getData();

    dgesv(&n, &nrhs, &a[0], &lda, &ipiv[0], &c[0], &ldb, &info);

    for (size_t e=0; e<numElements; e++)
        x[e] = c[e];
    
    if (info != 0)
        std::cerr << "ERROR: Failed to find solution to linear system!" << std::endl;
    return 0;
  }

  int checkConvergedReason()
  {
    /*KSPConvergedReason reason;
    KSPGetConvergedReason(m_kspSolver, &reason);

    if (reason == KSP_DIVERGED_NULL)
      std::cout << "Diverged null" << std::endl;
    else if (reason == KSP_DIVERGED_ITS)
      std::cout << "Diverged its" << std::endl;
    else if (reason == KSP_DIVERGED_NAN)
      std::cout << "Diverged nan" << std::endl;
    else if (reason == KSP_DIVERGED_BREAKDOWN_BICG)
      std::cout << "Diverged breakdown bicg" << std::endl;
    else if (reason == KSP_DIVERGED_DTOL)
      std::cout << "Diverged dtol" << std::endl;
    else if (reason == KSP_DIVERGED_INDEFINITE_PC)
      std::cout << "Diverged indefinite pc" << std::endl;
    else if (reason == KSP_DIVERGED_INDEFINITE_MAT)
      std::cout << "Diverged indefinite mat" << std::endl;
    else if (reason == KSP_DIVERGED_NONSYMMETRIC)
      std::cout << "Diverged nonsymmetric" << std::endl;
    else if (reason == KSP_DIVERGED_BREAKDOWN)
      std::cout << "Diverged breakdown" << std::endl;

    if (reason < 0) return -1;*/
    return 0;
  }

protected:

};

} // namespace BASim

#endif // USING_INTEL_COMPILER

#endif // MKLLINEARSOLVER_HH
