/**
 * \file SolverUtils.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 11/17/2009
 */

#include "BASim/Math"
#include "SolverUtils.hh"

namespace BASim {

SolverUtils* SolverUtils::m_instance = NULL;

SolverUtils* SolverUtils::instance()
{
  if (m_instance == NULL) m_instance = new SolverUtils();
  return m_instance;
}

SolverUtils::SolverType SolverUtils::getSolverType() const
{
  return solverType;
}

void SolverUtils::setSolverType(SolverType t) { solverType = t; }

SolverUtils::MatrixType SolverUtils::getMatrixType() const
{
  return matrixType;
}

void SolverUtils::setMatrixType(MatrixType t) {matrixType = t; }

MatrixBase*
SolverUtils::createSparseMatrix(int rows, int cols, int nnzPerRow) const
{
#ifdef HAVE_PETSC
  return new PetscMatrix(rows, cols, nnzPerRow);
#endif // HAVE_PETSC

  std::cerr << "createSparseMatrix failure" << std::endl;
  exit(-1);
  return NULL;
}

MatrixBase*
SolverUtils::createBandMatrix(int rows, int cols, int kl, int ku) const
{
  if (matrixType == BAND_MATRIX) {
    return new BandMatrix(rows, cols, kl, ku);
  }

#ifdef HAVE_PETSC
  if (matrixType == PETSC_MATRIX) {
    return new PetscMatrix(rows, cols, kl + ku + 1);
  }
#endif // HAVE_PETSC

  assert(matrixType == AUTO_MATRIX);
  return new BandMatrix(rows, cols, kl, ku);
}

LinearSolverBase* SolverUtils::createLinearSolver(MatrixBase* A) const
{
  if (solverType == CONJUGATE_GRADIENT)
    return new ConjugateGradient(*A);

#ifdef HAVE_PETSC
  if (solverType == PETSC_SOLVER)
    return new PetscLinearSolver(*A);
#endif // HAVE_PETSC

#ifdef HAVE_MKL
  if ((solverType == MKL_LINEAR_SOLVER || solverType == AUTO_SOLVER)
      && dynamic_cast<BandMatrix*>(A) != NULL) {
    return new MKLLinearSolver(dynamic_cast<BandMatrix&>(*A));
  }
#endif // HAVE_MKL

  std::cout << "cg" << std::endl;
  return new ConjugateGradient(*A);
}

} // namespace BASim
