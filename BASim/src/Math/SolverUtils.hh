/**
 * \file SolverUtils.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 11/17/2009
 */

#ifndef SOLVERUTILS_HH
#define SOLVERUTILS_HH

#include "BandMatrix.hh"
#include "ConjugateGradient.hh"

#ifdef HAVE_PETSC
#include "BASim/src/Math/Petsc/PetscMatrix.hh"
#include "BASim/src/Math/Petsc/PetscLinearSolver.hh"
#endif // HAVE_PETSC

#ifdef HAVE_MKL
#include "MKL/MKLLinearSolver.hh"
#endif // HAVE_MKL

namespace BASim {

class SolverUtils
{
public:

  static SolverUtils* instance();

  // define solver types
  enum SolverType {
#ifdef HAVE_PETSC
    PETSC_SOLVER,
#endif // HAVE_PETSC
#ifdef HAVE_MKL
    MKL_LINEAR_SOLVER,
#endif // HAVE_MKL
    CONJUGATE_GRADIENT,
    AUTO_SOLVER
  };

  // define matrix types
  enum MatrixType {
#ifdef HAVE_PETSC
    PETSC_MATRIX,
#endif // HAVE_PETSC
    BAND_MATRIX,
    AUTO_MATRIX
  };

  SolverType getSolverType() const;
  void setSolverType(SolverType t);

  MatrixType getMatrixType() const;
  void setMatrixType(MatrixType t);

  MatrixBase* createSparseMatrix(int rows, int cols, int nnzPerRow = 1) const;
  MatrixBase* createBandMatrix(int rows, int cols, int kl, int ku) const;

  LinearSolverBase* createLinearSolver(MatrixBase* A) const;

private:

  SolverUtils()
    : solverType(AUTO_SOLVER)
    , matrixType(AUTO_MATRIX)
  {}

  SolverUtils(const SolverUtils&) {}
  SolverUtils& operator=(const SolverUtils&) { return *this; }

  static SolverUtils* m_instance;

  SolverType solverType;
  MatrixType matrixType;
};

} // namespace BASim

#endif // SOLVERUTILS_HH
