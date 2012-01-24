/**
 * \file SolverUtils.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 11/17/2009
 */

#include "SolverUtils.hh"
#ifdef HAVE_LAPACK
#ifdef WETA
#include "MKL/MKLLinearSolver.hh"
#else
#include "BASim/src/Math/MKL/MKLLinearSolver.hh"
#endif
#endif // HAVE_LAPACK

#include "BASim/src/Math/SimpleSparseMatrix.hh"
#include "BASim/src/Math/EigenSparseMatrix.hh"
#include "BASim/src/Math/Eigen/EigenLinearSolver.hh"
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

std::string SolverUtils::getSolverName() const
{
  if (solverType == CONJUGATE_GRADIENT)
    return "CONJUGATE_GRADIENT";
  
  if(solverType == EIGEN_LDLT) 
     return "EIGEN_LDLT";
#ifdef HAVE_PARDISO
  if (solverType == PARDISO_SOLVER)
    return "PARDISO_SOLVER";
#endif

#ifdef HAVE_PETSC
  if (solverType == PETSC_SOLVER)
    return "PETSC_SOLVER";
#endif // HAVE_PETSC

#ifdef HAVE_LAPACK
  if ((solverType == MKL_LINEAR_SOLVER || solverType == AUTO_SOLVER) )
    //&& dynamic_cast<BandMatrix*>(A) != NULL)
    {
#ifdef HAVE_MKL
      return "MKL_LINEAR_SOLVER";
#else
      return "LAPACK_LINEAR_SOLVER";
#endif // HAVE_MKL
    }
#endif // HAVE_LAPACK

  return "CONJUGATE_GRADIENT";
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
#ifdef HAVE_PARDISO
  if (matrixType == PARDISO_MATRIX)
    return new PardisoMatrix(rows,cols);
#endif // HAVE_PARDISO  
  
#ifdef HAVE_PETSC
  if (matrixType == PETSC_MATRIX)
    return new PetscMatrix(rows, cols, nnzPerRow);
#endif // HAVE_PETSC

  return new EigenSparseMatrix(rows, cols, nnzPerRow);
  //return new SimpleSparseMatrix(rows, cols, nnzPerRow);

  std::cerr << "\033[31;1mWARNING IN SOLVERUTILS:\033[m Failed to create sparse matrix. " << std::endl;
  exit(-1);
  return NULL;
}

MatrixBase*
SolverUtils::createBandMatrix(int rows, int cols, int kl, int ku) const
{
  if (matrixType == BAND_MATRIX) {
    return new BandMatrix(rows, cols, kl, ku);
  }

#ifdef HAVE_PARDISO
  if (matrixType == PARDISO_MATRIX)
    return new PardisoMatrix(rows,cols);
#endif // HAVE_PARDISO

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

  if(solverType == EIGEN_LDLT)
     return new EigenLinearSolver(*A);

#ifdef HAVE_PARDISO
  if (solverType == PARDISO_SOLVER)
    return new PardisoLinearSolver(dynamic_cast<PardisoMatrix&>(*A));
#endif // HAVE_PARDISO
  
#ifdef HAVE_PETSC
  if (solverType == PETSC_SOLVER) {
    return new PetscLinearSolver(*A);
  }
#endif // HAVE_PETSC

#ifdef HAVE_LAPACK
  if ((solverType == MKL_LINEAR_SOLVER || solverType == AUTO_SOLVER)
      && dynamic_cast<BandMatrix*>(A) != NULL) {
    return new MKLLinearSolver(dynamic_cast<BandMatrix&>(*A));
  }
#endif // HAVE_LAPACK



  return new ConjugateGradient(*A);
}

} // namespace BASim
