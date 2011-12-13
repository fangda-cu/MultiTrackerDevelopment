/**
 * \file ConjugateGradient.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 11/16/2009
 */

#ifndef CONJUGATEGRADIENT_HH
#define CONJUGATEGRADIENT_HH

#include "DiagonalPreconditioner.hh"
#include "LinearSolverBase.hh"

namespace BASim {

/** Implements the conjugate gradient method for solving a linear
    system using a diagonal preconditioner. */
class ConjugateGradient : public LinearSolverBase
{
public:

  explicit ConjugateGradient(MatrixBase& A)
    : LinearSolverBase(A)
    , m_maxIterations(std::max(A.rows(), A.cols()))
    , m_currentIterations(0)
    , m_rnorm(1.0e-10)
    , m_preconditioner(NULL)
  {}

  ~ConjugateGradient()
  {
    if (m_preconditioner != NULL) delete m_preconditioner;
  }

  int getMaxIterations() const { return m_maxIterations; }
  void setMaxIterations(int m) { m_maxIterations = m; }
  int getCurrentIterations() const { return m_currentIterations; }
  Scalar getRNorm() const { return m_rnorm; }
  void setRNorm(Scalar r) { m_rnorm = r; }

  /**
   * Solves the equation \f$Ax=b\f$ for \f$x\f$, given \f$A\f$ and
   * \f$b\f$.
   */
  int solve(VecXd& x, const VecXd& b);
 

protected:

  int m_maxIterations;
  int m_currentIterations;

  Scalar m_rnorm;

  Preconditioner* m_preconditioner;

};

} // namespace BASim

#endif // CONJUGATEGRADIENT_HH
