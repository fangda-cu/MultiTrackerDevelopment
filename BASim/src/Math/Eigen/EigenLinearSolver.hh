/**
 * \file EigenLinearSolver.hh
 *
 * \author batty@cs.columbia.edu
 * \date Jan 24, 2012
 */

#ifndef EIGENLINEARSOLVER_HH
#define EIGENLINEARSOLVER_HH

#include "BASim/src/Math/LinearSolverBase.hh"

namespace BASim {

/** Linear solver that uses PETSc. */
class EigenLinearSolver : public LinearSolverBase
{
public:

  EigenLinearSolver (MatrixBase& A);

  ~EigenLinearSolver ()
  {
  
  }

  
  int solve(VecXd& x, const VecXd& b);

};

} // namespace BASim

#endif // EIGENLINEARSOLVER_HH
