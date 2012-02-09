/**
 * \file EigenCGSolver.hh
 *
 * \author batty@cs.columbia.edu
 * \date Feb 9, 2012
 */

#ifndef EIGENCGSOLVER_HH
#define EIGENCGSOLVER_HH

#include "BASim/src/Math/LinearSolverBase.hh"

namespace BASim {

/** Linear solver that uses PETSc. */
class EigenCGSolver : public LinearSolverBase
{
public:

  EigenCGSolver (MatrixBase& A);

  ~EigenCGSolver ()
  {
  
  }

  
  int solve(VecXd& x, const VecXd& b);

};

} // namespace BASim

#endif // EIGENCGSOLVER_HH
