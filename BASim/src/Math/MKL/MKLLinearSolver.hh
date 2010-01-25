/**
 * \file MKLLinearSolver.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 11/29/2009
 */

#ifndef MKLLINEARSOLVER_HH
#define MKLLINEARSOLVER_HH

#include "../LinearSolverBase.hh"
#include "../BandMatrix.hh"
#include "../../Core/Util.hh"

namespace BASim {

/** LU-based linear solver for band matrices. */
class MKLLinearSolver : public LinearSolverBase
{
public:

  MKLLinearSolver(BandMatrix& A);
  virtual ~MKLLinearSolver();

  virtual int solve(VecXd& x, const VecXd& b);

private:

  void* ipiv;
  double* ab;
};

} // namespace BASim

#endif // MKLLINEARSOLVER_HH
