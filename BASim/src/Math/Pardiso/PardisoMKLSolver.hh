/**
 * \file PardisoLinearSolver.hh
 *
 * \author batty@cs.columbia.edu
 * \date Jan 27, 2012
 */
 

#ifndef PARDISOMKLSOLVER_HH
#define PARDISOMKLSOLVER_HH

#include "BASim/src/Math/LinearSolverBase.hh"
#include "BASim/src/Math/EigenSparseMatrix.hh"

namespace BASim
{

/** 
 * Solves a sparse linear system using MKL's Pardiso API.
 */
class PardisoMKLSolver : public LinearSolverBase
{
public:

  explicit PardisoMKLSolver( EigenSparseMatrix& A );
  virtual ~PardisoMKLSolver();

  void parsePardisoError( int error ) const;

  virtual int solve( VecXd& x, const VecXd& b );

private:
  EigenSparseMatrix& prdsomat;
};

} // namespace BASim

#endif // PARDISOMKLSOLVER_HH
