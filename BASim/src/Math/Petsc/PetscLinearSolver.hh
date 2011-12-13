/**
 * \file PetscLinearSolver.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/08/2009
 */

#ifndef PETSCLINEARSOLVER_HH
#define PETSCLINEARSOLVER_HH

#include <petscksp.h>
#include "BASim/src/Math/LinearSolverBase.hh"
#include "PetscMatrix.hh"

namespace BASim {

/** Linear solver that uses PETSc. */
class PetscLinearSolver : public LinearSolverBase
{
public:

  PetscLinearSolver(MatrixBase& A);

  ~PetscLinearSolver()
  {
    KSPDestroy(m_kspSolver);
    VecDestroy(m_x);
    VecDestroy(m_b);
  }

  
  int solve(VecXd& x, const VecXd& b);


  int checkConvergedReason()
  {
    KSPConvergedReason reason;
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

    if (reason < 0) return -1;
    return 0;
  }

protected:

  KSP m_kspSolver;
  Vec m_x;
  Vec m_b;
};

} // namespace BASim

#endif // PETSCLINEARSOLVER_HH
