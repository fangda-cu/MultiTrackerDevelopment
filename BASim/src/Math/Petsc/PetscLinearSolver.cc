  #include "PetscLinearSolver.hh"
#include "BASim/src/Core/Util.hh"

  namespace BASim {
  
PetscLinearSolver::PetscLinearSolver(MatrixBase& A)
    : LinearSolverBase(A)
    , m_x(NULL)
    , m_b(NULL)
  {
    assert(m_A.rows() == m_A.cols());

    // create solver context and associate system matrix with it
    KSPCreate(PETSC_COMM_SELF, &m_kspSolver);
    Mat& pA = smart_cast<PetscMatrix&>(m_A).getPetscMatrix();
    KSPSetOperators(m_kspSolver, pA, pA, SAME_NONZERO_PATTERN);
    KSPSetFromOptions(m_kspSolver);
    /*
    const KSPType kspType;
    KSPGetType(m_kspSolver, &kspType);
    std::cout << "KSP type is " << kspType << std::endl;

    PC pc;
    KSPGetPC(m_kspSolver, &pc);
    const PCType pcType;
    PCGetType(pc, &pcType);
    std::cout << "PC type is " << pcType << std::endl;
    */
    // create PETSc solution and right-hand-side vectors
    VecCreateSeq(PETSC_COMM_SELF, m_A.rows(), &m_x);
    VecCreateSeq(PETSC_COMM_SELF, m_A.rows(), &m_b);
  }
  

  int PetscLinearSolver::solve(VecXd& x, const VecXd& b)
  {
     std::cout << "Copying to petscvectors\n";
    PetscUtils::copyToPetscVector(m_b, b);
    PetscUtils::copyToPetscVector(m_x, x);
    std::cout << "Performing solve\n";
    KSPSolve(m_kspSolver, m_b, m_x);
    /*PetscInt its;
    KSPGetIterationNumber(m_kspSolver, &its);
    std::cout << "converged in " << its << " iterations" << std::endl;
    */
    std::cout << "Done solve, copying vectors back\n";
    PetscUtils::copyFromPetscVector(x, m_x);

    std::cout << "Checking convergence\n";
    if (checkConvergedReason() == -1) return -1;
    std::cout << "Returning\n";
    return 0;
  }

  }