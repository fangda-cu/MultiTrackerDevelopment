#include "EigenLinearSolver.hh"
#include "BASim/src/Core/Util.hh"
#include "BASim/src/Math/EigenSparseMatrix.hh"
#include <Eigen/SparseExtra>

namespace BASim {

   EigenLinearSolver::EigenLinearSolver(MatrixBase& A)
      : LinearSolverBase(A)

   {
      assert(m_A.rows() == m_A.cols());
   }


   int EigenLinearSolver::solve(VecXd& x, const VecXd& b)
   {
      const Eigen::SparseMatrix<Scalar>& matrix = smart_cast<EigenSparseMatrix&>(m_A).getEigenMatrix();

      Eigen::SparseLDLT< Eigen::SparseMatrix<Scalar> > llt_of_A(matrix);
      
      if(!llt_of_A.succeeded()) {
         // decomposition failed
         std::cout << "Eigen's LDLT Decomposition failed\n";
         return -1;
      }
      x = llt_of_A.solve(b);
      //if(!llt_of_A.solve(b,&x)) {
      //   // solving failed
      //   return -1;
      //}

      std::cout << "Solve successful\n";
      return 0;
   }

}