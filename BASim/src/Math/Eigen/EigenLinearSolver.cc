#include "EigenLinearSolver.hh"
#include "BASim/src/Core/Util.hh"
#include "BASim/src/Math/EigenSparseMatrix.hh"
//#include <Eigen/SparseExtra>

namespace BASim {

   EigenLinearSolver::EigenLinearSolver(MatrixBase& A)
      : LinearSolverBase(A)

   {
      assert(m_A.rows() == m_A.cols());
   }


   int EigenLinearSolver::solve(VecXd& x, const VecXd& b)
   {
      const Eigen::SparseMatrix<Scalar>& matrix = smart_cast<EigenSparseMatrix&>(m_A).getEigenMatrix();

      Eigen::SimplicialLDLT< Eigen::SparseMatrix<Scalar> > llt_of_A;
      llt_of_A.compute(matrix);

      if(llt_of_A.info() != Eigen::Success) {
         // decomposition failed
         std::cout << "Eigen's LDLT Decomposition failed\n";
         return -1;
      }
      x = llt_of_A.solve(b);
      if(!llt_of_A.info() != Eigen::Success) {
         // solving failed
         std::cout << "Eigen's LDLT solve failed\n";
         return -1;
      }

      std::cout << "Solve successful\n";
      return 0;
   }

}
