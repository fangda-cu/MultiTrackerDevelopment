#include "EigenCGSolver.hh"
#include "BASim/src/Core/Util.hh"
#include "BASim/src/Math/EigenSparseMatrix.hh"


namespace BASim {

  EigenCGSolver::EigenCGSolver(MatrixBase& A)
      : LinearSolverBase(A)

   {
      assert(m_A.rows() == m_A.cols());
   }

   
   int EigenCGSolver::solve(VecXd& x, const VecXd& b)
   {
      const Eigen::SparseMatrix<Scalar, Eigen::RowMajor>& matrix = smart_cast<EigenSparseMatrix&>(m_A).getEigenMatrix();

      Eigen::ConjugateGradient< Eigen::SparseMatrix<Scalar, Eigen::RowMajor>, Eigen::Upper > solver;
      //Eigen::BiCGSTAB< Eigen::SparseMatrix<Scalar, Eigen::RowMajor> > solver;
      //defaults to using the lower triangle, with a diagonal preconditioner
      
      solver.setTolerance(1e-7); //is this absolute? relative?
      solver.setMaxIterations(10000);

      solver.compute(matrix);

      if(solver.info() != Eigen::Success) {
        // solving failed
        std::cout << "Eigen's CG failed to compute\n";
        return -1;
      }
      
      x = solver.solve(b);
      
      if(solver.info() != Eigen::Success) {
         // solving failed
         std::cout << "Eigen's CG failed to solve\n";
         return -1;
      }

      std::cout << "Solve successful\n";
      return 0;
   }

}