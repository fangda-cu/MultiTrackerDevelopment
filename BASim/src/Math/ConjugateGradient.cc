#include "ConjugateGradient.hh"

namespace BASim {

 /**
   * Solves the equation \f$Ax=b\f$ for \f$x\f$, given \f$A\f$ and
   * \f$b\f$.
   */
  int ConjugateGradient::solve(VecXd& x, const VecXd& b)
  {
    std::cout << "Performing conjugate gradient solve!\n";
    m_currentIterations = 0;
    m_preconditioner = new DiagonalPreconditioner(m_A);

    VecXd r = b;
    m_A.multiply(r, -1, x);

    VecXd z(m_A.rows());
    m_preconditioner->apply(z, r);

    Scalar rk_dot_zk = r.dot(z);
    VecXd p = z;
    VecXd w(m_A.rows());

    while (r.norm() >= m_rnorm && m_currentIterations < m_maxIterations) {

      w.setZero();
      m_A.multiply(w, 1, p);

      Scalar alpha = rk_dot_zk / p.dot(w);
      x += alpha * p;
      r -= alpha * w;

      m_preconditioner->apply(z, r);

      Scalar rk1_dot_zk1 = r.dot(z);
      Scalar beta = rk1_dot_zk1 / rk_dot_zk;
      p = z + beta * p;

      ++m_currentIterations;
      rk_dot_zk = rk1_dot_zk1;
    }

    delete m_preconditioner;
    m_preconditioner = NULL;

    // Check the inf norm of the residual
    #ifdef DEBUG
      VecXd residual(x.size());
      residual.setZero();
      m_A.multiply(residual,1.0,x);
      residual -= b;
      double infnorm = fabs(residual.maxCoeff());
      if( infnorm > 1.0e-6 )
      {
        std::cout << "\033[31;1mWARNING IN ConjugateGradient:\033[m Large residual detected. ||residual||_{inf} = " << infnorm << std::endl;
        return -1;
      }
    #endif

    return 0;
  }





}