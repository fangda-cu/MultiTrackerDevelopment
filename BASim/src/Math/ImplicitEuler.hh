/**
 * \file ImplicitEuler.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/04/2009
 */

#ifndef IMPLICITEULER_HH
#define IMPLICITEULER_HH

namespace BASim {

/** This class implements the implicit Euler time-stepping
    method. This assumes a pair of equations of the form
    \f{eqnarray*}\frac{dx}{dt} &=& v \\ M\frac{dv}{dt} &=& f(t,x) \
    . \f} The mass matrix \f$M\f$ is assumed to be diagonal.
*/
template <class ODE>
class ImplicitEuler : public DiffEqSolver
{
public:

  ImplicitEuler(ODE& ode, ObjectControllerBase::SolverLibrary solver)
    : m_diffEq(ode)
  {
    if ( solver == ObjectControllerBase::PETSC_SOLVER )
    {
      cerr << "Creating new Petsc solver\n";
      m_A = new PetscMatrix(m_diffEq.ndof(), m_diffEq.ndof(), 11);
      m_solver = new PetscLinearSolver(*m_A);
    }
    else if ( solver == ObjectControllerBase::MKL_SOLVER )
    {  
      cerr << "Creating new MKL solver\n";
      m_A = new MKLMatrix(m_diffEq.ndof(), m_diffEq.ndof(), 11);
      m_solver = new MKLLinearSolver(*m_A);
    }
  }

  ~ImplicitEuler()
  {
    delete m_A;
    delete m_solver;
  }

  void execute()
  {
    resize();
    setZero();
    const IntArray& fixed = m_diffEq.getFixedDofs();

    // copy start of step positions and velocities
    for (int i = 0; i < x0.size(); ++i) {
      v0(i) = m_diffEq.getV(i);
      x0(i) = m_diffEq.getX(i);
    }

    m_diffEq.evaluatePDot(m_rhs);
    m_rhs *= m_dt;
    int iter;
    double residual = 0;

    for (iter = 0; iter < m_maxIterations; ++iter) {
      m_diffEq.evaluatePDotDX(*m_A);
      m_A->finalize();
      m_A->scale(m_dt);
      
      if (iter == 0) {
        m_A->finalize();
        m_A->multiply(m_rhs, m_dt, v0);
      
      } else {
        for (int i = 0; i < m_diffEq.ndof(); ++i) {
          m_rhs(i) -= m_diffEq.getMass(i) * m_deltaV(i);
        }
      }

      m_diffEq.evaluatePDotDV(*m_A);
      m_A->finalize();
      m_A->scale(-m_dt);
      
      
      for (int i = 0; i < m_diffEq.ndof(); ++i) {
        m_A->add(i, i, m_diffEq.getMass(i));
      }
      m_A->finalize();
      
      for (size_t i = 0; i < fixed.size(); ++i) m_rhs(fixed[i]) = 0;
      m_A->zeroRows(fixed);
      
      m_solver->solve(m_increment, m_rhs);
      m_deltaV += m_increment;

/*
      if(!m_solver->solve(m_increment, m_rhs)) //solver returns 0 if it succeeds
	  m_deltaV += m_increment;
      else // using unconverged results but clamp them first
	  for(int i=0; i<m_deltaV.size(); i++)
	  {
	      Scalar inc, clamp;
	      clamp = 1.0;
	      if(m_increment[i] > clamp)
		  inc = clamp;
	      else if(m_increment[i] < -clamp)
		  inc = -clamp;
	      else
		  inc = m_increment[i];
	      m_deltaV[i] += inc;
	  }
*/	  
    

      for (int i = 0; i < m_deltaV.size(); ++i) {
        m_diffEq.setV(i, v0(i) + m_deltaV(i));
        m_diffEq.setX(i, x0(i) + m_dt * m_diffEq.getV(i));
      }

      m_diffEq.flush();

//      if (iter == m_maxIterations - 1) break;

      // check for convergence
      m_rhs.setZero();
      m_diffEq.evaluatePDot(m_rhs);
      m_rhs *= m_dt;
      for (int i = 0; i < m_increment.size(); ++i)
        m_increment(i) = m_diffEq.getMass(i) * m_deltaV(i);
      for (size_t i = 0; i < fixed.size(); ++i) m_rhs(fixed[i]) = 0;
      residual = (m_increment - m_rhs).norm();
      if (residual< 1.0e-10) {
        break;
      }
      m_increment.setZero();
      m_A->setZero();      
    }
    if(iter == m_maxIterations)
	std::cout<<"Newton failed to converge after "<<m_maxIterations<<" iterations with residual "<< residual<<std::endl;

  }

  void resize()
  {
    x0.resize(m_diffEq.ndof());
    v0.resize(m_diffEq.ndof());
    m_rhs.resize(m_diffEq.ndof());
    m_deltaV.resize(m_diffEq.ndof());
    m_increment.resize(m_diffEq.ndof());
    if (m_A->rows() != m_diffEq.ndof()) {
      delete m_A;
      //m_A = new PetscMatrix(m_diffEq.ndof(), m_diffEq.ndof(), 11);
      m_A = new MKLMatrix(m_diffEq.ndof(), m_diffEq.ndof(), 11);
      delete m_solver;
      //m_solver = new PetscLinearSolver(*m_A);
      m_solver = new MKLLinearSolver(*m_A);
    }
  }

  void setZero()
  {
    x0.setZero();
    v0.setZero();
    m_rhs.setZero();
    m_deltaV.setZero();
    m_increment.setZero();
    m_A->setZero();    
  }

protected:

  ODE& m_diffEq;

  VecXd x0;
  VecXd v0;
  VecXd m_rhs;
  VecXd m_deltaV;
  VecXd m_increment;

  MatrixBase* m_A;
  LinearSolverBase* m_solver;
};

} // namespace BASim

#endif // IMPLICITEULER_HH
