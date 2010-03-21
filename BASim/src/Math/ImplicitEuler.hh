/**
 * \file ImplicitEuler.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/04/2009
 */

#ifndef IMPLICITEULER_HH
#define IMPLICITEULER_HH

#include "TimeSteppingBase.hh"
#include "LinearSolverBase.hh"
#include "SolverUtils.hh"
#include "../Core/Timer.hh"

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

  ImplicitEuler(ODE& ode) : m_diffEq(ode)
  {
    m_A = m_diffEq.createMatrix();
    m_solver = SolverUtils::instance()->createLinearSolver(m_A);
  }

  ~ImplicitEuler()
  {
    delete m_A;
    delete m_solver;
  }

  void execute()
  {
    m_diffEq.startStep();

    resize();
    setZero();
    IntArray fixed;
    std::vector<Scalar> desired;
    m_diffEq.getScriptedDofs(fixed, desired);

    // copy start of step positions and velocities
    for (int i = 0; i < x0.size(); ++i) {
      v0(i) = m_diffEq.getV(i);
      x0(i) = m_diffEq.getX(i);
    }

    m_diffEq.evaluatePDot(m_rhs);
    m_rhs *= m_dt;

    for (int iter = 0; iter < m_maxIterations; ++iter) {
      START_TIMER("setup");
      m_diffEq.startIteration();
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

      for (size_t i = 0; i < fixed.size(); ++i) {
        int idx = fixed[i];
        m_rhs(idx) = (desired[i] - x0[idx]) / m_dt - (v0[idx] + m_deltaV[idx]);
      }
      m_A->zeroRows(fixed);
      STOP_TIMER("setup");

      START_TIMER("solver");
      m_solver->solve(m_increment, m_rhs);
      STOP_TIMER("solver");

      START_TIMER("setup");
      m_deltaV += m_increment;

      for (int i = 0; i < m_deltaV.size(); ++i) {
        m_diffEq.setV(i, v0(i) + m_deltaV(i));
        m_diffEq.setX(i, x0(i) + m_dt * m_diffEq.getV(i));
      }

      m_diffEq.endIteration();

      if (iter == m_maxIterations - 1) break;

      // check for convergence
      m_rhs.setZero();
      m_diffEq.evaluatePDot(m_rhs);
      m_rhs *= m_dt;
      for (int i = 0; i < m_increment.size(); ++i)
        m_increment(i) = m_diffEq.getMass(i) * m_deltaV(i);
      for (size_t i = 0; i < fixed.size(); ++i) {
        int idx = fixed[i];
        m_increment[idx] = m_deltaV[idx] + v0[idx];
        m_rhs[idx] = (desired[i] - x0[idx]) / m_dt;
      }
      if ((m_increment - m_rhs).norm() < 1.0e-10) {
        break;
      }
      m_increment.setZero();
      m_A->setZero();
      STOP_TIMER("setup");
    }

    m_diffEq.endStep();
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
      m_A = m_diffEq.createMatrix();
      delete m_solver;
      m_solver = SolverUtils::instance()->createLinearSolver(m_A);
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
