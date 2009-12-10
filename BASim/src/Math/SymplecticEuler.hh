/**
 * \file SymplecticEuler.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef SYMPLECTICEULER_HH
#define SYMPLECTICEULER_HH

namespace BASim {

/** This class implements the symplectic Euler time-stepping
    method. This assumes a pair of equations of the form
    \f{eqnarray*}\frac{dx}{dt} &=& v \\ M\frac{dv}{dt} &=& f(t,x) \
    . \f} Further, \f$M\f$ is assumed to be diagonal. */
template <class ODE>
class SymplecticEuler : public DiffEqSolver
{
public:

  SymplecticEuler(ODE& ode)
    : m_diffEq(ode)
  {}

  void execute()
  {
    m_pDot.resize(m_diffEq.ndof());
    m_pDot.setZero();

    m_diffEq.evaluatePDot(m_pDot);
    for (int i = 0; i < m_diffEq.ndof(); ++i) {
      Scalar v = m_diffEq.getV(i) + m_dt * m_pDot(i) / m_diffEq.getMass(i);
      m_diffEq.setV(i, v);
      m_diffEq.setX(i, m_diffEq.getX(i) + m_dt * v);
    }

    m_diffEq.flush();
  }

protected:

  ODE& m_diffEq;
  VecXd m_pDot;
};

} // namespace BASim

#endif // SYMPLECTICEULER_HH
