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

  explicit SymplecticEuler(ODE& ode)
    : m_diffEq(ode)
    , m_pDot()
    , m_x()
    , m_v()
    , m_m()
  {}

  bool execute()
  {
    std::vector<int> fixed;
    std::vector<Scalar> desired;
    std::vector<Scalar> desiredv;
    m_diffEq.getScriptedDofs(fixed, desired); // fixed are DOF indices, desired are corresponding desired values
    assert(fixed.size() == desired.size());
    desiredv.resize(fixed.size());
    
    m_pDot.resize(m_diffEq.ndof());
    m_pDot.setZero();
    m_diffEq.evaluatePDot(m_pDot);

    m_x.resize(m_diffEq.ndof());
    m_diffEq.getX(m_x);

    m_v.resize(m_diffEq.ndof());
    m_diffEq.getV(m_v);

    m_m.resize(m_diffEq.ndof());
    m_diffEq.getMass(m_m);
    
    for (size_t i = 0; i < fixed.size(); i++)
      desiredv[fixed[i]] = (desired[fixed[i]] - m_x[fixed[i]]) / m_dt;

    m_v.array() += m_dt*(m_pDot.array()/m_m.array());
    m_x += m_dt*m_v;
    
    for (size_t i = 0 ;i < fixed.size(); i++)
    {
      m_x[fixed[i]] = desired[fixed[i]];
      m_v[fixed[i]] = desiredv[fixed[i]];
    }
    std::cout << "v = " << m_v << std::endl;
    
    
    m_diffEq.set_qdot(m_v);
//    m_diffEq.set_q(m_x);

//    for (int i = 0; i < m_diffEq.ndof(); ++i) {
//      Scalar v = m_diffEq.getV(i) + m_dt * m_pDot(i) / m_diffEq.getMass(i);
//      m_diffEq.setV(i, v);
//      m_diffEq.setX(i, m_diffEq.getX(i) + m_dt * v);
//    }

    m_diffEq.endIteration();
    
    m_diffEq.endStep();
    
    return true;
  }

  std::string getName() const
  {
    return "Symplectic Euler";
  }

protected:

  ODE& m_diffEq;
  VecXd m_pDot;
  VecXd m_x;
  VecXd m_v;
  VecXd m_m;
};

} // namespace BASim

#endif // SYMPLECTICEULER_HH
