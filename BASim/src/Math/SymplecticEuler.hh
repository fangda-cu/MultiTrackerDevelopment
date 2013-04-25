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
    // adaptive time stepping
    Scalar dt = m_dt;
    Scalar dt_left = dt;
    Scalar dt_largest_possible = 0;

    int substep_count = 0;
    while (dt_left > 0)
    {
      m_diffEq.startStep();
      
      m_diffEq.startIteration();
      
      std::vector<int> fixed;
      std::vector<Scalar> desired;
      std::vector<Scalar> desiredv;
      m_diffEq.getScriptedDofs(fixed, desired); // fixed are DOF indices, desired are corresponding desired values
      assert(fixed.size() == desired.size());
      desiredv.resize(fixed.size());
      
      m_pDot.resize(m_diffEq.ndof());
      m_pDot.setZero();
      m_diffEq.evaluatePDot(m_pDot);
      
      // determine max dt
      dt_largest_possible = m_diffEq.determineMaxDt(m_pDot);
      std::cout << "dt_left = " << dt_left << " max dt = " << dt_largest_possible << std::endl;
      if (dt_largest_possible >= dt_left)
        dt_largest_possible = dt_left;

      m_x.resize(m_diffEq.ndof());
      m_diffEq.getX(m_x);

      m_v.resize(m_diffEq.ndof());
      m_diffEq.getV(m_v);

      m_m.resize(m_diffEq.ndof());
      m_diffEq.getMass(m_m);
      
      for (size_t i = 0; i < fixed.size(); i++)
        desiredv[i] = (desired[i] - m_x[fixed[i]]) / dt_largest_possible;

      m_v = m_pDot;
      m_x += dt_largest_possible * m_v;
      dt_left -= dt_largest_possible;
      
      for (size_t i = 0 ;i < fixed.size(); i++)
      {
        m_x[fixed[i]] = desired[i];
        m_v[fixed[i]] = desiredv[i];
      }
      
      m_diffEq.set_qdot(m_v);
      m_diffEq.set_q(m_x);

      m_diffEq.endIteration();
      
      m_diffEq.endStep();
      
      substep_count++;
    }
    
    std::cout << "Symplectic Euler time step finished after " << substep_count << " substeps." << std::endl;
    
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
