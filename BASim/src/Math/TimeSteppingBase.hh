/**
 * \file DiffEqSolver.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef DIFFEQSOLVER_HH
#define DIFFEQSOLVER_HH

namespace BASim {

/** Base class for all time steppers. */
class DiffEqSolver
{
public:

  virtual ~DiffEqSolver() {}

  Scalar getTime() const { return m_time; }

  void setTime(Scalar time) { m_time = time; }

  Scalar getTimeStep() const { return m_dt; }

  void setTimeStep(Scalar dt) { m_dt = dt;} 

  virtual void execute() = 0;

  int getMaxIterations() const { return m_maxIterations; }

  void setMaxIterations(int iterations)
  {
    m_maxIterations = std::max(1, iterations);
  }

protected:

  DiffEqSolver(Scalar time = 0, Scalar dt = 0.1)
    : m_time(time)
    , m_dt(dt)
    , m_maxIterations(1)
  {}

  Scalar m_time; ///< the current time
  Scalar m_dt; ///< size of the time step

  int m_maxIterations; ///< maximum number of iterations (only for implicit methods)
};

} // namespace BASim

#endif // DIFFEQSOLVER_HH
