/**
 * \file SymmetricImplicitEuler.hh
 *
 * \author smith@cs.columbia.edu
 * \date 06/26/2010
 */

#ifndef SYMMETRICIMPLICITEULER_HH
#define SYMMETRICIMPLICITEULER_HH

#include "TimeSteppingBase.hh"
#include "MatrixBase.hh"
#include "LinearSolverBase.hh"
#include "SolverUtils.hh"
#include "../Core/Timer.hh"
//#include "../Physics/ElasticRods/MinimalRodStateBackup.hh"
#include "../Core/StatTracker.hh"
#include "../Util/TextLog.hh"

namespace BASim
{

/** This class implements the implicit Euler time-stepping
 method. This assumes a pair of equations of the form
 \f{eqnarray*}\frac{dx}{dt} &=& v \\ M\frac{dv}{dt} &=& f(t,x) \
   . \f} The mass matrix \f$M\f$ is assumed to be diagonal.
 */
template<class ODE>
class SymmetricImplicitEuler: public DiffEqSolver
{
public:

    explicit SymmetricImplicitEuler(ODE& ode) :
        m_diffEq(ode), m_ndof(-1), m_mass(), m_mass_set(false), x0(), v0(), m_rhs(), m_deltaX(), m_deltaX_save(),
                m_increment(), m_fixed(), m_desired(), m_initial_residual(0), m_residual(0), m_A(NULL), m_solver(NULL)
    {
        m_A = m_diffEq.createMatrix();
        m_solver = SolverUtils::instance()->createLinearSolver(m_A);

#ifdef TIMING_ON
        IntStatTracker::getIntTracker("INITIAL_ITERATE_1_SUCCESSES",0);
        IntStatTracker::getIntTracker("INITIAL_ITERATE_2_SUCCESSES",0);
        IntStatTracker::getIntTracker("INITIAL_ITERATE_3_SUCCESSES",0);
        IntStatTracker::getIntTracker("INITIAL_ITERATE_4_SUCCESSES",0);
        IntStatTracker::getIntTracker("INITIAL_ITERATE_5_SUCCESSES",0);
#endif
    }

    ~SymmetricImplicitEuler()
    {
        assert(m_A != NULL);
        assert(m_solver != NULL);

        if (m_A != NULL)
        {
            delete m_A;
            m_A = NULL;
        }

        if (m_solver != NULL)
        {
            delete m_solver;
            m_solver = NULL;
        }
    }

    bool execute();

    std::string getName() const
    {
        return "Symmetric Implicit Euler";
    }

   void resize();

    void setZero();

    // Computes the inf norm of the residual
    Scalar computeResidual()
    {
        // Sanity checks for NANs
        //        assert((x0.cwise() == x0).all());
        //        assert((m_deltaX.cwise() == m_deltaX).all());

        // rhs == h*h*forces
       
        TraceStream(g_log, "") << "SymmetricImplicitEuler::computeResidual: evaluating PDot...\n";
        m_rhs.setZero();
        m_diffEq.evaluatePDot(m_rhs);
      
        // FD 20121217: First order physics
//        m_rhs *= m_dt * m_dt;
        m_rhs *= m_dt;

        // lhs == M*deltaV == M*(deltax-h*v_n)
//        m_rhs.array() -= m_mass.array() * (m_deltaX - m_dt * v0).array();
        // FD 20121217: No mass
        m_rhs -= m_deltaX;
      
        for (int i = 0; i < (int) m_fixed.size(); ++i)
            m_rhs(m_fixed[i]) = 0.0;

        // Save the infinity norm
        m_infnorm = m_rhs.lpNorm<Eigen::Infinity> ();

        // TODO: EG: Be consistent. Write all norms to member fields, or return all norms, but not half and half

        // Return the L2 norm
        return m_rhs.norm();
    }

    bool isConverged()
    {
        TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged: residual = " << m_residual << " infnorm = "
                << m_infnorm << " rel residual = " << m_residual / m_initial_residual << " inc norm = " << m_increment.norm()
                << " atol = " << m_atol << " inftol = " << m_inftol << " rtol = " << m_rtol << " stol = " << m_stol << '\n';

        // L2 norm of the residual is less than tolerance
        if (m_residual < m_atol)
        {
            TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged(): converged atol: residual = " << m_residual
                    << " < " << m_atol << " = atol " << '\n';
            return true;
        }
        // Infinity norm of residual is less than tolerance
        if (m_infnorm < m_inftol)
        {
            TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged(): converged inftol: |residual|_inf = " << m_infnorm
                    << " < " << m_inftol << " = inftol" << '\n';
            return true;
        }
        if (m_residual <= m_rtol * m_initial_residual)
        {
            TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged(): converged rtol: residual = " << m_residual
                    << " <= " << " (rtol = " << m_rtol << ") * (init. residual = " << m_initial_residual << ") = " << m_rtol
                    * m_initial_residual << '\n';
            return true;
        }
        // L2 norm of change in solution at last step of solve is less than tolerance
       /* if (m_alpha * m_increment.norm() < m_stol)
        {
            TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged(): converged stol: " << " |increment|_L2 < "
                    << m_stol << " = stol " << '\n';
            return true;
        }*/
        TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged(): convergence test fails" << '\n';
        //std::cout << "Convergence test failed\n";
        return false;
    }

protected:

    // Initial guess based on rigid motion of the first two vertices
    bool generateInitialIterate0(VecXd& dx);

    void generateInitialIterate1(VecXd& dx)
    {
        dx = m_dt * v0; // explicit inertial step
    }

    void generateInitialIterate2(VecXd& dx)
    {
        TraceStream(g_log, "") << "Initial guess = initial position\n";
        dx.setZero(); // zero motion (not moving)
    }

    void generateInitialIterate3(VecXd& dx)
    {
        dx = 0.5 * m_dt * v0; // midway through explicit inertial step
    }

    void generateInitialIterate4(VecXd& dx)
    {
        dx = 1.5 * m_dt * v0; // 150% of explicit intertial step
    }

    void generateInitialIterate5(VecXd& dx)
    {
        dx = 2.0 * m_dt * v0; // 200% of explicit intertial step
    }

    bool position_solve(int guess_to_use); // Implementation moved to SymmetricImplicitEuler.cc to save compilation time

    ODE& m_diffEq;

    int m_ndof;

    VecXd m_mass;
    bool m_mass_set;
    VecXd x0;
    VecXd v0;
    VecXd m_rhs;
    VecXd m_deltaX;
    VecXd m_deltaX_save;
    //VecXd m_deltaV;
    VecXd m_increment;

    Scalar m_alpha;

    IntArray m_fixed;
    std::vector<Scalar> m_desired;

    Scalar m_initial_residual;
    Scalar m_residual;

    MatrixBase* m_A;
    LinearSolverBase* m_solver;

};

} // namespace BASim

#endif // SYMMETRICIMPLICITEULER_HH
