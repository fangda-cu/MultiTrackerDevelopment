/**
 * \file StaticSolver.hh
 *
 * \author Samantha Ainsley and Eitan Grinspun
 * \date 06/7/2011
 */

#ifndef STATICSOLVER_HH
#define STATICSOLVER_HH

#include "TimeSteppingBase.hh"
#include "MatrixBase.hh"
#include "LinearSolverBase.hh"
#include "SolverUtils.hh"
#include "../Core/Timer.hh"
#include "../Physics/ElasticRods/MinimalRodStateBackup.hh"
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
class StaticSolver: public DiffEqSolver
{
public:

    explicit StaticSolver(ODE& ode) :
        m_diffEq(ode), m_ndof(-1), x0(), m_rhs(), m_deltaX(), m_prevDeltaX(),
                m_increment(), m_fixed(), m_desired(), m_initl2norm(0), m_l2norm(0), m_energy(0), m_A(NULL), m_solver(NULL)
    {
        m_A = m_diffEq.createMatrix();
        m_solver = SolverUtils::instance()->createLinearSolver(m_A);

        m_maxlsit = 5;
    }

    ~StaticSolver()
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

    bool execute()
    {
        START_TIMER("Staticsolver::execute");
	START_TIMER("Staticsolver::execute/backup");
        m_diffEq.backupResize();
        m_diffEq.backup();
        STOP_TIMER("Staticsolver::execute/backup");
        if (position_solve())
        {
            START_TIMER("Staticsolver::execute/backup");
            m_diffEq.backupClear();
            STOP_TIMER("Staticsolver::execute/backup");
	    STOP_TIMER("Staticsolver::execute");
            return true;
        }
	else 
	{          
	  STOP_TIMER("Staticsolver::execute");
	  return false;
	}
    }

    std::string getName() const
    {
        return "StaticSolver";
    }

    void resize()
    {
        m_ndof = m_diffEq.ndof();
        if (x0.size() != m_ndof)
        {
            x0.resize(m_ndof);
            m_rhs.resize(m_ndof);
            m_deltaX.resize(m_ndof);
            m_prevDeltaX.resize(m_ndof);
            m_increment.resize(m_ndof);
        }
        assert(m_A->rows() == m_A->cols());
        if (m_A->rows() != m_ndof)
        {
            assert(m_A != NULL);
            delete m_A;
            m_A = m_diffEq.createMatrix();
            assert(m_solver != NULL);
            delete m_solver;
            m_solver = SolverUtils::instance()->createLinearSolver(m_A);
        }
    }

    void setZero()
    {
        x0.setZero();
        m_rhs.setZero();
        m_deltaX.setZero();
        m_increment.setZero();
        m_A->setZero();
	m_energy = 0;
    }

    // update computation of: forces (the RHS), residual, energy
    void updatePositionBasedQuantities()
    {
        // rhs == forces
        TraceStream(g_log, "") << "Staticsolver::computeResidual: evaluating PDot...\n";
        m_rhs.setZero();
	m_energy = 0;
        m_diffEq.evaluateConservativeForcesEnergy(m_rhs, m_energy);

        // For prescribed (fixed) DOFS, overwrite heuristic initial guess
        // Set deltaX for fixed DOFs
        for (int i = 0; i < (int) m_fixed.size(); ++i)
        {
            int dof = m_fixed[i];
            m_rhs(dof) = 0; // set zero desired increment
        }

        // Save the infinity norm
        m_infnorm = m_rhs.lpNorm<Eigen::Infinity> ();
        m_l2norm  = m_rhs.norm();
    }

    bool isConverged()
    {
        TraceStream(g_log, "") << "Staticsolver::isConverged: residual = " << m_l2norm << " infnorm = "
                << m_infnorm << " rel residual = " << m_l2norm / m_initl2norm << " inc norm = " << m_increment.norm()
                << " atol = " << m_atol << " inftol = " << m_inftol << " rtol = " << m_rtol << " stol = " << m_stol << '\n';

	// Energy was halved... good enough!
        if (m_energy < m_initEnergy)
        {
            TraceStream(g_log, "") << "Staticsolver::isConverged(): converged: energy halved from " << m_initEnergy << " to " << m_energy << "\n";
            return true;
        }

        // L2 norm of the residual is less than tolerance
        if (m_l2norm < m_atol)
        {
            TraceStream(g_log, "") << "Staticsolver::isConverged(): converged atol: residual = " << m_l2norm
                    << " < " << m_atol << " = atol " << '\n';
            return true;
        }
        // Infinity norm of residual is less than tolerance
        if (m_infnorm < m_inftol)
        {
            TraceStream(g_log, "") << "Staticsolver::isConverged(): converged inftol: |residual|_inf = " << m_infnorm
                    << " < " << m_inftol << " = inftol" << '\n';
            return true;
        }
        if (m_l2norm <= m_rtol * m_initl2norm)
        {
            TraceStream(g_log, "") << "Staticsolver::isConverged(): converged rtol: residual = " << m_l2norm
                    << " <= " << " (rtol = " << m_rtol << ") * (init. residual = " << m_initl2norm << ") = " << m_rtol
                    * m_initl2norm << '\n';
            return true;
        }
        // L2 norm of change in solution at last step of solve is less than tolerance
        if (m_alpha * m_increment.norm() < m_stol)
        {
            TraceStream(g_log, "") << "Staticsolver::isConverged(): converged stol: " << " |increment|_L2 < "
                    << m_stol << " = stol " << '\n';
            return true;
        }
        TraceStream(g_log, "") << "Staticsolver::isConverged(): convergence test fails" << '\n';
        return false;
    }

protected:

    bool position_solve()
    {
        START_TIMER("StaticSolver::position_solve/setup");

        // Chapter 0: Basic housekeeping
        ////////////////////////////////////////////////////

        bool successful_solve = true;

        m_diffEq.startStep();

        resize();
        setZero();
        m_diffEq.getScriptedDofs(m_fixed, m_desired); // m_fixed are DOF indices, m_desired are corresponding desired values
        assert(m_fixed.size() == m_desired.size());

#ifndef NDEBUG // ensure indices in valid range
        for (int i = 0; i < (int) m_fixed.size(); ++i)
            assert(m_fixed[i] >= 0);
        for (int i = 0; i < (int) m_fixed.size(); ++i)
            assert(m_fixed[i] < m_ndof);
#endif


        // Chapter 1: Set up initial guess for Newton Solver
        ////////////////////////////////////////////////////

        // m_deltaX is the difference between end of step and start of step positions

        // Copy start of step positions
        m_diffEq.getX(x0);
	  
	// Zero the initial guess
        m_deltaX.setZero();

        // Set deltaX for prescribed DOFs
        for (int i = 0; i < (int) m_fixed.size(); ++i)
        {
            int dof = m_fixed[i];
            m_rhs(dof) = m_desired[i] - x0(dof); // set desired position
        }

        // Signal the differential equation that it should get 
        // ready for the first iteration. In practice this is where the ODE
        // can precompute some reusable quantities that depend on the state (position & velocity)
        m_diffEq.endIteration();

        // Based on the initial guess, 
        //   1) set up right hand side (RHS = F).
        //   2) compute the potential energy
        //   3) compute the residual
	updatePositionBasedQuantities();

	// save the initial residual
        m_initl2norm = m_l2norm;

	m_initEnergy = m_energy;

        TraceStream(g_log, "")
                << "StaticSolver::position_solve: starting Newton solver. Initial energy = " << m_energy << " and residual = "
                << m_l2norm << ", convergence test will use thresholds atol = " << m_atol << " inftol = " << m_inftol
                << " rtol = " << m_initl2norm * m_rtol << " stol = " << m_stol << '\n';

        STOP_TIMER("StaticSolver::position_solve/setup");


        // Chapter 2: Iterate using Newton's Method
        ////////////////////////////////////////////////////

	Scalar lambda = 1e-8;

        int curit = 0;
        for (curit = 0; curit < m_maxit; ++curit)
        {
	    lambda *= 2.;

            TraceStream(g_log, "") << "\nStaticSolver::position_solve: Newton iteration = " << curit << " lambda = " << lambda << "\n";

            // TODO: Assert m_A, increment are zero
 
            START_TIMER("StaticSolver::position_solve/setup");

            // Set up LHS Matrix
            ////////////////////////

            // TODO: make the finalize() not virtual

            // The LHS is the Hessian of the potential energy
            // m_A = -dF/dx
            m_diffEq.evaluatePDotDX(-1, *m_A);
            m_A->finalize();
            assert(m_A->isApproxSymmetric(1.0e-6));

	    // Spectral shift
            for (int i = 0; i < m_ndof; ++i)
                m_A->add(i, i, lambda);
            m_A->finalize();
            assert(m_A->isApproxSymmetric(1.0e-6));

            // Set the rows and columns corresponding to fixed degrees of freedom to 0
            m_A->zeroRows(m_fixed, 1.0);
            m_A->finalize();

            m_A->zeroCols(m_fixed, 1.0);
            m_A->finalize();

            // Finalize the nonzero structure before the linear solve (for sparse matrices only)
            m_A->finalizeNonzeros();
            STOP_TIMER("StaticSolver::position_solve/setup");

            assert(m_A->isApproxSymmetric(1.0e-6));

            // Solve the linear system for the "Newton direction" m_increment
            //
            // Later, we will take the Newton step m_deltaX += m_increment
            // (or some scaled multiple of m_increment, as per line search)
            //////////////////////////////////////////////////////////////////

            START_TIMER("StaticSolver::position_solve/solver");
            int status = m_solver->solve(m_increment, m_rhs);
            STOP_TIMER("StaticSolver::position_solve/solver");
            if (status < 0)
            {
                DebugStream(g_log, "") << "\033[31;1mWARNING IN StaticSolver:\033[m Problem during linear solve detected. "
                        << '\n';
                return false;
            }

            START_TIMER("StaticSolver::position_solve/ls");

            m_alpha = 1.; // actual step will m_deltaX += alpha * m_increment

            // Save m_deltaX and residual for later
            m_prevDeltaX = m_deltaX;

	    Scalar prevEnergy = m_energy;

            // Attempt a full Newton step (alpha = 1)
            m_deltaX = m_prevDeltaX + m_alpha * m_increment;

	    // line search loop
            for (int i = 0;; i++)
            {
                // Evaluate residual for attempted increment
                //////////////////////////////////////////////

                // Update the differential equation with the current guess
                m_diffEq.set_qdot(m_deltaX / m_dt);
                m_diffEq.set_q(x0 + m_deltaX);

                // Signal the differential equation that it should recompute cached quantities
                m_diffEq.endIteration();

                updatePositionBasedQuantities();

                bool converged = isConverged();

                TraceStream(g_log, "") << "Staticsolver::position_solve: summary of line search i " << i
                        << ": increment " << m_increment.norm() << " energy " << m_energy << ", previous "
                        << prevEnergy << " " << " initial " << m_initEnergy << "\n";

                // Is this step acceptable?
                /////////////////////////////////////////////////////////////

                if (m_energy < prevEnergy || converged)
                {
                    TraceStream(g_log, "") << "Succeeded (done).\n\n";
                    break;
                }
                else if (i >= m_maxlsit)
                {
                    TraceStream(g_log, "")
                            << "Exceeded max iterations.\nStaticsolver::position_solve/line search: \033[31;1mWARNING IN IMPLICITEULER:\033[m Line search failed. Proceeding anyway.\n\n";
                    break;
                }
                else
                {
                    TraceStream(g_log, "") << "cutting increment and iterating.\n\n";
                }

                // Attempt a smaller step
                //////////////////////////////

                m_alpha *= .5;
                m_deltaX = m_prevDeltaX + m_alpha * m_increment;
            }

            STOP_TIMER("Staticsolver::position_solve/ls");

            // After the line search...
            ///////////////////////////////

            // Check for convergence.
            if (isConverged())
                break;

            // Check for exceeding limit on number of Newton iterations
            if (curit == m_maxit - 1)
            {
                DebugStream(g_log, "") << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Newton solver reached max iterations: "
                        << m_maxit << "\n";

		if (m_energy > m_initEnergy)
		{
		    // Solver failed -- undo the changes
		    /////////////////////////////////////////////////

		    // Update the differential equation with the current guess
		    m_diffEq.set_q(x0);

		    // Signal the differential equation that it should recompute cached quantities
		    m_diffEq.endIteration();

		    updatePositionBasedQuantities();
		    return true;
		}

		return false;
            }

            START_TIMER("Staticsolver::position_solve/setup");

            m_increment.setZero();
            m_A->setZero();

            // Allow the nonzero structure to be modified again (for sparse matrices only)
            m_A->resetNonzeros();

            STOP_TIMER("Staticsolver::position_solve/setup");

            // Now go back and begin next Newton iteration...
        }

        TraceStream(g_log, "") << "Staticsolver::position_solve: completed " << curit + 1 << " Newton iterations."
                << '\n';

        m_diffEq.endStep();

#ifndef NDEBUG      // Ensure that fixed DOFs are at their desired values
        if (successful_solve)
        {
            VecXd xf(m_ndof);
            m_diffEq.getX(xf);
            for (int i = 0; i < (int) m_fixed.size(); ++i)
                assert(approxEq(m_desired[i], xf(m_fixed[i]), 1.0e-6));
        }
#endif

        return true;
    }

    ODE& m_diffEq;

    int m_ndof;

    VecXd x0;
    VecXd m_rhs;
    VecXd m_deltaX;
    VecXd m_prevDeltaX;
    VecXd m_increment;

    Scalar m_alpha;

    IntArray m_fixed;
    std::vector<Scalar> m_desired;

    MatrixBase* m_A;
    LinearSolverBase* m_solver;

    Scalar m_energy; // EG: move to DiffEqSolver
    Scalar m_initEnergy;
    Scalar m_l2norm;
    Scalar m_initl2norm;
};

} // namespace BASim

#endif // STATICSOLVER_HH
