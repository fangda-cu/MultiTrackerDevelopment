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

static int static_solve_counter = 0;

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
        m_diffEq(ode), m_ndof(-1), x0(), m_rhs(), m_deltaX(),
                m_fixed(), m_desired(), m_l2norm(0), m_energy(0), m_A(NULL), m_solver(NULL)
    {
        m_A = m_diffEq.createMatrix();
        m_solver = SolverUtils::instance()->createLinearSolver(m_A);

        m_maxlsit = 5;

	m_lambdamin = 1e-8;
	m_lambdamax = 1e+10;
	m_lambda    = m_lambdamin;
	m_gearup    = 2.00; // above 1.0
	m_geardown  = 0.80; // below 1.0
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
	m_diffEq.updateCachedQuantities();
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


protected:

    inline Scalar clipvalue( Scalar minvalue, Scalar variable, Scalar maxvalue )
    {
	return fmin( fmax( variable, minvalue), maxvalue);
    }


    bool isConverged()
    {
        TraceStream(g_log, "StaticSolver::isConverged") << "residual = " << m_l2norm << " infnorm = "
                << m_infnorm << " atol = " << m_atol << " inftol = " << m_inftol << '\n';

        // L2 norm of the residual is less than tolerance
        if (m_l2norm < m_atol)
        {
            TraceStream(g_log, "StaticSolver::isConverged") << "converged atol: residual = " << m_l2norm
                    << " < " << m_atol << " = atol " << '\n';
            return true;
        }
        // Infinity norm of residual is less than tolerance
        if (m_infnorm < m_inftol)
        {
            TraceStream(g_log, "StaticSolver::isConverged") << "converged inftol: |residual|_inf = " << m_infnorm
                    << " < " << m_inftol << " = inftol" << '\n';
            return true;
        }
        // if (m_deltaX.norm() < m_stol)
        // {
        //     TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged(): converged stol: " << " |increment|_L2 < "
        //             << m_stol << " = stol " << '\n';
        //     return true;
        // }
        TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged(): convergence test fails" << '\n';
        return false;
    }


    bool position_solve()
    {
        START_TIMER("StaticSolver::newton_step/setup");

        // Chapter 0: Basic housekeeping
        ////////////////////////////////////////////////////

        bool successful_solve = true;        

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

        // Based on the initial guess, 
        //   1) set up right hand side (RHS = F).
        //   2) compute the potential energy
        //   3) compute the residual
	updatePositionBasedQuantities();

	// save the initial residual
        Scalar initl2norm = m_l2norm;

	Scalar initEnergy = m_energy;
	Scalar initLambda = m_lambda;

        STOP_TIMER("StaticSolver::position_solve/setup");



        // Chapter 2: Iterate using Newton's Method
        ////////////////////////////////////////////////////

	TraceStream(g_log, "StaticSolver::position_solve") << "call #" << ++static_solve_counter << " initial energy = " << initEnergy << " lambda = " << m_lambda << "\n";

        // TODO: Assert m_A, increment are zero
 
        START_TIMER("StaticSolver::position_solve/setup");

        // Set up LHS Matrix
        ////////////////////////

        // TODO: make the finalize() not virtual

        // The LHS is the minus the conservative force Jacobian
        // m_A = -dF/dx
        m_diffEq.evaluatePDotDX(-1, *m_A);
        m_A->finalize();
        assert(m_A->isApproxSymmetric(1.0e-6));

        for (int i = 0; i < m_ndof; ++i)
        {
	    // Spectral shift (Tikhonov regularization)
            //m_A->add(i, i, m_lambda);

	    // Levenberg-Marquardt diagonal shift
	    Scalar d = (*m_A)(i,i);
	    //TraceStream(g_log, "StaticSolver::position_solve") << "lambda = " << m_lambda << " d[" << i << "] = " << d << " (1.+m_lambda) * d = " << (1.+m_lambda) * d << "\n";	    
            m_A->set(i, i, (1.+m_lambda) * d);
	}
        m_A->finalize();
        assert(m_A->isApproxSymmetric(1.0e-6));

        // Boundary conditions: Set the rows and columns corresponding to fixed degrees of freedom to 0
        m_A->zeroRows(m_fixed, 1.0);
        m_A->finalize();

        m_A->zeroCols(m_fixed, 1.0);
        m_A->finalize();

        // Finalize the nonzero structure before the linear solve (for sparse matrices only)
        m_A->finalizeNonzeros();
        STOP_TIMER("StaticSolver::position_solve/setup");
	assert(m_A->isApproxSymmetric(1.0e-6));


        // Solve the linear system for the Newton step m_deltaX
        //
        //////////////////////////////////////////////////////////////////

        START_TIMER("StaticSolver::position_solve/solver");
        int status = m_solver->solve(m_deltaX, m_rhs);
        STOP_TIMER("StaticSolver::position_solve/solver");
        if (status < 0)
        {
	    // shrink trust region (increase regularization)
	    m_lambda = clipvalue( m_lambdamin, m_lambda * m_gearup, m_lambdamax );

            DebugStream(g_log, "StaticSolver::position_solve") << "\033[31;1mWARNING IN StaticSolver:\033[m Problem during linear solve detected. "
				   << " new lambda = " << m_lambda << "\n";

	    m_diffEq.set_q(x0);
	    m_diffEq.updateCachedQuantities();
	    
            return false;
        }


	// Evaluate residual for attempted increment
        //////////////////////////////////////////////

        // Update the differential equation with the current guess
        m_diffEq.set_q(x0 + m_deltaX);
        m_diffEq.updateCachedQuantities();

        updatePositionBasedQuantities();

        // Update the Levenberg-Marquardt trust region size
        /////////////////////////////////////////////////////////////

        if (isConverged())
	{
	    // Leave lambda alone

	    TraceStream(g_log, "StaticSolver::position_solve") << "new energy = " << m_energy << "; new residual = " << m_l2norm << "; retaining step; converged! new lambda = " << m_lambda << "\n";
	}
        else if (m_energy < initEnergy) // || m_l2norm < initl2norm)
        {
	    // Decrease lambda (= increase trust region size = decrease regularization)
	    m_lambda = clipvalue( m_lambdamin, m_lambda * m_geardown, m_lambdamax );
            TraceStream(g_log, "StaticSolver::position_solve") << "new energy = " << m_energy << "; new residual = " << m_l2norm << "; retaining step; new lambda = " << m_lambda << "\n";
        }
        else
        {
	    // Increase lambda (= decrease trust region size = increase regularization)
	    m_lambda = clipvalue( m_lambdamin, m_lambda * m_gearup, m_lambdamax );
	    TraceStream(g_log, "StaticSolver::position_solve") << "new energy " << m_energy << "; new residual = " << m_l2norm << "; discarding step; new lambda = " << m_lambda << "\n";
	    
	    // Solver failed -- undo the changes
	    /////////////////////////////////////////////////
	    
	    m_diffEq.set_q(x0);
	    m_diffEq.updateCachedQuantities();
        }

	START_TIMER("Staticsolver::position_solve/setup");

	// Allow the nonzero structure to be modified again (for sparse matrices only)
	m_A->resetNonzeros();

	STOP_TIMER("Staticsolver::position_solve/setup");


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

    Scalar m_lambda;
    Scalar m_lambdamin, m_lambdamax;
    Scalar m_gearup, m_geardown;

    IntArray m_fixed;
    std::vector<Scalar> m_desired;

    MatrixBase* m_A;
    LinearSolverBase* m_solver;

    Scalar m_energy; // EG: move to DiffEqSolver
    Scalar m_l2norm;
};

} // namespace BASim

#endif // STATICSOLVER_HH






