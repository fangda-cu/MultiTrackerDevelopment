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

    bool execute()
    {
        START_TIMER("SymmetricImplicitEuler::execute");
        START_TIMER("SymmetricImplicitEuler::execute/backup");
        m_diffEq.backupResize();
        m_diffEq.backup();
        STOP_TIMER("SymmetricImplicitEuler::execute/backup");
        if (position_solve(0))
        {
#ifdef TIMING_ON
            IntStatTracker::getIntTracker("INITIAL_ITERATE_1_SUCCESSES") += 1;
#endif
            START_TIMER("SymmetricImplicitEuler::execute/backup");
            m_diffEq.backupClear();
            STOP_TIMER("SymmetricImplicitEuler::execute/backup");
            STOP_TIMER("SymmetricImplicitEuler::execute");
            return true;
        }

        //DebugStream(g_log, "") << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Newton solver failed to converge in max iterations: " << m_maxit << ". Attempting alternate initial iterate 2." << '\n';
        START_TIMER("SymmetricImplicitEuler::execute/backup");
        m_diffEq.backupRestore();
        STOP_TIMER("SymmetricImplicitEuler::execute/backup");
        if (position_solve(1))
        {
#ifdef TIMING_ON
            IntStatTracker::getIntTracker("INITIAL_ITERATE_2_SUCCESSES") += 1;
#endif
            START_TIMER("SymmetricImplicitEuler::execute/backup");
            m_diffEq.backupClear();
            STOP_TIMER("SymmetricImplicitEuler::execute/backup");
            STOP_TIMER("SymmetricImplicitEuler::execute");
            return true;
        }

        //DebugStream(g_log, "") << "                          Newton solver failed to converge in max iterations: " << m_maxit << ". Attempting alternate initial iterate 3." << '\n';
        START_TIMER("SymmetricImplicitEuler::execute/backup");
        m_diffEq.backupRestore();
        STOP_TIMER("SymmetricImplicitEuler::execute/backup");
        if (position_solve(2))
        {
#ifdef TIMING_ON
            IntStatTracker::getIntTracker("INITIAL_ITERATE_3_SUCCESSES") += 1;
#endif
            START_TIMER("SymmetricImplicitEuler::execute/backup");
            m_diffEq.backupClear();
            STOP_TIMER("SymmetricImplicitEuler::execute/backup");
            STOP_TIMER("SymmetricImplicitEuler::execute");
            return true;
        }

        //DebugStream(g_log, "") << "                          Newton solver failed to converge in max iterations: " << m_maxit << ". Attempting alternate initial iterate 4." << '\n';
        START_TIMER("SymmetricImplicitEuler::execute/backup");
        m_diffEq.backupRestore();
        STOP_TIMER("SymmetricImplicitEuler::execute/backup");
        if (position_solve(3))
        {
#ifdef TIMING_ON
            IntStatTracker::getIntTracker("INITIAL_ITERATE_4_SUCCESSES") += 1;
#endif
            START_TIMER("SymmetricImplicitEuler::execute/backup");
            m_diffEq.backupClear();
            STOP_TIMER("SymmetricImplicitEuler::execute/backup");
	    STOP_TIMER("SymmetricImplicitEuler::execute");
            return true;
        }

        //DebugStream(g_log, "") << "                          Newton solver failed to converge in max iterations: " << m_maxit << ". Attempting alternate initial iterate 4." << '\n';
        START_TIMER("SymmetricImplicitEuler::execute/backup");
        m_diffEq.backupRestore();
        STOP_TIMER("SymmetricImplicitEuler::execute/backup");
        if (position_solve(4))
        {
#ifdef TIMING_ON
            IntStatTracker::getIntTracker("INITIAL_ITERATE_5_SUCCESSES") += 1;
#endif
            START_TIMER("SymmetricImplicitEuler::execute/backup");
            m_diffEq.backupClear();
            STOP_TIMER("SymmetricImplicitEuler::execute/backup");
	    STOP_TIMER("SymmetricImplicitEuler::execute");
            return true;
        }

      //  DebugStream(g_log, "") << "\033[31;1mWARNING IN SYM IMPLICITEULER:\033[m Newton solver failed to converge in max iterations: "
      //          << m_maxit << "." << '\n';
        START_TIMER("SymmetricImplicitEuler::execute/backup");
        m_diffEq.backupClear();
        STOP_TIMER("SymmetricImplicitEuler::execute/backup");

        STOP_TIMER("SymmetricImplicitEuler::execute");
        return false;
    }

    std::string getName() const
    {
        return "Symmetric Implicit Euler";
    }

    void resize()
    {
        m_ndof = m_diffEq.ndof();
        if (m_mass.size() != m_ndof)
        {
            m_mass.resize(m_ndof);
            x0.resize(m_ndof);
            v0.resize(m_ndof);
            m_rhs.resize(m_ndof);
            m_deltaX.resize(m_ndof);
            m_deltaX_save.resize(m_ndof);
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
        v0.setZero();
        m_rhs.setZero();
        m_deltaX.setZero();
        m_increment.setZero();
        m_A->setZero();
    }

    // Computes the inf norm of the residual
    Scalar computeResidual()
    {
        // Sanity checks for NANs
//        assert((x0.cwise() == x0).all());
//        assert((m_deltaX.cwise() == m_deltaX).all());

        // rhs == h*h*forces
        m_rhs.setZero();
        m_diffEq.evaluatePDot(m_rhs);
        m_rhs *= m_dt * m_dt;

        // lhs == M*deltaV == M*(deltax-h*v_n)
        m_rhs.array() -= m_mass.array() * (m_deltaX - m_dt * v0).array();

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
         TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged: residual = " << m_residual
	   << " infnorm = " << m_infnorm
	   << " rel residual = " << m_residual / m_initial_residual
	   << " inc norm = " << m_increment.norm() 
	   << " atol = " << m_atol
	   << " inftol = " << m_inftol
	   << " rtol = " << m_rtol
	   << " stol = " << m_stol << '\n';
         
        // L2 norm of the residual is less than tolerance
        if (m_residual < m_atol)
        {
	    TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged(): converged atol: residual = " << m_residual << " < " << m_atol << " = atol " << '\n';
            return true;
        }
        // Infinity norm of residual is less than tolerance
        if (m_infnorm < m_inftol)
        {
	    TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged(): converged inftol: |residual|_inf = " << m_infnorm << " < " << m_inftol << " = inftol" << '\n';
            return true;
        }
        if (m_residual <= m_rtol * m_initial_residual)
        {
	    TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged(): converged rtol: residual = " << m_residual << " <= " << " (rtol = " << m_rtol << ") * (init. residual = " << m_initial_residual << ") = " << m_rtol * m_initial_residual << '\n';
            return true;
        }
        // L2 norm of change in solution at last step of solve is less than tolerance
        if (m_increment.norm() < m_stol)
        {
 	    TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged(): converged stol: " << " |increment|_L2 < " << m_stol << " = stol "<< '\n';
            return true;
        }
	TraceStream(g_log, "") << "SymmetricImplicitEuler::isConverged(): convergence test fails" << '\n';
        return false;
    }

protected:

    void generateInitialIterate1(VecXd& dx)
    {
        dx = m_dt * v0; // explicit inertial step
    }

    void generateInitialIterate2(VecXd& dx)
    {
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

    bool position_solve(int guess_to_use)
    {
        START_TIMER("SymmetricImplicitEuler::position_solve/setup");

	// Chapter 0: Basic housekeeping
	////////////////////////////////////////////////////

        bool successful_solve = true;

        m_diffEq.startStep();

        resize();
        setZero();
        m_diffEq.getScriptedDofs(m_fixed, m_desired); // m_fixed are DOF indices, m_desired are corresponding desired values
        assert(m_fixed.size() == m_desired.size());

#ifndef NDEBUG // ensure indices in valid range
        for( int i = 0; i < (int) m_fixed.size(); ++i ) assert( m_fixed[i] >= 0 );
        for( int i = 0; i < (int) m_fixed.size(); ++i ) assert( m_fixed[i] < m_ndof );
#endif

	// TODO: EG: Is this copying of data REALLY NEEDED? Probably not.
        // Copy masses.
        if (!m_mass_set) // Assuming masses are constant
        {
            m_diffEq.getMass(m_mass);
            m_mass_set = true;
        }



	// Chapter 1: Set up initial guess for Newton Solver
	////////////////////////////////////////////////////

	// m_deltaX is the difference between end of step and start of step positions

        // Copy start of step positions and velocities.
        m_diffEq.getX(x0);
        m_diffEq.getV(v0);

        // Initialize guess for the root.
        switch (guess_to_use)
        {
        case 0:
        {
            generateInitialIterate1(m_deltaX);
            break;
        }
        case 1:
        {
            generateInitialIterate2(m_deltaX);
            break;
        }
        case 2:
        {
            generateInitialIterate3(m_deltaX);
            break;
        }
        case 3:
        {
            generateInitialIterate4(m_deltaX);
            break;
        }
        case 4:
        {
            generateInitialIterate5(m_deltaX);
            break;
        }
        default:
        {
            DebugStream(g_log, "") << "\033[31;1mERROR IN IMPLICITEULER:\033[m Invalid initial iterate requested, exiting." << '\n';
	    STOP_TIMER("SymmetricImplicitEuler::position_solve/setup");
            exit(1);
            break;
        }
        }

	// For prescribed (fixed) DOFS, overwrite heuristic initial guess
        // Set deltaX and v0 for fixed DOFs
        for (int i = 0; i < (int) m_fixed.size(); ++i)
        {
            int dof = m_fixed[i];
            m_deltaX(dof) = m_desired[i] - x0(dof);     // set desired position
            v0(dof) = (m_desired[i] - x0(dof)) / m_dt;  // reverse engineer desired velocity
        }

     
        #ifndef NDEBUG // check that the desired velocity was correctly computed by taking virtual forward step
            for(int i = 0; i < (int) m_fixed.size(); ++i) 
	      assert( approxEq(m_desired[i], /* =approx= */  x0(m_fixed[i]) + v0[m_fixed[i]]*m_dt, 1.0e-6) );
        #endif

        // Update the differential equation with the current guess
	m_diffEq.set_qdot(m_deltaX / m_dt);  // set velocity
        m_diffEq.set_q   (x0 + m_deltaX  );  // set position

        // Signal the differential equation that it should get 
        // ready for the first iteration. In practice this is where the ODE
        // can precompute some reusable quantities that depend on the state (position & velocity)
        m_diffEq.endIteration();


        // Based on the initial guess, 
        //   1) set up right hand side (RHS) of implicit Euler: RHS = M(m_dt*v_n-m_deltaX) + h^2*F.
        //   2) compute the residual of the ODE
        //   3) cache the residual as m_initial_residual, for convergence test later
        m_initial_residual = m_residual = computeResidual();

        TraceStream(g_log, "") << "SymmetricImplicitEuler::position_solve: starting Newton solver. Initial guess has residual = " << m_residual
	     << ", convergence test will use thresholds atol = " << atol << " infnorm = " << m_infnorm
	     << " rtol = " << m_residual / m_initial_residual
	     << " stol = " << m_increment.norm() << '\n';


        STOP_TIMER("SymmetricImplicitEuler::position_solve/setup");



	// Chapter 2: Iterate using Newton's Method
	////////////////////////////////////////////////////

        int curit = 0;
        for (curit = 0; curit < m_maxit; ++curit)
        {
	  // TraceStream(g_log, "") << "\nSymmetricImplicitEuler::position_solve: iteration = " << m_curit << "\n\n" << '\n';

            // TODO: Assert m_A, increment are zero

            START_TIMER("SymmetricImplicitEuler::position_solve/setup");


	    // Set up RHS
	    ///////////////////////

	    // Update the RHS for the fixed DOFs 
	    // TODO: Note, this really should be done at the end of the Newton iteration, since it was already done once in Chapter 1
            for (int i = 0; i < (int) m_fixed.size(); ++i)
                m_rhs(m_fixed[i]) = m_dt * v0(m_fixed[i]) - m_deltaX(m_fixed[i]);

	    

	    // Set up LHS Matrix
	    ////////////////////////

	    // TODO: make the finalize() not virtual

            // Consider LHS arising from potential forces (function of position)
            // m_A = -h^2*dF/dx
            m_diffEq.evaluatePDotDX(-m_dt * m_dt, *m_A);
            m_A->finalize();
	    assert( m_A->isApproxSymmetric(1.0e-6) );

	    // Consider LHS arising from dissipative forces (function of velocity)
            // m_A = -h*dF/dv -h^2*dF/dx
            m_diffEq.evaluatePDotDV(-m_dt, *m_A);
            m_A->finalize();
            assert( m_A->isApproxSymmetric(1.0e-6) );

	    // Consider inertial contribution from mass matrix
            // m_A = M -h*dF/dv -h^2*dF/dx
            for (int i = 0; i < m_ndof; ++i)
                m_A->add(i, i, m_mass(i));
            m_A->finalize();
            assert( m_A->isApproxSymmetric(1.0e-6) );

            // Set the rows and columns corresponding to fixed degrees of freedom to 0
            m_A->zeroRows(m_fixed, 1.0);
            m_A->finalize();

            m_A->zeroCols(m_fixed, 1.0);
            m_A->finalize();

            // Finalize the nonzero structure before the linear solve (for sparse matrices only)
            m_A->finalizeNonzeros();
            STOP_TIMER("SymmetricImplicitEuler::position_solve/setup");

            assert(m_A->isApproxSymmetric(1.0e-6));



	    // Solve the linear system for the "Newton direction" m_increment
	    //
	    // Later, we will take the Newton step m_deltaX += m_increment
	    // (or some scaled multiple of m_increment, as per line search)
	    //////////////////////////////////////////////////////////////////

            START_TIMER("SymmetricImplicitEuler::position_solve/solver");
            int status = m_solver->solve(m_increment, m_rhs);
            STOP_TIMER("SymmetricImplicitEuler::position_solve/solver");
            if (status < 0)
            {
                DebugStream(g_log, "") << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. " << '\n';
                return false;
            }



	    START_TIMER("SymmetricImplicitEuler::position_solve/ls");

            //double alpha = 1.; // actual step will m_deltaX += alpha * m_increment

            // Save m_deltaX and residual for later
            m_deltaX_save = m_deltaX;
            double previous_residual = m_residual;

	    // Attempt a full Newton step (alpha = 1)
            m_deltaX += m_increment;

            for (int i = 0;; i++)
            {
		// Evaluate residual for attempted increment
		//////////////////////////////////////////////

                // Update the differential equation with the current guess
                m_diffEq.set_qdot(m_deltaX / m_dt);
                m_diffEq.set_q   (x0 + m_deltaX  );

	        // Signal the differential equation that it should recompute cached quantities
                m_diffEq.endIteration();

                // Calling computeResidual also sets m_rhs = M(m_dt*v_n-m_deltaX) + h^2*F.
                m_residual = computeResidual();

                TraceStream(g_log, "") << "\nSymmetricImplicitEuler::position_solve/line search: i "<<i<<", increment " << m_increment.norm() << " previous "<<previous_residual<<", residual "<<m_residual<<'\n';

		
		// Is this residual (hence the increment) acceptable?
		/////////////////////////////////////////////////////////////

                if (m_residual < .9 * previous_residual || isConverged())
                {
                    TraceStream(g_log, "") << "Line search succeeded." << '\n';
                    break;
	        }
                else if (i >= m_maxlsit)
                {
                    TraceStream(g_log, "") << "SymmetricImplicitEuler::position_solve/line search: \033[31;1mWARNING IN IMPLICITEULER:\033[m Line search failed. Proceeding anyway." << '\n';
                    //return false;
		    break;
                }
		else 
		{
		    TraceStream(g_log, "") << "SymmetricImplicitEuler::position_solve/line search: cutting increment and iterating." << '\n';
		}

		// Attempt a smaller step
		//////////////////////////////

                //alpha       *= .5;
		m_increment *= .5;
                m_deltaX = m_deltaX_save + m_increment;
            }
	    STOP_TIMER("SymmetricImplicitEuler::position_solve/ls");


	    // After the line search...
	    ///////////////////////////////

            // Check for convergence.
            if (isConverged())
                break;

	    // Check for exceeding limit on number of Newton iterations
            if ( curit == m_maxit - 1)
	    {
                DebugStream(g_log, "") << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Newton solver reached max iterations: " << m_maxit << '\n';
		return false;
            }

           
	    START_TIMER("SymmetricImplicitEuler::position_solve/setup");

            m_increment.setZero();
            m_A->setZero();

            // Allow the nonzero structure to be modified again (for sparse matrices only)
            m_A->resetNonzeros();

	    STOP_TIMER("SymmetricImplicitEuler::position_solve/setup");

	    // Now go back and begin next Newton iteration...
        } 

        TraceStream(g_log, "") << "SymmetricImplicitEuler::position_solve: completed " << curit+1 << " Newton iterations." << '\n';
          
        m_diffEq.endStep();

#ifndef NDEBUG      // Ensure that fixed DOFs are at their desired values
        if( successful_solve )
        {
            VecXd xf(m_ndof);
            m_diffEq.getX(xf);
            for( int i = 0; i < (int) m_fixed.size(); ++i ) assert( approxEq(m_desired[i],xf(m_fixed[i]),1.0e-6) );
        }
#endif

        return true;
    }

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

    IntArray m_fixed;
    std::vector<Scalar> m_desired;

    Scalar m_initial_residual;
    Scalar m_residual;

    MatrixBase* m_A;
    LinearSolverBase* m_solver;

};

} // namespace BASim

#endif // SYMMETRICIMPLICITEULER_HH
