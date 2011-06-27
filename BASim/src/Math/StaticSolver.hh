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
        m_diffEq(ode), m_ndof(-1), m_x0(), m_rhs(), m_deltaX(),
                m_fixed(), m_desired(), m_l2norm(0), m_energy(0), m_A(NULL), m_solver(NULL)
    {
        m_A = m_diffEq.createMatrix();
        m_solver = SolverUtils::instance()->createLinearSolver(m_A);

        m_maxlsit = 5;

        // NOTE (sainsley) : this may be something we want to configure
        // on the maya side
        m_lambdamin = 1e-8;
        m_lambdamax = 1e+10;
        m_lambda    = 1e-3;
        m_gearup    = 3.00; // above 1.0
        m_geardown  = 0.50; // below 1.0
        m_failurecount = 0;
        m_successcount = 1;
        m_keepUpdating = true;
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
    	if ( m_keepUpdating )
    	{
    		position_solve();
    	}
    	return true;
    }

    std::string getName() const
    {
        return "StaticSolver";
    }

    void resize()
    {
        m_ndof = m_diffEq.ndof();
        if (m_x0.size() != m_ndof)
        {
            m_x0.resize(m_ndof);
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
        m_x0.setZero();
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

    inline Scalar funnyclipvalue( Scalar minvalue, Scalar variable, Scalar maxvalue )
    {
    	// funny wrap-around behavior ensures we don't get "stuck" at lambda=maxvalue
    	if (variable > maxvalue) variable = minvalue;

    	// clip
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


    bool examine_solution()
    {
        bool keepSolution = false;

        // Evaluate residual for attempted increment
        //////////////////////////////////////////////

        // Update the differential equation with the current guess
        m_diffEq.set_q(m_x0 + m_deltaX);

        updatePositionBasedQuantities();

        // Update the Levenberg-Marquardt trust region size
        /////////////////////////////////////////////////////////////

        assert(approxEq(m_energy*m_dt, m_l2norm, 1e-8));
        TraceStream(g_log, "StaticSolver::position_solve energy-force check") << m_energy << " " << m_l2norm << " " << m_diffEq.getTimeStep() << " " << m_energy<< "\n";
        if (isConverged())
        {
        	// Leave lambda alone
        	TraceStream(g_log, "StaticSolver::position_solve") << "prev / new energy = " << m_initEnergy << " / " << m_energy << "; new residual = " << m_l2norm << "; retaining step; converged! keeping same lambda = " << m_lambda << "\n";
        	m_keepUpdating = false;
        	keepSolution = true;
        }
        // we've reached a point of lower energy
        else if (m_energy <= m_initEnergy)
        {
        	// Decrease lambda (= increase trust region size = decrease regularization)
        	m_lambda = funnyclipvalue( m_lambdamin, m_lambda * m_geardown / m_successcount, m_lambdamax );
            TraceStream(g_log, "StaticSolver::position_solve") << "prev / new energy = " << m_initEnergy << " / " << m_energy << "; new residual = " << m_l2norm << "; retaining step; growing trust region: new lambda = " << m_lambda << "\n";
            keepSolution = true;
        }
        else
        {	    
        	// Solver failed
        	/////////////////////////////////////////////////
        	// discard the step
        	TraceStream(g_log, "StaticSolver::position_solve") << "prev / new energy = " << m_initEnergy << " / " << m_energy << "; new residual = " << m_l2norm << "; discarding step \n";
        }


	#ifndef NDEBUG      // Ensure that fixed DOFs are at their desired values
        VecXd xf(m_ndof);
        m_diffEq.getX(xf);
        for (int i = 0; i < (int) m_fixed.size(); ++i)
        {
        	assert(approxEq(m_desired[i], xf(m_fixed[i]), 1.0e-6));
        }
	#endif

        if (keepSolution)
        {
        	m_failurecount = 0;
        	m_successcount++;
        }
        else
        {
        	m_failurecount++;
        	m_successcount = 1;
        }
        return keepSolution;
    }


    bool position_solve()
    {
        START_TIMER("StaticSolver::newton_step/setup");

        TraceStream(g_log, "StaticSolver::position_solve") << "call #" << ++static_solve_counter << "\n";

        //for (int vidx = 0; vidx < (int) m_diffEq.getRod()->nv(); vidx++)
         //{
        //	std::cout << m_diffEq.getRod()->getVertex(vidx) << std::endl;
         //}

        // Chapter 0: Basic housekeeping
        ////////////////////////////////////////////////////

        bool successful_solve = true;        


        resize();
        setZero();
        m_diffEq.getScriptedDofs(m_fixed, m_desired); // m_fixed are DOF indices, m_desired are corresponding desired values
        assert(m_fixed.size() == m_desired.size());

        m_diffEq.set_qdot( m_deltaX ); // used to visualize increments (as velocities), initially zero

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
        m_diffEq.getX(m_x0);
	  
        // Zero the initial guess
        m_deltaX.setZero();

        // Based on the initial guess, 
        //   1) set up right hand side (RHS = F).
        //   2) compute the potential energy
        //   3) compute the residual
        updatePositionBasedQuantities();

        // save the initial energy
        m_initEnergy = m_energy;

        STOP_TIMER("StaticSolver::position_solve/setup");

        // Chapter 2: Iterate using Newton's Method
        ////////////////////////////////////////////////////

        TraceStream(g_log, "StaticSolver::position_solve") <<
        		"Initial energy = " << m_initEnergy << " lambda = " << m_lambda << "\n";

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
            m_A->add(i, i, m_lambda);

            // Levenberg-Marquardt diagonal shift
            //Scalar d = (*m_A)(i,i);
            //TraceStream(g_log, "StaticSolver::position_solve") << "lambda = " << m_lambda << " d[" << i << "] = " << d << " (1.+m_lambda) * d = " << (1.+m_lambda) * d << "\n";
            //m_A->set(i, i, (1.+m_lambda) * d);
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

        START_TIMER("Staticsolver::position_solve/setup");
        // Allow the nonzero structure to be modified again (for sparse matrices only)
        m_A->resetNonzeros();
        STOP_TIMER("Staticsolver::position_solve/setup");

        if (status < 0)
        {
        	m_failurecount++;
        	// shrink trust region (increase regularization)
        	m_lambda = funnyclipvalue( m_lambdamin, m_lambda * m_gearup * m_failurecount, m_lambdamax );

            DebugStream(g_log, "StaticSolver::position_solve") << "\033[31;1mWARNING IN StaticSolver:\033[m Problem during linear solve detected. "
				   << " new lambda = " << m_lambda << "\n";

            m_diffEq.set_q(m_x0);
            m_diffEq.updateCachedQuantities();
	    
            return false;
        }

        // visualize pre-filtering velocities
        m_diffEq.set_qdot( m_deltaX ); // used to visualize increments (as velocities)

        bool done = examine_solution();

        if (!done)
        {
        	TraceStream(g_log, "StaticSolver::position_solve") << "filtering delta_X and trying again...\n";

        	filterDeltaX();

        	// visualize post-filtering velocities
        	m_diffEq.set_qdot( m_deltaX ); // used to visualize increments (as velocities)

        	done = examine_solution();

        	if (!done)
        	{
                TraceStream(g_log, "StaticSolver::position solve") << "before step positions : " << m_x0 << "\n";
                TraceStream(g_log, "StaticSolver::position solve") << "a step positions : " << m_x0+m_deltaX << "\n";
        		// Increase lambda (= decrease trust region size = increase regularization)
        		m_lambda = funnyclipvalue( m_lambdamin, m_lambda * m_gearup * m_failurecount, m_lambdamax );
        		TraceStream(g_log, "StaticSolver::position_solve") << "shrinking trust region: new lambda = " << m_lambda << "\n";
        		ElasticRod &r = *m_diffEq.getRod();
        		int   ia  = r.globalRodIndex;
        		InfoStream(g_log, "StaticSolver::position_solve") << "rod_idx " << ia << " failure count " << m_failurecount
        			<< " shrinking trust region: new lambda = " << m_lambda << "\n";

        		m_diffEq.set_q(m_x0);
        		m_diffEq.updateCachedQuantities();
        	}
        }

        return true;
    }


  void filterDeltaX() // this code assumes HAIR: the first two vertices are prescribed, the rest are free
  { 

    // ////////////////////// EG: I'm not sure that this code helps, so let's skip the strain rate filter until it's better evaluated ///////////////
    // WarningStream(g_log, "SymmetricImplicitEuler::filterDeltaX/skipped", MsgInfo::kOncePerId) << "WARNING: Skipping the strain rate limiting in Newton solver line search!\n";
    // return;
    // return;
    // return;
    // ////////////////////// EG: I'm not sure that this code helps, so let's skip the strain rate filter until it's better evaluated ///////////////


    // using l_ to denote local variables, to avoid confusion with similarly-named member variables m_

    TraceStream(g_log, "StaticSolver::filterDeltaX") << "before filtering: " << m_deltaX << "\n";

    assert(m_diffEq.getRod());
    ElasticRod &r = *m_diffEq.getRod();

    int   ia  = r.vertIdx(1,0);
    Vec3d xaP = m_x0.segment<3>     (ia);  // start-of-step position of second vertex
    Vec3d xaD = m_deltaX.segment<3> (ia);
    Vec3d xaN = xaP + xaD;                 // end-of-step position of second vertex

    Vec3d disp(0,0,0);
    
    for (int i = 2; i < r.nv(); ++i)
    {
      int   ib  = r.vertIdx(i,0);
      Vec3d xbP = m_x0.segment<3>     (ib);
      Vec3d xbD = m_deltaX.segment<3> (ib) + disp;
      Vec3d xbN = xbP + xbD;                 

      double lP = (xbP - xaP).norm();
      double lN = (xbN - xaN).norm();

      double lNrev = lN;

      double maxStrainRate = 0.001;

      lNrev = lP;

      // if (lNrev > (1. + maxStrainRate) * lP) 
      // {
      // 	lNrev = (1. + maxStrainRate) * lP;
      // }
      // else if (lNrev < (1. - maxStrainRate) * lP) 
      // {
      // 	lNrev = (1. - maxStrainRate) * lP;
      // }

      // compute and store revised delta
      Vec3d xbNrev = xaN + (xbN - xaN) * lNrev / lN;
      Vec3d xbDrev = xbNrev - xbP;
      m_deltaX.segment<3> (ib) = xbDrev;

      disp += xbNrev - xbN;

      //double lNrevised = (x1Nrevised - x0N).norm();

      //Vec3d v1revised = (x1Nrevised - x1) / m_dt;

      ia  = ib;
      xaP = xbP;
      xaN = xbNrev;
    }

    TraceStream(g_log, "StaticSolver::filterDeltaX") << "after filtering: " << m_deltaX << "\n";

    // TraceStream(g_log, "") << "Edge lengths after inextensibility: ";
    // for (int i = 0; i < m_rods[rodidx]->nv() - 1; ++i)
    // {
    //     Vec3d x0 = m_xn.segment<3> (rodbase + 3 * i);
    //     Vec3d v0 = m_vnphalf.segment<3> (rodbase + 3 * i);
    //     Vec3d x1 = m_xn.segment<3> (rodbase + 3 * i + 3);
    //     Vec3d v1 = m_vnphalf.segment<3> (rodbase + 3 * i + 3);
    //  Vec3d x0N = x0 + m_dt * v0;
    //  Vec3d x1N = x1 + m_dt * v1;
    //  Vec3d eN = (x1N - x0N);
    //  TraceStream(g_log, "") << " " << eN.norm();
    // }
    // TraceStream(g_log, "") << '\n';

    // DebugStream(g_log, "") << "Velocity Filter end: rod " << rodidx << " vertex 0 = " << m_xn.segment<3>(rodbase) << '\n';
  }



    ODE& m_diffEq;

    int m_ndof;

    VecXd m_x0;
    VecXd m_rhs;
    VecXd m_deltaX;

    Scalar m_lambda;
    Scalar m_lambdamin, m_lambdamax;
    Scalar m_gearup, m_geardown;
    Scalar m_failurecount, m_successcount;
    bool m_keepUpdating;

    IntArray m_fixed;
    std::vector<Scalar> m_desired;

    MatrixBase* m_A;
    LinearSolverBase* m_solver;

    Scalar m_energy; // EG: move to DiffEqSolver
    Scalar m_l2norm;

    Scalar m_initEnergy;
};

} // namespace BASim

#endif // STATICSOLVER_HH






