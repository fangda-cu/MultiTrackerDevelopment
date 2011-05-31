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

        //std::cerr << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Newton solver failed to converge in max iterations: " << m_maxit << ". Attempting alternate initial iterate 2." << std::endl;
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

        //std::cerr << "                          Newton solver failed to converge in max iterations: " << m_maxit << ". Attempting alternate initial iterate 3." << std::endl;
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

        //std::cerr << "                          Newton solver failed to converge in max iterations: " << m_maxit << ". Attempting alternate initial iterate 4." << std::endl;
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

        //std::cerr << "                          Newton solver failed to converge in max iterations: " << m_maxit << ". Attempting alternate initial iterate 4." << std::endl;
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

      //  std::cerr << "\033[31;1mWARNING IN SYM IMPLICITEULER:\033[m Newton solver failed to converge in max iterations: "
      //          << m_maxit << "." << std::endl;
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

        // Return the L2 norm
        return m_rhs.norm();
    }

    bool isConverged()
    {
        /*
         std::cout << "atol " << m_residual
         << ", infnorm " << m_infnorm
         << ", rtol " << m_residual / m_initial_residual
         << ", stol " << m_increment.norm() << std::endl;
         */
        // L2 norm of the residual is less than tolerance
        if (m_residual < m_atol)
        {
            //std::cout << "converged atol" << std::endl;
            return true;
        }
        // Infinity norm of residual is less than tolerance
        if (m_infnorm < m_inftol)
        {
            //std::cout << "converged inftol" << std::endl;
            return true;
        }
        if (m_residual <= m_rtol * m_initial_residual)
        {
            //std::cout << "converged rtol" << std::endl;
            return true;
        }
        // L2 norm of change in solution at last step of solve is less than tolerance
        if (m_increment.norm() < m_stol)
        {
            //std::cout << "converged stol" << std::endl;
            return true;
        }
        return false;
    }

protected:

    void generateInitialIterate1(VecXd& dx)
    {
        dx = m_dt * v0;
    }

    void generateInitialIterate2(VecXd& dx)
    {
        dx.setZero();
    }

    void generateInitialIterate3(VecXd& dx)
    {
        dx = 0.5 * m_dt * v0;
    }

    void generateInitialIterate4(VecXd& dx)
    {
        dx = 1.5 * m_dt * v0;
    }

    void generateInitialIterate5(VecXd& dx)
    {
        dx = 2.0 * m_dt * v0;
    }

    bool position_solve(int guess_to_use)
    {
        START_TIMER("SymmetricImplicitEuler::position_solve/setup");

        bool successfull_solve = true;

        m_diffEq.startStep();

        resize();
        setZero();
        m_diffEq.getScriptedDofs(m_fixed, m_desired);
        assert(m_fixed.size() == m_desired.size());
#ifdef DEBUG
        for( int i = 0; i < (int) m_fixed.size(); ++i ) assert( m_fixed[i] >= 0 );
        for( int i = 0; i < (int) m_fixed.size(); ++i ) assert( m_fixed[i] < m_ndof );
#endif

#ifdef TIMING_ON
        m_rhs.setZero();
        m_diffEq.evaluatePDot(m_rhs);
        double timing_force_at_start = m_rhs.norm();
#endif

        // Copy masses.
        if (!m_mass_set) // Assuming masses are constant
        {
            m_diffEq.getMass(m_mass);
            m_mass_set = true;
        }

        // Copy start of step positions and velocities.
        m_diffEq.getX(x0);
        m_diffEq.getV(v0);

        // Set velocites for fixed degrees of freedom.
        //for( int i = 0; i < (int) m_fixed.size(); ++i )
        //{
        //  int dof = m_fixed[i];
        //  v0(dof) = (m_desired[i]-x0(dof))/m_dt;
        //}

        //#ifdef DEBUG
        //  for( int i = 0; i < (int) m_fixed.size(); ++i ) assert( approxEq(m_desired[i],x0(m_fixed[i])+v0[m_fixed[i]]*m_dt,1.0e-6) );
        //#endif

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
            std::cerr << "\033[31;1mERROR IN IMPLICITEULER:\033[m Invalid initial iterate requested, exiting." << std::endl;
	    STOP_TIMER("SymmetricImplicitEuler::position_solve/setup");
            exit(1);
            break;
        }
        }

        // Set deltaX and v0 for fixed degrees of freedom.
        for (int i = 0; i < (int) m_fixed.size(); ++i)
        {
            int dof = m_fixed[i];
            m_deltaX(dof) = m_desired[i] - x0(dof);
            v0(dof) = (m_desired[i] - x0(dof)) / m_dt;
        }

        // Update the differential equation with the current guess
        m_diffEq.set_qdot(m_deltaX / m_dt);
        m_diffEq.set_q(x0 + m_deltaX);

        m_diffEq.endIteration();

        // Calling computeResidual also sets m_rhs = M(m_dt*v_n-m_deltaX) + h^2*F.
        m_initial_residual = m_residual = computeResidual();
        /*
         std::cout << "atol " << m_residual
         << ", infnorm " << m_infnorm
         << ", rtol " << m_residual / m_initial_residual
         << ", stol " << m_increment.norm() << std::endl;
         */

#ifdef DEBUG
        for( int i = 0; i < (int) m_fixed.size(); ++i )
        {
            //if( !approxEq(m_desired[i],x0(m_fixed[i])+m_deltaX(m_fixed[i]),1.0e-6) )
            //{
            //  std::cout << "Residual: " << -m_desired[i]+x0(m_fixed[i])+m_deltaX(m_fixed[i]) << std::endl;
            //  std::cout << "Initial guess: " << guess_to_use << std::endl;
            //}
            assert( approxEq(m_desired[i],x0(m_fixed[i])+m_deltaX(m_fixed[i]),1.0e-6) );
        }
#endif

        STOP_TIMER("SymmetricImplicitEuler::position_solve/setup");

        int m_curit = 0;
        for (; m_curit < m_maxit; ++m_curit)
        {
	  std::cout << "\nSymmetricImplicitEuler::position_solve: iteration = " << m_curit << "\n\n" << std::endl; 

            // TODO: Assert m_A, increment are zero
            START_TIMER("SymmetricImplicitEuler::position_solve/setup");

            for (int i = 0; i < (int) m_fixed.size(); ++i)
                m_rhs(m_fixed[i]) = m_dt * v0(m_fixed[i]) - m_deltaX(m_fixed[i]);

            // m_A = -h^2*dF/dx
            m_diffEq.evaluatePDotDX(-m_dt * m_dt, *m_A);
            m_A->finalize();

            //std::cout << "m_A = -h^2*dF/dx" << std::endl;
            //m_A->print();

            // m_A = -h*dF/dv -h^2*dF/dx
            m_diffEq.evaluatePDotDV(-m_dt, *m_A);
            m_A->finalize();
            assert( m_A->isApproxSymmetric(1.0e-6) );

            //std::cout << "m_A = -h*dF/dv -h^2*dF/dx" << std::endl;
            //m_A->print();

            // m_A = M -h*dF/dv -h^2*dF/dx
            for (int i = 0; i < m_ndof; ++i)
                m_A->add(i, i, m_mass(i));
            m_A->finalize();
            assert( m_A->isApproxSymmetric(1.0e-6) );

            //std::cout << "m_A = -h*dF/dv -h^2*dF/dx" << std::endl;
            //m_A->print();

            // Set the rows and columns corresponding to fixed degrees of freedom to 0
            m_A->zeroRows(m_fixed, 1.0);
            m_A->finalize();

            //std::cout << "m_A rows cleared" << std::endl;
            //m_A->print();

            m_A->zeroCols(m_fixed, 1.0);
            m_A->finalize();

            //std::cout << "m_A cols cleared" << std::endl;
            //m_A->print();

            // Finalize the nonzero structure before the linear solve (for sparse matrices only)
            m_A->finalizeNonzeros();
            STOP_TIMER("SymmetricImplicitEuler::position_solve/setup");

            assert(m_A->isApproxSymmetric(1.0e-6));

            START_TIMER("SymmetricImplicitEuler::position_solve/solver");
            int status = m_solver->solve(m_increment, m_rhs);
            STOP_TIMER("SymmetricImplicitEuler::position_solve/solver");
            if (status < 0)
            {
                successfull_solve = false;
                std::cerr << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. " << std::endl;
                return successfull_solve;
            }

            // Save m_deltaX for line search purposes
            m_deltaX_save = m_deltaX;

            // Verify that we have the correct linearization
            if (0)
            {
                VecXd residual0(m_rhs);
                Scalar residual0_norm = residual0.norm();
                Scalar steps[7] = { 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1 };
                std::cout << "check:\n  x0 " << x0.norm() << ", deltaX " << m_deltaX.norm() << ", increment "
                        << m_increment.norm() << "\n";
                for (int ss = 0; ss < 7; ss++)
                {
                    Scalar s = steps[ss];
                    m_deltaX = m_deltaX_save + s * m_increment;
                    m_diffEq.set_qdot(m_deltaX / m_dt);
                    m_diffEq.set_q(x0 + m_deltaX);
                    m_diffEq.endIteration();
                    computeResidual();
                    VecXd change(m_rhs - residual0), predicted(-s * residual0);
                    Scalar error = (change - predicted).norm() / predicted.norm();
                    std::cout << "  step " << s << ", change " << change.norm() << ", predicted " << predicted.norm()
                            << ", error " << error << std::endl;
                }
            }

	    START_TIMER("SymmetricImplicitEuler::position_solve/ls");
            double alpha = 1;
            double previous_residual = m_residual;
            m_deltaX += m_increment;
            for (int i = 0;; i++)
            {
                m_diffEq.set_qdot(m_deltaX / m_dt);
                m_diffEq.set_q(x0 + m_deltaX);
                m_diffEq.endIteration();
                // Calling computeResidual also sets m_rhs = M(m_dt*v_n-m_deltaX) + h^2*F.
                m_residual = computeResidual();
                std::cout << "\nSymmetricImplicitEuler::position_solve: line search: i "<<i<<", alpha "<<alpha<<", previous "<<previous_residual<<", residual "<<m_residual<<std::endl;
                if (m_residual < previous_residual)
                    break;
                else if (i >= 20)
                {
                    std::cerr << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Line search failed" << std::endl;
                    successfull_solve = false;
                    break;
                }
                alpha *= .5;
                m_deltaX = m_deltaX_save + alpha * m_increment;
            }
	    STOP_TIMER("SymmetricImplicitEuler::position_solve/ls");

            if (!successfull_solve)
                break;
            
            if (m_curit == m_maxit - 1)
                break;

            // Check for convergence.
            if (isConverged())
                break;

	    START_TIMER("SymmetricImplicitEuler::position_solve/setup");

            m_increment.setZero();
            m_A->setZero();

            // Allow the nonzero structure to be modified again (for sparse matrices only)
            m_A->resetNonzeros();

	    STOP_TIMER("SymmetricImplicitEuler::position_solve/setup");
        }

        //std::cout << "Iterations: " << m_curit << std::endl;
        if (m_curit == m_maxit - 1)
        {
            successfull_solve = false;
            //std::cerr << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Newton solver failed to converge in max iterations: " << m_maxit << std::endl;
        }

#ifdef TIMING_ON_DISABLE // FIXME
        std::string itrstrng = toString(m_curit+1);
        while( itrstrng.size() < 3 ) itrstrng = "0" + itrstrng;
        itrstrng = "NEWTON_SOLVES_OF_LENGTH_" + itrstrng;
        IntStatTracker::getIntTracker(itrstrng) += 1;

        PairVectorBase::insertPair( "ForcesVsIterations", std::pair<double,int>(timing_force_at_start,(m_curit+1)/m_dt), PairVectorBase::StringPair("ForceNorm","ImplicitIterationsPerSec") );

        ObjPropHandle<double> ophndl;
        if( m_diffEq.getRod()->property_exists(ophndl,"collision_induced_strain") )
        {
            m_diffEq.getRod()->property_handle(ophndl,"collision_induced_strain");
            double strain = m_diffEq.getRod()->property(ophndl);
            PairVectorBase::insertPair( "StrainVsIterations", std::pair<double,int>(strain,(m_curit+1)/m_dt), PairVectorBase::StringPair("Strain","ImplicitIterationsPerSec") );
        }
        //m_rods[i]->add_property(ophndl,"collision_induced_strain");
        //m_rods[i]->property(ophndl) = total_strain;

        if( successfull_solve )
        {
            ObjPropHandle<std::pair<double,int> > ophndl2;
            if( m_diffEq.getRod()->property_exists(ophndl2,"collision_induced_force_change") )
            {
                m_diffEq.getRod()->property_handle(ophndl2,"collision_induced_force_change");
                m_diffEq.getRod()->property(ophndl2).second += m_curit+1;
                //std::cout << m_diffEq.getRod()->property(ophndl2).first << "  -  " << m_diffEq.getRod()->property(ophndl2).second << std::endl;
            }
        }
#endif

        m_diffEq.endStep();

#ifdef DEBUG
        if( successfull_solve )
        {
            VecXd xf(m_ndof);
            m_diffEq.getX(xf);
            for( int i = 0; i < (int) m_fixed.size(); ++i ) assert( approxEq(m_desired[i],xf(m_fixed[i]),1.0e-6) );
        }
#endif

        //std::cout << m_curit << std::endl;
        return successfull_solve;
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
