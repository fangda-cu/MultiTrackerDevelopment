/*
 * SymmetricImplicitEuler.cc
 *
 *  Created on: 6/06/2011
 *      Author: jaubry
 */

typedef double Scalar;
#include <string>
#include "SymmetricImplicitEuler.hh"
#include "../Physics/ElasticRods/MultipleRodTimeStepper.hh"
#include "../Physics/ElasticRods/RodTimeStepper.hh"
#include "../Physics/DeformableObjects/DefoObjTimeStepper.hh"
//#include "../Math/DummyMatrix.hh"

namespace BASim
{

template<class ODE>
bool SymmetricImplicitEuler<ODE>::execute()
{
    START_TIMER("SymmetricImplicitEuler::execute");

    START_TIMER("SymmetricImplicitEuler::execute/backup");
    
    m_diffEq.backupResize();
    m_diffEq.backup();
    STOP_TIMER("SymmetricImplicitEuler::execute/backup");
    for (int guess = 0; guess <= 5; ++guess)
    {
        if (position_solve(guess))
        {
            START_TIMER("SymmetricImplicitEuler::execute/backup");
            m_diffEq.backupClear();
            STOP_TIMER("SymmetricImplicitEuler::execute/backup");

            STOP_TIMER("SymmetricImplicitEuler::execute");
            return true;
        }
        START_TIMER("SymmetricImplicitEuler::execute/backup");
        m_diffEq.backupRestore();
        STOP_TIMER("SymmetricImplicitEuler::execute/backup");
    }
    
    
    // Falling back to rigid motion if solver failed, so at least the rod doesn't mess up collision detection for this time step
    generateInitialIterate0(m_deltaX);
    // Force the scripted vertices to their place. NB: if the first two vertices are inextensibly scripted this shouldn't be necessary.
    for (unsigned int i = 0; i < m_fixed.size(); ++i)
    {
        const int dof = m_fixed[i];
        m_deltaX(dof) = m_desired[i] - x0(dof); // set desired position
    }
    m_diffEq.set_q(x0 + m_deltaX);

    START_TIMER("SymmetricImplicitEuler::execute/backup");
    m_diffEq.backupClear();
    STOP_TIMER("SymmetricImplicitEuler::execute/backup");

    STOP_TIMER("SymmetricImplicitEuler::execute");

    return false;
}

template<class ODE>
bool SymmetricImplicitEuler<ODE>::position_solve(int guess_to_use)
{
    START_TIMER("SymmetricImplicitEuler::position_solve/setup");

    // Chapter 0: Basic housekeeping
    ////////////////////////////////////////////////////
    
    bool successful_solve = true;

    m_diffEq.startStep();

    resize();
    setZero();

    //clear the list of constraints to start!
    m_fixed.clear();
    m_desired.clear();
    m_diffEq.getScriptedDofs(m_fixed, m_desired); // m_fixed are DOF indices, m_desired are corresponding desired values
    assert(m_fixed.size() == m_desired.size());
    
#ifndef NDEBUG // ensure indices in valid range
    for (int i = 0; i < (int) m_fixed.size(); ++i)
        assert(m_fixed[i] >= 0);
    for (int i = 0; i < (int) m_fixed.size(); ++i)
        assert(m_fixed[i] < m_ndof);
#endif

    // TODO: EG: Is this copying of data REALLY NEEDED? Probably not.
    // Copy masses.
    
    //if (!m_mass_set) // Assuming masses are constant
    {
      //printf("Actual copying masses:\n");
        m_diffEq.getMass(m_mass);
        //for(int i = 0; i < m_mass.size(); ++i) {
        //   printf("Mass: %e\n", m_mass(i));
        //}
        //m_mass_set = true;
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
        if (!generateInitialIterate0(m_deltaX))
            return false;
        break;
    }
    case 1:
    {
        generateInitialIterate1(m_deltaX);
        break;
    }
    case 2:
    {
        generateInitialIterate2(m_deltaX);
        break;
    }
    case 3:
    {
        generateInitialIterate3(m_deltaX);
        break;
    }
    case 4:
    {
        generateInitialIterate4(m_deltaX);
        break;
    }
    case 5:
    {
        generateInitialIterate5(m_deltaX);
        break;
    }
    default:
    {
        ErrorStream(g_log, "") << "\033[31;1mERROR IN IMPLICITEULER:\033[m Invalid initial iterate requested, exiting." << '\n';
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
        m_deltaX(dof) = m_desired[i] - x0(dof); // set desired position
        v0(dof) = (m_desired[i] - x0(dof)) / m_dt; // reverse engineer desired velocity
    }

#ifndef NDEBUG // check that the desired velocity was correctly computed by taking virtual forward step
    for (int i = 0; i < (int) m_fixed.size(); ++i)
        assert(approxEq(m_desired[i], /* =approx= */x0(m_fixed[i]) + v0[m_fixed[i]] * m_dt, 1.0e-6));
#endif

    
    // Update the differential equation with the current guess
    m_diffEq.set_qdot(m_deltaX / m_dt); // set velocity
    m_diffEq.set_q(x0 + m_deltaX); // set position

   
    // Signal the differential equation that it should get
    // ready for the first iteration. In practice this is where the ODE
    // can precompute some reusable quantities that depend on the state (position & velocity)
    m_diffEq.endIteration();

    // Based on the initial guess,
    //   1) set up right hand side (RHS) of implicit Euler: RHS = M(m_dt*v_n-m_deltaX) + h^2*F.
    //   2) compute the residual of the ODE
    //   3) cache the residual as m_initial_residual, for convergence test later
    m_initial_residual = m_residual = computeResidual();


    TraceStream(g_log, "") << "SymmetricImplicitEuler::position_solve: starting Newton solver. Initial guess has residual = "
            << m_residual << ", convergence test will use thresholds atol = " << m_atol << " inftol = " << m_inftol
            << " rtol = " << m_initial_residual * m_rtol << " stol = " << m_stol << '\n';

    STOP_TIMER("SymmetricImplicitEuler::position_solve/setup");

    // Chapter 2: Iterate using Newton's Method
    ////////////////////////////////////////////////////
    
    int curit = 0;
    for (curit = 0; curit < m_maxit; ++curit)
    {

        TraceStream(g_log, "") << "\nSymmetricImplicitEuler::position_solve: Newton iteration = " << curit << "\n";

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
        m_diffEq.evaluatePDotDX(-m_dt * m_dt, *m_A); // NB m_A is set to zero at construction time and at the end of this loop.
        m_A->finalize();

        // Consider LHS arising from dissipative forces (function of velocity)
        // m_A = -h*dF/dv -h^2*dF/dx
        m_diffEq.evaluatePDotDV(-m_dt, *m_A);
        m_A->finalize();

        // Consider inertial contribution from mass matrix
        // m_A = M -h*dF/dv -h^2*dF/dx
        for (int i = 0; i < m_ndof; ++i) {
            m_A->add(i, i, m_mass(i));
        }
        m_A->finalize();

        // Set the rows and columns corresponding to fixed degrees of freedom to 0
        m_A->zeroRows(m_fixed, 1.0);
        m_A->finalize();

        m_A->zeroCols(m_fixed, 1.0);
        m_A->finalize();

        // Finalize the nonzero structure before the linear solve (for sparse matrices only)
        m_A->finalizeNonzeros();
        assert(isSymmetric(*m_A));
        STOP_TIMER("SymmetricImplicitEuler::position_solve/setup");

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
          std::cout << "***Solver failed.***\n";
            DebugStream(g_log, "") << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Problem during linear solve detected. "
                    << '\n';
            return false;
        }
        START_TIMER("SymmetricImplicitEuler::position_solve/ls");

        m_alpha = 1.; // actual step will m_deltaX += alpha * m_increment

        // Save m_deltaX and residual for later
        m_deltaX_save = m_deltaX;
        double previous_residual = m_residual;

        // Attempt a full Newton step (alpha = 1)
        m_deltaX = m_deltaX_save + m_alpha * m_increment;

        for (int i = 0;; i++)
        {
            // Evaluate residual for attempted increment
            //////////////////////////////////////////////

            // Update the differential equation with the current guess
            m_diffEq.set_qdot(m_deltaX / m_dt);
            m_diffEq.set_q(x0 + m_deltaX);

            // Signal the differential equation that it should recompute cached quantities
            m_diffEq.endIteration();

            // Calling computeResidual also sets m_rhs = M(m_dt*v_n-m_deltaX) + h^2*F.
            m_residual = computeResidual();

            bool converged = isConverged();
            TraceStream(g_log, "") << "SymmetricImplicitEuler::position_solve: summary of line search i " << i
                    << ": increment " << m_increment.norm() << " previous " << previous_residual << ", residual " << m_residual
                    << '\n';

            // Is this residual (hence the increment) acceptable?
            /////////////////////////////////////////////////////////////

            if (m_residual < .9 * previous_residual || converged)
            {
                TraceStream(g_log, "") << "Succeeded (done).\n\n";
                break;
            }
            else if (i >= m_maxlsit)
            {
               //std::cout << "Exceeded max iterations.\n\n";
                TraceStream(g_log, "")
                        << "Exceeded max iterations.\nSymmetricImplicitEuler::position_solve/line search: \033[31;1mWARNING IN IMPLICITEULER:\033[m Line search failed. Proceeding anyway.\n\n";
                break;
            }
            else
            {
                TraceStream(g_log, "") << "cutting increment and iterating.\n\n";
            }

            // Attempt a smaller step
            //////////////////////////////

            m_alpha *= .5;
            m_deltaX = m_deltaX_save + m_alpha * m_increment;
        }

        STOP_TIMER("SymmetricImplicitEuler::position_solve/ls");

        // After the line search...
        ///////////////////////////////

        // Check for convergence.
        if (isConverged())
            break;

        // Check for exceeding limit on number of Newton iterations
        if (curit == m_maxit - 1)
        {
            TraceStream(g_log, "") << "\033[31;1mWARNING IN IMPLICITEULER:\033[m Newton solver reached max iterations: "
                    << m_maxit << " with initial guess " << guess_to_use << '\n';
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

    //std::cout << "SymmetricImplicitEuler solve completed.\n";
    TraceStream(g_log, "") << "SymmetricImplicitEuler::position_solve: completed " << curit + 1 << " Newton iterations."
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

inline static Vec3d RigidMotion(const Vec3d& x, const Vec3d& p0, const Eigen::AngleAxis<double>& rotation, const Vec3d& w0,
        const double dt)
{
    return dt * w0 + p0 - x + rotation._transformVector(x - p0);
}

template<class ODE>
void SymmetricImplicitEuler<ODE>::setZero()
{
   x0.setZero();
   v0.setZero();
   m_rhs.setZero();
   m_deltaX.setZero();
   m_increment.setZero();
   m_A->setZero();
   m_A->resetNonzeros();
}

template<class ODE>
void SymmetricImplicitEuler<ODE>::resize()
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

// Initial guess based on rigid motion of the first two vertices, assuming their distance remains constant.
template<>
bool SymmetricImplicitEuler<RodTimeStepper>::generateInitialIterate0(VecXd& dx)
{
    const Vec3d p0 = x0.segment<3> (0);
    const Vec3d p1 = x0.segment<3> (4);
    const Vec3d w0 = v0.segment<3> (0);
    const Vec3d w1 = v0.segment<3> (4);
    const Vec3d q0 = p0 + m_dt * w0;
    const Vec3d q1 = p1 + m_dt * w1;

    const double cosAngle = (p1 - p0).dot(q1 - q0) / ((p1 - p0).norm() * (q1 - q0).norm());
    const double angle = cosAngle >= 1.0 ? 0.0 : cosAngle <= -1.0 ? M_PI : acos(cosAngle);
    assert(!isnan(angle));

    Vec3d normal = (p1 - p0).cross(w1 - w0);
    const double normalNorm = normal.norm();
    if (normalNorm < std::numeric_limits<double>::epsilon())
        normal = Vec3d(1, 0, 0); // The rotation is either identity or central symmetry, any axis will do...
    else
        normal = normal / normalNorm;
    assert(approxEq(normal.norm(), 1.0));
    // DebugStream(g_log, "") << "Initial guess by rigid motion: p0 = " << p0 << " q0 = " << q0 << " w0 = " << w0 << " p1 = "
    //         << p1 << " q1 = " << q1 << " w1 = " << w1 << " normal = " << normal << " angle = " << angle << '\n';

    // Set initial vertex coordinates by rigid motion and initial twist angles to zero. TODO: surely there is a better way than explicitly computing the rotation matrix!
    const Eigen::AngleAxis<double> rotation(angle, normal);
    for (int i = 0; i < m_ndof - 3; i += 4)
    {
        dx.segment<3> (i) = RigidMotion(x0.segment<3> (i), p0, rotation, w0, m_dt);
        // DebugStream(g_log, "") << "Predicted vertex " << i / 4 << " = " << x0.segment<3> (i) + dx.segment<3> (i) << '\n';
        dx(i + 3) = 0;
    }
    dx.segment<3> (m_ndof - 3) = RigidMotion(x0.segment<3> (m_ndof - 3), p0, rotation, w0, m_dt);

    return true;
}

template<class ODE>
bool SymmetricImplicitEuler<ODE>::generateInitialIterate0(VecXd& dx)
{
    WarningStream(g_log, "")
            << "SymmetricImplicitEuler<ODE>::generateInitialIterate0 not implemented yet for this ODE, ignoring.\n";

    return false;
}

// Explicit template instantiations
template class SymmetricImplicitEuler<RodTimeStepper> ;
template class SymmetricImplicitEuler<MultipleRodTimeStepper> ;
template class SymmetricImplicitEuler<DefoObjTimeStepper> ;

}
