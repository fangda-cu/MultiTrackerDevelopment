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
template<typename ODE>
class StaticSolver: public DiffEqSolver
{
public:

    explicit StaticSolver(ODE& ode);

    ~StaticSolver();

    bool execute();

    std::string getName() const
    {
        return "StaticSolver";
    }

    void resize();
    void setZero();
    // update computation of: forces (the RHS), residual, energy
    void updatePositionBasedQuantities();

    static int solveCounter;

protected:
    static Scalar funnyclipvalue(Scalar minvalue, Scalar variable, Scalar maxvalue);
    bool isConverged();
    bool examine_solution();
    bool position_solve();
    void filterDeltaX();

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
