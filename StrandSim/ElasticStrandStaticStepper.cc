/*
 * ElasticStrandStaticStepper.cc
 *
 *  Created on: 8/07/2011
 *      Author: jaubry
 */

#include "ElasticStrandStaticStepper.hh"

namespace strandsim
{

ElasticStrandStaticStepper::ElasticStrandStaticStepper() :
    m_lambda( 1.0e-3 ), m_previousLambda( 0.0 ), m_failurecount( 0 ), m_successcount( 1 )
{
}

ElasticStrandStaticStepper::~ElasticStrandStaticStepper()
{
}

inline Scalar clipValue( Scalar minvalue, Scalar variable, Scalar maxvalue )
{
    // funny wrap-around behavior ensures we don't get "stuck" at lambda=maxvalue
    if ( variable > maxvalue )
        variable = minvalue;

    return fmin( fmax( variable, minvalue ), maxvalue );
}

void ElasticStrandStaticStepper::execute( ElasticStrand& strand )
{
    assert( strand.readyForSolving() );

    static const int numberOfFixedDOFs = 7;

    const Scalar E = strand.getTotalEnergy();
    std::cout << "Energy before = " << E << '\n';

    VecXd& F = strand.getTotalForces();
    F.segment<numberOfFixedDOFs> ( 0 ).setZero(); // Enforce fixed DOFs
    std::cout << "Forces norm before = " << F.norm() << '\n';

    ElasticStrand::JacobianMatrixType& J = strand.getTotalJacobian();
    // Change the sign!!!
    J *= -1.0;
    // Add a constant diagonal to J to enforce a trust region.
    std::cout << "Regularising m_lambda = " << m_lambda << '\n';
    J.addConstantDiagonal( m_lambda - m_previousLambda );
    // Enforce fixed DOFs
    J.fixFirstDOFs<numberOfFixedDOFs> ();

    BandMatrixLinearSolver<10, 10> linearSolver;
    VecXd& newDOFs = strand.getNewDegreesOfFreedom();

    linearSolver.solve( newDOFs, J, F ); // X = J^{-1} F
    newDOFs += strand.getDegreesOfFreedom(); // X = X_0 + J^{-1} F
    std::cout << "Proposed new DOFS = " << newDOFs << '\n';

    // Update the new position's frames and forces
    strand.prepareForExamining();

    const Scalar newE = strand.getNewTotalEnergy();
    std::cout << "Energy after = " << newE << '\n';

    // Set the new force on the first two vertices to zero, for comparison purposes.
    strand.getNewTotalForces().segment<numberOfFixedDOFs> ( 0 ).setZero();

    std::cout << "Forces norm after = " << strand.getNewTotalForces().norm() << '\n';

    if ( newE <= E )
    {
        m_lambda = clipValue( m_lambdamin, m_lambda * m_geardown / m_successcount, m_lambdamax );
        std::cout << "Accepting position, m_lambda = " << m_lambda << '\n';

        // Accept the new position. This also updates the strand for the next solve.
        strand.acceptNewPositions();
        m_previousLambda = 0;

        m_failurecount = 0;
        m_successcount++;
    }
    else
    {
        std::cout << "Filtering...\n";
        strand.filterNewGeometryLength();
        std::cout << "Proposed new DOFS = " << newDOFs << '\n';
        strand.prepareForExamining();
        const Scalar newFilteredE = strand.getNewTotalEnergy();
        std::cout << "Energy after filtering = " << newFilteredE << '\n';
        // Set the new force on the first two vertices to zero, for comparison purposes.
        strand.getNewTotalForces().segment<numberOfFixedDOFs> ( 0 ).setZero();

        std::cout << "Forces norm after filtering = " << strand.getNewTotalForces().norm() << '\n';

        if ( newFilteredE <= E )
        {
            m_lambda = clipValue( m_lambdamin, m_lambda * m_geardown / m_successcount, m_lambdamax );
            std::cout << "Accepting position, m_lambda = " << m_lambda << '\n';

            // Accept the new position. This also updates the strand for the next solve.
            strand.acceptNewPositions();
            m_previousLambda = 0;

            m_failurecount = 0;
            m_successcount++;
        }
        else
        {
            m_failurecount++;
            m_successcount = 1;
            m_previousLambda = m_lambda; // Because we already added m_lambda to the Jacobian, that we want to keep we'll need only to add the increment
            m_lambda = clipValue( m_lambdamin, m_lambda * m_gearup * m_failurecount, m_lambdamax );
            std::cout << "Rejecting position, m_lambda = " << m_lambda << '\n';
        }
    }

}

}
