/*
 * ElasticStrandStaticStepper.cc
 *
 *  Created on: 8/07/2011
 *      Author: jaubry
 */

#include "ElasticStrandStaticStepper.hh"

namespace strandsim
{

ElasticStrandStaticStepper::ElasticStrandStaticStepper()
{
    // TODO Auto-generated constructor stub

}

ElasticStrandStaticStepper::~ElasticStrandStaticStepper()
{
    // TODO Auto-generated destructor stub
}

void ElasticStrandStaticStepper::execute(ElasticStrand& strand) const
{
    assert(strand.readyForSolving());
    static const int numberOfFixedDOFs = 7;

    const Scalar E = strand.getTotalEnergy();
    ElasticStrand::ForceVectorType& F = strand.getTotalForces();
    ElasticStrand::JacobianMatrixType& J = strand.getTotalJacobian();
    VecXd& newDOFs = strand.getNewDegreesOfFreedom();

    // Enforce fixed first two vertices by setting the relevant portions of force and Jacobian to zero
    F.segment<numberOfFixedDOFs> (0).setZero();
    J.localStencilReplace<numberOfFixedDOFs> (0, Eigen::Matrix<Scalar, numberOfFixedDOFs, numberOfFixedDOFs>());

    // Add a constant diagonal to J to enforce a trust region.
    J.addConstantDiagonal(1.0);

    LinearSolver<ElasticStrand::JacobianMatrixType> linearSolver;
    linearSolver.solve(newDOFs, J, F); // X = J^{-1} F
    newDOFs += strand.getDegreesOfFreedom(); // X = X_0 + J^{-1} F
    // AND THE MASSES????

    // Update the new position's frames and forces
    strand.prepareForExamining();
    // Set the new force on the first two vertices to zero, for comparison purposes.
    strand.getNewTotalForces().segment<numberOfFixedDOFs> (0).setZero();

    // TODO: Examine here the new position.

    // Accept the new position. This also updates the strand for the next solve.
    strand.acceptNewPositions();
}


}
