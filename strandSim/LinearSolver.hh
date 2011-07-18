/*
 * LinearSolver.hh
 *
 *  Created on: 18/07/2011
 *      Author: jaubry
 */

#ifndef LINEARSOLVER_HH_
#define LINEARSOLVER_HH_

#include "Definitions.hh"

namespace strandsim
{

template<typename MatrixT>
class LinearSolver
{
public:
    LinearSolver();
    virtual ~LinearSolver();

    int solve(VecXd& newDOFs, const MatrixT& J, const VecXd& F) const;
};

}

#endif /* LINEARSOLVER_HH_ */
