/*
 * BandMatrixLinearSolver.hh
 *
 *  Created on: 18/07/2011
 *      Author: jaubry
 */

#ifndef LINEARSOLVER_HH_
#define LINEARSOLVER_HH_

#include "Definitions.hh"
#include "BandMatrix.hh"

namespace strandsim
{

template<int kl, int ku>
class BandMatrixLinearSolver
{
public:
    BandMatrixLinearSolver();
    virtual ~BandMatrixLinearSolver();

    int solve( VecXd& x, const BandMatrix<Scalar, kl, ku>& A, const VecXd& b );

private:
    std::vector<Scalar> m_ab;
    std::vector<int> m_ipiv;
};

}

#endif /* LINEARSOLVER_HH_ */
