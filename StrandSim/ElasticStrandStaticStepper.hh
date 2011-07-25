/*
 * ElasticStrandStaticStepper.hh
 *
 *  Created on: 8/07/2011
 *      Author: jaubry
 */

#ifndef ELASTICSTRANDSTATICSTEPPER_HH_
#define ELASTICSTRANDSTATICSTEPPER_HH_

#include "StepperBase.hh"
#include "ElasticStrand.hh"
#include "LinearSolver.hh"

namespace strandsim
{

class ElasticStrandStaticStepper: public StepperBase
{
public:
    ElasticStrandStaticStepper();
    virtual ~ElasticStrandStaticStepper();

    void execute( ElasticStrand& strand );

private:
    static const int m_maxlsit = 5;
    static const Scalar m_lambdamin = 1.0e-8;
    static const Scalar m_lambdamax = 1.0e+10;
    static const Scalar m_gearup = 2.00; // above 1.0
    static const Scalar m_geardown = 0.80; // below 1.0

    int m_successcount;
    int m_failurecount;
    Scalar m_lambda;
    Scalar m_previousLambda;

    BandMatrixLinearSolver<10, 10> linearSolver;
};

}

#endif /* ELASTICSTRANDSTATICSTEPPER_HH_ */
