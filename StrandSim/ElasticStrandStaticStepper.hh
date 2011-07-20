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

    void execute(ElasticStrand& strand) const;
};

}

#endif /* ELASTICSTRANDSTATICSTEPPER_HH_ */
