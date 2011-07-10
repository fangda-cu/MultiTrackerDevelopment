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

namespace strandsim
{

class ElasticStrandStaticStepper: public StepperBase
{
public:
    ElasticStrandStaticStepper();
    virtual ~ElasticStrandStaticStepper();

    void execute(ElasticStrand& strand) const;

    template<typename MatrixT>
    void solve(VecXd& newDOFs, const MatrixT J, const VecXd& F) const;
};

}

#endif /* ELASTICSTRANDSTATICSTEPPER_HH_ */
