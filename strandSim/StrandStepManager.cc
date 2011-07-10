/*
 * StrandStepManager.cc
 *
 *  Created on: 8/07/2011
 *      Author: jaubry
 */

#include "StrandStepManager.hh"

namespace strandsim
{

template<typename StepperT>
StrandStepManager<StepperT>::StrandStepManager(const std::vector<ElasticStrand*>& strands) :
    m_strands(strands)
{
    // TODO Auto-generated constructor stub

}

template<typename StepperT>
StrandStepManager<StepperT>::~StrandStepManager()
{
    // TODO Auto-generated destructor stub
}

}

// Explicit template instantiations

#include "ElasticStrand.hh"
#include "ElasticStrandStaticStepper.hh"

namespace strandsim
{

template class StrandStepManager<ElasticStrandStaticStepper> ;

}
