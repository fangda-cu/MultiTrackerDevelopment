/*
 * StrandStepManager.hh
 *
 *  Created on: 8/07/2011
 *      Author: jaubry
 */

#ifndef STRANDSTEPMANAGER_HH_
#define STRANDSTEPMANAGER_HH_

#include <vector>
#include "ElasticStrand.hh"

namespace strandsim
{

template<typename StepperT>
class StrandStepManager
{
public:
    explicit StrandStepManager(const std::vector<ElasticStrand*>& strands);

    virtual ~StrandStepManager();

    void execute()
    {
        StepperT stepper;

        for (std::vector<ElasticStrand*>::const_iterator s = m_strands.begin(); s != m_strands.end(); ++s)
            stepper.execute(**s);
    }

private:
    const std::vector<ElasticStrand*>& m_strands;
};

}

#endif /* STRANDSTEPMANAGER_HH_ */
