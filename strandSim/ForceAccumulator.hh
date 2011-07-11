/*
 * ForceAccumulator.hh
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#ifndef FORCEACCUMULATOR_HH_
#define FORCEACCUMULATOR_HH_

namespace strandsim
{

class ElasticStrand;

template<typename ForceT>
class ForceAccumulator
{
public:
    ForceAccumulator();
    virtual ~ForceAccumulator();

    static void accumulate(ElasticStrand& strand);
};

}

#endif /* FORCEACCUMULATOR_HH_ */
