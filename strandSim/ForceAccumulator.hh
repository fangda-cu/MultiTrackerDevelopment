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
class StrandGeometry;

template<typename ForceT>
class ForceAccumulator
{
public:
    ForceAccumulator();
    virtual ~ForceAccumulator();

    // Compute global energy, force, Jacobian
    static void accumulate(ElasticStrand& strand, const StrandGeometry& geometry);
};

}

#endif /* FORCEACCUMULATOR_HH_ */
