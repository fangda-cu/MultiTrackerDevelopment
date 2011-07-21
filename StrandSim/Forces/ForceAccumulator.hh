/*
 * ForceAccumulator.hh
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#ifndef FORCEACCUMULATOR_HH_
#define FORCEACCUMULATOR_HH_

#include "../Definitions.hh"
#include "../BandMatrix.hh"

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

    static void accumulate(Scalar& globalEnergy, VecXd& globalForce, const ElasticStrand& strand,
            const StrandGeometry& geometry);

    static void accumulate(Scalar& globalEnergy, VecXd& globalForce, strandsim::BandMatrix<Scalar, 10, 10>& globalJacobian,
            const ElasticStrand& strand, const StrandGeometry& geometry);

    static void accumulate(strandsim::BandMatrix<Scalar, 10, 10>& globalJacobian,
            const ElasticStrand& strand, const StrandGeometry& geometry);

};

}

#endif /* FORCEACCUMULATOR_HH_ */
