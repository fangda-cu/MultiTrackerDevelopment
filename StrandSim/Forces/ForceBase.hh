/*
 * ForceBase.hh
 *
 *  Created on: 12/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef FORCE_HH_
#define FORCE_HH_

#include "../Definitions.hh"
#include "../BandMatrix.hh"

namespace strandsim
{

class ElasticStrand;
class StrandGeometry;

/**
 * This is the interface that any external force should expose. Each accumulate* method accumulates energy, force and/or Jacobian in geometry
 */
class ForceBase
{
public:
    ForceBase();
    virtual ~ForceBase();

    virtual std::string getName() = 0;

    virtual void accumulateEFJ( StrandGeometry& geometry, const ElasticStrand& strand ) const = 0;

    virtual void accumulateEF( StrandGeometry& geometry, const ElasticStrand& strand ) const = 0;

    virtual void accumulateJ( StrandGeometry& geometry, const ElasticStrand& strand ) const = 0;

};

}

#endif /* FORCE_HH_ */
