/*
 * StrandBase.hh
 *
 *  Created on: 7/07/2011
 *      Author: Jean-Marie Aubry <jaubry@wetafx.co.nz>
 */

#ifndef STRANDBASE_H_
#define STRANDBASE_H_

#include "Definitions.hh"

namespace strandsim
{

class StrandBaseParameters
{

};

class StrandBase
{
    typedef StrandBaseParameters ParametersType;

public:
    StrandBase();
    virtual ~StrandBase();
};

}

#endif /* STRANDBASE_H_ */
