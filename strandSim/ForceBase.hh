/*
 * ForceBase.hh
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#ifndef FORCE_HH_
#define FORCE_HH_

#include "Definitions.hh"

namespace strandsim
{

template<typename StrandT>
class ForceBase
{
public:
    ForceBase();
    virtual ~ForceBase();
};

}

#endif /* FORCE_HH_ */
