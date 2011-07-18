/*
 * ForceBase.cc
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#include "ForceBase.hh"
#include "ElasticStrand.hh"

namespace strandsim
{

template<typename StrandT>
ForceBase<StrandT>::ForceBase()
{
    // TODO Auto-generated constructor stub

}

template<typename StrandT>
ForceBase<StrandT>::~ForceBase()
{
    // TODO Auto-generated destructor stub
}

template class ForceBase<ElasticStrand>;

}
