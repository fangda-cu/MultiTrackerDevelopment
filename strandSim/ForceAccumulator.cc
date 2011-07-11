/*
 * ForceAccumulator.cc
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#include "ForceAccumulator.hh"
#include "ElasticStrand.hh"
#include "StretchingForce.hh"
#include "BendingForce.hh"
#include "TwistingForce.hh"

namespace strandsim
{

template<typename ForceT>
ForceAccumulator<ForceT>::ForceAccumulator()
{
    // TODO Auto-generated constructor stub

}

template<typename ForceT>
ForceAccumulator<ForceT>::~ForceAccumulator()
{
    // TODO Auto-generated destructor stub
}

template<typename ForceT>
void ForceAccumulator<ForceT>::accumulate(ElasticStrand& strand)
{
    for (IndexType vtx = 0; vtx < strand.m_numVertices; ++vtx)
    {
        strand.m_totalEnergy += ForceT::localEnergy(strand, vtx);
        ForceT::addInPosition(strand.m_totalForces, vtx, ForceT::localForce(strand, vtx));
        ForceT::addInPosition(strand.m_totalJacobian, vtx, ForceT::localJacobian(strand, vtx));
    }
}


// Explicit template instantiations
template class ForceAccumulator<StretchingForce>;
template class ForceAccumulator<BendingForce>;
template class ForceAccumulator<TwistingForce>;

}
