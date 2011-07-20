/*
 * ForceAccumulator.cc
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#include "ForceAccumulator.hh"
#include "../ElasticStrand.hh"
#include "StretchingForce.hh"
#include "BendingForce.hh"
#include "TwistingForce.hh"
#include "GravitationForce.hh"

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
void ForceAccumulator<ForceT>::accumulate(Scalar& globalEnergy, VecXd& globalForce, const ElasticStrand& strand,
        const StrandGeometry& geometry)
{
    assert(geometry.m_framesUpToDate);

    for (IndexType vtx = ForceT::s_first; vtx < strand.m_numVertices - ForceT::s_last; ++vtx)
    {
        globalEnergy += ForceT::localEnergy(strand, geometry, vtx);
        ForceT::addInPosition(globalForce, vtx, ForceT::localForce(strand, geometry, vtx));
    }
}

template<typename ForceT>
void ForceAccumulator<ForceT>::accumulate(Scalar& globalEnergy, VecXd& globalForce,
        strandsim::BandMatrix<Scalar, 10, 10> globalJacobian, const ElasticStrand& strand, const StrandGeometry& geometry)
{
    assert(geometry.m_framesUpToDate);

    for (IndexType vtx = ForceT::s_first; vtx < strand.m_numVertices - ForceT::s_last; ++vtx)
    {
        globalEnergy += ForceT::localEnergy(strand, geometry, vtx);
        ForceT::addInPosition(globalForce, vtx, ForceT::localForce(strand, geometry, vtx));
        ForceT::addInPosition(globalJacobian, vtx, ForceT::localJacobian(strand, geometry, vtx));
    }
}

template<typename ForceT>
void ForceAccumulator<ForceT>::accumulate(strandsim::BandMatrix<Scalar, 10, 10>& globalJacobian, const ElasticStrand& strand,
        const StrandGeometry& geometry)
{
    assert(geometry.m_framesUpToDate);

    for (IndexType vtx = ForceT::s_first; vtx < strand.m_numVertices - ForceT::s_last; ++vtx)
    {
        ForceT::addInPosition(globalJacobian, vtx, ForceT::localJacobian(strand, geometry, vtx));
    }
}

// Explicit template instantiations
template class ForceAccumulator<StretchingForce> ;
template class ForceAccumulator<BendingForce> ;
template class ForceAccumulator<TwistingForce> ;
template class ForceAccumulator<GravitationForce> ;

}
