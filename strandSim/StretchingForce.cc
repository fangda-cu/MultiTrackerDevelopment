/*
 * StretchingForce.cc
 *
 *  Created on: 12/07/2011
 *      Author: jaubry
 */

#include "StretchingForce.hh"
#include "ElasticStrand.hh"

namespace strandsim
{

StretchingForce::StretchingForce()
{
    // TODO Auto-generated constructor stub

}

StretchingForce::~StretchingForce()
{
    // TODO Auto-generated destructor stub
}

Scalar StretchingForce::localEnergy(const ElasticStrand& strand, const IndexType vtx)
{
    return 0;
}

}
