/*
 * StaticSolver.cc
 *
 *  Created on: 6/06/2011
 *      Author: jaubry
 */

#include "../Core/Definitions.hh"
#include "StaticSolver.hh"
#include "../Physics/ElasticRods/GroomingTimeStepper.hh"

namespace BASim {

// Static variable definition
template<class ODE> int StaticSolver<ODE>::solveCounter = 0;

// Explicit template instantiation.
template class StaticSolver<GroomingTimeStepper>;

}
