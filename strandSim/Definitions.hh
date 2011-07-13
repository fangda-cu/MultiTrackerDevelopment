/*
 * Definitions.hh
 *
 *  Created on: 8/07/2011
 *      Author: jaubry
 */

#ifndef DEFINITIONS_HH_
#define DEFINITIONS_HH_

#include <stdint.h>
#include "../BASim/src/Core/Definitions.hh"

namespace strandsim
{
using namespace BASim;

typedef uint16_t IndexType;

typedef Eigen::Matrix<Scalar, 11, 1> Vec11d;
typedef Eigen::Matrix<Scalar, 11, 11> Mat11d;
typedef std::pair<Mat11d, Mat11d> Mat11dPair;

}
#endif /* DEFINITIONS_HH_ */
