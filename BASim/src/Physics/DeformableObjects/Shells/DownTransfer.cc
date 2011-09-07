#include <cmath>

#include "DownTransfer.hh"

namespace BASim {

Scalar DownTransfer::computeFunction(Scalar& x)
{
	return (pow(x, 1.0 / 6.0) * pow(M_PI, 5.0 / 6.0));
}



Scalar DownTransfer::computeFirstDerivative(Scalar& x)
{
	return (pow(M_PI, 5.0 / 6.0) / (6.0 * pow(x, 5.0 / 6.0)));
}



Scalar DownTransfer::computeSecondDerivative(Scalar& x)
{
	return (-5.0 * pow(M_PI, 5.0 / 6.0) / (36.0 * pow(x, 11.0 / 6.0)));
}

}