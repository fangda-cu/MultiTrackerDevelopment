#include <cmath>

#include "UpTransfer.hh"

namespace BASim {

Scalar UpTransfer::computeFunction(Scalar& x)
{
	return tan(x / 2.0);
	//return (pow(x, 6.0) / pow(M_PI, 5.0));
}



Scalar UpTransfer::computeFirstDerivative(Scalar & x)
{
	return (1.0 / (2.0 * pow(cos(x / 2.0), 2.0)));
	//return (6.0 * pow(x, 5.0) / pow(M_PI, 5.0));
}



Scalar UpTransfer::computeSecondDerivative(Scalar & x)
{
	return (sin(x / 2.0) / (2.0 * pow(cos(x / 2.0), 3.0)));
	//return (30.0 * pow(x, 4.0) / pow(M_PI, 5.0));
}

}