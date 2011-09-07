#include "LinearTransfer.hh"

namespace BASim {

Scalar LinearTransfer::computeFunction(Scalar & x)
{
  return x;
}


Scalar LinearTransfer::computeFirstDerivative(Scalar & x)
{
  return 1.0;
}


Scalar LinearTransfer::computeSecondDerivative(Scalar & x)
{
  return 0.0;
}

}