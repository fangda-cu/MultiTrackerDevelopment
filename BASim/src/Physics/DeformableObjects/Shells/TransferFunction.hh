#ifndef _TRANSFER_FUNCTION_H_
#define _TRANSFER_FUNCTION_H_

#include "BASim/src/Core/Definitions.hh"

namespace BASim {

class TransferFunction
{
public:
  enum TYPE {
    LINEAR
  };

  TransferFunction() {}
  virtual ~TransferFunction() {}

  virtual Scalar computeFunction(Scalar& x) = 0;
  virtual Scalar computeFirstDerivative(Scalar& x) = 0;
  virtual Scalar computeSecondDerivative(Scalar& x) = 0;
};

}

#endif // _TRANSFER_FUNCTION_H_
