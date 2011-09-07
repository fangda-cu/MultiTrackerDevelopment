#ifndef _UP_TRANSFER_H_
#define _UP_TRANSFER_H_

#include "TransferFunction.hh"

namespace BASim {

class UpTransfer : public TransferFunction
{
public:
  UpTransfer() {}
  virtual ~UpTransfer() {}

  virtual Scalar computeFunction(Scalar& x);
  virtual Scalar computeFirstDerivative(Scalar& x);
  virtual Scalar computeSecondDerivative(Scalar& x);
};

}

#endif // _UP_TRANSFER_H_
