#ifndef _DOWN_TRANSFER_H_
#define _DOWN_TRANSFER_H_

#include "TransferFunction.hh"

namespace BASim {

class DownTransfer : public TransferFunction
{
public:
  DownTransfer() {}
  virtual ~DownTransfer() {}

  virtual Scalar computeFunction(Scalar& x);
  virtual Scalar computeFirstDerivative(Scalar& x);
  virtual Scalar computeSecondDerivative(Scalar& x);
};

}

#endif // _DOWN_TRANSFER_H_
