#ifndef _LINEAR_TRANSFER_H_
#define _LINEAR_TRANSFER_H_

#include "TransferFunction.hh"

namespace BASim {

class LinearTransfer : public TransferFunction
{
public:
  LinearTransfer() {}
  virtual ~LinearTransfer() {}

  virtual Scalar computeFunction(Scalar& x);
  virtual Scalar computeFirstDerivative(Scalar& x);
  virtual Scalar computeSecondDerivative(Scalar& x);
};

}

#endif // _LINEAR_TRANSFER_H_
