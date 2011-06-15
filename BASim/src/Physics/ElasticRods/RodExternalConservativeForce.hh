/**
 * \file RodExternalConservativeForce.hh
 *
 * \author eitan@cs.columbia.edu
 * \date 06/06/2011
 */

#include "RodExternalForce.hh"

#ifndef RODEXTERNALCONSERVATIVEFORCE_HH
#define RODEXTERNALCONSERVATIVEFORCE_HH

namespace BASim
{

class ElasticRod;
class MatrixBase;

/** Base class for external forces applied to a rod. */
class RodExternalConservativeForce : public RodExternalForce
{
public:

    explicit RodExternalConservativeForce(bool implicit = true)
    {
      m_conservative = true;
    }

    virtual ~RodExternalConservativeForce()
    {
    }

    virtual Scalar computeEnergy(const ElasticRod& rod) const = 0;
 
    virtual void computeForceEnergy(const ElasticRod& rod, VecXd& force, Scalar& energy) const = 0;
};

} // namespace BASim

#endif // RODEXTERNALCONSERVATIVEFORCE_HH
