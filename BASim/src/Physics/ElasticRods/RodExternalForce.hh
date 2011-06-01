/**
 * \file RodExternalForce.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/04/2009
 */

#ifndef RODEXTERNALFORCE_HH
#define RODEXTERNALFORCE_HH

namespace BASim
{

class ElasticRod;
class MatrixBase;

/** Base class for external forces applied to a rod. */
class RodExternalForce
{
public:

    explicit RodExternalForce(bool implicit = true) :
        m_implicit(implicit), m_name("name not given!")
    {
    }

    virtual ~RodExternalForce()
    {
    }

    virtual void computeForce(const ElasticRod& rod, VecXd& force) const = 0;

    virtual void computeForceDX(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) const = 0;

    virtual void computeForceDV(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J) const = 0;

    bool isImplicit() const
    {
        return m_implicit;
    }

    void setImplicit(bool implicit)
    {
        m_implicit = implicit;
    }

    std::string getName() const
    {
        return m_name;
    }

protected:
    bool m_implicit;

    std::string m_name;
};

} // namespace BASim

#endif // RODEXTERNALFORCE_HH
