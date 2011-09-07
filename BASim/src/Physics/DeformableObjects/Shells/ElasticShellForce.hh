/**
 * \file ElasticShellForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date April 18, 2011
 */

#ifndef ELASTICSHELLFORCE_H
#define ELASTICSHELLFORCE_H

#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShell.hh"

namespace BASim {

/** Base class for a force that acts on shells. */
class ElasticShellForce
{
public:
  
  ElasticShellForce(ElasticShell& shell, const std::string& name = "ElasticShellForce") : m_shell(shell), m_name(name) {}
  virtual ~ElasticShellForce() {}

  std::string getName() const { return m_name; }

  virtual Scalar globalEnergy() const = 0;
  virtual void globalForce(VecXd& force) const = 0;
  virtual void globalJacobian(Scalar scale, MatrixBase& Jacobian) const = 0;

protected:

  ElasticShell& m_shell;
  std::string m_name;

};


} // namespace BASim

#endif // ELASTICSHELLFORCE_H
