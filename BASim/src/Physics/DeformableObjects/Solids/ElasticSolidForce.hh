/**
 * \file ElasticSolidForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date May 22, 2011
 */

#ifndef ELASTICSOLIDFORCE_H
#define ELASTICSOLIDFORCE_H

#include "BASim/src/Physics/DeformableObjects/Solids/ElasticSolid.hh"

namespace BASim {

/** Base class for a force that acts on shells. */
class ElasticSolidForce
{
public:
  
  ElasticSolidForce(ElasticSolid& solid, const std::string& name = "ElasticSolidForce") : m_solid(solid), m_name(name) {}
  virtual ~ElasticSolidForce() {}

  std::string getName() const { return m_name; }

  virtual Scalar globalEnergy() const = 0;
  virtual void globalForce(VecXd& force) const = 0;
  virtual void globalJacobian(Scalar scale, MatrixBase& Jacobian) const = 0;

  virtual void update() {};

protected:

  ElasticSolid& m_solid;
  std::string m_name;

};


} // namespace BASim

#endif // ELASTICSHELLFORCE_H
