/**
 * \file DrainingBubblePressureForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date Jan 4, 2012
 */

#ifndef DRAININGBUBBLEPRESSUREFORCE_H
#define DRAININGBUBBLEPRESSUREFORCE_H

#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShellForce.hh"

namespace BASim {

class DrainingBubblePressureForce : public ElasticShellForce {

public:

  DrainingBubblePressureForce  (ElasticShell& shell, const std::string& name = "DrainingBubblePressureForce", std::vector<EdgeHandle> holeEdges = std::vector<EdgeHandle>(), std::vector<EdgeHandle> baseEdges = std::vector<EdgeHandle>(), Scalar gas_density = 0, Scalar timestep = 1);
  virtual ~DrainingBubblePressureForce  () {}

  std::string getName() const;

  Scalar globalEnergy() const;
  void globalForce(VecXd& force) const;
  void globalJacobian(Scalar scale, MatrixBase& Jacobian) const;
  
  void update();

protected:
  
  std::vector<EdgeHandle> m_hole_edges; //list of edges representing the hole in the bubble
  std::vector<EdgeHandle> m_base_edges; //list of edges representing the base hole of the bubble

  Scalar m_old_volume;
  Scalar m_gas_density;
  Scalar m_timestep;
  Scalar m_current_pressure;

};




}


#endif //DRAININGBUBBLEPRESSUREFORCE_H
