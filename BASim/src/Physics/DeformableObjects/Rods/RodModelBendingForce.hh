/**
 * \file RodModelBendingForce.h
 *
 * \author fang@cs.columbia.edu
 * \date May 12, 2011
 */

#ifndef RODMODELBENDINGFORCE_HH
#define RODMODELBENDINGFORCE_HH

#include "BASim/src/Physics/DeformableObjects/Rods/RodModelForce.hh"

namespace BASim 
{
  class RodModelBendingForce : public RodModelForce
  {
  public:
    typedef Eigen::Matrix<Scalar, 11, 1> ElementForce;
    typedef Eigen::Matrix<Scalar, 11, 2> ElementBiForce;
    typedef Eigen::Matrix<Scalar, 11, 11> ElementJacobian;
    typedef std::pair<ElementJacobian, ElementJacobian> ElementBiJacobian;
    
    struct Stencil
    {
      VertexHandle v;
      EdgeHandle e1;
      EdgeHandle e2;
      IntArray dofindices;
    };
    
  public:
    void addStencil(Stencil & s) { m_stencils.push_back(s); }
    std::vector<Stencil> & stencils() { return m_stencils; }
    const std::vector<Stencil> & stencils() const { return m_stencils; }
    
  public:
    RodModelBendingForce(ElasticRodModel & rod, Scalar youngs_modulus, Scalar youngs_modulus_damping, Scalar timestep);
    virtual ~RodModelBendingForce();
    
  public:
    Scalar globalEnergy();
    void globalForce(VecXd & force);
    void globalJacobian(Scalar scale, MatrixBase & Jacobian);
    
  protected:
    Mat2d computeStiffness(Stencil & s, bool viscous);

    Scalar localEnergy(Stencil & s, bool viscous);
    void localForce(ElementForce & force, Stencil & s, bool viscous);
    void localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous);
    
  protected:
    std::vector<Stencil> m_stencils;
    
    Scalar m_youngs_modulus;
    Scalar m_youngs_modulus_damping;
    
    Scalar m_timestep;
    
  };
  
}


#endif // RODMODELBENDINGFORCE_HH
