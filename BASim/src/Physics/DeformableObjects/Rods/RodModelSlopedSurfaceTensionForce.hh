/**
 * \file RodModelSlopedSurfaceTensionForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date Sept 18, 2012
 */

#ifndef RODMODELSLOPEDSURFACETENSIONFORCE_HH
#define RODMODELSLOPEDSURFACETENSIONFORCE_HH

#include "BASim/src/Physics/DeformableObjects/Rods/RodModelForce.hh"

namespace BASim 
{

  class RodModelSlopedSurfaceTensionForce : public RodModelForce
  {
  public:
    typedef Eigen::Matrix<Scalar, 9, 1> ElementForce;
    typedef Eigen::Matrix<Scalar, 9, 9> ElementJacobian;

    struct Stencil : public ElasticRodModel::JointStencil
    {
      Stencil(const ElasticRodModel::JointStencil & s) : ElasticRodModel::JointStencil(s) { }
      
    };

  public:
    RodModelSlopedSurfaceTensionForce   (ElasticRodModel & rod, const std::vector<ElasticRodModel::JointStencil> & stencils, Scalar surface_tension_coeff, Scalar timestep);
    virtual ~RodModelSlopedSurfaceTensionForce   ();

  public:
    void addStencil(Stencil & s) { m_stencils.push_back(s); }
    std::vector<Stencil> & stencils() { return m_stencils; }
    const std::vector<Stencil> & stencils() const { return m_stencils; }
    
  public:
   
    Scalar globalEnergy();
    void globalForce(VecXd & force);
    void globalJacobian(Scalar scale, MatrixBase & Jacobian);

  protected:
    Scalar localEnergy(Stencil & s);
    void localForce(ElementForce & force, Stencil & s);
    void localJacobian(ElementJacobian & jacobian, Stencil & s);
    
  protected:
    std::vector<Stencil> m_stencils;
    Scalar m_surface_tension_coeff;
  };
  
}


#endif // RODMODELSLOPEDSURFACETENSIONFORCE_HH
