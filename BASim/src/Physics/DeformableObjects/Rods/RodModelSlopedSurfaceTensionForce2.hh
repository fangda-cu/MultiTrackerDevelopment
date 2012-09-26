/**
 * \file RodModelSlopedSurfaceTensionForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date Sept 18, 2012
 */

#ifndef RODMODELSLOPEDSURFACETENSIONFORCE2_HH
#define RODMODELSLOPEDSURFACETENSIONFORCE2_HH

#include "BASim/src/Physics/DeformableObjects/Rods/RodModelForce.hh"

namespace BASim 
{

  class RodModelSlopedSurfaceTensionForce2 : public RodModelForce
  {
  public:
    typedef Eigen::Matrix<Scalar, 12, 1> ElementForce;
    typedef Eigen::Matrix<Scalar, 12, 12> ElementJacobian;

    struct Stencil : public ElasticRodModel::ThreeEdgeStencil
    {
      Stencil(const ElasticRodModel::ThreeEdgeStencil & s) : ElasticRodModel::ThreeEdgeStencil(s) { }
      
    };

  public:
    RodModelSlopedSurfaceTensionForce2(ElasticRodModel & rod, const std::vector<ElasticRodModel::ThreeEdgeStencil> & stencils, Scalar surface_tension_coeff);
    virtual ~RodModelSlopedSurfaceTensionForce2   ();

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
    
    void collectData(const Stencil& s, std::vector<Scalar>& deformed_data, Vec3d& vols, bool& prevValid, bool& nextValid);
  

  protected:
    std::vector<Stencil> m_stencils;
    Scalar m_surface_tension_coeff;
  };
  
}


#endif // RODMODELSLOPEDSURFACETENSIONFORCE2_HH
