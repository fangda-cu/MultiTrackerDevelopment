/**
 * \file RodModelSlopedSurfaceTensionForce.hh
 *
 * \author batty@cs.columbia.edu
 * \date Sept 18, 2012
 */

#ifndef RODMODELSLOPEDSURFACETENSIONFORCE_HH
#define RODMODELSLOPEDSURFACETENSIONFORCE_HH

#include "BASim/src/Physics/DeformableObjects/Rods/RodModelForce.hh"

//This sets up a surface tension force as described in the discrete viscous threads paper.
//i.e. it sets up truncated cones connecting midpoints of edges, using the radii of each
//edge as the end-radii of the cones. (The radii in turn are dictated by the volume and
//length of the edge, assuming a straight cylinder for each edge.)

namespace BASim 
{

  class RodModelSlopedSurfaceTensionForce : public RodModelForce
  {
  public:
    typedef Eigen::Matrix<Scalar, 12, 1> ElementForce;
    typedef Eigen::Matrix<Scalar, 12, 12> ElementJacobian;

  public:
    RodModelSlopedSurfaceTensionForce   (ElasticRodModel & rod, Scalar surface_tension_coeff);
    virtual ~RodModelSlopedSurfaceTensionForce   ();

  public:
    //void addStencil(Stencil & s) { m_stencils.push_back(s); }
    //std::vector<Stencil> & stencils() { return m_stencils; }
    //const std::vector<Stencil> & stencils() const { return m_stencils; }
    
  public:
   
    Scalar globalEnergy();
    void globalForce(VecXd & force);
    void globalJacobian(Scalar scale, MatrixBase & Jacobian);

  protected:
    //For interior vertices.
    Scalar localEnergy(const ElasticRodModel::JointStencil & s);
    void localForce(ElementForce & force, const ElasticRodModel::JointStencil & s);
    void localJacobian(ElementJacobian & jacobian, const ElasticRodModel::JointStencil & s);

    //For end vertices.
    Scalar localEndEnergy(VertexHandle& vh);
    void localEndForce(VertexHandle& vh, VertexHandle& vh2, ElementForce & force);
    void localEndJacobian(VertexHandle& vh, VertexHandle& vh2, ElementJacobian & jacobian);
    
  protected:
    //std::vector<Stencil> m_stencils;
    Scalar m_surface_tension_coeff;
  };
  
}


#endif // RODMODELSLOPEDSURFACETENSIONFORCE_HH
