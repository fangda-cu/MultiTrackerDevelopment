/**
 * \file RodModelStretchingForce.h
 *
 * \author fang@cs.columbia.edu
 * \date May 12, 2011
 */

#ifndef RODMODELSTRETCHINGFORCE_HH
#define RODMODELSTRETCHINGFORCE_HH

#include "BASim/src/Physics/DeformableObjects/Rods/RodModelForce.hh"
//#include "BASim/src/Physics/ElasticRods/RodStretchingForce.hh"

namespace BASim 
{
  class RodModelStretchingForce : public RodModelForce
  {
  public:
    struct Stencil
    {
      EdgeHandle e;
      VertexHandle v1;
      VertexHandle v2;
    };
    
  public:
    Scalar globalEnergy();
    void globalForce(VecXd & force);
    void globalJacobian(Scalar scale, MatrixBase & Jacobian);
  
  protected:
    
    
  };
  
}


#endif // RODMODELSTRETCHINGFORCE_HH
