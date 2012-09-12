/**
 * \file RodTwistShellFaceCouplingForce.h
 *
 * \author fang@cs.columbia.edu
 * \date 09/07/2012
 */

#ifndef RODTWISTSHELLFACECOUPLINGFORCE_HH
#define RODTWISTSHELLFACECOUPLINGFORCE_HH

#include "BASim/src/Physics/DeformableObjects/DefoObjForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Physics/DeformableObjects/Rods/ElasticRodModel.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShell.hh"

#include "BASim/src/Math/ADT/adreal.h"
#include "BASim/src/Math/ADT/advec.h"
#include "BASim/src/Math/ADT/mat3t.h"

namespace BASim 
{
  class RodTwistShellFaceCouplingForce : public DefoObjForce
  {
  public:
    const static int NumDof = 10;

    typedef Eigen::Matrix<Scalar, NumDof, 1>      ElementForce;
    typedef Eigen::Matrix<Scalar, NumDof, NumDof> ElementJacobian;

    typedef CVec3T<Scalar> Vector3d;
    
    struct Stencil
    {
      Stencil() { }
      
      IntArray dofindices;

      // stencil coverage
      EdgeHandle e;
      VertexHandle v;
      
      // cached stiffness
      Scalar stiffness;
      Scalar viscous_stiffness;
      
      // reference strain
      Scalar undeformed_delta;
      Scalar damping_undeformed_delta;
      
      // cached properties (none)
    };

  public:
    RodTwistShellFaceCouplingForce(ElasticRodModel & rod, ElasticShell & shell, const std::vector<Stencil> & stencils, Scalar stiffness, Scalar stiffness_damp, Scalar timestep);
    virtual ~RodTwistShellFaceCouplingForce();

  public:
    void addStencil(Stencil & s) { m_stencils.push_back(s); }
    std::vector<Stencil> & stencils() { return m_stencils; }
    const std::vector<Stencil> & stencils() const { return m_stencils; }
    
  public:
    void updateStiffness();
    void updateViscousReferenceStrain();
    void updateProperties();
    
  public:
    Scalar globalEnergy();
    void globalForce(VecXd & force);
    void globalJacobian(Scalar scale, MatrixBase & Jacobian);

  public:
    ElasticRodModel & rod() { assert(m_rod); return *m_rod; }
    ElasticShell & shell() { assert(m_shell); return *m_shell; }
    
    DeformableObject & defoObj() { assert(m_rod); return m_rod->getDefoObj(); }
    
  protected:
    template <int DO_HESS>
//    adreal<NumDof, DO_HESS, Scalar> adEnergy(const RodTwistShellFaceCouplingForce & mn, const std::vector<Scalar> & deformed, const std::vector<Scalar> & undeformed, Scalar stiffness);
    adreal<NumDof, DO_HESS, Scalar> adEnergy(const RodTwistShellFaceCouplingForce & mn, const Vec3d & A, const Vec3d & B, const Vec3d & C, Scalar theta, const Vec3d & ref1, const Vec3d & ref2, Scalar undeformed_delta, Scalar stiffness);

  protected:
    Scalar localEnergy(Stencil & s, bool viscous);
    void localForce(ElementForce & force, Stencil & s, bool viscous);
    void localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous);

    void computeReferenceStrain();
    
    Vector3d vec2vector(const Vec3d & v);

  protected:
    ElasticRodModel * m_rod;
    ElasticShell * m_shell;
    
    std::vector<Stencil> m_stencils;
    
    Scalar m_stiffness;
    Scalar m_stiffness_damp;
  };
  
}


#endif // RODMODELSTRETCHINGFORCE_HH
