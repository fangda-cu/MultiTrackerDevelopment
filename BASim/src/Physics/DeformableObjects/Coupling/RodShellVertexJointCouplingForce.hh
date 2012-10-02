/**
 * \file RodShellVertexJointCouplingForce.h
 *
 * \author fang@cs.columbia.edu
 * \date 10/01/2012
 */

#ifndef RODSHELLVERTEXJOINTCOUPLINGFORCE_HH
#define RODSHELLVERTEXJOINTCOUPLINGFORCE_HH

#include "BASim/src/Physics/DeformableObjects/DefoObjForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Physics/DeformableObjects/Rods/ElasticRodModel.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/ElasticShell.hh"

#include "BASim/src/Math/ADT/adreal.h"
#include "BASim/src/Math/ADT/advec.h"
#include "BASim/src/Math/ADT/mat3t.h"

namespace BASim 
{
  class RodShellVertexJointCouplingForce : public DefoObjForce
  {
  public:
    const static int NumDof = 13;

    typedef Eigen::Matrix<Scalar, NumDof, 1>      ElementForce;
    typedef Eigen::Matrix<Scalar, NumDof, NumDof> ElementJacobian;

    typedef CVec3T<Scalar> Vector3d;
    
    struct Stencil
    {
      Stencil() { }
      
      IntArray dofindices;

      // stencil coverage
      EdgeHandle e;   // rod edge AD
      FaceHandle f;   // shell face ABC
      
//      // cached stiffness
//      Scalar stiffness;
//      Scalar viscous_stiffness;
      
      // reference strain
      Vec3d undeformed_AB;
      Vec3d undeformed_AC;
      Vec3d damping_undeformed_AB;
      Vec3d damping_undeformed_AC;
      
      // cached properties
      Vec3d AB;  // shell face edge AB in rod's material frame
      Vec3d AC;  // shell face edge AC in rod's material frame
    };

  public:
    RodShellVertexJointCouplingForce(ElasticRodModel & rod, ElasticShell & shell, const std::vector<Stencil> & stencils, Scalar stiffness, Scalar stiffness_damp, Scalar timestep);
    virtual ~RodShellVertexJointCouplingForce();

  public:
    void addStencil(Stencil & s) { m_stencils.push_back(s); }
    std::vector<Stencil> & stencils() { return m_stencils; }
    const std::vector<Stencil> & stencils() const { return m_stencils; }
    
  public:
    void updateStiffness();
    void updateViscousReferenceStrain();
    void updateProperties();
    
  public:
    void startStep(Scalar time, Scalar timestep);
    void endStep(Scalar time, Scalar timestep);
    void startIteration(Scalar time, Scalar timestep);
    void endIteration(Scalar time, Scalar timestep);
    
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
    adreal<NumDof, DO_HESS, Scalar> adEnergy(const RodShellVertexJointCouplingForce & mn, const Vec3d & A, const Vec3d & B, const Vec3d & C, const Vec3d & D, Scalar theta, const Vec3d & ref1, const Vec3d & ref2, const Vec3d & undeformed_AB, const Vec3d & undeformed_AC, Scalar stiffness);

  protected:
    Scalar localEnergy(Stencil & s, bool viscous);
    void localForce(ElementForce & force, Stencil & s, bool viscous);
    void localJacobian(ElementJacobian & jacobian, Stencil & s, bool viscous);

    void computeReferenceStrain();
    
    Vector3d vec2vector(const Vec3d & v);
    std::vector<VertexHandle> getVertices(const Stencil & s);
    
  protected:
    ElasticRodModel * m_rod;
    ElasticShell * m_shell;
    
    std::vector<Stencil> m_stencils;
    
    Scalar m_stiffness;
    Scalar m_stiffness_damp;
  };
  
}


#endif // RODSHELLVERTEXJOINTCOUPLINGFORCE_HH
