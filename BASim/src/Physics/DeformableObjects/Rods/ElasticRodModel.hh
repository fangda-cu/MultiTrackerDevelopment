/**
 * \file ElasticShell.h
 *
 * \author fang@cs.columbia.edu
 * \date May 11, 2012
 */
#ifndef ELASTICRODMODEL_H
#define ELASTICRODMODEL_H

#include "BASim/src/Physics/DeformableObjects/PhysicalModel.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"

namespace BASim 
{
  class RodModelForce;
  class RodModelStretchingForce;
  class RodModelBendingForce;
  class RodModelTwistingForce;
  
  const int ELASTIC_ROD_DOFS_PER_VERTEX = 0;  //nodal position vectors
  const int ELASTIC_ROD_DOFS_PER_EDGE = 1;    //edge twist
    
  class ElasticRodModel : public PhysicalModel
  {
  public:
    struct Stencil
    {
      Stencil() : id(0), dofindices() { }
      Stencil(const Stencil & s) : id(s.id), dofindices(s.dofindices) { }
      
      int id;
      IntArray dofindices;
    };
    
    // two common types of stencils (reusable by forces)
    struct EdgeStencil : Stencil  // one rod edge
    {
      EdgeStencil() { }
      EdgeStencil(const EdgeStencil & s) : Stencil(s), e(s.e), v1(s.v1), v2(s.v2) { copyData(s); }
      
      // stencil coverage
      EdgeHandle e;
      VertexHandle v1;
      VertexHandle v2;
      
      // per-stencil data (cached properties, updated by updateProperties() automatically)
      
      // copying per-stencil data (used by ElasticRodModel to propagate its computation to forces' stencils)
      void copyData(const EdgeStencil & s) { }
    };
    
    struct JointStencil : Stencil // a pair of incident rod edges
    {      
      JointStencil() { }
      JointStencil(const JointStencil & s) : Stencil(s), e1(s.e1), e2(s.e2), v1(s.v1), v2(s.v2), v3(s.v3), e1flip(s.e1flip), e2flip(s.e2flip) { copyData(s); }
      
      // stencil coverage
      EdgeHandle e1;
      EdgeHandle e2;
      VertexHandle v1;
      VertexHandle v2;
      VertexHandle v3;
      bool e1flip;
      bool e2flip;
      
      // per-stencil data (cached properties, updated by updateProperties() automatically)
      Vec3d curvatureBinormal;
      Scalar referenceTwist;
      Scalar voronoiLength;
      
      // copying per-stencil data (used by ElasticRodModel to propagate its computation to forces' stencils)
      void copyData(const JointStencil & s)
      {
        curvatureBinormal = s.curvatureBinormal;
        referenceTwist = s.referenceTwist;
        voronoiLength = s.voronoiLength;
      }
    };
    
  public:
    ElasticRodModel(DeformableObject* object, const std::vector<EdgeHandle> & rodedges, Scalar timestep); ////////////////
    
    ~ElasticRodModel();
    
    void setup(Scalar youngs_modulus, Scalar youngs_modulus_damping, Scalar shear_modulus, Scalar shear_modulus_damping, Scalar timestep);
    
    //*Inherited from PhysicalModel
    void computeForces(VecXd& force);
    void computeJacobian(Scalar scale, MatrixBase& J);
    
    const Scalar& getDof(const DofHandle& hnd) const;
    void setDof(const DofHandle& hnd, const Scalar& dof);
    const Scalar& getVel(const DofHandle& hnd) const;
    void setVel(const DofHandle& hnd, const Scalar& vel);
    const Scalar& getMass(const DofHandle& hnd) const;
    
    int numVertexDofs() const { return ELASTIC_ROD_DOFS_PER_VERTEX; }
    int numEdgeDofs() const { return ELASTIC_ROD_DOFS_PER_EDGE; }
    int numFaceDofs() const { return 0; }
    int numTetDofs() const { return 0; }
    
    bool isVertexActive(const VertexHandle& v) const { return false; }
    bool isEdgeActive(const EdgeHandle& e) const { return m_edge_active[e]; }
    bool isFaceActive(const FaceHandle& f) const { return false; }
    bool isTetActive(const TetHandle& t) const { return false; }
    
    void startStep(Scalar time, Scalar timestep);
    void endStep(Scalar time, Scalar timestep);
    
    void startIteration(Scalar time, Scalar timestep);
    void endIteration(Scalar time, Scalar timestep);
    
    const std::vector<RodModelForce *> & getForces() const { return m_forces; }
    void addForce(RodModelForce * force) { assert(force); m_forces.push_back(force); }
    
    void setEdgeThetas(const EdgeProperty<Scalar>& theta);
    void setEdgeThetaVelocities(const EdgeProperty<Scalar>& vels);
    void setEdgeUndeformedThetas(const EdgeProperty<Scalar>& undef);

    Scalar getEdgeTheta(const EdgeHandle& eh) const { return m_theta[eh]; }
    Scalar getEdgeThetaVelocity(const EdgeHandle& eh) const { return m_theta_vel[eh]; }
    Scalar getEdgeUndeformedTheta(const EdgeHandle& eh) const { return m_undef_theta[eh]; }
    Scalar getEdgeDampingUndeformedTheta(const EdgeHandle& eh) const { return m_damping_undef_theta[eh]; }
    
    const VertexProperty<Scalar> & getVertexMasses() const { return m_vertex_masses; }  // inherited from PhysicalModel
    void computeMasses();
    
    void setDensity(Scalar density);
    void setRadii(Scalar ra, Scalar rb);
    
    Scalar getMass(const VertexHandle& v) const { return getDefoObj().getVertexMass(v); }
    Scalar getMass(const EdgeHandle& e) const { return m_edge_masses[e]; }
    
    Vec2d getRadii(const EdgeHandle& e) const { return m_radii[e]; }
    void setRadii(const EdgeHandle& e, const Vec2d & r) { m_radii[e] = r; }

    Scalar getVolume(const EdgeHandle& e) const { return m_volumes[e]; }
    
    // Cached properties, in support for the forces (similar to BASim)
    void updateProperties();
    
    // edge properties
    Vec3d & getEdge(const EdgeHandle & e) { return m_properties_edge[e]; }
    Vec3d & getEdgeTangent(const EdgeHandle & e) { return m_properties_edge_tangent[e]; }
    Scalar & getEdgeLength(const EdgeHandle & e) { return m_properties_edge_length[e]; }
    Vec3d & getReferenceDirector1(const EdgeHandle & e) { return m_properties_reference_director1[e]; }
    Vec3d & getReferenceDirector2(const EdgeHandle & e) { return m_properties_reference_director2[e]; }
    Vec3d & getMaterialDirector1(const EdgeHandle & e) { return m_properties_material_director1[e]; }
    Vec3d & getMaterialDirector2(const EdgeHandle & e) { return m_properties_material_director2[e]; }
        
    // reference director 1 is the only derived property that needs initialization. this will only be used in setup().
    // if not specified, the reference directors will be generated randomly in setup().
    void setUndeformedReferenceDirector1(const EdgeProperty<Vec3d> & undeformed_reference_director1);
    
    // set undeformed positions, used in setup(). this overrides the one in position dofs model. this allows the rod
    //  to have a different undeformed configuration than the undef config in position dofs model that serves as a 
    //  common undef config for all models.
    // if not specified, the undef config in position dofs model will be used in setup().
    void setUndeformedPositions(const VertexProperty<Vec3d> & undeformed_positions);
    
    // dof scripting interface inherited from PhysicalModel
    void getScriptedDofs(IntArray & dofIndices, std::vector<Scalar> & dofValues, Scalar time) const;
    
    // scripting on edge dofs
    void constrainEdge(const EdgeHandle & e, Scalar t);
    void constrainEdgeVel(const EdgeHandle & e, Scalar init_value, Scalar velocity, Scalar start_time); //time varying constraint
    void releaseEdge(const EdgeHandle & e);
    bool isConstrained(const EdgeHandle & e) const;
    
  protected:
    void updateRadii();
            
    // list of stencils in the mesh that are part of the rod
    std::vector<EdgeStencil> m_edge_stencils;
    std::vector<JointStencil> m_joint_stencils;
    
    // active flag for each edge, used for isEdgeActive() query which is required by defo obj's dof indexing
    EdgeProperty<char> m_edge_active;
    
    // rod dofs
    EdgeProperty<Scalar> m_theta;
    EdgeProperty<Scalar> m_theta_vel;
    EdgeProperty<Scalar> m_undef_theta;    
    EdgeProperty<Scalar> m_damping_undef_theta;
    
    VertexProperty<Scalar> m_vertex_masses;
    EdgeProperty<Scalar> m_edge_masses;
        
    EdgeProperty<Vec2d> m_radii;
    EdgeProperty<Scalar> m_volumes;
        
    Scalar m_density;

    // forces, including the three internal forces
    std::vector<RodModelForce *> m_forces;
    RodModelStretchingForce * m_stretching_force;
    RodModelBendingForce * m_bending_force;
    RodModelTwistingForce * m_twisting_force;
        
    // cached properties
    EdgeProperty<Vec3d>   m_properties_edge;
    EdgeProperty<Vec3d>   m_properties_edge_tangent;
    EdgeProperty<Scalar>  m_properties_edge_length;
    EdgeProperty<Vec3d>   m_properties_reference_director1;
    EdgeProperty<Vec3d>   m_properties_reference_director2;
    EdgeProperty<Vec3d>   m_properties_material_director1;
    EdgeProperty<Vec3d>   m_properties_material_director2;
    
    // initialization of the cached properties
    const EdgeProperty<Vec3d> * m_undeformed_reference_director1;
    
    // undeformed configuration storage: this is used in setup. this allows specifying for the rod a different 
    //  undeformed configuration than the one in the position dofs model, which will serve as the common undeformed
    //  configuration for all models. if this member is not assigned, the one in position dofs model will be used
    //  in setup().
    const VertexProperty<Vec3d> * m_undeformed_positions;
    
    // dof scripting
    std::vector<std::pair<EdgeHandle, Scalar> > m_edge_constraints;
    std::vector<std::pair<EdgeHandle, Vec3d> > m_edge_vel_constraints;
    
  };
  
}


#endif //ELASTICRODMODEL_H
