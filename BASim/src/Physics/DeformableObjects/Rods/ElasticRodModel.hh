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
//#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"

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
      IntArray dofindices;
    };
    
    // two common types of stencils (reusable by forces)
    struct EdgeStencil : Stencil  // one rod edge
    {
      EdgeHandle e;
      VertexHandle v1;
      VertexHandle v2;
    };
    
    struct JointStencil : Stencil // a pair of incident rod edges
    {
      EdgeHandle e1;
      EdgeHandle e2;
      VertexHandle v1;
      VertexHandle v2;
      VertexHandle v3;
    };
    
  public:
    ElasticRodModel(
      DeformableObject* object, 
      const std::vector<EdgeHandle> & rodedges, 
      Scalar timestep, 
      Scalar youngs_modulus, 
      Scalar youngs_modulus_damping, 
      Scalar shearing_modulus, 
      Scalar shearing_modulus_damping); ////////////////
    
    ~ElasticRodModel();
    
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
    
//    void setEdgeActive(const EdgeHandle & e)/////////////////////
    
    const std::vector<RodModelForce *> & getForces() const { return m_forces; } /////////////////////
    void addForce(RodModelForce * force) { assert(force); m_forces.push_back(force); }  /////////////////////
    
    /////////////////////
    void setEdgeThetas(const EdgeProperty<Scalar>& theta);
    void setEdgeThetaVelocities(const EdgeProperty<Scalar>& vels);
    void setEdgeUndeformedThetas(const EdgeProperty<Scalar>& undef);

    /////////////////////
    Scalar getEdgeTheta(const EdgeHandle& eh) const { return m_theta[eh]; }
    Scalar getEdgeThetaVelocity(const EdgeHandle& eh) const { return m_theta_vel[eh]; }
    Scalar getEdgeUndeformedTheta(const EdgeHandle& eh) const { return m_undef_theta[eh]; }
    Scalar getEdgeDampingUndeformedTheta(const EdgeHandle& eh) const { return m_damping_undef_theta[eh]; }
    
    const VertexProperty<Scalar> & getVertexMasses() const { return m_vertex_masses; }  // inherited from PhysicalModel
    void computeMasses();/////////////////////
    
    void setDensity(Scalar density);
    void setRadii(Scalar ra, Scalar rb);/////////////////////
    
    Scalar getMass(const VertexHandle& v) const { return getDefoObj().getVertexMass(v); }
    Scalar getMass(const EdgeHandle& e) const { return m_edge_masses[e]; }
    
    Vec2d getRadii(const EdgeHandle& e) const { return m_radii[e]; } /////////////////////
    void setRadii(const EdgeHandle& e, const Vec2d & r) { m_radii[e] = r; } /////////////////////
    
//    Scalar getThickness(const VertexHandle& vh) const;/////////////////////
//    Scalar getMaxThickness () const;/////////////////////
//    Scalar getMinThickness () const;/////////////////////

    Scalar getVolume(const EdgeHandle& e) const { return m_volumes[e]; }/////////////////////
//    Scalar getArea(const FaceHandle& f, bool current = true) const;/////////////////////
    
//    void getThickness(VertexProperty<Scalar> & vThickness) const;/////////////////////
          
    // Cached properties, in support for the forces (similar to BASim)
    void updateProperties();
    
    // edge properties
    Vec3d & getEdge(const EdgeHandle & e) { return m_properties_edge[e]; }
    Vec3d & getEdgeTangent(const EdgeHandle & e) { return m_properties_edge_tangent[e]; }
    Scalar & getEdgeLength(const EdgeHandle & e) { return m_properties_edge_length[e]; }
    Vec3d & getReferenceDirector1(const EdgeHandle & e) { return m_properties_reference_director1[e]; }
    Vec3d & getReferenceDirector2(const EdgeHandle & e) { return m_properties_reference_director2[e]; }
    Vec3d & getMaterialDirector1(const EdgeHandle & e) { return m_properties_reference_director1[e]; }
    Vec3d & getMaterialDirector2(const EdgeHandle & e) { return m_properties_reference_director2[e]; }
    
    // vertex properties
    Scalar & getVoronoiLength(const VertexHandle & v) { return m_properties_voronoi_length[v]; }
    Scalar & getReferenceTwist(const VertexHandle & v) { return m_properties_reference_twist[v]; }
    Vec3d & getCurvatureBinormal(const VertexHandle & v) { return m_properties_curvature_binormal[v]; }
    
    
  protected:
    void updateRadii();/////////////////////
            
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
    
    VertexProperty<Scalar>  m_properties_voronoi_length;
    VertexProperty<Scalar>  m_properties_reference_twist;
    VertexProperty<Vec3d>   m_properties_curvature_binormal;    
    
  };
  
}


#endif //ELASTICRODMODEL_H
