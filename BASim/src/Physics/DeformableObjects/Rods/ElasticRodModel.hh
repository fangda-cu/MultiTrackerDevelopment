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
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"

namespace BASim 
{
  
  const int ELASTIC_ROD_DOFS_PER_VERTEX = 0;  //nodal position vectors
  const int ELASTIC_ROD_DOFS_PER_EDGE = 1;    //edge twist
    
  class ElasticRodModel : public PhysicalModel
  {
    
  public:
    ElasticRodModel(DeformableObject* object, const std::vector<EdgeHandle> & rodedges, Scalar timestep); ////////////////
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
    
    bool isVertexActive(const VertexHandle& v) const;  /////////////////////
    bool isEdgeActive(const EdgeHandle& e) const;   /////////////////////
    bool isFaceActive(const FaceHandle& f) const { return false; }
    bool isTetActive(const TetHandle& t) const { return false; }
    
    void startStep(Scalar time, Scalar timestep);
    void endStep(Scalar time, Scalar timestep);
    
    //*Elastic Shell-specific
//    void setFaceActive(const FaceHandle& f) {m_active_faces[f] = true; }/////////////////////
    
//    const std::vector<ElasticShellForce*>& getForces() const;/////////////////////
//    void addForce(ElasticShellForce* force);/////////////////////
    
    /////////////////////
    void setEdgeXis(const EdgeProperty<Scalar>& xi);
    void setEdgeVelocities(const EdgeProperty<Scalar>& vels);
    void setEdgeUndeformed(const EdgeProperty<Scalar>& undef);

    /////////////////////
    Scalar getEdgeXi(const EdgeHandle& eh) const { return m_xi[eh]; }
    Scalar getEdgeVelocity(const EdgeHandle& eh) const { return m_xi_vel[eh]; }
    Scalar getEdgeUndeformedXi(const EdgeHandle& eh) const { return m_undef_xi[eh]; }
    Scalar getDampingUndeformedXi(const EdgeHandle& eh) const { return m_damping_undef_xi[eh]; }
    
    const VertexProperty<Scalar> & getVertexMasses() const { return m_vertex_masses; }
    void computeMasses();/////////////////////
    
    void setDensity(Scalar density);
    void setRadii(Scalar ra, Scalar rb);/////////////////////
    
    Scalar getMass(const VertexHandle& v) const { return m_obj->getVertexMass(v); }
    Scalar getMass(const EdgeHandle& e) const { return m_edge_masses[e]; }/////////////////////
    Scalar getThickness(const FaceHandle& f) const { }// return m_thicknesses[f]; }/////////////////////
    void setThickness(const FaceHandle& f, Scalar thick) { } //m_thicknesses[f] = thick; }/////////////////////
    Scalar getThickness(const VertexHandle& vh) const;/////////////////////
    Scalar getMaxThickness () const;/////////////////////
    Scalar getMinThickness () const;/////////////////////
    Scalar getVolume(const FaceHandle& f) const {return m_volumes[f]; }/////////////////////
    Scalar getArea(const FaceHandle& f, bool current = true) const;/////////////////////
    
    void getThickness(VertexProperty<Scalar> & vThickness) const;/////////////////////
            
  protected:
    void updateThickness();/////////////////////
            
    //Various shell data
    EdgeProperty<Scalar> m_xi;
    EdgeProperty<Scalar> m_xi_vel;
    EdgeProperty<Scalar> m_undef_xi;    
    EdgeProperty<Scalar> m_damping_undef_xi;
    
    VertexProperty<Scalar> m_vertex_masses;
    EdgeProperty<Scalar> m_edge_masses;
        
    EdgeProperty<Vec2d> m_radii;
    EdgeProperty<Scalar> m_volumes;
        
    Scalar m_density;
    
//    FaceProperty<char> m_active_faces; //list of faces to which this model is applied  /////////////////////
    //Note: this should ideally use booleans, but std::vector<bool> doesn't support references, which we need. (vector<bool> isn't technically a container)
    
    //The base object, and the list of forces
    DeformableObject* m_obj;
//    std::vector<ElasticShellForce*> m_shell_forces;/////////////////////
        
  };
  
}


#endif //ELASTICRODMODEL_H
