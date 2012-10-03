//
//  PositionDofsModel.hh
//  BASim
//
//  Created by Fang Da (fang@cs.columbia.edu) on 5/11/12.
//  Copyright (c) 2012 Columbia. All rights reserved.
//

#ifndef BASim_PositionDofsModel_hh
#define BASim_PositionDofsModel_hh

#include "BASim/src/Physics/DeformableObjects/PhysicalModel.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"



namespace BASim 
{

  class DefoObjForce;
  //
  // This is a default model that stores position data. Its sole purpose is to carry
  // vertex position Dofs, which are shared between common models.
  // It also supports generic external forces that apply only to vertex positions, such as gravity.
  //
  // The deformable object is specifically aware of this model and will
  // query it when asked for position Dofs.
  //
  // Other than sharing position Dofs through this dummy model, there is no
  // Dof sharing mechanism in the current design. Each model has to allocate
  // its own Dofs, and duplication cannot be avoided if multiple models have
  // the common Dofs.
  //
  class PositionDofsModel : public PhysicalModel
  {
  public:
    PositionDofsModel(DeformableObject * obj) : 
      PhysicalModel(*obj), 
      m_positions(obj), 
      m_velocities(obj), 
      m_vertex_masses(obj), 
      m_undeformed_positions(obj),
      m_damping_undeformed_positions(obj)
    { }
      
    virtual ~PositionDofsModel();

  public:
    //The number of DOF's of each type that the object needs.
    virtual int numVertexDofs() const { return 3; }
    virtual int numEdgeDofs() const { return 0; }
    virtual int numFaceDofs() const { return 0; }
    virtual int numTetDofs() const { return 0; }
    
  public:
    //Determines whether a particular simplex is active in the given model.
    virtual bool isVertexActive(const VertexHandle& v) const { return true; }
    virtual bool isEdgeActive(const EdgeHandle& e) const { return false; }
    virtual bool isFaceActive(const FaceHandle& f) const { return false; }
    virtual bool isTetActive(const TetHandle& t) const { return false; }

  public:
    //Functions to compute force and Jacobians for the specific model
    virtual void computeForces(VecXd& force);
    virtual void computeJacobian(Scalar scale, MatrixBase& J);
    virtual void computeConservativeForcesEnergy(VecXd& f, Scalar& energy) { }
    
  public:
    //Accessors for DOFs, velocities, and masses
    virtual const Scalar& getDof(const DofHandle& hnd) const
    {
      assert(hnd.getType() == DofHandle::VERTEX_DOF);
      const VertexHandle& vh = static_cast<const VertexHandle&>(hnd.getHandle());
      return m_positions[vh][hnd.getNum()];
    }
    
    virtual void setDof(const DofHandle& hnd, const Scalar& dof)
    {
      assert(hnd.getType() == DofHandle::VERTEX_DOF);
      const VertexHandle& vh = static_cast<const VertexHandle&>(hnd.getHandle());
      m_positions[vh][hnd.getNum()] = dof;
    }
    
    virtual const Scalar& getVel(const DofHandle& hnd) const
    {
      assert(hnd.getType() == DofHandle::VERTEX_DOF);
      const VertexHandle& vh = static_cast<const VertexHandle&>(hnd.getHandle());
      return m_velocities[vh][hnd.getNum()];
    }
    
    virtual void setVel(const DofHandle& hnd, const Scalar& vel)
    {
      assert(hnd.getType() == DofHandle::VERTEX_DOF);
      const VertexHandle& vh = static_cast<const VertexHandle&>(hnd.getHandle());
      m_velocities[vh][hnd.getNum()] = vel;
    }
    
    virtual const Scalar& getMass(const DofHandle& hnd) const
    {
      assert(hnd.getType() == DofHandle::VERTEX_DOF);
      const VertexHandle& vh = static_cast<const VertexHandle&>(hnd.getHandle());
      return m_vertex_masses[vh];
    }
    
    virtual void setMass(const DofHandle& hnd, const Scalar& mass)
    {
      assert(hnd.getType() == DofHandle::VERTEX_DOF);
      const VertexHandle& vh = static_cast<const VertexHandle&>(hnd.getHandle());
      m_vertex_masses[vh] = mass;
    }
    
  public:
    // Accessors based on VertexHandle
    //All DOFS at once
    const VertexProperty<Vec3d>& getPositions() const                   { return m_positions; }
    const VertexProperty<Vec3d>& getVelocities() const                  { return m_velocities; }
    const VertexProperty<Vec3d>& getUndeformedPositions() const         { return m_undeformed_positions; }
    const VertexProperty<Vec3d>& getDampingUndeformedPositions() const  { return m_damping_undeformed_positions; }
    
    void setPositions                 (const VertexProperty<Vec3d>& pos) { m_positions = pos; }
    void setVelocities                (const VertexProperty<Vec3d>& vel) { m_velocities = vel; }
    void setUndeformedPositions       (const VertexProperty<Vec3d>& pos) { m_undeformed_positions = pos; }
    void setDampingUndeformedPositions(const VertexProperty<Vec3d>& pos) { m_damping_undeformed_positions = pos; }
    
    //Individual DOFs
    Vec3d getPosition                 (const VertexHandle& v) const { return m_positions[v]; }
    Vec3d getVelocity                 (const VertexHandle& v) const { return m_velocities[v]; }
    Scalar getMass                    (const VertexHandle& v) const { return m_vertex_masses[v]; }
    Vec3d getUndeformedPosition       (const VertexHandle& v) const { return m_undeformed_positions[v]; }
    Vec3d getDampingUndeformedPosition(const VertexHandle& v) const { return m_damping_undeformed_positions[v]; }
    
    void setPosition                  (const VertexHandle& v, const Vec3d& pos) { m_positions[v] = pos; }
    void setVelocity                  (const VertexHandle& v, const Vec3d& vel) { m_velocities[v] = vel; }
    void setMass                      (const VertexHandle& v, Scalar m)         { m_vertex_masses[v] = m; }
    void setUndeformedPosition        (const VertexHandle& v, const Vec3d& pos) { m_undeformed_positions[v] = pos; }
    void setDampingUndeformedPosition (const VertexHandle& v, const Vec3d& pos) { m_damping_undeformed_positions[v] = pos; }
    
    // Masses computation
    // Masses are computed by summing the masses computed by all models
    void clearMasses() { m_vertex_masses.assign(0); }
    void accumulateMasses(const VertexProperty<Scalar> & masses);
    void accumulateMass(const VertexHandle & v, Scalar mass) { m_vertex_masses[v] += mass; }

    // in compliance with the PhysicalModel interface, these methods have to be implemented; but are
    // actually useless here due to the special role of this model.
    virtual const VertexProperty<Scalar> & getVertexMasses() const { return m_vertex_masses; }
    virtual const Scalar getModelVertexMass(const BASim::VertexHandle & vh) const { return m_vertex_masses[vh]; }

    // dof scripting interface inherited from PhysicalModel
    void getScriptedDofs(IntArray & dofIndices, std::vector<Scalar> & dofValues, Scalar time) const;

    // scripting on position dofs
    void constrainVertex(const VertexHandle & v, const Vec3d & pos);
    void constrainVertex(const VertexHandle & v, PositionConstraint * p); //time varying constraint
    void releaseVertex(const VertexHandle & v);
    bool isConstrained(const VertexHandle & v) const;
    
    void addForce(DefoObjForce* force) { m_position_forces.push_back(force); }

  public:
    virtual void startStep(Scalar time, Scalar timestep);
    virtual void endStep(Scalar time, Scalar timestep);
    
  protected:
    // Dof, Dof dot, and Dof mass
    VertexProperty<Vec3d> m_positions;
    VertexProperty<Vec3d> m_velocities;
    VertexProperty<Scalar> m_vertex_masses;
    
    // Dof bar
    VertexProperty<Vec3d> m_undeformed_positions;
    VertexProperty<Vec3d> m_damping_undeformed_positions; 

    // Position dof constraints 
    std::vector<VertexHandle> m_constrained_vertices;
    std::vector<PositionConstraint*> m_constraint_positions;    

    // Generic position dof forces (gravity, etc)
    std::vector<DefoObjForce*> m_position_forces;
  };
  
}


#endif
