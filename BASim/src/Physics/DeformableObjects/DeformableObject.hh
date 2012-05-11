/**
 * \file DeformableObject.h
 *
 * \author batty@cs.columbia.edu
 * \date April 18, 2011
 */

#ifndef DEFORMABLEOBJECT_H
#define DEFORMABLEOBJECT_H

#include "BASim/src/Physics/PhysObject.hh"
#include "BASim/src/Physics/DeformableObjects/PhysicalModel.hh"
#include "BASim/src/Core/TopologicalObject/TopObjProperty.hh"

namespace BASim {

class MatrixBase;
class PositionDofsModel;

//A class to manage position constraints.
class PositionConstraint {
public:
  virtual Vec3d operator()(Scalar time) = 0;
};

class FixedPositionConstraint : public PositionConstraint {
  Vec3d m_position;

public:
  FixedPositionConstraint(const Vec3d& pos) {m_position = pos;}
  Vec3d operator()(Scalar time){return m_position;}

};

class FixedVelocityConstraint : public PositionConstraint {
  Vec3d m_initial_position;
  Vec3d m_velocity;
  Scalar m_start_time;

public:
  FixedVelocityConstraint(const Vec3d& pos, const Vec3d& velocity, Scalar creationTime)
    : m_initial_position(pos), m_velocity(velocity), m_start_time(creationTime) {}

  Vec3d operator()(Scalar time){return m_initial_position + (time-m_start_time)*m_velocity;}

};


/** Base class for dynamic objects comprised of simplices.
*/
class DeformableObject : public PhysObject
{
public:

  DeformableObject();
  virtual ~DeformableObject();

  /** \name Inherited from PhysObject */
  //@{
  virtual void computeForces(VecXd& force);
  virtual void computeJacobian(Scalar scale, MatrixBase& J);

  virtual const Scalar& getDof(int i) const;
  virtual void setDof(int i, const Scalar& dof);

  virtual const Scalar& getVel(int i) const;
  virtual void setVel(int i, const Scalar& vel);

  virtual const Scalar& getMass(int i) const;
  //@}

  // position dofs access (accessible by all models because position dofs are shared)
  //All DOFS at once
  const VertexProperty<Vec3d>& getVertexPositions() const;
  const VertexProperty<Vec3d>& getVertexVelocities() const;
  const VertexProperty<Vec3d>& getVertexUndeformedPositions() const;
  const VertexProperty<Vec3d>& getVertexDampingUndeformedPositions() const;
  
  void setVertexPositions                 (const VertexProperty<Vec3d>& pos);
  void setVertexVelocities                (const VertexProperty<Vec3d>& vel);
  void setVertexUndeformedPositions       (const VertexProperty<Vec3d>& pos);
  void setVertexDampingUndeformedPositions(const VertexProperty<Vec3d>& pos);
  
  //Individual DOFs
  Vec3d getVertexPosition                 (const VertexHandle& v) const;
  Vec3d getVertexVelocity                 (const VertexHandle& v) const;
  Vec3d getVertexUndeformedPosition       (const VertexHandle& v) const;
  Vec3d getVertexDampingUndeformedPosition(const VertexHandle& v) const;
  
  void setVertexPosition                  (const VertexHandle& v, const Vec3d& pos);
  void setVertexVelocity                  (const VertexHandle& v, const Vec3d& vel);
  void setVertexUndeformedPosition        (const VertexHandle& v, const Vec3d& pos);
  void setVertexDampingUndeformedPosition (const VertexHandle& v, const Vec3d& pos);
  
  void clearMasses();
  void accumulateMasses(const VertexProperty<Scalar>& masses);  
  void accumulateMass(const VertexHandle&v, Scalar mass);
  
  void getScriptedDofs(IntArray& dofIndices, std::vector<Scalar>& dofValues, Scalar time) const;

  void setTimeStep(Scalar dt);
  Scalar getTimeStep();
  void setTime(Scalar time);
  Scalar getTime();

  void startStep() { for(unsigned int i = 0; i < m_models.size(); ++i) m_models[i]->startStep(m_time, m_dt); }
  void endStep() { for(unsigned int i = 0; i < m_models.size(); ++i) m_models[i]->endStep(m_time, m_dt); }

  void addModel(PhysicalModel* model);
  PhysicalModel* getModel(int i) const { assert(i < (int)m_models.size()); return m_models[i]; }
  int numModels() const { return m_models.size(); }

  //Sets up the mapping from a linear list of DOFs to whatever internal DOFs that the associated models have requested.
  void computeDofIndexing();

protected:

  std::vector<PhysicalModel*> m_models; ///< physical models layered on this object (each with its own forces)
  PositionDofsModel * m_posdofsmodel;
  Scalar m_dt;  ///< size of time step
  Scalar m_time; //current time 

  std::vector<int> m_dofModels; //for each dof, which model does it belong to 
  std::vector<DofHandle> m_dofHandles; //for each dof, the information to look it up in the model (handle, type, DOF number).

};

} // namespace BASim

#endif // DEFORMABLEOBJECT_H
