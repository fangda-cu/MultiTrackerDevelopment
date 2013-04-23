/**
 * \file PhysicalModel.h
 *
 * \author batty@cs.columbia.edu
 * \date April 18, 2011
 */
#ifndef PHYSICALMODEL_H
#define PHYSICALMODEL_H

#include "BASim/src/Core/Definitions.hh"
#include "BASim/src/Core/TopologicalObject/TopologicalObject.hh"
#include "BASim/src/Physics/DegreeOfFreedom.hh"
#include "BASim/src/Core/TopologicalObject/TopObjProperty.hh"

namespace BASim {

//The base class for different kinds of objects that rely on a DeformableObject for connectivity/geometry information.
//eg. Rods, Shells, Solids
//Each will have different parameters, properties, forces, and behaviours.
//A DeformableObject will have a list of all the models that are applied to it.
//This class also manages the reverse DOF mapping from the internal model's simplices to the global picture.

class DeformableObject;
class MatrixBase;

class PhysicalModel {

public:
  
  PhysicalModel(DeformableObject& obj);
  virtual ~PhysicalModel() {}

  //The number of DOF's of each type that the object needs.
  virtual int numVertexDofs() const = 0;
  virtual int numEdgeDofs() const = 0;
  virtual int numFaceDofs() const = 0;
  virtual int numTetDofs() const = 0;
  
  //Determines whether a particular simplex is active in the given model.
  virtual bool isVertexActive(const VertexHandle& v) const = 0;
  virtual bool isEdgeActive(const EdgeHandle& e) const = 0;
  virtual bool isFaceActive(const FaceHandle& f) const = 0;
  virtual bool isTetActive(const TetHandle& t) const = 0;

  //Returns the starting dof index for a given simplex
  int getVertexDofBase(const VertexHandle& vh) const { return m_vertexDofIdxs[vh]; }
  int getEdgeDofBase(const EdgeHandle& eh) const { return m_edgeDofIdxs[eh]; }
  int getFaceDofBase(const FaceHandle& fh) const { return m_faceDofIdxs[fh]; }
  int getTetDofBase(const TetHandle& th) const { return m_tetDofIdxs[th]; }

  //Sets the initial dof index for each simplex type. This is called by the DeformableModel when the indexing is set up.
  void setVertexDofBase(const VertexHandle&v, int index) { m_vertexDofIdxs[v] = index; }
  void setEdgeDofBase(const EdgeHandle&e, int index) { m_edgeDofIdxs[e] = index; }
  void setFaceDofBase(const FaceHandle&f, int index) { m_faceDofIdxs[f] = index; }
  void setTetDofBase(const TetHandle&t, int index) { m_tetDofIdxs[t] = index; }

  //Functions to compute force and Jacobians for the specific model
  virtual void computeForces(VecXd& force) = 0;
  virtual void computeJacobian(Scalar scale, MatrixBase& J) = 0;
  
  virtual void computeConservativeForcesEnergy(VecXd& f, Scalar& energy) = 0;

  //Accessors for DOFs, velocities, and masses
  virtual const Scalar& getDof(const DofHandle& hnd) const = 0;
  virtual void setDof(const DofHandle& hnd, const Scalar& dof) = 0;

  virtual const Scalar& getVel(const DofHandle& hnd) const = 0;
  virtual void setVel(const DofHandle& hnd, const Scalar& vel) = 0;

//  virtual const Scalar& getMass(const DofHandle& hnd) const = 0;
  
  virtual void startStep(Scalar time, Scalar timestep) = 0;
  virtual void endStep(Scalar time, Scalar timestep) = 0;
  virtual void startIteration(Scalar time, Scalar timestep) { } // these two are not required (so that this addition does not break existing shell code)
  virtual void endIteration(Scalar time, Scalar timestep) { }

  //Masses for the shared position DOFS (return the mass specific to that model, and we'll add them up in posdofsmodel)
//  virtual const VertexProperty<Scalar> & getVertexMasses() const = 0;
//  virtual const Scalar getModelVertexMass(const VertexHandle& vh) const = 0;

  //For constraining particular DOFs
  virtual void getScriptedDofs(IntArray& dofIndices, std::vector<Scalar>& dofValues, Scalar time) const {}
  virtual bool isDofScripted(const DofHandle & hnd) const { return false; }
  
  //Accessor for the deformable object that this model is attached to.
  DeformableObject& getDefoObj() const { return m_obj; }

protected:
  DeformableObject& m_obj;
  
  //For each simplex, store the DOF index of the first DOF for that simplex. 
  //eg. for a given vertex, the first DOF associated to it. The remaining DOFs for that vertex follow in order.
  VertexProperty<int> m_vertexDofIdxs;
  EdgeProperty<int> m_edgeDofIdxs;
  FaceProperty<int> m_faceDofIdxs;
  TetProperty<int> m_tetDofIdxs;
};


}

#endif PHYSICALMODEL