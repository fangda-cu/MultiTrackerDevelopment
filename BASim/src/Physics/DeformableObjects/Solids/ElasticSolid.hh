/**
 * \file ElasticSolid.h
 *
 * \author batty@cs.columbia.edu
 * \date May 22, 2012
 */
#ifndef ELASTICSOLID_H
#define ELASTICSOLID_H

#include "BASim/src/Physics/DeformableObjects/PhysicalModel.hh"
#include "BASim/src/Core/TopologicalObject/TopObjProperty.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"

namespace BASim {

const int ELASTIC_SOLID_DOFS_PER_VERTEX = 3; //nodal position vectors

class DeformableObject;
class ElasticSolidForce;

class ElasticSolid : public PhysicalModel {

public:
  ElasticSolid(DeformableObject* object, const TetProperty<char>& solidTets, Scalar timestep);
  ~ElasticSolid();

  //*Inherited from PhysicalModel
  void computeForces(VecXd& force);
  void computeJacobian(Scalar scale, MatrixBase& J);
  
  const Scalar& getDof(const DofHandle& hnd) const;
  void setDof(const DofHandle& hnd, const Scalar& dof);
  const Scalar& getVel(const DofHandle& hnd) const;
  void setVel(const DofHandle& hnd, const Scalar& vel);
  const Scalar& getMass(const DofHandle& hnd) const;
 
  int numVertexDofs() const { return ELASTIC_SOLID_DOFS_PER_VERTEX; }
  int numEdgeDofs() const { return 0; }
  int numFaceDofs() const { return 0; }
  int numTetDofs() const { return 0; }

  bool isVertexActive(const VertexHandle& v) const;
  bool isEdgeActive(const EdgeHandle& e) const; 
  bool isFaceActive(const FaceHandle& f) const;
  bool isTetActive(const TetHandle& t) const { return m_active_tets[t] == 1; }

  void getScriptedDofs(IntArray& dofIndices, std::vector<Scalar>& dofValues, Scalar time) const;

  void startStep(Scalar time, Scalar timestep);
  void endStep(Scalar time, Scalar timestep);

  //*Elastic solid-specific
  void setTetActive(const TetHandle& t) {m_active_tets[t] = true; }

  const std::vector<ElasticSolidForce*>& getForces() const;
  void addForce(ElasticSolidForce* force);

  //All DOFs at once
  // these methods should have be removed because position access is now provided by DeformableObject; but 
  // too much code in other parts of the codebase need to change because they depend on this, so these
  // methods are kept and implemented to redirect the calls
  void setVertexPositions(const VertexProperty<Vec3d>& positions) { m_obj->setVertexPositions(positions); }
  void setVertexVelocities(const VertexProperty<Vec3d>& velocities) { m_obj->setVertexVelocities(velocities); }
  void setVertexUndeformed(const VertexProperty<Vec3d>& undef) { m_obj->setVertexUndeformedPositions(undef); }
  
  //Individual DOFs
  Vec3d getVertexUndeformed(const VertexHandle& v) const { return m_obj->getVertexUndeformedPosition(v); }
  Vec3d getVertexPosition(const VertexHandle& v) const { return m_obj->getVertexPosition(v); }
  Vec3d getVertexVelocity(const VertexHandle& v) const { return m_obj->getVertexVelocity(v); }
  Vec3d getVertexDampingUndeformed(const VertexHandle& v) const { return m_obj->getVertexDampingUndeformedPosition(v); }

  const VertexProperty<Vec3d>& getVertexPositions() const{ return m_obj->getVertexPositions(); }

  void setUndeformedVertexPosition(const VertexHandle& v, const Vec3d& pos) { m_obj->setVertexUndeformedPosition(v, pos); }
  void setVertexPosition(const VertexHandle& v, const Vec3d& pos) { m_obj->setVertexPosition(v, pos); }
  void setVertexVelocity(const VertexHandle& v, const Vec3d& vel) { m_obj->setVertexVelocity(v, vel); }

  const VertexProperty<Scalar> & getVertexMasses() const { return m_vertex_masses; }
  void computeMasses();

  void setDensity(Scalar density);
  
  Scalar getMass(const VertexHandle& v) const { return m_obj->getVertexMass(v); }
  
protected:

  VertexProperty<Scalar> m_vertex_masses;
  
  //"undeformed" configuration that is updated at each step to support Rayleigh damping/viscosity
  //This is also used as the "start of step" configuration for eltopo collision resolution

  Scalar m_density;

  TetProperty<char> m_active_tets; //list of tets to which this model is applied
  //Note: this should ideally use booleans, but std::vector<bool> doesn't support references, which we need. (vector<bool> isn't technically a container)

  //The base object, and the list of forces
  DeformableObject* m_obj;
  std::vector<ElasticSolidForce*> m_solid_forces;

  
};

}


#endif //ELASTICSOLID_H
