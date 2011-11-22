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

  void getScriptedDofs(IntArray& dofIndices, std::vector<Scalar>& dofValues) const;

  void setTimeStep(Scalar dt);
  Scalar getTimeStep();
  void startStep() { for(unsigned int i = 0; i < m_models.size(); ++i) m_models[i]->startStep(); }
  void endStep() { for(unsigned int i = 0; i < m_models.size(); ++i) m_models[i]->endStep(); }

  void addModel(PhysicalModel* model);
  PhysicalModel* getModel(int i) const { assert(i < (int)m_models.size()); return m_models[i]; }
  int numModels() const { return m_models.size(); }

  //Sets up the mapping from a linear list of DOFs to whatever internal DOFs that the associated models have requested.
  void computeDofIndexing();

protected:

  std::vector<PhysicalModel*> m_models; ///< physical models layered on this object (each with its own forces)
  Scalar m_dt;  ///< size of time step

  std::vector<int> m_dofModels; //for each dof, which model does it belong to 
  std::vector<DofHandle> m_dofHandles; //for each dof, the information to look it up in the model (handle, type, DOF number).

};

} // namespace BASim

#endif // DEFORMABLEOBJECT_H