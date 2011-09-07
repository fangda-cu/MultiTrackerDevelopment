/**
 * \file PhysObject.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef PHYSOBJECT_HH
#define PHYSOBJECT_HH

#include "BASim/src/Core/TopologicalObject/TopologicalObject.hh"
#include "BASim/src/Physics/DegreeOfFreedom.hh"

namespace BASim {

class MatrixBase;

/** Base class for physics objects. Provides an array-like interface
    to the degrees of freedom of the object. Derived classes and/or
    other classes that use a class derived from this one are
    responsible for setting up the indexing of the degrees of freedom
    (how to map indices to the degrees of freedom and vice-versa).
*/
class PhysObject : public TopologicalObject
{
public:

  PhysObject();

  virtual ~PhysObject() {}

  /** The number of degrees of freedom this object has. */
  int ndof() const;

  /** This function must be called after all initialization to the
      object has been made. */
  virtual void setup() {}

  /** \name Accessors

      Provide access to degrees of freedom, velocities, and masses. */

  //@{

  virtual const Scalar& getDof(int i) const = 0;
  virtual void setDof(int i, const Scalar& dof) = 0;

  virtual const Scalar& getVel(int i) const = 0;
  virtual void setVel(int i, const Scalar& vel) = 0;

  virtual const Scalar& getMass(int i) const = 0;
  //virtual void setMass(int i, const Scalar& mass) = 0;


  //@}

  /** \name Internal force and Jacobian computation */

  //@{

  virtual void computeForces(VecXd& force) {}
  virtual void computeJacobian(MatrixBase& J) {}
  
  //@}

  void addMapping(const DofHandle& dof, int index) { m_map.addMapping(dof, index); ++m_ndof;}
  
protected:
  
  int m_ndof; ///< number of degrees of freedom
  DOFMap m_map; ///< mapping from indices to degrees of freedom
  
};


/** Handle for referring to a PhysObject. */
class PhysObjectHandle : public HandleBase
{
public:

  explicit PhysObjectHandle(int idx = -1) : HandleBase(idx) {}
};

#include "PhysObject.inl"

} // namespace BASim

#endif // PHYSOBJECT_HH
