/**
 * \file ObjectControllerBase.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#ifndef OBJECTCONTROLLERBASE_HH
#define OBJECTCONTROLLERBASE_HH

#include "BASim/src/Core/ObjectBase.hh"

namespace BASim {

/** Base class for all object controllers. */
class ObjectControllerBase
{
public:
  enum SolverLibrary { PETSC_SOLVER, MKL_SOLVER }; 

  ObjectHandle addObject(ObjectBase& object);

  ObjectBase& getObject(const ObjectHandle& handle);

  const ObjectBase& getObject(const ObjectHandle& handle) const;

  std::vector<ObjectBase*>& getObjects();

  const std::vector<ObjectBase*>& getObjects() const;

  virtual void execute() = 0;
  
  bool isEnabled();
  void setEnabled( bool enabled );

  ObjectControllerBase() { m_enabled = true; }
  virtual ~ObjectControllerBase() {}
protected:

  std::vector<ObjectBase*> m_objects;
  bool m_enabled;
};

/**
 * Object controller handle.
 */

class ObjectControllerHandle : public HandleBase
{
public:

  explicit ObjectControllerHandle(int idx = -1) : HandleBase(idx) {}
};

#include "ObjectControllerBase.inl"

} // namespace BASim

#endif // OBJECTCONTROLLERBASE_HH
