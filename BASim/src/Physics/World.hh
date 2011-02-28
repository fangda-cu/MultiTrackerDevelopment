/**
 * \file World.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/09/2009
 */

#ifndef WORLD_HH
#define WORLD_HH

#ifdef WETA
#include "../Core/ObjectBase.hh"
#include "../Core/ObjectControllerBase.hh"
#include "../Render/RenderBase.hh"
#else
#include "BASim/src/Core/ObjectBase.hh"
#include "BASim/src/Core/ObjectControllerBase.hh"
#include "BASim/src/Render/RenderBase.hh"
#endif

namespace BASim {

/** Class that collects all of the objects in the world. */
class World : public ObjectBase
{
public:
typedef std::vector<ObjectBase*> Objects;
typedef std::vector<ObjectControllerBase*> Controllers;
typedef std::vector<RenderBase*> Renderers;

  World(); //{}
  ~World(); //{}

  void initialize(int argc, char** argv);
//  {
//    #ifdef HAVE_PETSC
//      PetscUtils::initializePetsc(&argc, &argv);
//    #endif // HAVE_PETSC
//  }

  void finalize();
//  {
//    #ifdef HAVE_PETSC
//      PetscUtils::finalizePetsc();
//    #endif // HAVE_PETSC
//  }

  ObjectHandle addObject(ObjectBase* object);
//  {
//    assert( object != NULL );
//
//    int idx = m_objects.size();
//    m_objects.push_back(object);
//
//    return ObjectHandle(idx);
//  }

  ObjectBase& getObject(const ObjectHandle& oh);
//  {
//    assert(oh.isValid());
//    assert(oh.idx() >= 0);
//    assert((size_t) oh.idx() < m_objects.size());
//
//    return *m_objects[oh.idx()];
//  }

  const Objects& getObjects() const;
//  {
//    return m_objects;
//  }

  Objects& getObjects();
//  {
//    return m_objects;
//  }

  ObjectControllerHandle addController(ObjectControllerBase* controller);
//  {
//    assert( controller != NULL );
//    
//    int idx = m_controllers.size();
//    m_controllers.push_back(controller);
//
//    return ObjectControllerHandle(idx);
//  }

  ObjectControllerBase& getController(const ObjectControllerHandle& och);
//  {
//    assert(och.isValid());
//    assert(och.idx() >= 0);
//    assert(och.idx() < (int) m_controllers.size());
//
//    return *m_controllers[och.idx()];
//  }

  Controllers& getControllers();
//  {
//    return m_controllers;
//  }

  void addRenderer(RenderBase* renderer);
//  {
//    assert( renderer != NULL );
//    m_renderers.push_back(renderer);
//  }
  
  Renderers& getRenderers();
//  {
//    return m_renderers;
//  }

  void execute();
//  {
//    Controllers::iterator it;
//    for (it = m_controllers.begin(); it != m_controllers.end(); ++it) {
//      (*it)->execute();
//    }
//  }

//  void serialize( std::ofstream& of )
//  {
//    assert( of.is_open() );
//    
//    // Serialize the object base this class inherits from
//    ObjectBase::serialize(of);
//
//    // Serialize this class
//    std::cout << "World::serialize()" << std::endl;
//  }

protected:

  Objects m_objects;
  Controllers m_controllers;
  Renderers m_renderers;
};

} // namespace BASim

#endif // WORLD_HH
