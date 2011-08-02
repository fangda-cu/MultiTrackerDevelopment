#include "World.hh"

namespace BASim 
{
  
World::World()
{}

World::~World()
{}

void World::initialize(int argc, char** argv)
{
  #ifdef HAVE_PETSC
    PetscUtils::initializePetsc(&argc, &argv);
  #endif // HAVE_PETSC
}

void World::finalize()
{
  #ifdef HAVE_PETSC
    PetscUtils::finalizePetsc();
  #endif // HAVE_PETSC
}

ObjectHandle World::addObject(ObjectBase* object)
{
  assert( object != NULL );

  int idx = (int)( m_objects.size() );
  m_objects.push_back(object);

  return ObjectHandle(idx);
}

ObjectBase& World::getObject(const ObjectHandle& oh)
{
  assert(oh.isValid());
  assert(oh.idx() >= 0);
  assert((size_t) oh.idx() < m_objects.size());

  return *m_objects[oh.idx()];
}

const World::Objects& World::getObjects() const
{
  return m_objects;
}

World::Objects& World::getObjects()
{
  return m_objects;
}

ObjectControllerHandle World::addController(ObjectControllerBase* controller)
{
  assert( controller != NULL );
  
  int idx = m_controllers.size();
  m_controllers.push_back(controller);

  return ObjectControllerHandle(idx);
}

ObjectControllerBase& World::getController(const ObjectControllerHandle& och)
{
  assert(och.isValid());
  assert(och.idx() >= 0);
  assert(och.idx() < (int) m_controllers.size());

  return *m_controllers[och.idx()];
}

World::Controllers& World::getControllers()
{
  return m_controllers;
}

void World::addRenderer(RenderBase* renderer)
{
  assert( renderer != NULL );
  m_renderers.push_back(renderer);
}

World::Renderers& World::getRenderers()
{
  return m_renderers;
}

void World::execute()
{

   Controllers::iterator it;
  for (it = m_controllers.begin(); it != m_controllers.end(); ++it) {
    (*it)->execute();
  }
}

} // namespace BASim
