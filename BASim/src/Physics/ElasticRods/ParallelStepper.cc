/**
 * \file ParallelStepper.cc
 *
 * \author smith@cs.columbia.edu
 * \date 05/23/2010
 */

#include "ParallelStepper.hh"


namespace BASim 
{
  
ParallelStepper::ParallelStepper()
: m_controllers()
, m_scripting_controllers()
{}

ParallelStepper::~ParallelStepper()
{}
  
  
void ParallelStepper::addController( ObjectControllerBase* controller )
{
  assert( controller != NULL );
  m_controllers.push_back( controller );
}

void ParallelStepper::addScriptingController( ObjectControllerBase* controller )
{
  assert( controller != NULL );
  m_scripting_controllers.push_back( controller );
}
  
bool ParallelStepper::execute()
{
  #ifdef HAVE_OPENMP
  #pragma omp parallel for
  #endif
  for( int i = 0; i < (int) m_scripting_controllers.size(); ++i ) 
  {
    assert( m_scripting_controllers[i] != NULL );
    m_scripting_controllers[i]->execute();
  }

  #ifdef HAVE_OPENMP
  #pragma omp parallel for
  #endif
  for( int i = 0; i < (int) m_controllers.size(); ++i ) 
  {
    assert( m_controllers[i] != NULL );
    m_controllers[i]->execute();
  }
  
  // TODO: actually check return values to ensure operation was cool.
  return true;
}


}
