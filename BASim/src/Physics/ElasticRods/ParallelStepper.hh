/**
 * \file ParallelStepper.hh
 *
 * \author smith@cs.columbia.edu
 * \date 05/23/2010
 */


#ifndef PARALLELSTEPPER_HH
#define PARALLELSTEPPER_HH

#ifdef WETA
#include "../..//Core/ObjectControllerBase.hh"
#include "ElasticRod.hh"
#include "RodTimeStepper.hh"
#include "RodMassDamping.hh"
#include "RodGravity.hh"
#include "../../Math/Math.hh"
#else
#include "BASim/src/Core/ObjectControllerBase.hh"
#include "BASim/src/Physics/ElasticRods/ElasticRod.hh"
#include "BASim/src/Physics/ElasticRods/RodTimeStepper.hh"
#include "BASim/src/Physics/ElasticRods/RodMassDamping.hh"
#include "BASim/src/Physics/ElasticRods/RodGravity.hh"
#include "BASim/src/Math/Math.hh"
#endif

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include <limits>

namespace BASim {

  /**
   * Class to execute a number of controllers in parallel. This assumes
   * that the controllers can safely execute in parallel. 
   */
  class ParallelStepper : public ObjectControllerBase
  {

  public:
    /**
     * Default constructor.
     */
    ParallelStepper();

    /**
     * Destructor.
     */
    virtual ~ParallelStepper();

    /**
     * Adds a rod that will be evolved in time using this BridsonStepper.
     */
    void addController( ObjectControllerBase* controller );

    /**
     * Adds a controller that, for example, sets boundary conditions and MUST
     * be executed before other controllers. 
     */
    void addScriptingController( ObjectControllerBase* controller );

    /**
     * Executes all inserted controllers in parallel. 
     */
    bool execute();
    
  private:
    std::vector<ObjectControllerBase*> m_controllers;
    std::vector<ObjectControllerBase*> m_scripting_controllers;

  };
  
} // namespace BASim

#endif // PARALLELSTEPPER_HH

