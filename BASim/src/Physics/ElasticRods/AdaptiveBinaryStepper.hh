/**
 * \file AdaptiveBinaryStepper.hh
 *
 * \author smith@cs.columbia.edu
 * \date 06/21/2010
 */


#ifndef ADAPTIVEBINARYSTEPPER_HH
#define ADAPTIVEBINARYSTEPPER_HH

#include "RodTimeStepper.hh"
#include "../../Core/ObjectControllerBase.hh"
#include "MinimalRodStateBackup.hh"

namespace BASim 
{

class AdaptiveBinaryStepper : public RodTimeStepper
{
  
public:

  AdaptiveBinaryStepper( ElasticRod* rod, RodTimeStepper* rod_stepper, double minstep = 1.0e-9 )
  : RodTimeStepper(*rod)
  , m_rod(rod)
  , m_stepper(rod_stepper)
  , m_min_step(minstep)
  {
    assert( m_rod != NULL );
    assert( m_stepper != NULL );
    
  #ifdef TIMING_ON
    IntStatTracker::getIntTracker("ADAPTIVE_STEPS_NEEDED",0);
    DoubleStatTracker::getDoubleTracker("MIN_TIME_STEP_ENCOUNTERED",std::numeric_limits<double>::infinity());
  #endif    
  }

  virtual ~AdaptiveBinaryStepper()
  {
    if( m_stepper != NULL )
    {
      delete m_stepper;
      m_stepper = NULL;
    }
  }

  virtual void setTimeStep( double dt )
  {
    assert( m_rod != NULL );
    assert( m_stepper != NULL );

//std::cout << "ADAP dt " << dt << "\n";

    // Sets timestep of rod too
    m_stepper->setTimeStep(dt);
    RodTimeStepper::setTimeStep(dt);
  }

  Scalar getTimeStep() const
  {
    return m_stepper->getTimeStep();
  }


  Scalar getTime() const
  {
    double t = m_stepper->getTime();
    std::cout << "AdaptiveBinaryStepper::getTime() = " << t << std::endl;
    return t;
  }

  virtual void setTime(const double time)
  {
    std::cout << "Setting time in AdaptiveBinaryStepper to be " << time << std::endl;
    m_stepper->setTime(time);
    RodTimeStepper::setTime(time);
  }
  
  bool execute()
  {
    assert( m_rod != NULL );
    assert( m_stepper != NULL );

    // Clean this up in a minute
    double dt = m_stepper->getTimeStep();
    bool returnstatus = recursiveSolve( m_stepper->getTimeStep() );
    assert( dt == m_stepper->getTimeStep() );
    return returnstatus;
  }
  
  bool recursiveSolve( double dt )
  {
    #ifdef TIMING_ON
      if( dt < DoubleStatTracker::getDoubleTracker("MIN_TIME_STEP_ENCOUNTERED").getVal() ) DoubleStatTracker::getDoubleTracker("MIN_TIME_STEP_ENCOUNTERED") = dt;
    #endif
    
    if( dt < m_min_step ) 
    {
      std::cerr << "\033[31;1mERROR IN ADAPTIVE BINARY STEPPER:\033[m Timestep fell below minimum, aborting recursive solve." << std::endl;
      exit(1);
    }

    this->setTimeStep(dt);

    // Backup the rod in case the solve fails.
    m_stepper->backupResize();
    m_stepper->backup();

    // TMP DBG
    //RodState prestepstate;
    //prestepstate.copyState(*m_rod);
    //prestepstate.print(*m_rod);
    // END TMP DBG    
    
    // If solve was successfull, we are done.
    if( m_stepper->execute() ) return true;
    std::cerr << "\033[31;1mWARNING IN ADAPTIVE BINARY STEPPER:\033[m Solve failed for timestep: " << dt << ". Taking two timesteps of half length." << std::endl;

    #ifdef TIMING_ON
      IntStatTracker::getIntTracker("ADAPTIVE_STEPS_NEEDED") += 1;
    #endif

    // The solve failed, restore the rod to its previous state.
    m_stepper->backupRestore();
    m_stepper->backupClear();

    // TMP DBG
    //prestepstate.restoreState(*m_rod);
    //prestepstate.compareProperties(*m_rod);
    //std::cout << "FAILED STATE" << std::endl;
    //backupstate2.print(*m_rod);
    // END TMP DBG    
    
    // Otherwise, we must run two solves with half the timestep.
    bool firstsucceeded = recursiveSolve( 0.5*dt );
    if( !firstsucceeded ) { setTimeStep(dt); return false; }
    bool secondsucceeded = recursiveSolve( 0.5*dt );
    if( !secondsucceeded ) { setTimeStep(dt); return false; }
    
    this->setTimeStep(dt);
    return firstsucceeded && secondsucceeded;
  }

private:
  ElasticRod* m_rod;
  RodTimeStepper* m_stepper;
  double m_min_step;
};
  
} // namespace BASim

#endif // ADAPTIVEBINARYSTEPPER_HH

