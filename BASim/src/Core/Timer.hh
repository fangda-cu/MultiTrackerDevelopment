#ifndef TIMER_HH
#define TIMER_HH

#include "Definitions.hh"

// Include a high-res timer library
#include "Cycle.hh"

namespace BASim {

// macros to compile timing in/out
#define TIMING_ON

#ifdef TIMING_ON

#define START_TIMER(name)                               \
  Timer::getTimer(name).beginBlock();			
 
#define STOP_TIMER(name)                               \
  Timer::getTimer(name).endBlock();                   

#define PRINT_TIMER(name)                                         \
  {                                                               \
    std::cout << name << ": " << Timer::getTimer(name).getTotal() \
              << std::endl;                                       \
  }

#define CLEAR_TIMER(name) \
  {                                             \
    Timer::getTimer(name).clear();              \
  }                                             \

#else

#define START_TIMER(name)                               
 
#define STOP_TIMER(name)                               

#define PRINT_TIMER(name)                                         

#define CLEAR_TIMER(name) 

#endif

// #else // timing off
// #define PRINT_TIMER(name) {}
// #define CLEAR_TIMER(name) {}
// #endif

/** A very simple timer class. */
class Timer
{
public:

  typedef std::map<std::string, Timer*> TimerMap;

  /// returns the timer for a given name; reuses timers if they exist
  static Timer& getTimer(std::string name);

  /// deletes all timers created via getTimer() function
  static void deleteAllTimers();

  /// dumps timing information to stdout
  static void report();

  /// reset any previous timings
  void clear();

  /// start the timer by increasing the nesting level
  void beginBlock();

  /// decreases the nesting level, and stops the timer if the outermost block is reached
  /// when timer is stopped, the elapsed time is added the total elapsed time
  void endBlock();

  /// returns the total elapsed time
  double getElapsed();

  /// returns the total elapsed time
  double getCumulativeElapsed();

  /// returns the number of timings elapsed time represents
  int getCount()
  {
    return count;
  }

  /// returns the number of timings elapsed time represents
  int getCumulativeCount()
  {
    return cumulativeCount;
  }

  int getNestingLevel() 
  {
    return nestingLevel;
  }

  int getDeepestNestingLevel() 
  {
    return deepestNestingLevel;
  }

  std::vector<double> getTimes()
  {
    return times;
  }

  // Timer& operator-(const Timer& t)
  // {
  //   total = total - t.total;
  //   return *this;
  // }

  friend std::ostream& operator<<(std::ostream& os, Timer& timer);

  void lap(); ///< adds "elapsed" and "count" to their cumulative counterpart, then resets them


private:
  static TimerMap timerMap; ///< static map of name->Timer

  /// private timer constructor; create timers via static getTimer() function
  Timer();

  std::vector<double> times;

  ticks  startTicks;  ///< last time the timer's start function was called
  double elapsedT;    ///< time measured since last report
  int count;          ///< the number of outer blocks encountered since last report

  int nestingLevel;  ///< if nestingLevel is greater than zero, the timer is on

  double cumulativeElapsedT; ///< total time over all reports
  double cumulativeCount;    ///< total count over all reports
  int deepestNestingLevel;   ///< deepest nesting level encountered (helps spot unbalanced begin/end-block)
};

} // namespace BASim

#endif // TIMER_HH
