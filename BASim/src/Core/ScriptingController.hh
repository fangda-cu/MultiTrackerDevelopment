/**
 * \file ScriptingController.hh
 *
 * \author smith@cs.columbia.edu
 * \date 07/21/2010
 */

#ifndef SCRIPTINGCONTROLLER_HH
#define SCRIPTINGCONTROLLER_HH

#include "BASim/src/Core/ObjectControllerBase.hh"

namespace BASim {

  class ScriptingController : public ObjectControllerBase
  {
  public:    
    ScriptingController( double time, double dt );

    void setTime( double time );
    void setDt( double dt );

    double getTime();
    double getDt();

  protected:
    double m_time;
    double m_dt;
  };
  
} // namespace BASim

#endif // STATTRACKER_HH
