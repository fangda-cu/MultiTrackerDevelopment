#ifndef FENV_H
#define FENV_H
#pragma fenv_access (on)
#include <float.h>

//define the functions we need. I hate Microsoft! Just implement C99 already, c'mon.
#define FE_DOWNWARD   _RC_DOWN
#define FE_UPWARD     _RC_UP
#define FE_TONEAREST  _RC_NEAR
#define FE_TOWARDZERO _RC_CHOP

inline void fesetround(unsigned int choice) {
  _controlfp(choice, _MCW_RC);
}

#endif