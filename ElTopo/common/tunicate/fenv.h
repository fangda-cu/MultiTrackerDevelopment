#ifndef FENV_H
#define FENV_H
#include <float.h>
#include <stdio.h>
#include <errno.h>

#pragma fenv_access(on)


//define the functions we need. I hate Microsoft! Just implement C99 already, c'mon.
#define FE_DOWNWARD   _RC_DOWN
#define FE_UPWARD     _RC_UP
#define FE_TONEAREST  _RC_NEAR
#define FE_TOWARDZERO _RC_CHOP

inline void fesetround(unsigned int choice) {
  _controlfp(choice, _MCW_RC);
}

#endif