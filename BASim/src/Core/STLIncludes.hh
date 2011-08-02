/**
 * \file STLIncludes.hh
 *
 * \author miklos@cs.columbia.edu
 * \date 09/15/2009
 */

#ifndef STLINCLUDES_HH
#define STLINCLUDES_HH

#include <algorithm>
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#ifndef _MSC_VER
#include <sys/resource.h>
#include <sys/time.h>
#else
typedef unsigned int uint;
typedef unsigned long int uint32_t;
//int isnan(double x) { return _isnan(x); }
//int isinf(double x) { return !_finite(x); }
#define isnan(x) _isnan(x)
#define isinf(x) !_finite(x)
#endif

#endif // STLINCLUDES_HH
