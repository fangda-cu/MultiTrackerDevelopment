cmake_minimum_required (VERSION 2.6)

project (BASim)
set (BASIM_MAJOR_VERSION "0")
set (BASIM_MINOR_VERSION "1")
set (BASIM_BUGFIX_VERSION "0")
set (BASIM_VERSION_NUMBER "${BASIM_MAJOR_VERSION}.${BASIM_MINOR_VERSION}.${BASIM_BUGFIX_VERSION}")
set (BASIM_VERSION "${BASIM_VERSION_NUMBER}")
add_definitions (-DBASIM_VERSION="${BASIM_VERSION}")

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Debug CACHE STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# add directory with macros
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# add cmake file
include (BASimCommon)

# set options
option (BUILD_APPS "Build applications" ON)
option (TEST_ROD_STRETCHING_FORCE "Tests RodStretchingForce" OFF)
option (TEST_ROD_TWISTING_FORCE "Tests RodTwistingForce" OFF)
option (TEST_ROD_BENDING_FORCE "Tests RodBendingForce" OFF)
option (BUILD_MAYA_NODES "Build Maya nodes" ON)

# add some extra warnings
if (CMAKE_COMPILER_IS_GNUCXX)
  if (CMAKE_SYSTEM_NAME MATCHES Linux)
    set (CMAKE_C_COMPILER gcc-4.1)
    set ( CMAKE_CXX_COMPILER  g++-4.1)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnon-virtual-dtor -ansi -Wcast-align -Wall -Wextra -Wpointer-arith -Wwrite-strings -Wformat-security -fno-common -Wno-unused")
    if (CMAKE_BUILD_TYPE MATCHES Release)
      add_definitions (-DNDEBUG)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mfpmath=sse -msse2 -fno-check-new -fstrict-aliasing") # -fno-exceptions
    endif (CMAKE_BUILD_TYPE MATCHES Release)
  endif (CMAKE_SYSTEM_NAME MATCHES Linux)
endif (CMAKE_COMPILER_IS_GNUCXX)

if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions (-DDEBUG)
endif (CMAKE_BUILD_TYPE MATCHES Debug)

# Eigen library is required
find_package (Eigen REQUIRED)
include_directories (${EIGEN_INCLUDE_DIR})

#find_package (FLENS REQUIRED)
#include_directories (${FLENS_INCLUDE_DIR})

find_package (PETSc REQUIRED)
include_directories (${PETSC_INCLUDE_PATH})

find_package (PNG)
if (PNG_FOUND)
  add_definitions (${PNG_DEFINITIONS})
  include_directories (${PNG_INCLUDE_DIR})
else (PNG_FOUND)
  add_definitions (-DNO_LIBPNG)
endif (PNG_FOUND)

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory (BASim)

if (BUILD_APPS)
  add_subdirectory (Apps)
endif (BUILD_APPS)

add_subdirectory (doc)
