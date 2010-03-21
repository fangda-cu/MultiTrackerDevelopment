cmake_minimum_required (VERSION 2.6)

project (BASim)
set (BASIM_MAJOR_VERSION "0")
set (BASIM_MINOR_VERSION "1")
set (BASIM_BUGFIX_VERSION "0")
set (BASIM_VERSION_NUMBER "${BASIM_MAJOR_VERSION}.${BASIM_MINOR_VERSION}.${BASIM_BUGFIX_VERSION}")
set (BASIM_VERSION "${BASIM_VERSION_NUMBER}")
add_definitions (-DBASIM_VERSION="${BASIM_VERSION}")

if (NOT CMAKE_BUILD_TYPE)
  set (CMAKE_BUILD_TYPE Release CACHE STRING
    "Choose the type of build, options are: None(CMAKE_CXX_FLAGS or CMAKE_C_FLAGS used) Debug Release RelWithDebInfo MinSizeRel."
    FORCE)
endif (NOT CMAKE_BUILD_TYPE)

# add directory with macros
set (CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${PROJECT_SOURCE_DIR}/cmake)

# add cmake file
include (BASimCommon)

# set options
option (USE_INTEL_COMPILER "Use Intel Compiler" OFF)
option (BUILD_APPS "Build applications" ON)
option (TEST_ROD_STRETCHING_FORCE "Tests RodStretchingForce" OFF)
option (TEST_ROD_TWISTING_FORCE "Tests RodTwistingForce" OFF)
option (TEST_ROD_BENDING_FORCE "Tests RodBendingForce" OFF)
option (TIMING_ON "Enables timing macros" OFF)
if (TIMING_ON)
  add_definitions (-DTIMING_ON)
endif (TIMING_ON)
option (BUILD_MAYA_NODES "Build Maya nodes" ON)

# Deal with the intel compiler vs gcc
if (CMAKE_COMPILER_IS_GNUCXX)
  if (CMAKE_SYSTEM_NAME MATCHES Linux)
    if (USE_INTEL_COMPILER)
      set (CMAKE_C_COMPILER icc)
      set (CMAKE_CXX_COMPILER  icpc)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSING_INTEL_COMPILER -openmp -pthread")
      # setup some things and disable a bunch of warnings because Intel is painfully strict
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -xhost -ipo -parallel -openmp -DCONST=const -w1 -diag-disable vec -diag-disable ipo -vec-report0 -opt-report 0 -wd383 -wd981 -wd1572 -wd1684 -wd2012 -wd2017 -wd2304 -std=c++0x") 
    else (USE_INTEL_COMPILER)
      #set (CMAKE_C_COMPILER gcc-4.1)
      #set ( CMAKE_CXX_COMPILER  g++)
    endif (USE_INTEL_COMPILER)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnon-virtual-dtor -ansi -Wcast-align -Wall -Wextra -Wpointer-arith -Wwrite-strings -Wformat-security -fno-common -Wno-unused -Wno-empty-body -Wno-deprecated -fopenmp -pthread")
    if (CMAKE_BUILD_TYPE MATCHES Release)
      add_definitions (-DNDEBUG)
      set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -mfpmath=sse -msse2 -fno-check-new -fstrict-aliasing") # -fno-exceptions
    endif (CMAKE_BUILD_TYPE MATCHES Release)
  endif (CMAKE_SYSTEM_NAME MATCHES Linux)
endif (CMAKE_COMPILER_IS_GNUCXX)

if (CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions (-DDEBUG)
endif (CMAKE_BUILD_TYPE MATCHES Debug)

# Eigen library is required
find_package (Eigen REQUIRED)
include_directories (${EIGEN_INCLUDE_DIR})

option (USE_PNG "Builds in support for exporting to png" ON)
if (USE_PNG)
  find_package (PNG)
  if (PNG_FOUND)
    add_definitions (-DHAVE_PNG)
    add_definitions (${PNG_DEFINITIONS})
    include_directories (${PNG_INCLUDE_DIR})
    set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${PNG_LIBRARIES})
  endif (PNG_FOUND)
endif (USE_PNG)

option (USE_PETSc "Builds in support for PETSc" OFF)
if (USE_PETSc)
  find_package (PETSc)
  if (PETSC_FOUND)
    add_definitions (-DHAVE_PETSC)
    include_directories (${PETSC_INCLUDE_PATH})
    set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${PETSC_LIBRARIES})
  endif (PETSC_FOUND)
endif (USE_PETSc)

#find_package (FLENS REQUIRED)
#include_directories (${FLENS_INCLUDE_DIR})

option (USE_MKL "Builds in support for Intel's MKL" OFF)
if (USE_MKL)
  find_package (MKL)
  if (MKL_FOUND)
    add_definitions (-DHAVE_LAPACK)
    add_definitions (-DHAVE_MKL)
    include_directories (${MKL_INCLUDES})
    set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${MKL_LIBRARIES})
  endif (MKL_FOUND)
endif (USE_MKL)

option (USE_LAPACK "Builds in support for generic LAPACK" ON)
if (USE_LAPACK)
  find_package (LAPACK)
  if (LAPACK_FOUND)
    add_definitions (-DHAVE_LAPACK)
    set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${LAPACK_LIBRARIES})
  endif(LAPACK_FOUND)
endif (USE_LAPACK)

include_directories (${CMAKE_CURRENT_SOURCE_DIR})
add_subdirectory (BASim)

if (BUILD_APPS)
  add_subdirectory (Apps)
endif (BUILD_APPS)

add_subdirectory (doc)
