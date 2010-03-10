append_files (Headers "hh" . Problems)
append_files (Templates "tcc" . Problems)
append_files (Sources "cc" . Problems)

find_package (TCLAP REQUIRED)
include_directories (${TCLAP_INCLUDE_DIR})

#find_package (wxWidgets REQUIRED base core gl)
#include (${wxWidgets_USE_FILE})

#find_package (FLENS REQUIRED)
#include_directories (${FLENS_INCLUDE_DIR})

#find_package (PETSc REQUIRED)

find_package (Maya REQUIRED)
include_directories (${MAYA_INCLUDE_DIR})

#find_package (CUDA REQUIRED)

#$find_package (OpenCurrent REQUIRED)
#$include_directories (${OCU_INCLUDE_DIR})

# bad, should have a findMKL.cmake
#include_directories ($ENV{MKLROOT}/include)

###########################################################################################
#
# I'm still working on the dependencies for MKL etc, I'll tidy this file up as I figure
# out how to get it running on my workstation and my laptop.
#
###########################################################################################
if (CMAKE_SYSTEM_NAME MATCHES Linux)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MAYA_CXX_FLAGS}")
endif (CMAKE_SYSTEM_NAME MATCHES Linux)

ADD_LIBRARY(wmBunsen SHARED ${Headers} ${Templates} ${Sources})

#######
## We absolutely must link with -Bsymbolic or Maya will crash when the plugin is running!
#if (USE_INTEL_COMPILER)
  SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -shared -fPIC -Wl,-Bsymbolic -pthread -lGLU -openmp -pthread")
#else (USE_INTEL_COMPILER)
#  SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -shared -fPIC -Wl,-Bsymbolic -pthread -lGLU -openmp")
#endif (USE_INTEL_COMPILER)

if (USE_INTEL_COMPILER)
	target_link_libraries (wmBunsen ${MAYA_LIBRARIES} BASim)
else (USE_INTEL_COMPILER)
	target_link_libraries (wmBunsen ${MAYA_LIBRARIES} -fopenmp -lpthread BASim)
endif (USE_INTEL_COMPILER)
