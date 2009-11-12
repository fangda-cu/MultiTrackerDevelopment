append_files (Headers "hh" . Problems)
append_files (Templates "tcc" . Problems)
append_files (Sources "cc" . Problems)

find_package (OpenGL REQUIRED)
include_directories (${OPENGL_INCLUDE_DIR})

find_package (GLUT REQUIRED glut)
include_directories (${GLUT_INCLUDE_DIR})

find_package (TCLAP REQUIRED)
include_directories (${TCLAP_INCLUDE_DIR})

#find_package (wxWidgets REQUIRED base core gl)
#include (${wxWidgets_USE_FILE})

#find_package (FLENS REQUIRED)
#include_directories (${FLENS_INCLUDE_DIR})

find_package (PETSc REQUIRED)

find_package (Maya REQUIRED)
include_directories (${MAYA_INCLUDE_DIR})

if (USE_INTEL_COMPILER)
    find_package (MKL REQUIRED)
    include_directories (${MKL_INCLUDE_DIR})
endif (USE_INTEL_COMPILER)

find_package (CUDA REQUIRED)

find_package (OpenCurrent REQUIRED)
include_directories (${OCU_INCLUDE_DIR})

#include_directories ("/usr/include/c++/4.1.2/")

# bad, should have a findMKL.cmake
include_directories (${MKLROOT}/include)
include_directories (${TBB_INSTALL_DIR}/include)

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
SET (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -shared -fPIC -Wl,-Bsymbolic -pthread -lGLU -openmp")

target_link_libraries (wmBunsen BASim ${MAYA_LIBRARIES} ${GLUT_glut_LIBRARY} ${OPENMESH_LIBRARY} ${SMV_LIBRARY} ${wxWidgets_LIBRARIES} ${FLENS_LIBRARIES} ${PETSC_LIBRARIES}  $ENV{MKLROOT}/lib/em64t/libmkl_solver_lp64.a -Wl,--start-group $ENV{MKLROOT}/lib/em64t/libmkl_intel_lp64.a $ENV{MKLROOT}/lib/em64t/libmkl_intel_thread.a $ENV{MKLROOT}/lib/em64t/libmkl_core.a -Wl,--end-group -openmp -lpthread)

