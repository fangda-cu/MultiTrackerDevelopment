add_definitions(-D__BIG_ENDIAN__)
add_definitions(-DUSE_FORTRAN_BLAS)
add_definitions(-DNO_GUI)

set (Directories
 .
 common
 common/tunicate
 common/newsparse
 eltopo3d
)

append_files (Headers "h" ${Directories})
append_files (Sources "cpp" ${Directories})

find_package (OpenGL REQUIRED)
include_directories (${OPENGL_INCLUDE_DIR})

find_package (GLUT REQUIRED glut)
include_directories (${GLUT_INCLUDE_DIR})

add_library(ElTopo SHARED ${Headers} ${Sources})
target_link_libraries (ElTopo ${DEFAULT_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY})


