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

add_executable (BASimulator ${Headers} ${Templates} ${Sources})
target_link_libraries (BASimulator BASim ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY} ${OPENMESH_LIBRARY} ${SMV_LIBRARY} ${wxWidgets_LIBRARIES} ${FLENS_LIBRARIES} ${PETSC_LIBRARIES})
