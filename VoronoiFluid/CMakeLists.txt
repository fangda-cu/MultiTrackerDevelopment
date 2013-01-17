append_files (Headers "h" . sparse)
append_files (Sources "cpp" . sparse)

find_package (OpenGL REQUIRED)
include_directories (${OPENGL_INCLUDE_DIR})

find_package (CGAL REQUIRED)
include_directories (${CGAL_INCLUDE_DIR})

find_package (GLUT REQUIRED glut)
include_directories (${GLUT_INCLUDE_DIR})

include_directories (${EIGEN_INCLUDE_DIR})

add_executable (VoronoiFluid ${Headers} ${Templates} ${Sources})
target_link_libraries (VoronoiFluid ${CGAL_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY})

