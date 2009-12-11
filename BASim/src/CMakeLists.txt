find_package (OpenGL REQUIRED)
include_directories (${OPENGL_INCLUDE_DIR})

find_package (GLUT REQUIRED glut)
include_directories (${GLUT_INCLUDE_DIR})

set (Directories
  .
  Core
  Core/TopologicalObject
  Math
  Physics
  Physics/ElasticRods
  Render
  Util
  Collisions
)

include_directories (${CMAKE_CURRENT_SOURCE_DIR})

#find_package (FLENS REQUIRED)
#include_directories (${FLENS_INCLUDE_DIR})

append_files (Headers "hh" ${Directories})
append_files (Inlines "inl" ${Directories})
append_files (Templates "tcc" ${Directories})
append_files (Sources "cc" ${Directories})

if (PNG_FOUND)
  set (Headers ${Headers} Render/YImage/YImage.hh)
  set (Sources ${Sources} Render/YImage/YImage.cc)
endif (PNG_FOUND)

if (TEST_ROD_STRETCHING_FORCE)
  add_definitions (-DTEST_ROD_STRETCHING)
  set (Headers ${Headers} Physics/ElasticRods/Tests/RodStretchingTest.hh)
  set (Inlines ${Inlines} Physics/ElasticRods/Tests/stretchingEnergy.inl)
  set (Inlines ${Inlines} Physics/ElasticRods/Tests/stretchingForce.inl)
  set (Inlines ${Inlines} Physics/ElasticRods/Tests/stretchingJacobian.inl)
endif (TEST_ROD_STRETCHING_FORCE)

if (TEST_ROD_TWISTING_FORCE)
  add_definitions (-DTEST_ROD_TWISTING)
  set (Headers ${Headers} Physics/ElasticRods/Tests/RodTwistingTest.hh)
  set (Inlines ${Inlines} Physics/ElasticRods/Tests/twistingEnergy.inl)
  set (Inlines ${Inlines} Physics/ElasticRods/Tests/twistingForce.inl)
  set (Inlines ${Inlines} Physics/ElasticRods/Tests/twistingJacobian.inl)
endif (TEST_ROD_TWISTING_FORCE)

if (TEST_ROD_BENDING_FORCE)
  add_definitions (-DTEST_ROD_BENDING)
  set (Headers ${Headers} Physics/ElasticRods/Tests/RodBendingTest.hh)
  set (Inlines ${Inlines} Physics/ElasticRods/Tests/bendingEnergy.inl)
  set (Inlines ${Inlines} Physics/ElasticRods/Tests/bendingForce.inl)
  set (Inlines ${Inlines} Physics/ElasticRods/Tests/bendingJacobian.inl)
  set (Sources ${Sources} Physics/ElasticRods/Tests/RodBendingTest.cc)
endif (TEST_ROD_BENDING_FORCE)

add_library (BASim SHARED ${Headers} ${Inlines} ${Templates} ${Sources})
target_link_libraries (BASim ${PETSC_LIBRARIES} ${PNG_LIBRARIES} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY})
