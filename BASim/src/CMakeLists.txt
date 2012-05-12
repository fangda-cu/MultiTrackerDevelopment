find_package (OpenGL REQUIRED)
include_directories (${OPENGL_INCLUDE_DIR})
set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${OPENGL_LIBRARIES})

find_package (GLUT REQUIRED glut)
include_directories (${GLUT_INCLUDE_DIR})
set (DEFAULT_LIBRARIES ${DEFAULT_LIBRARIES} ${GLUT_glut_LIBRARY})

set (Directories
  .
  Core
  Core/TopologicalObject
  Collisions
  Collisions/ElTopo
  IO
  Math
  Math/Eigen
  Physics
  Physics/ElasticRods
  Physics/DeformableObjects
  Physics/DeformableObjects/Shells
  Physics/DeformableObjects/Rods
  Render
  Util
  Threads
)

include_directories (${CMAKE_CURRENT_SOURCE_DIR})

append_files (Headers "hh" ${Directories})
append_files (Inlines "inl" ${Directories})
append_files (Templates "tcc" ${Directories})
append_files (Sources "cc" ${Directories})

if (PNG_FOUND)
  set (Headers ${Headers} Render/YImage/YImage.hh)
  set (Sources ${Sources} Render/YImage/YImage.cc)
endif (PNG_FOUND)

if (PETSC_FOUND)
  append_files (Headers "hh" Math/Petsc)
  append_files (Inlines "inl" Math/Petsc)
  append_files (Templates "tcc" Math/Petsc)
  append_files (Sources "cc" Math/Petsc)
endif (PETSC_FOUND)

if (MKL_FOUND)
  append_files (Headers "hh" Math/MKL)
  append_files (Inlines "inl" Math/MKL)
  append_files (Templates "tcc" Math/MKL)
  append_files (Sources "cc" Math/MKL)
endif (MKL_FOUND)

if (LAPACK_FOUND)
  append_files (Headers "hh" Math/MKL)
  append_files (Inlines "inl" Math/MKL)
  append_files (Templates "tcc" Math/MKL)
  append_files (Sources "cc" Math/MKL)
endif (LAPACK_FOUND)

if (USE_PARDISO)
  append_files (Headers "hh" Math/Pardiso)
  append_files (Inlines "inl" Math/Pardiso)
  append_files (Templates "tcc" Math/Pardiso)
  append_files (Sources "cc" Math/Pardiso)
endif (USE_PARDISO)

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

# create source groups for the IDE support
set (CREATE_SOURCE_GROUPS true)
if (CREATE_SOURCE_GROUPS)
  source_group ("Others" FILES ${Headers} ${Inlines} ${Templates} ${Sources})

  foreach (d ${Directories})
    string (COMPARE EQUAL ${d} "." currentdir)
    file (GLOB files "${d}/*")
    if (${currentdir})
      set (groupname "Sources")
    else ()
      string (REPLACE / \\ groupname ${d})
    endif ()
    source_group (${groupname} FILES ${files})
  endforeach ()

  # top level headers
  file (GLOB TopLevelHeaders "../*")
  foreach (d ${TopLevelHeaders})
    string (REGEX MATCH ^.*/src$ regexmatch1 ${d})
    string (REGEX MATCH ^.*/CMakeLists.txt$ regexmatch2 ${d})
    if (regexmatch1 OR regexmatch2)
      list (REMOVE_ITEM TopLevelHeaders ${d})
    endif ()
  endforeach ()
  source_group ("TopLevel" FILES ${TopLevelHeaders})

endif (CREATE_SOURCE_GROUPS)

# add files in alphabetical order helps when generating IDE projects
set (Files ${Headers} ${Inlines} ${Templates} ${Sources} ${TopLevelHeaders})
list (SORT Files)

add_library (BASim SHARED ${Files})
target_link_libraries (BASim ElTopo ${DEFAULT_LIBRARIES})
