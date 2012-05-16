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

add_library(ElTopo SHARED ${Headers} ${Sources})
target_link_libraries (ElTopo ${DEFAULT_LIBRARIES} "-framework vecLib")

