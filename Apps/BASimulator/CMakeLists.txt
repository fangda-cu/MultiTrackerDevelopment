append_files (Headers "hh" . Problems)
append_files (Templates "tcc" . Problems)
append_files (Sources "cc" . Problems)

find_package (TCLAP REQUIRED)
include_directories (${TCLAP_INCLUDE_DIR})

#find_package (wxWidgets REQUIRED base core gl)
#include (${wxWidgets_USE_FILE})

add_executable (BASimulator ${Headers} ${Templates} ${Sources})
target_link_libraries (BASimulator BASim)
