append_files (Headers "h" . drivers curlnoise)
append_files (Sources "cpp" . drivers curlnoise)

add_executable (talpa ${Headers} ${Sources})
target_link_libraries (talpa ElTopo)

