#ifndef TETMESHIO_H
#define TETMESHIO_H

#include <vec.h>

bool load_ele_file( const char* ele_file_name, std::vector<ElTopo::Vec4st>& tets );
bool write_ele_file( const char* ele_file_name, std::vector<ElTopo::Vec4st>& tets );
bool load_node_file( const char* node_file_name, std::vector<ElTopo::Vec3f>& xs );
bool write_node_file( const char* node_file_name, std::vector<ElTopo::Vec3f>& xs );

#endif
