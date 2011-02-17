/**
 * \file MinimalTriangleMeshBackup.cc
 *
 * \author smith@cs.columbia.edu
 * \date 07/21/2010
 */

#include "MinimalTriangleMeshBackup.hh"

namespace BASim 
{
  
void MinimalTriangleMeshBackup::resize( const TriangleMesh& mesh )
{
  m_vp.resize(mesh.nv());
}

void MinimalTriangleMeshBackup::backupMesh( TriangleMesh& mesh )
{
  for( TriangleMesh::vertex_iter i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i )
  {
    m_vp[i->idx()] = mesh.getVertex(*i);
  }
}

void MinimalTriangleMeshBackup::restoreMesh( TriangleMesh& mesh )
{
  for( TriangleMesh::vertex_iter i = mesh.vertices_begin(); i != mesh.vertices_end(); ++i )
  {
    mesh.getVertex(*i) = m_vp[i->idx()];
  }
}

void MinimalTriangleMeshBackup::clear()
{
  m_vp.clear();
}

}