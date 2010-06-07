/**
 * \file TriangleMesh.hh
 *
 * \author smith@cs.columbia.edu
 * \date 03/15/2010
 */

#ifndef TRIANGLEMESH_HH
#define TRIANGLEMESH_HH

#ifdef WETA
#include "TopologicalObject/TopologicalObject.hh"
#else
#include "BASim/src/Core/TopologicalObject/TopologicalObject.hh"
#endif

namespace BASim {

/**
 * A triangle mesh.
 */
class TriangleMesh : public TopologicalObject
{
public:
  
  TriangleMesh()
  {
    add_property(m_vertex_positions, "vertex_positions", Vec3d(0,0,0));
  }
  
  /**
   * Access and modify the mesh's vertices.
   */
  Vec3d& getVertex( const TopologicalObject::vertex_handle& vh )
  {
    return property(m_vertex_positions)[vh];
  }
  
  /**
   * Const version.
   */
  const Vec3d& getVertex( const TopologicalObject::vertex_handle& vh ) const
  {
    return property(m_vertex_positions)[vh];
  }

private:
  VPropHandle<Vec3d> m_vertex_positions;
};

} // namespace BASim

#endif // TRIANGLEMESH_HH
