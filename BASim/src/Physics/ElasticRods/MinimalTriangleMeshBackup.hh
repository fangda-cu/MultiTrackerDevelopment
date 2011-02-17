/**
 * \file MinimalTriangleMeshBackup.hh
 *
 * \author smith@cs.columbia.edu
 * \date 07/21/2010
 */


#ifndef MINIMALTRIANGLEMESHBACKUP_HH
#define MINIMALTRIANGLEMESHBACKUP_HH

#include <fstream>

#ifdef WETA
#include "../../Core/TriangleMesh.hh"
#else
#include "BASim/src/Core/TriangleMesh.hh"
#endif

namespace BASim 
{
  
  class MinimalTriangleMeshBackup
  {
  public:
    void resize( const TriangleMesh& mesh );

    void backupMesh( TriangleMesh& mesh );

    void restoreMesh( TriangleMesh& mesh );

    void clear();

  private:
    std::vector<Vec3d> m_vp;
  };
  
} // namespace BASim

#endif // MINIMALRODSTATEBACKUP_HH

