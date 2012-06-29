#include "BASim/src/Physics/DeformableObjects/Shells/ShellPointForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjIterators.hh"

namespace BASim {

ShellPointForce ::ShellPointForce( ElasticShell& shell, const std::string& name, const std::vector<VertexHandle>& verts, const std::vector<Vec3d>& forces) : 
    ElasticShellForce(shell, name), m_vertices(verts), m_forces(forces)
{

}

Scalar ShellPointForce  ::globalEnergy() const
{
  return 0;
}

void ShellPointForce  ::globalForce( VecXd& force ) const
{
  
  DeformableObject& obj = m_shell.getDefoObj();
  for(unsigned int i= 0; i < m_forces.size(); ++i) {
    VertexHandle vh = m_vertices[i];
    int dofIdx = m_shell.getDefoObj().getPositionDofBase(vh);
    force[dofIdx]   += m_forces[i][0];
    force[dofIdx+1] += m_forces[i][1];
    force[dofIdx+2] += m_forces[i][2];
  }
  
}

void ShellPointForce  ::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  return;
}

} //namespace BASim