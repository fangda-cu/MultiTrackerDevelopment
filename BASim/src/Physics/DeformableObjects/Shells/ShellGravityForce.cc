#include "BASim/src/Physics/DeformableObjects/Shells/ShellGravityForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjIterators.hh"

namespace BASim {

ShellGravityForce::ShellGravityForce( ElasticShell& shell, const std::string& name, const Vec3d& gravity_vector ) : 
    ElasticShellForce(shell, name), m_gravity(gravity_vector)
{

}

Scalar ShellGravityForce::globalEnergy() const
{
  Scalar energy = 0;
  DeformableObject& obj = m_shell.getDefoObj();
  for(VertexIterator vit = obj.vertices_begin(); vit != obj.vertices_end(); ++vit) {
    VertexHandle& vh = *vit;
    int dofIdx = m_shell.getDefoObj().getPositionDofBase(vh);
    Vec3d pos = m_shell.getDefoObj().getVertexPosition(vh);
    energy -= m_shell.getMass(vh)*m_gravity.dot(pos);
  }
  return energy;
}

void ShellGravityForce::globalForce( VecXd& force ) const
{
  
  DeformableObject& obj = m_shell.getDefoObj();
  for(VertexIterator vit = obj.vertices_begin(); vit != obj.vertices_end(); ++vit) {
    VertexHandle& vh = *vit;

    int dofIdx = m_shell.getDefoObj().getPositionDofBase(vh);
    force[dofIdx] += m_gravity[0] * m_shell.getMass(vh);
    force[dofIdx+1] += m_gravity[1] * m_shell.getMass(vh);
    force[dofIdx+2] += m_gravity[2] * m_shell.getMass(vh);
    
  }
  
}

void ShellGravityForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  //Jacobian is constant zero
  return;
}

} //namespace BASim