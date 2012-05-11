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
  //TODO: Given a plane as input, we can use the distance to the plane to determine the gravitational potential. Yes?
  return 0;
}

void ShellGravityForce::globalForce( VecXd& force ) const
{
  
  DeformableObject& obj = m_shell.getDefoObj();
  for(VertexIterator vit = obj.vertices_begin(); vit != obj.vertices_end(); ++vit) {
    VertexHandle& vh = *vit;
////////////////////    int dofIdx = m_shell.getVertexDofBase(vh);
    int dofIdx = m_shell.getDefoObj().getPositionDofBase(vh);
    force[dofIdx] += m_gravity[0] * m_shell.getMass(vh);
    force[dofIdx+1] += m_gravity[1] * m_shell.getMass(vh);
    force[dofIdx+2] += m_gravity[2] * m_shell.getMass(vh);
    
  }
  
}

void ShellGravityForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  //Assuming explicit integration of gravity for now, which I believe avoids the need for the Jacobian.
  return;
}

} //namespace BASim