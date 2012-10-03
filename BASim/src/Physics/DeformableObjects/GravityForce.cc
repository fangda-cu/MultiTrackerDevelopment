#include "BASim/src/Physics/DeformableObjects/GravityForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjIterators.hh"
#include "BASim/src/Physics/DeformableObjects/DefoObjForce.hh"

namespace BASim {

GravityForce::GravityForce( DeformableObject& obj, Scalar timestep, const std::string& name, const Vec3d& gravity_vector ) : 
    DefoObjForce(obj, timestep, name), m_gravity(gravity_vector)
{

}

Scalar GravityForce::globalEnergy() 
{
  Scalar energy = 0;
  
  for(VertexIterator vit = m_obj.vertices_begin(); vit != m_obj.vertices_end(); ++vit) {
    VertexHandle& vh = *vit;
    int dofIdx = m_obj.getPositionDofBase(vh);
    Vec3d pos = m_obj.getVertexPosition(vh);
    energy -= m_obj.getVertexMass(vh)*m_gravity.dot(pos);
  }
  return energy;
}

void GravityForce::globalForce( VecXd& force ) 
{
  
  for(VertexIterator vit = m_obj.vertices_begin(); vit != m_obj.vertices_end(); ++vit) {
    VertexHandle& vh = *vit;
    int dofIdx = m_obj.getPositionDofBase(vh);
    force[dofIdx] += m_gravity[0] * m_obj.getVertexMass(vh);
    force[dofIdx+1] += m_gravity[1] * m_obj.getVertexMass(vh);
    force[dofIdx+2] += m_gravity[2] * m_obj.getVertexMass(vh);
    
  }
  
}

void GravityForce::globalJacobian( Scalar scale, MatrixBase& Jacobian )
{
  //Jacobian is constant zero
  return;
}

} //namespace BASim