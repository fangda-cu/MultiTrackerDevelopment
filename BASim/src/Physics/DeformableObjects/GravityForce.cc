#include "BASim/src/Physics/DeformableObjects/GravityForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjIterators.hh"

namespace BASim {

GravityForce::GravityForce( PositionDofsModel& obj, const std::string& name, const Vec3d& gravity_vector ) : 
    DefoObjForce(obj, name), m_gravity(gravity_vector)
{

}

Scalar GravityForce::globalEnergy() const
{
  Scalar energy = 0;
  DeformableObject& obj = m_model.getDefoObj();
  for(VertexIterator vit = obj.vertices_begin(); vit != obj.vertices_end(); ++vit) {
    VertexHandle& vh = *vit;
    int dofIdx = obj.getPositionDofBase(vh);
    Vec3d pos = obj.getVertexPosition(vh);
    energy -= m_model.getMass(vh)*m_gravity.dot(pos);
  }
  return energy;
}

void GravityForce::globalForce( VecXd& force ) const
{
  
  DeformableObject& obj = m_model.getDefoObj();
  for(VertexIterator vit = obj.vertices_begin(); vit != obj.vertices_end(); ++vit) {
    VertexHandle& vh = *vit;
    int dofIdx = obj.getPositionDofBase(vh);
    force[dofIdx] += m_gravity[0] * m_model.getMass(vh);
    force[dofIdx+1] += m_gravity[1] * m_model.getMass(vh);
    force[dofIdx+2] += m_gravity[2] * m_model.getMass(vh);
    
  }
  
}

void GravityForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  //Jacobian is constant zero
  return;
}

} //namespace BASim