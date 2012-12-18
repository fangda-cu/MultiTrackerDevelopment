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
//    Scalar m = m_obj.getVertexMass(vh);
    Scalar m = 1;
    energy -= m*m_gravity.dot(pos);
  }
  return energy;
}

void GravityForce::globalForce( VecXd& force ) 
{
  
  for(VertexIterator vit = m_obj.vertices_begin(); vit != m_obj.vertices_end(); ++vit) {
    VertexHandle& vh = *vit;
    int dofIdx = m_obj.getPositionDofBase(vh);
//    Scalar m = m_obj.getVertexMass(vh);
    Scalar m = 1;
    force[dofIdx] += m_gravity[0] * m;
    force[dofIdx+1] += m_gravity[1] * m;
    force[dofIdx+2] += m_gravity[2] * m;
    
  }
  
}

void GravityForce::globalJacobian( Scalar scale, MatrixBase& Jacobian )
{
  //Jacobian is constant zero
  return;
}

} //namespace BASim