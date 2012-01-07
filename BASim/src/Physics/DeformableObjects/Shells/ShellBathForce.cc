#include "BASim/src/Physics/DeformableObjects/Shells/ShellBathForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjIterators.hh"

namespace BASim {

ShellBathForce::ShellBathForce( ElasticShell& shell, const std::string& name, const Vec3d& gravity_vector, Scalar density, Scalar bathHeight ) : 
    ElasticShellForce(shell, name), m_gravity(gravity_vector), m_density(density), m_bath_height(bathHeight)
{

}

Scalar ShellBathForce::globalEnergy() const
{
  return 0;
}

void ShellBathForce::globalForce( VecXd& force ) const
{
  
  DeformableObject& obj = m_shell.getDefoObj();
  for(FaceIterator fit = obj.faces_begin(); fit != obj.faces_end(); ++fit) {
    FaceHandle& fh = *fit;
    FaceVertexIterator fvit = obj.fv_iter(fh);
    VertexHandle v0 = *fvit; ++fvit;
    VertexHandle v1 = *fvit; ++fvit;
    VertexHandle v2 = *fvit;


    int dofIdx0 = m_shell.getVertexDofBase(v0);
    int dofIdx1 = m_shell.getVertexDofBase(v1);
    int dofIdx2 = m_shell.getVertexDofBase(v2);

    Vec3d pos0 = m_shell.getVertexPosition(v0);
    Vec3d pos1 = m_shell.getVertexPosition(v1);
    Vec3d pos2 = m_shell.getVertexPosition(v2);
    
    Vec3d barycentre = (pos0+pos1+pos2)/3.0;
    Scalar area = 0.5*((pos1-pos0).cross(pos2-pos0)).norm();

    Scalar pressure = (barycentre[1] - m_bath_height);
    force[dofIdx0+1] += m_density * m_gravity[1] * pressure * area / 3;
    force[dofIdx1+1] += m_density * m_gravity[1] * pressure * area / 3;
    force[dofIdx2+1] += m_density * m_gravity[1] * pressure * area / 3;
    
  }
  
}

void ShellBathForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  return;
}

} //namespace BASim