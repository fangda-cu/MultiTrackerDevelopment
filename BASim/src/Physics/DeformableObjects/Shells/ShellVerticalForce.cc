#include "BASim/src/Physics/DeformableObjects/Shells/ShellVerticalForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjIterators.hh"

namespace BASim {

ShellVerticalForce ::ShellVerticalForce ( ElasticShell& shell, const std::string& name, const Vec3d& gravity_vector, Scalar forcePerUnitArea) : 
    ElasticShellForce(shell, name), m_gravity(gravity_vector), m_strength(forcePerUnitArea)
{

}

Scalar ShellVerticalForce ::globalEnergy() const
{
  Scalar energy = 0;
  DeformableObject& obj = m_shell.getDefoObj();
  Vec3d vertical = m_gravity.normalized();
  for(FaceIterator fit = obj.faces_begin(); fit != obj.faces_end(); ++fit) {
    FaceHandle fh = *fit;
    Scalar area = m_shell.getArea(fh, false);
    Scalar totalForce = m_strength * area;

    //find the barycentre
    Vec3d barycentre(0,0,0);
    for(FaceVertexIterator fvit = obj.fv_iter(fh); fvit; ++fvit) {
      VertexHandle vh = *fvit;
      barycentre += m_shell.getDefoObj().getVertexPosition(vh);
    }
    barycentre /= 3;
    
    //add the resulting energy for the load.
    energy -= totalForce*barycentre.dot(vertical);
  }
  
  return energy;
}

void ShellVerticalForce ::globalForce( VecXd& force ) const
{
  
  DeformableObject& obj = m_shell.getDefoObj();
  Vec3d vertical = m_gravity.normalized();
  for(FaceIterator fit = obj.faces_begin(); fit != obj.faces_end(); ++fit) {
    FaceHandle fh = *fit;
    Scalar area = m_shell.getArea(fh, false); //use rest area? load shouldn't change if area changes...
    Scalar totalForce = m_strength * area;
    
    //divide up force to the vertices
    for(FaceVertexIterator fvit = obj.fv_iter(fh); fvit; ++fvit) {
      VertexHandle vh = *fvit;
      int dofIdx = m_shell.getDefoObj().getPositionDofBase(vh);
      force[dofIdx]   += vertical[0] * totalForce / 3;
      force[dofIdx+1] += vertical[1] * totalForce / 3;
      force[dofIdx+2] += vertical[2] * totalForce / 3;
    }
  }
  
}

void ShellVerticalForce ::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  return;
}

} //namespace BASim