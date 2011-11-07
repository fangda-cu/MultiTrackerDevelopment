#include "BASim/src/Physics/DeformableObjects/Shells/ShellRadialForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjIterators.hh"

namespace BASim {

ShellRadialForce::ShellRadialForce( ElasticShell& shell, const std::string& name, const Vec3d& centre, Scalar strength ) : 
    ElasticShellForce(shell, name), m_centre(centre), m_strength(strength)
{

}

Scalar ShellRadialForce::globalEnergy() const
{
  
  return 0;
}

void ShellRadialForce::globalForce( VecXd& force ) const
{
  
  DeformableObject& obj = m_shell.getDefoObj();
  /*
  for(FaceIterator fit = obj.faces_begin(); fit != obj.faces_end(); ++fit) {
    FaceHandle fh = *fit;
    FaceVertexIterator fvit = obj.fv_iter(fh);
    VertexHandle vhs[3];
    Vec3d vpos[3];
    for(int i = 0;fvit;++fvit,++i) {
      vhs[i] = *fvit; 
      vpos[i] = m_shell.getVertexPosition(vhs[i]);
    }
    
    //compute a constant outward force scaled by area
    
    //Use triangle normals as the force direction
    Vec3d triForce = -0.5 * m_strength * (vpos[1]-vpos[0]).cross(vpos[2]-vpos[0]) / 3.0;
    
    //Use barycentre to sphere centre as the force direction.
    //Vec3d direction = (vpos[0]+vpos[1]+vpos[2])/3-m_centre;
    //direction.normalize();
    //Scalar area = 0.5*((vpos[1]-vpos[0]).cross(vpos[2]-vpos[0])).norm();
    //Vec3d triForce = m_strength*direction*area / 3;

    for(int i = 0; i < 3; ++i) {
      int dofIdx = m_shell.getVertexDofBase(vhs[i]);
      force[dofIdx] += triForce[0];
      force[dofIdx+1] += triForce[1];
      force[dofIdx+2] += triForce[2];
    }

  }
  */
  int c = 0;
  for(VertexIterator vit = obj.vertices_begin(); vit != obj.vertices_end(); ++vit) {
    VertexHandle& vh = *vit;
    int dofIdx = m_shell.getVertexDofBase(vh);
    
    //compute r^2
    Vec3d direction = m_shell.getVertexPosition(vh) - m_centre;
    Scalar r_len2 = direction.squaredNorm();
    
    direction.normalize();
    direction *= m_strength/r_len2;
    
    //visit all triangles and accumulate: ref_area / area 
    Scalar curArea = 0, refArea = 0;
    for(VertexFaceIterator vfit = obj.vf_iter(vh); vfit; ++vfit) {
      FaceHandle fh = *vfit;
      curArea += m_shell.getArea(fh, true)/3.0;
      refArea += m_shell.getArea(fh, false)/3.0;
    }
    
    direction *= refArea;
    force[dofIdx] += direction[0];
    force[dofIdx+1] += direction[1];
    force[dofIdx+2] += direction[2];
  }
  
}

void ShellRadialForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  //Assuming explicit integration of gravity for now, which I believe avoids the need for the Jacobian.
  return;
}

} //namespace BASim