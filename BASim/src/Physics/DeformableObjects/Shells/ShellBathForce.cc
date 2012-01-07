#include "BASim/src/Physics/DeformableObjects/Shells/ShellBathForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjIterators.hh"
#include "BASim/src/Math/MatrixBase.hh"

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

    Scalar pressure = m_density * m_gravity[1] * (m_bath_height - barycentre[1]);
    Vec3d areaNormal = 0.5*(pos1-pos0).cross(pos2-pos0) / 3;
          
    force[dofIdx0]   += pressure*areaNormal[0];
    force[dofIdx0+1] += pressure*areaNormal[1];
    force[dofIdx0+2] += pressure*areaNormal[2];
    
    force[dofIdx1]   += pressure*areaNormal[0];
    force[dofIdx1+1] += pressure*areaNormal[1];
    force[dofIdx1+2] += pressure*areaNormal[2];

    force[dofIdx2]   += pressure*areaNormal[0];
    force[dofIdx2+1] += pressure*areaNormal[1];
    force[dofIdx2+2] += pressure*areaNormal[2];

  }
  
}

void ShellBathForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
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
    
    std::vector<int> indices;
    indices.push_back(dofIdx0); indices.push_back(dofIdx0+1); indices.push_back(dofIdx0+2);
    indices.push_back(dofIdx1); indices.push_back(dofIdx1+1); indices.push_back(dofIdx1+2);
    indices.push_back(dofIdx2); indices.push_back(dofIdx2+1); indices.push_back(dofIdx2+2);

    Vec3d pos0 = m_shell.getVertexPosition(v0);
    Vec3d pos1 = m_shell.getVertexPosition(v1);
    Vec3d pos2 = m_shell.getVertexPosition(v2);

    //Values
    Scalar pressure = m_density * m_gravity[1] * (m_bath_height - (pos0[1]+pos1[1]+pos2[1])/3.0);
    Vec3d areaNormal = 0.5*(pos1-pos0).cross(pos2-pos0) / 3;
 
    //jac = dF/dx = d(pAn)/dx = (dp/dx)*An + p*d(An)/dx = part1 + part2
    Eigen::Matrix<Scalar,9,9> part1, part2, jac;

    //Part 1
    Eigen::Matrix<Scalar,9,1> dpdx;
    Eigen::Matrix<Scalar,1,9> areaNormals;
    Scalar dpdh = -m_density*m_gravity[1]/3.0;
    
    dpdx << 0, dpdh, 0, 
            0, dpdh, 0, 
            0, dpdh, 0;
    areaNormals << areaNormal[0], areaNormal[1], areaNormal[2], 
                   areaNormal[0], areaNormal[1], areaNormal[2], 
                   areaNormal[0], areaNormal[1], areaNormal[2];

    part1 = dpdx*areaNormals; 
    
    //Part 2
    //These derivatives leave out the 1/6 factor in the areaNormal term above. 
    //I multiply it back in below.
    Eigen::Matrix<Scalar,3,3> dAn_db, dAn_dc, dAn_da;
    dAn_db <<                0, +pos2[2]-pos0[2],    +pos2[1]-pos0[1],
      -pos2[2]+pos0[2],                0,    -pos2[0]+pos0[0],
      -pos2[1]+pos0[1],  pos2[0]-pos0[0],                   0;

    dAn_dc <<        0, +pos1[2]-pos0[2],    +pos1[1]-pos0[1],
      -pos1[2]+pos0[2],                0,    -pos1[0]+pos0[0],
      -pos1[1]+pos0[1], +pos1[0]-pos0[0],                   0;

    dAn_da = -(dAn_db+dAn_dc);
    
    //Build this 3x9 chunk of the matrix
    part2.block(0,0, 3,3) += pressure * dAn_da / 6;
    part2.block(0,3, 3,3) += pressure * dAn_db / 6;
    part2.block(0,6, 3,3) += pressure * dAn_dc / 6;

    //copy it downwards, since the effect on each particle is the same
    part2.block(3,0, 3,9) = part2.block(0,0,3,9);
    part2.block(6,0, 3,9) = part2.block(0,0,3,9);

    //Combine it
    jac = part1 + part2; 

    //Scatter back into the main matrix
    for(int i = 0; i < 9; ++i) {
      int ind0 = indices[i];
      for(int j = 0; j < 9; ++j) {
        int ind1 = indices[j];
        Jacobian.add(ind0,ind1, scale*jac(i,j));
      }
    }
    
  }
  
}

} //namespace BASim