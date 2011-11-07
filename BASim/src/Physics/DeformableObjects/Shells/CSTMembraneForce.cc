#include "BASim/src/Physics/DeformableObjects/Shells/CSTMembraneForce.hh"
#include "BASim/src/Core/EigenIncludes.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/MatrixBase.hh"

using Eigen::Matrix;

namespace BASim {

//some convenience functions
static Scalar lenSq(const Vec3d& vec) {
  return vec.dot(vec);
}
static Scalar len(const Vec3d& vec) {
  return sqrt(vec.dot(vec));
}
static Scalar dot(const Vec3d& v1, const Vec3d& v2) {
  return v1.dot(v2);
}

CSTMembraneForce::CSTMembraneForce( ElasticShell& shell, const std::string& name, Scalar Youngs, Scalar Poisson, Scalar Youngs_damping, Scalar Poisson_damping, Scalar timestep )
: ElasticShellForce(shell, name), m_Youngs(Youngs), m_Poisson(Poisson), m_Youngs_damp(Youngs_damping), m_Poisson_damp(Poisson_damping), m_timestep(timestep)
{
  
}

bool CSTMembraneForce::gatherDOFs(const FaceHandle& fh, std::vector<Vec3d>& undeformed, std::vector<Vec3d>& damp_undeformed, std::vector<Vec3d>& deformed, std::vector<int>& indices) const {
  assert(undeformed.size() == 3);
  assert(deformed.size() == 3);
  assert(indices.size() == 9);
  
  if(!m_shell.isFaceActive(fh)) return false;

  //extract the relevant data for the local element
  FaceVertexIterator fv_it = m_shell.getDefoObj().fv_iter(fh);
  int i = 0;
  for(;fv_it; ++fv_it) {
    const VertexHandle& vh = *fv_it;
    undeformed[i] = m_shell.getUndeformedPosition(vh);
    deformed[i] = m_shell.getVertexPosition(vh);
    damp_undeformed[i] = m_shell.getDampingUndeformedPosition(vh);
    int dofBase = m_shell.getVertexDofBase(vh);
    indices[i*3] = dofBase;
    indices[i*3+1] = dofBase+1;
    indices[i*3+2] = dofBase+2;
    ++i;
  }
  
  return true;
}

Scalar CSTMembraneForce::globalEnergy() const
{
  Scalar energy = 0;
  std::vector<int> indices(9);
  std::vector<Vec3d> undeformed(3), deformed(3), undeformed_damp(3);

  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;
    
    gatherDOFs(fh, undeformed, undeformed_damp, deformed, indices);
    
    Scalar thickness = m_shell.getThickness(fh);
    
    //determine the energy for this element
    energy += elementEnergy(undeformed, deformed, m_Youngs, m_Poisson, thickness);

  }
  return energy;
}

void CSTMembraneForce::globalForce( VecXd& force )  const
{

  std::vector<int> indices(9);
  std::vector<Vec3d> undeformed(3), undeformed_damp(3), deformed(3);
  Eigen::Matrix<Scalar, 9, 1> localForce;

  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;
   
    if(_debugFlag) {
      //check if the face is above the central axis, and skip it if not
      Vec3d centrePoint(0,0,0);
      for(FaceVertexIterator fvit = m_shell.getDefoObj().fv_iter(fh); fvit; ++fvit) {
        centrePoint += m_shell.getVertexPosition(*fvit);
      }
      centrePoint /= 3.0;
      if(centrePoint[2] > 0) continue;
    }

    bool valid = gatherDOFs(fh, undeformed, undeformed_damp, deformed, indices);
    if(!valid) continue;

    Scalar thickness = m_shell.getThickness(fh);

    //determine the elastic forces for this element
    if(m_Youngs != 0) {
      elementForce(undeformed, deformed, localForce, m_Youngs, m_Poisson, thickness);
      for (unsigned int i = 0; i < indices.size(); ++i)
        force(indices[i]) += localForce(i);
    }

    ////determine the damping forces for this element
    if(m_Youngs_damp != 0) {
      elementForce(undeformed_damp, deformed, localForce, m_Youngs_damp, m_Poisson_damp, thickness);
      for (unsigned int i = 0; i < indices.size(); ++i)
        force(indices[i]) += localForce(i) / m_timestep;  //division by timestep to do Viscous Threads-style viscosity/damping
    }

  }
}

void CSTMembraneForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  std::vector<int> indices(9);
  std::vector<Vec3d> undeformed(3), undeformed_damp(3), deformed(3);
  Eigen::Matrix<Scalar, 9, 9> localMatrix;
  
  FaceIterator fit = m_shell.getDefoObj().faces_begin();
  for (;fit != m_shell.getDefoObj().faces_end(); ++fit) {
    const FaceHandle& fh = *fit;

    bool valid = gatherDOFs(fh, undeformed, undeformed_damp, deformed, indices);
    if(!valid) continue;
    
    Scalar thickness = m_shell.getThickness(fh);

    //determine the elastic forces for this element
    if(m_Youngs != 0) {
      elementJacobian(undeformed, deformed, localMatrix, m_Youngs, m_Poisson, thickness);
      for (unsigned int i = 0; i < indices.size(); ++i)
        for(unsigned int j = 0; j < indices.size(); ++j)
          Jacobian.add(indices[i], indices[j], scale * localMatrix(i,j));
    }

    //determine the viscous forces for this element
    if(m_Youngs_damp != 0) {
      elementJacobian(undeformed_damp, deformed, localMatrix, m_Youngs_damp, m_Poisson_damp, thickness);
      
      //divide by timestep to do Viscous Threads-style symmetric viscous force
      for (unsigned int i = 0; i < indices.size(); ++i)
        for(unsigned int j = 0; j < indices.size(); ++j)
          Jacobian.add(indices[i], indices[j], scale / m_timestep * localMatrix(i,j));

    }

  }
}


Scalar CSTMembraneForce::elementEnergy(const std::vector<Vec3d>& undeformed,
                                       const std::vector<Vec3d>& deformed,
                                      Scalar Young, Scalar Poisson, Scalar thickness) const
{
  
  Eigen::Matrix<Scalar, 3, 3> Tm;
  Scalar lenSqv[3];

  computeHash(undeformed, Tm, lenSqv, Young, Poisson, thickness);

  // difference of squared lengths of edges
  Scalar s[3] = {lenSq(deformed[2] - deformed[1]) - lenSqv[0],
    lenSq(deformed[0] - deformed[2]) - lenSqv[1],
    lenSq(deformed[1] - deformed[0]) - lenSqv[2]};

  // energy = transpose(s) * Tm * s
  Scalar energy = 0;

  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++)
      energy += Tm(i, j) * s[i] * s[j];

  return energy;
}

void CSTMembraneForce::elementForce(const std::vector<Vec3d>& undeformed,
                                  const std::vector<Vec3d>& deformed,
                                  Eigen::Matrix<Scalar, 9, 1>& force,
                                  Scalar Young, Scalar Poisson, Scalar thickness) const
{
  
  force.setZero();

  Eigen::Matrix<Scalar, 3, 3> Tm;
  Scalar lenSqv[3];
  computeHash(undeformed, Tm, lenSqv, Young, Poisson, thickness);

  Vec3d v[3] = {deformed[2] - deformed[1], 
                deformed[0] - deformed[2], 
                deformed[1] - deformed[0]};

  // difference of squared lengths of edges
  Scalar  s[3] = {lenSq(v[0]) - lenSqv[0],
                  lenSq(v[1]) - lenSqv[1],
                  lenSq(v[2]) - lenSqv[2]};

  // taking derivative of the energy with respect to degree of freedom pkl
  // (l-th component of vertex k) gives the formula below for the force on
  // degree of freedom pkl
  for (int k = 0; k < 3; k++) {
    int k2 = (k + 3 - 2) % 3;
    int k1 = (k + 3 - 1) % 3;

    for (int l = 0; l < 3; l++) {
      for (int i = 0; i < 3; i++) {
        force[3 * k + l] -= 4 * (Tm(i, k2) * s[i] * v[k2](l) - Tm(i, k1) * s[i] * v[k1](l));
      }
    }
  }
}

void CSTMembraneForce::elementJacobian(const std::vector<Vec3d>& undeformed,
                                  const std::vector<Vec3d>& deformed,
                                  Eigen::Matrix<Scalar,9,9>&J,
                                  Scalar Young, Scalar Poisson, Scalar thickness) const
{
  
  J.setZero();

  Eigen::Matrix<Scalar, 3, 3> Tm;
  Scalar lenSqv[3];
  computeHash(undeformed, Tm, lenSqv, Young, Poisson, thickness);

  Vec3d v[3] = {deformed[2] - deformed[1], 
                deformed[0] - deformed[2], 
                deformed[1] - deformed[0]};

  // difference of squared lengths of edges
  Scalar s[3] = {lenSq(v[0]) - lenSqv[0],
                  lenSq(v[1]) - lenSqv[1],
                  lenSq(v[2]) - lenSqv[2]};

  // taking derivative of the force on degree of freedom pij (j-th component
  // of vertex i) with respect to degree of freedom pkl (l-th component of
  // vertex k) gives the formula below for the jacobian
  for (int k = 0; k < 3; k++) {
    int k2 = (k + 3 - 2) % 3;
    int k1 = (k + 3 - 1) % 3;

    for (int m = 0; m < 3; m++) {
      int m2 = (m + 3 - 2) % 3;
      int m1 = (m + 3 - 1) % 3;

      for (int l = 0; l < 3; l++) {
        for (int n = 0; n < 3; n++) {
          J(3 * k + l, 3 * m + n) -= 8 * (Tm(m2, k2) * v[k2](l) * v[m2](n) -
            Tm(m1, k2) * v[k2](l) * v[m1](n) -
            Tm(m2, k1) * v[k1](l) * v[m2](n) +
            Tm(m1, k1) * v[k1](l) * v[m1](n));
        }
      }
    }

    for (int l = 0; l < 3; l++) {
      for (int i = 0; i < 3; i++) {
        J(3 * k + l, 3 * k + l) -= 4 * s[i] * (Tm(i, k2) + Tm(i, k1));
        J(3 * k + l, 3 * k1 + l) += 4 * s[i] * Tm(i, k2);
        J(3 * k + l, 3 * k2 + l) += 4 * s[i] * Tm(i, k1);
      }
    }
  }
 
}


void
CSTMembraneForce::computeHash(const std::vector<Vec3d>& undeformed, Eigen::Matrix<Scalar, 3, 3>& Tm, Scalar* lenSqv, 
                                 Scalar Youngs, Scalar Poisson, Scalar thickness) const
{

   Vec3d v[3] = {undeformed[2] - undeformed[1],
    undeformed[0] - undeformed[2],
    undeformed[1] - undeformed[0]};

  for (int i = 0; i < 3; i++)
    lenSqv[i] = lenSq(v[i]);

  Scalar A = 0.5 * len((undeformed[1] - undeformed[0]).cross(undeformed[2] - undeformed[0]));
  
  Eigen::Matrix<Scalar, 3, 3> Tm1, Tm2;
  
  int i, j, k, l, m, n;
  // T_{il1} = (1/8/A^4)*(dot(v_k,v_m)*dot(v_j, v_n)+dot(v_k,v_n)*dot(v_j,v_m))
  for (i = 0; i < 3; i++) {
    j = (i + 1) % 3;
    k = (i + 2) % 3;

    for (l = i; l < 3; l++) {
      m = (l + 1) % 3;
      n = (l + 2) % 3;

      Tm1(i, l) = Tm1(l, i) = dot(v[k], v[m]) * dot(v[j], v[n]) + dot(v[k], v[n]) * dot(v[j], v[m]);
      Tm2(i, l) = Tm2(l, i) = dot(v[k], v[j]) * dot(v[m], v[n]);
    }
  }

  //Dropping the 4th multiplication by A accounts for the fact that we are integrating the energy/force/etc
  //over the area of the triangle. This makes the forces/energies/jacobians consistent as you refine the mesh.
  //The factors of 1/2 from the s terms (strain definition) have been folded in here as well.

  Tm1 = 1.0 / (128.0 * (A * A * A)) * Tm1;
  Tm2 = 1.0 / (64.0 * (A * A * A)) * Tm2;
  Tm = (Youngs*thickness)/2/(1.0-Poisson*Poisson)*((1.0-Poisson)*Tm1 + Poisson*Tm2);

}



} //namespace BASim