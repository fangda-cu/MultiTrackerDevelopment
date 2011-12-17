#include "BASim/src/Physics/DeformableObjects/Shells/DSBendingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/LinearTransfer.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/UpTransfer.hh"
#include "BASim/src/Physics/DeformableObjects/Shells/DownTransfer.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/Math.hh"
#include "BASim/src/Math/MatrixBase.hh"

namespace BASim {

DSBendingForce::DSBendingForce( ElasticShell& shell, const std::string& name, Scalar stiffness, Scalar damping, Scalar timestep) : 
    ElasticShellForce(shell, name), m_stiffness(stiffness), m_damping(damping), m_timestep(timestep), m_func(new LinearTransfer())
{
  
}

bool DSBendingForce::gatherDOFs(const EdgeHandle& eh, 
                                std::vector<Vec3d>& undeformed, 
                                std::vector<Vec3d>& undeformed_damp,
                                std::vector<Vec3d>& deformed, 
                                std::vector<int>& indices ) const {
  assert(undeformed.size() == 4);
  assert(deformed.size() == 4);
  assert(indices.size() == 12);

  if(!m_shell.isEdgeActive(eh)) return false;

  const DeformableObject& defo = m_shell.getDefoObj();

  VertexHandle from_vh, to_vh;
  from_vh = defo.fromVertex(eh);
  to_vh = defo.toVertex(eh);

  //Just use the first two faces we hit. TODO: Non-manifold case...
  EdgeFaceIterator ef_it = defo.ef_iter(eh);
  //if(!ef_it) return false; //if there is no face at all, we can't apply bending.
  const FaceHandle& fh = *ef_it;
  ++ef_it;
  if(!ef_it) return false; //there is no 2nd face, we also can't apply bending.
  const FaceHandle& fh2 = *ef_it;
  ++ef_it;
  assert(fh != fh2); //we should never hit the same face
  assert(!ef_it); //this edge should have no more faces. If it does we need to write a non-manifold extension.

  //find the 3rd vertex in the first face
  FaceVertexIterator fv_it = defo.fv_iter(fh);
  while((*fv_it == from_vh) || (*fv_it == to_vh)) ++fv_it;
  VertexHandle f1_vh = *fv_it;

  //find the 3rd vertex in the second face
  FaceVertexIterator fv_it2 = defo.fv_iter(fh2);
  while((*fv_it2 == from_vh) || (*fv_it2 == to_vh)) ++fv_it2;
  VertexHandle f2_vh = *fv_it2;

  //push back the edge's two vertices, then the opposing vertices of the two face flaps
  undeformed[0] = m_shell.getVertexUndeformed(from_vh);
  undeformed_damp[0] = m_shell.getVertexDampingUndeformed(from_vh);
  deformed[0] = m_shell.getVertexPosition(from_vh);
  int baseIdx = m_shell.getVertexDofBase(from_vh);
  indices[0] = baseIdx;
  indices[1] = baseIdx+1;
  indices[2] = baseIdx+2;

  undeformed[1] = m_shell.getVertexUndeformed(to_vh);
  undeformed_damp[1] = m_shell.getVertexDampingUndeformed(to_vh);
  deformed[1] = m_shell.getVertexPosition(to_vh);
  baseIdx = m_shell.getVertexDofBase(to_vh);
  indices[3] = baseIdx;
  indices[4] = baseIdx+1;
  indices[5] = baseIdx+2;

  undeformed[2] = m_shell.getVertexUndeformed(f1_vh);
  undeformed_damp[2] = m_shell.getVertexDampingUndeformed(f1_vh);
  deformed[2] = m_shell.getVertexPosition(f1_vh);
  baseIdx = m_shell.getVertexDofBase(f1_vh);
  indices[6] = baseIdx;
  indices[7] = baseIdx+1;
  indices[8] = baseIdx+2;

  undeformed[3] = m_shell.getVertexUndeformed(f2_vh);
  undeformed_damp[3] = m_shell.getVertexDampingUndeformed(f2_vh);
  deformed[3] = m_shell.getVertexPosition(f2_vh);
  baseIdx = m_shell.getVertexDofBase(f2_vh);
  indices[9] = baseIdx;
  indices[10] = baseIdx+1;
  indices[11] = baseIdx+2;
  
  return true;
}

void DSBendingForce::gatherDOFs(const EdgeHandle& eh, 
                                const FaceHandle& fh,
                                const FaceHandle& fh2,
                                std::vector<Vec3d>& undeformed, 
                                std::vector<Vec3d>& undeformed_damp,
                                std::vector<Vec3d>& deformed, 
                                std::vector<int>& indices ) const 
{
  assert(undeformed.size() == 4);
  assert(deformed.size() == 4);
  assert(indices.size() == 12);

  const DeformableObject& defo = m_shell.getDefoObj();

  VertexHandle from_vh, to_vh;
  from_vh = defo.fromVertex(eh);
  to_vh = defo.toVertex(eh);

  //find the 3rd vertex in the first face
  FaceVertexIterator fv_it = defo.fv_iter(fh);
  while((*fv_it == from_vh) || (*fv_it == to_vh)) ++fv_it;
  VertexHandle f1_vh = *fv_it;

  //find the 3rd vertex in the second face
  FaceVertexIterator fv_it2 = defo.fv_iter(fh2);
  while((*fv_it2 == from_vh) || (*fv_it2 == to_vh)) ++fv_it2;
  VertexHandle f2_vh = *fv_it2;

  //push back the edge's two vertices, then the opposing vertices of the two face flaps
  undeformed[0] = m_shell.getVertexUndeformed(from_vh);
  undeformed_damp[0] = m_shell.getVertexDampingUndeformed(from_vh);
  deformed[0] = m_shell.getVertexPosition(from_vh);
  int baseIdx = m_shell.getVertexDofBase(from_vh);
  indices[0] = baseIdx;
  indices[1] = baseIdx+1;
  indices[2] = baseIdx+2;

  undeformed[1] = m_shell.getVertexUndeformed(to_vh);
  undeformed_damp[1] = m_shell.getVertexDampingUndeformed(to_vh);
  deformed[1] = m_shell.getVertexPosition(to_vh);
  baseIdx = m_shell.getVertexDofBase(to_vh);
  indices[3] = baseIdx;
  indices[4] = baseIdx+1;
  indices[5] = baseIdx+2;

  undeformed[2] = m_shell.getVertexUndeformed(f1_vh);
  undeformed_damp[2] = m_shell.getVertexDampingUndeformed(f1_vh);
  deformed[2] = m_shell.getVertexPosition(f1_vh);
  baseIdx = m_shell.getVertexDofBase(f1_vh);
  indices[6] = baseIdx;
  indices[7] = baseIdx+1;
  indices[8] = baseIdx+2;

  undeformed[3] = m_shell.getVertexUndeformed(f2_vh);
  undeformed_damp[3] = m_shell.getVertexDampingUndeformed(f2_vh);
  deformed[3] = m_shell.getVertexPosition(f2_vh);
  baseIdx = m_shell.getVertexDofBase(f2_vh);
  indices[9] = baseIdx;
  indices[10] = baseIdx+1;
  indices[11] = baseIdx+2;

}


void DSBendingForce::getEdgeFacePairs(EdgeHandle eh, std::vector< std::pair<FaceHandle,FaceHandle> >& facePairs) const {
  DeformableObject& defo = m_shell.getDefoObj();

  std::vector<FaceHandle> faces;;
  for(EdgeFaceIterator efit = defo.ef_iter(eh); efit; ++efit) {
    faces.push_back(*efit);
  }
  
  if(faces.size() == 2)
    facePairs.push_back(std::make_pair(faces[0], faces[1])); //just one pair
  else if(faces.size() == 3) {
    facePairs.push_back(std::make_pair(faces[0], faces[1])); //all pairs
    facePairs.push_back(std::make_pair(faces[0], faces[2]));
    facePairs.push_back(std::make_pair(faces[1], faces[2]));
  }
  else if(faces.size() >= 4) {
    std::cerr << "Non-manifold edges with > 3 faces not yet supported.\n";
    //TODO Need to work out which faces are closest to each other.
  }

}

Scalar DSBendingForce::globalEnergy() const
{
  Scalar energy = 0;
  DeformableObject& defo = m_shell.getDefoObj();
  std::vector<Vec3d> undeformed(4), undeformed_damp(4), deformed(4);
  std::vector<int> indices(12);

  EdgeIterator eit = defo.edges_begin();
  for (;eit != defo.edges_end(); ++eit) {
    const EdgeHandle& eh = *eit;

    std::vector< std::pair<FaceHandle, FaceHandle> > pairs;
    getEdgeFacePairs(eh, pairs);

    for(unsigned int i = 0; i < pairs.size(); ++i) {
      gatherDOFs(eh, pairs[i].first, pairs[i].second, undeformed, undeformed_damp, deformed, indices);

      //determine the energy for this element
      energy += elementEnergy(undeformed, deformed);
    }
  }

  return energy;
}

void DSBendingForce::globalForce( VecXd& force ) const
{

  Eigen::Matrix<Scalar, 12, 1> localForce;
  std::vector<Vec3d> undeformed(4), undeformed_damp(4), deformed(4);
  std::vector<int> indices(12);

  DeformableObject& defo = m_shell.getDefoObj();

  EdgeIterator eit = defo.edges_begin();
  for (;eit != defo.edges_end(); ++eit) {
    const EdgeHandle& eh = *eit;
    
    //gather the data
    std::vector< std::pair<FaceHandle, FaceHandle> > pairs;
    getEdgeFacePairs(eh, pairs);
    
    for(unsigned int p = 0; p < pairs.size(); ++p) {
      
      gatherDOFs(eh, pairs[p].first, pairs[p].second, undeformed, undeformed_damp, deformed, indices);

      //account for double-counted faces in the non-manifold case
      Scalar local_stiffness = (pairs.size() == 2? m_stiffness : 0.5*m_stiffness);
      Scalar local_damping = (pairs.size() == 2? m_damping : 0.5*m_damping);
      
      //determine the elastic forces for this element
      if(m_stiffness != 0) {
        elementForce(undeformed, deformed, localForce);
        for (unsigned int i = 0; i < indices.size(); ++i)
          force(indices[i]) += local_stiffness * localForce(i);
      }

      //////determine the (Rayleigh) damping / viscous forces for this element
      if(m_damping != 0) {
        elementForce(undeformed_damp, deformed, localForce);
        
        //Thickness dependent damping, coefficient mu*h^3/3, div by timestep for symmetric dissipative potential viscosity
        Scalar thickness = getEdgeThickness(eh);
        Scalar h3 = thickness*thickness*thickness;
        Scalar scale_factor = local_damping * h3 / 3.0 / m_timestep;
        for (unsigned int i = 0; i < indices.size(); ++i) {
           force(indices[i]) += scale_factor * localForce(i);
        }
      }
    }
  

  }
}

void DSBendingForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{
  std::vector<Vec3d> undeformed(4), undeformed_damp(4), deformed(4);
  std::vector<int> indices(12);
  Eigen::Matrix<Scalar, 12, 12> localJacobian;
  DeformableObject& defo = m_shell.getDefoObj();

  EdgeIterator eit = defo.edges_begin();
  for (;eit != defo.edges_end(); ++eit) {
    const EdgeHandle& eh = *eit;
    
    
    //gather the data
    std::vector< std::pair<FaceHandle, FaceHandle> > pairs;
    getEdgeFacePairs(eh, pairs);

    for(unsigned int p = 0; p < pairs.size(); ++p) {

      //gather the relevant data for the local element
      //and map the DOF's.
      gatherDOFs(eh, pairs[p].first, pairs[p].second, undeformed, undeformed_damp, deformed, indices);

      //account for double-counted faces in the non-manifold case
      Scalar local_stiffness = (pairs.size() == 2? m_stiffness : 0.5*m_stiffness);
      Scalar local_damping = (pairs.size() == 2? m_damping : 0.5*m_damping);

      //determine the forces for this element
      if(m_stiffness != 0) {
        elementJacobian(undeformed, deformed, localJacobian);
        for (unsigned int i = 0; i < indices.size(); ++i)
          for(unsigned int j = 0; j < indices.size(); ++j) 
            Jacobian.add(indices[i],indices[j], scale * local_stiffness * localJacobian(i,j));
      }
      
      if(m_damping != 0) {
        elementJacobian(undeformed_damp, deformed, localJacobian);
        Scalar thickness = getEdgeThickness(eh);
        Scalar h3 = thickness*thickness*thickness;
        Scalar scale_factor = scale * local_damping * h3 / 3.0 / m_timestep;

        for (unsigned int i = 0; i < indices.size(); ++i)
          for(unsigned int j = 0; j < indices.size(); ++j) 
            Jacobian.add(indices[i],indices[j], scale_factor * localJacobian(i,j));
      }
    }

  }
}

static Scalar len(const Vec3d& vec) { return sqrt(vec.dot(vec)); }
static Scalar lenSq(const Vec3d& vec) { return vec.dot(vec); }
static Vec3d dir(const Vec3d& c) { Scalar l = len(c); assert(l != 0); return c/l; }

//does Eigen not have a 2-parameter version of these? it's way more readable.
static Vec3d cross(const Vec3d& a, const Vec3d& b) { return a.cross(b); }
static Scalar dot(const Vec3d& v1, const Vec3d& v2) { return v1.dot(v2); }



static void computeAngle(Scalar& theta, const std::vector<Vec3d>& positions)
{
  Vec3d p1 = positions[0], 
    p2 = positions[1], 
    q1 = positions[2], 
    q2 = positions[3];

  Vec3d ev = p2 - p1;

  // compute normals of the two deformed triangles and the angle between them
  //
  Vec3d n1 = cross(ev, q1 - p1), 
        n2 = cross(q2 - p1, ev);
  theta = atan2(dot(cross(n1, n2), dir(ev)), dot(n1, n2));

}


static void computeEdgeLength(Scalar& e, const std::vector<Vec3d>& positions)
{
  Vec3d p1 = positions[0], 
    p2 = positions[1];
    //q1 = positions[2], 
    //q2 = positions[3];

  e = len(p2 - p1);

}


static void computeH(Scalar& h, const std::vector<Vec3d>& positions)
{
  Vec3d p1 = positions[0], 
    p2 = positions[1], 
    q1 = positions[2], 
    q2 = positions[3];

  Scalar A1 = 0.5 * len(cross(p2 - p1, q1 - p1));
  Scalar A2 = 0.5 * len(cross(p2 - p1, q2 - p1));
  Scalar e;
  computeEdgeLength(e, positions);
  h = 2.0 * (A1 + A2) / (3.0 * e);

}

static void computeConstant(Scalar& c, const std::vector<Vec3d>& positions)
{
  Vec3d p1 = positions[0], 
    p2 = positions[1], 
    q1 = positions[2], 
    q2 = positions[3];

  // compute normals of the two triangles and their areas
  //
  Vec3d ev = p2 - p1;
  Vec3d n1 = cross(ev, q1 - p1), n2 = cross(q2 - p1, ev);
  Scalar A1doubled = len(n1);
  Scalar A2doubled = len(n2);

  c = 3.0 * lenSq(ev) / (A1doubled + A2doubled);

}

Scalar DSBendingForce::getEdgeThickness(const EdgeHandle& eh) const {
  EdgeFaceIterator efit = m_shell.getDefoObj().ef_iter(eh);
  Scalar thickness = 0;
  int count = 0;
  for(;efit; ++efit) {
    FaceHandle fh = *efit;
    thickness += m_shell.getThickness(fh);
    ++count;
  }
  assert(count == 2);
  thickness /= 2.0; //Thickness is just the average of incident face thicknesses.
  
  return thickness;
}

Scalar DSBendingForce::elementEnergy(const std::vector<Vec3d>& undeformed,
                                   const std::vector<Vec3d>& deformed) const 
{
  Scalar f, theta;
  Scalar energy;
  computeAngle(theta, deformed);
  f = m_func->computeFunction(theta);

  Scalar f0, e0, h0, theta0;
  
  //TODO Could be precomputed
  computeEdgeLength(e0, undeformed);
  computeH(h0, undeformed);
  computeAngle(theta0, undeformed);
  f0 = m_func->computeFunction(theta0);
  Scalar strength = 1;

  //hard edge
  if (fabs(theta0) > 2*M_PI) strength = 1.0e-03;

  Scalar x = (f - f0);
  energy = e0 / h0 * x * x * strength;
  return energy;
}

void DSBendingForce::elementForce(const std::vector<Vec3d>& undeformed,
                             const std::vector<Vec3d>& deformed,
                             Eigen::Matrix<Scalar,12,1>& force) const
{
  assert(undeformed.size() == deformed.size());

  force.setZero();

  Scalar f, theta;

  computeAngle(theta, deformed);
  f = m_func->computeFunction(theta);

  //TODO These can be precomputed and stored in the elastic case.
  Scalar f0, e0, h0, theta0;
  computeAngle(theta0, undeformed);
  f0 = m_func->computeFunction(theta0);
  computeEdgeLength(e0, undeformed);
  computeH(h0, undeformed);
  Scalar strength = 1;

  //hard edge
  if (fabs(theta0) > 2*M_PI) strength = 1.0e-03;


  Scalar x = (f - f0);

  Vec3d p1 = deformed[0], 
    p2 = deformed[1], 
    q1 = deformed[2], 
    q2 = deformed[3];

  Vec3d dp1, dp2, dq1, dq2;
  ComputeDihedralAngleDerivatives(dp1, dp2, dq1, dq2, p1, p2, q1, q2);

  Scalar fp = m_func->computeFirstDerivative(theta);
  
  for (int i = 0; i < 3; i++)
  {
    force[i] = -2 * e0 / h0 * x * fp * dp1[i] * strength;
    force[3+i] = -2 * e0 / h0 * x * fp * dp2[i] * strength;
    force[6+i] = -2 * e0 / h0 * x * fp * dq1[i] * strength;
    force[9+i] = -2 * e0 / h0 * x * fp * dq2[i] * strength;
  }

}

void DSBendingForce::elementJacobian(const std::vector<Vec3d>& undeformed,
                                     const std::vector<Vec3d>& deformed,
                                     Eigen::Matrix<Scalar,12,12>& jac) const
{

  jac.setZero();

  Scalar f, theta;

  computeAngle(theta, deformed);
  f = m_func->computeFunction(theta);

  Scalar f0, e0, h0, theta0=0;
  computeAngle(theta0, undeformed);
  f0 = m_func->computeFunction(theta0);
  computeEdgeLength(e0, undeformed);
  computeH(h0, undeformed);
  Scalar strength = 1;

  //hard edge
  if (fabs(theta0) > 2*M_PI) strength = 1.0e-03;


  Vec3d p1 = deformed[0], 
    p2 = deformed[1], 
    q1 = deformed[2], 
    q2 = deformed[3];

  static Vec3d dpq[4];
  ComputeDihedralAngleDerivatives(dpq[0], dpq[1], dpq[2], dpq[3], p1, p2, q1, q2);

  EnergyHessian J(4);
  ComputeDihedralAngleSecondDerivatives(J, 1, p1, p2, q2, q1, 0, 1, 3, 2);

  Scalar x = (f - f0);
  Scalar fp = m_func->computeFirstDerivative(theta);
  Scalar fpp = m_func->computeSecondDerivative(theta);

  for (int i = 0; i < 12; i++) {
    for (int j = 0; j < 12; j++) {
      jac(i, j) =
        -2.0 * e0 / h0 * x * fp * strength *
        J.m_Hessian[i/3][j/3](i % 3, j % 3);
    }
  }

  for (int i = 0; i < 4; i++) {
    for (int j = 0; j < 4; j++) {
      Mat3d m = outerProd(dpq[i], dpq[j]);

      for (int k = 0; k < 3; k++) {
        for (int l = 0; l < 3; l++) {
          jac(3*i + k, 3*j + l) +=
            -2.0 * e0 / h0 * (fp * fp + x * fpp) * strength *
            m(k, l);
        }
      }
    }
  }


}

} //namespace BASim