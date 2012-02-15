#include "BASim/src/Physics/DeformableObjects/Shells/MNBendingForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/Math.hh"
#include "BASim/src/Math/MatrixBase.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"



namespace BASim {


struct MNBendDofStruct
{
  Vector3d* p;
  Vector3d* q;
  Real* xi;
};


MNBendingForce::MNBendingForce(ElasticShell& shell, const std::string& name, Scalar Youngs, Scalar Poisson, Scalar Youngs_damping, Scalar Poisson_damping, Scalar timestep) : 
ElasticShellForce(shell, name), m_Youngs(Youngs), m_Poisson(Poisson), m_Youngs_damp(Youngs_damping), m_Poisson_damp(Poisson_damping), m_timestep(timestep), m_precomputed(&(shell.getDefoObj())), m_initialized(false)
{
  
}

void fillVert(std::vector<Scalar>& data, int startInd, const Vec3d& vert) {
  data[startInd] = vert[0];
  data[startInd+1] = vert[1];
  data[startInd+2] = vert[2];
}

void fillVertNaN(std::vector<Scalar>& data, int startInd) {
  Scalar zero = 0.0;
  Scalar NaN = 0.0 / zero;
  data[startInd] = NaN;
  data[startInd+1] = NaN;
  data[startInd+2] = NaN;
}

void MNBendingForce::update() {
  
  if(m_initialized)
    return;
  
  std::cout << "Initial update() called.\n";

  DeformableObject& defo = m_shell.getDefoObj();
  std::vector<Scalar> undeformed(MNBendStencilSize), undeformed_damp(MNBendStencilSize), deformed(MNBendStencilSize);
  std::vector<int> indices(MNBendStencilSize);

  FaceIterator fit = defo.faces_begin();
  for (;fit != defo.faces_end(); ++fit) {
    const FaceHandle& fh = *fit;

    bool valid = gatherDOFs(fh, undeformed, undeformed_damp, deformed, indices);
    if(!valid) continue;

    updatePrecomp(fh, undeformed, deformed, &m_precomputed[fh]);
  }
}

bool MNBendingForce::gatherDOFs(const FaceHandle& fh, 
                                std::vector<Scalar>& undeformed, 
                                std::vector<Scalar>& undeformed_damp,
                                std::vector<Scalar>& deformed,
                                std::vector<int>& indices ) const 
{
  assert(undeformed.size() == MNBendStencilSize);
  assert(undeformed_damp.size() == MNBendStencilSize);
  assert(deformed.size() == MNBendStencilSize);
  assert(indices.size() == MNBendStencilSize);
  //data consists of 3 face vertices, then 3 (opposing, ordered) flap vertices,
  //then the xi values on the 3 edges.

  if(!m_shell.isFaceActive(fh)) return false;

  const DeformableObject& defo = m_shell.getDefoObj();

  FaceEdgeIterator feit = defo.fe_iter(fh);
  int v = 0;
  int e = 9;
  for(; feit; ++feit) {
    EdgeHandle eh = *feit;
    
    //get the edge DOF
    indices[e] = m_shell.getEdgeDofBase(eh);
    deformed[e] = m_shell.getEdgeXi(eh);
    undeformed[e] = m_shell.getEdgeUndeformedXi(eh);
    undeformed_damp[e] = m_shell.getDampingUndeformedXi(eh);

    //vertex opposite the edge in the central tri
    VertexHandle ph; 
    getFaceThirdVertex(defo, fh, eh, ph);
    int baseID = m_shell.getVertexDofBase(ph);
    indices[v] = baseID;
    indices[v+1] = baseID+1;
    indices[v+2] = baseID+2;

    fillVert(deformed, v, m_shell.getVertexPosition(ph));
    fillVert(undeformed, v, m_shell.getVertexUndeformed(ph));
    fillVert(undeformed_damp, v, m_shell.getVertexDampingUndeformed(ph));

    //vertex opposite the edge on the flap triangle for this edge
    FaceHandle flapFace = getEdgeOtherFace(defo, eh, fh);
    if(flapFace.isValid()) {
      VertexHandle qh; 
      getFaceThirdVertex(defo, flapFace, eh, qh);
      int baseID = m_shell.getVertexDofBase(qh);
      fillVert(deformed, v+9, m_shell.getVertexPosition(qh));
      fillVert(undeformed, v+9, m_shell.getVertexUndeformed(qh));
      fillVert(undeformed_damp, v+9, m_shell.getVertexDampingUndeformed(qh));
    }
    else {
      fillVertNaN(deformed, v+9);
      fillVertNaN(undeformed, v+9);
      fillVertNaN(undeformed_damp, v+9);
    }

    e++;
    v+=3;
  }
  return true;
}

void mapDof(const std::vector<Scalar>& src, MNBendDofStruct& dest) {
  dest.p = (Vector3d*)&src[0];
  dest.q = (Vector3d*)&src[9];
  dest.xi = (Real*)&src[18];
}

void MNBendingForce::doPrecomputation() const {
  
  if(m_initialized)
    return;

  DeformableObject& defo = m_shell.getDefoObj();
  std::vector<Scalar> undeformed(MNBendStencilSize), undeformed_damp(MNBendStencilSize), deformed(MNBendStencilSize);
  std::vector<int> indices(MNBendStencilSize);

  FaceIterator fit = defo.faces_begin();
  for (;fit != defo.faces_end(); ++fit) {
    const FaceHandle& fh = *fit;
    
    bool valid = gatherDOFs(fh, undeformed, undeformed_damp, deformed, indices);
    if(!valid) continue;

    //construct the "precomputed" data needed for this element
    
    initializePrecomp(fh, undeformed, deformed, &m_precomputed[fh]);
  }

  m_initialized = true;
}

Scalar MNBendingForce::globalEnergy() const
{
  doPrecomputation();

  Scalar energy = 0;
  DeformableObject& defo = m_shell.getDefoObj();
  std::vector<Scalar> undeformed(MNBendStencilSize), undeformed_damp(MNBendStencilSize), deformed(MNBendStencilSize);
  std::vector<int> indices(NumMNBendDof);

  FaceIterator fit = defo.faces_begin();
  for (;fit != defo.faces_end(); ++fit) {
    const FaceHandle& fh = *fit;

    bool valid = gatherDOFs(fh, undeformed, undeformed_damp, deformed, indices);
    if(!valid) continue;

    //get the precomputed data for this element
    MNPrecomputed* pre = &m_precomputed[fh];

    //determine the energy for this element
    energy += elementEnergy(undeformed, deformed, pre);

  }

  return energy;
}


void MNBendingForce::globalForce( VecXd& force ) const
{
  
  doPrecomputation();

  Eigen::Matrix<Scalar, NumMNBendDof, 1> localForce;
  std::vector<Scalar> undeformed(MNBendStencilSize), undeformed_damp(MNBendStencilSize), deformed(MNBendStencilSize);
  std::vector<int> indices(NumMNBendDof);

  DeformableObject& defo = m_shell.getDefoObj();

  FaceIterator fit = defo.faces_begin();
  for (;fit != defo.faces_end(); ++fit) {
    const FaceHandle& fh = *fit;
    //gather the data

    bool valid = gatherDOFs(fh, undeformed, undeformed_damp, deformed, indices);
    if(!valid) continue;

    //determine the elastic forces for this element
    if(m_Youngs != 0) {
      //get the precomputed data for this element
      MNPrecomputed* pre = &m_precomputed[fh];

      elementForce(undeformed, deformed, localForce, pre);
      for (unsigned int i = 0; i < indices.size(); ++i) {
        force(indices[i]) += localForce(i);
      }

    }

    //determine the (Rayleigh) damping / viscous forces for this element
    /*
    if(m_Youngs_damp != 0) {
      std::cout << "Warning, doing damping code\n";
      //TODO Precomp
      MNPrecomputed* pre = &m_precomputed[fh];

      elementForce(undeformed_damp, deformed, localForce, pre);
      for (unsigned int i = 0; i < indices.size(); ++i) {
        //divide by timestep for symmetric dissipative potential viscosity
        force(indices[i]) += localForce(i) / m_timestep;
      }
    }
    */

  }
}


void MNBendingForce::globalJacobian( Scalar scale, MatrixBase& Jacobian ) const
{

  doPrecomputation();

  std::vector<Scalar> undeformed(MNBendStencilSize), undeformed_damp(MNBendStencilSize), deformed(MNBendStencilSize);
  std::vector<int> indices(NumMNBendDof);
  Eigen::Matrix<Scalar, NumMNBendDof, NumMNBendDof> localJacobian;
  DeformableObject& defo = m_shell.getDefoObj();

  FaceIterator fit = defo.faces_begin();
  for (;fit != defo.faces_end(); ++fit) {
    const FaceHandle& fh = *fit;

    //gather the relevant data for the local element
    //and map the DOF's.
    bool valid = gatherDOFs(fh, undeformed, undeformed_damp, deformed, indices);
    if(!valid) continue;

    //determine the forces for this element
    if(m_Youngs != 0) {
      //get the precomputed data for this element
      MNPrecomputed* pre = &m_precomputed[fh];

      elementJacobian(undeformed, deformed, localJacobian, pre);
      for (unsigned int i = 0; i < indices.size(); ++i)
        for(unsigned int j = 0; j < indices.size(); ++j) 
          Jacobian.add(indices[i],indices[j], scale * localJacobian(i,j));
    }

    /*
    if(m_Youngs_damp != 0) {
      //TODO Precomp
      MNPrecomputed* pre = &m_precomputed[fh];
      elementJacobian(undeformed_damp, deformed, localJacobian, pre);

      for(unsigned int i = 0; i < indices.size(); ++i)
        for(unsigned int j = 0; j < indices.size(); ++j) 
          Jacobian.add(indices[i],indices[j], scale / m_timestep * localJacobian(i,j));
    }
    */

  }
}


//  compute a set of quantities associated with a triangle 
//  used in several functions
template <class RealType>
void ComputeTriangleAttrib(
  const CVec3T<RealType> p[NumTriPoints],
  CVec3T<RealType> v[NumTriPoints], // vectors along edges v[i] = p[k]-p[j]
  CVec3T<RealType> t[NumTriPoints], // side perpendiculars, pointing outwards, of length equal to v[i]
  CVec3T<RealType>& n,               // unit normal 
  RealType& A                        // area
  ) 
{ 

  v[0] = p[2]-p[1]; v[1] = p[0]-p[2]; v[2] = p[1]-p[0];
  n  = cross(v[1],v[2]);
  A  = 0.5*len(n); assert(A > 0.0 );
  n = n/(2.0*A);
  t[0]  = cross(v[0],n);  t[1] = cross(v[1],n);  t[2] = cross(v[2],n); 
}

// compute normals of the 3 triangles adjacent to the element triangle
template <class RealType>
void ComputeFlapNormals(const CVec3T<RealType> p[NumTriPoints], const CVec3T<RealType> q[NumTriPoints], const CVec3T<RealType> v[NumTriPoints],
  CVec3T<RealType> nopp[NumTriPoints] // normals to the adjacent triangles
) {
  CVec3T<RealType> n =    dir(cross(v[1],v[2]));
  if(!isnan(q[0](0))) {
    nopp[0] = cross(q[0]-p[1],v[0]); normalize(nopp[0]);
  }  else nopp[0] = n;

  if(!isnan(q[1](0))) {
    nopp[1] = cross(q[1]-p[2],v[1]); normalize(nopp[1]);
  }  else nopp[1] = n;

  if(!isnan(q[2](0))) {
    nopp[2] = cross(q[2]-p[0],v[2]);  normalize(nopp[2]);
  }  else nopp[2] = n;
}


// computes perpendiculars to average edge normals (tau) and c which are
// equivalent 1/projections of triangle normals to average normals

void ComputeEdgeFrameParams(
  const CVec3T<Real> v[NumTriPoints], const CVec3T<Real> t[NumTriPoints],  const CVec3T<Real> nopp[NumTriPoints],//in
  CVec3T<Real> tau[NumTriPoints], Real c[NumTriPoints])// out
{ 
  CVec3T<Real> topp[NumTriPoints], tunit[NumTriPoints];
  topp[0] = cross(nopp[0],v[0]);    topp[1] = cross(nopp[1],v[1]);    topp[2] = cross(nopp[2],v[2]); 
  tunit[0] = dir(t[0]);             tunit[1] = dir(t[1]);              tunit[2] = dir(t[2]); 

  tau[0] = dir(t[0]-topp[0]); tau[1] = dir(t[1]-topp[1]); tau[2] = dir(t[2]-topp[2]);

  c[0] = dot(tunit[0],tau[0]); assert(fabs(c[0]) > FLT_EPSILON); c[0] =   1.0/c[0];
  c[1] = dot(tunit[1],tau[1]); assert(fabs(c[1]) > FLT_EPSILON); c[1] =   1.0/c[1];
  c[2] = dot(tunit[2],tau[2]); assert(fabs(c[2]) > FLT_EPSILON); c[2] =   1.0/c[2];  
}



template <int DO_HESS>
adreal<NumMNBendDof,DO_HESS,Real> MNEnergy(const MNBendingForce& mn, const std::vector<Scalar>& undeformed, const std::vector<Scalar>& deformed, MNPrecomputed* pre) {  

  // typedefs to simplify code below
  typedef adreal<NumMNBendDof,DO_HESS,Real> adrealMN;
  typedef CVec3T<adrealMN> advecMN;

  MNPrecomputed* pc = pre;

  MNBendDofStruct  s_undeformed, s_deformed;
  mapDof( undeformed, s_undeformed );    
  mapDof( deformed, s_deformed );

  const Vector3d* qu  = s_undeformed.q;

  // independent variables
  advecMN   p[NumTriPoints]; // vertex positions
  set_independent( p[0], s_deformed.p[0], 0 );
  set_independent( p[1], s_deformed.p[1], 3 );
  set_independent( p[2], s_deformed.p[2], 6 );    

  adrealMN  xi[NumTriPoints]; // mid-edge normal variables
  xi[0].set_independent( s_deformed.xi[0], 9 );
  xi[1].set_independent( s_deformed.xi[1], 10 );
  xi[2].set_independent( s_deformed.xi[2], 11 );        

  // dependent variables
  advecMN v[NumTriPoints], t[NumTriPoints], n;   
  adrealMN   A;  
  adrealMN  w[NumTriPoints];

  ComputeTriangleAttrib(p,v,t,n,A); 

  bool nbrValid0 = !isnan(qu[0](0));
  bool nbrValid1 = !isnan(qu[1](0));
  bool nbrValid2 = !isnan(qu[2](0));

  w[0] = (-dot(n,pc->tau0[0])  + Real(pc->s[0])*xi[0])*pc->c0[0];
  w[1] = (-dot(n,pc->tau0[1])  + Real(pc->s[1])*xi[1])*pc->c0[1];
  w[2] = (-dot(n,pc->tau0[2])  + Real(pc->s[2])*xi[2])*pc->c0[2];

  if(!nbrValid0) w[0] = 0;
  if(!nbrValid1) w[1] = 0;
  if(!nbrValid2) w[2] = 0;

  adrealMN e(0);
  for(int i= 0; i < NumTriPoints; i++)
    for(int j= 0; j < NumTriPoints; j++)
      e += (w[i]-pc->w_undef[i])*(w[j]-pc->w_undef[j])*pc->T[NumTriPoints*i+j];  
  e *= 0.5;    

  return e;
}


// compute T and w_undef from undeformed configuration
void MNBendingForce::initializePrecomp( const FaceHandle& face, const std::vector<Scalar>& undeformed, const std::vector<Scalar>& deformed, MNPrecomputed* pc) const 
{ 
  // map generic variables to energy-specific structures
  std::cout << "init";
  MNBendDofStruct s_undeformed;
  mapDof( undeformed, s_undeformed );

  // aliases for brevity
  const Vector3d* pu  = s_undeformed.p;
  const Vector3d* qu  = s_undeformed.q;
  const Real* xi_undef = s_undeformed.xi;

  Vector3d vu[NumTriPoints], tu[NumTriPoints], noppu[NumTriPoints];
  Vector3d nu;
  Real Au;
  Vector3d tauu[NumTriPoints];
  Real cu[NumTriPoints];

  ComputeTriangleAttrib(pu,vu,tu,nu,Au);
  ComputeFlapNormals(pu,qu,vu,noppu);
  ComputeEdgeFrameParams(vu,tu,noppu,tauu,cu);  

  //determine ownership according to whichever face comes first in the relevant edge's iterator.
  FaceEdgeIterator feit = m_shell.getDefoObj().fe_iter(face);
  int i = 0;
  for(; feit; ++feit ) {
    EdgeHandle eh = *feit;
    EdgeFaceIterator efit = m_shell.getDefoObj().ef_iter(eh);
    pc->s[i] = (*efit == face) ? 1 : -1;
    ++i;
  }
  bool nbrValid0 = !isnan(qu[0](0));
  bool nbrValid1 = !isnan(qu[1](0));
  bool nbrValid2 = !isnan(qu[2](0));

  pc->w_undef[0] = (-dot(nu,tauu[0])  + Real(pc->s[0])*xi_undef[0])*cu[0];
  pc->w_undef[1] = (-dot(nu,tauu[1])  + Real(pc->s[1])*xi_undef[1])*cu[1];
  pc->w_undef[2] = (-dot(nu,tauu[2])  + Real(pc->s[2])*xi_undef[2])*cu[2];

  if(!nbrValid0) pc->w_undef[0] = 0;
  if(!nbrValid1) pc->w_undef[1] = 0;
  if(!nbrValid2) pc->w_undef[2] = 0;

  //Note that it would ordinarily be divided by area^2, but since we integrate over area to compute energy,
  //one of them goes away.
  Real Tc1[9], Tc2[9];
  for(int i=0; i <= 2; i++)
    for(int j=0; j <= 2; j++) {
      Tc1[NumTriPoints*j+i] = dot(vu[i],vu[j])*dot(vu[i],vu[j])/Au/len(vu[i])/len(vu[j]); //sqr(dot(vu[i],vu[j]))/Au/len(vu[i])/len(vu[j]);                 
      Tc2[NumTriPoints*j+i] = len(vu[i])*len(vu[j])/Au;
    }

    Real Y = m_Youngs;
    Real h = m_shell.getThickness(face);
    Real poisson = m_Poisson;

    for(int i=0; i<9; i++ )
      pc->T[i] = (Y*h*h*h/12.0)/(1.0-poisson*poisson)*((1.0-poisson)*Tc1[i] + poisson*Tc2[i]);
   
}

// use current deformed positions to compute reference coordinate system quantities 
void MNBendingForce::updatePrecomp(const FaceHandle& face, const std::vector<Scalar>& undeformed, const std::vector<Scalar>& deformed, MNPrecomputed* pre) const
{    
  std::cout << "upd";

  // map generic variables to energy-specific structures
  MNPrecomputed* pc = pre;
  MNBendDofStruct s_deformed;
  mapDof( deformed, s_deformed );

  const Vector3d* p0  = s_deformed.p;
  const Vector3d* q0  = s_deformed.q;

  Vector3d* tau0 = pc->tau0;
  Real* c0 = pc->c0;

  Vector3d v0[NumTriPoints], t0[NumTriPoints],  n0; 
  Real   A0;  
  Vector3d   nopp0[NumTriPoints]; 

  ComputeTriangleAttrib(p0,v0,t0,n0,A0);  
  ComputeFlapNormals(p0,q0,v0,nopp0);
  ComputeEdgeFrameParams(v0,t0,nopp0,tau0,c0);  

}


Scalar MNBendingForce::elementEnergy(const std::vector<Scalar>& undeformed,
  const std::vector<Scalar>& deformed, MNPrecomputed* pre) const 
{
  assert(pre);

  adreal<NumMNBendDof,0,Real> e = MNEnergy<0>( *this, undeformed, deformed, const_cast<MNPrecomputed*>(pre) );
  Scalar energy = e.value();

  return energy;
}

void MNBendingForce::elementForce(const std::vector<Scalar>& undeformed,
  const std::vector<Scalar>& deformed,
  Eigen::Matrix<Scalar,NumMNBendDof,1>& force, 
  MNPrecomputed* pre) const
{
  assert( pre );
  assert( undeformed.size() == deformed.size() );
//  assert( force.numDof() == undeformed.size() );

  adreal<NumMNBendDof,0,Real> e = MNEnergy<0>(*this, undeformed, deformed, const_cast<MNPrecomputed*>(pre) );     
  for( uint i = 0; i < NumMNBendDof; i++ )
  {
    force[i] = -e.gradient(i);
  }

}

void MNBendingForce::elementJacobian(const std::vector<Scalar>& undeformed,
  const std::vector<Scalar>& deformed,
  Eigen::Matrix<Scalar,NumMNBendDof,NumMNBendDof>& jac, 
  MNPrecomputed* pre) const
{

  jac.setZero();

  assert( pre );
  assert( undeformed.size() == NumMNBendDof);
  assert( undeformed.size() == deformed.size() );
//  assert( jac.numDof() == undeformed.size() );


  adreal<NumMNBendDof,1,Real> e = MNEnergy<1>(*this, undeformed, deformed, const_cast<MNPrecomputed*>(pre) );     
  // insert in the element Jacobian matrix
  for( uint i = 0; i < NumMNBendDof; i++ )
  {
    for( uint j = 0; j < NumMNBendDof; j++ )
    {
      jac(i,j) = -e.hessian(i,j);
    }
  }


}

} //namespace BASim
