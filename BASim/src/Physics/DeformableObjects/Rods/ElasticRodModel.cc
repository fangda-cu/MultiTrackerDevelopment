#include "BASim/src/Physics/DeformableObjects/Rods/ElasticRodModel.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/Math.hh"

namespace BASim 
{  
  ElasticRodModel::ElasticRodModel(DeformableObject* object, const std::vector<EdgeHandle> & rodedge, Scalar timestep) : 
  PhysicalModel(*object), m_obj(object), 
//  m_active_faces(shellFaces), /////////////////////
  m_undef_xi(object),
  m_damping_undef_xi(object),
  m_vertex_masses(object),
  m_edge_masses(object),
  m_radii(object),
  m_volumes(object),
  m_xi(object), /////////////////////
  m_xi_vel(object),/////////////////////
  m_density(1)
  {
  
  }
  
  ElasticRodModel::~ElasticRodModel() 
  {
  
  }
  
  void ElasticRodModel::computeForces( VecXd& force )
  {
/////////////////////
//    const std::vector<ElasticShellForce*>& forces = getForces();
//    std::vector<ElasticShellForce*>::const_iterator fIt;
//    
//    VecXd curr_force(force.size());
//    for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
//      curr_force.setZero();
//      (*fIt)->globalForce(curr_force);
//      
//      force += curr_force;
//    }
    
  }
  
  void ElasticRodModel::computeJacobian( Scalar scale, MatrixBase& J )
  {
/////////////////////
//    const std::vector<ElasticShellForce*>& forces = getForces();
//    std::vector<ElasticShellForce*>::const_iterator fIt;
//    
//    for (fIt = forces.begin(); fIt != forces.end(); ++fIt)
//      (*fIt)->globalJacobian(scale, J);
  }
  
/////////////////////
//  const std::vector<ElasticShellForce*>& ElasticShell::getForces() const
//  {
//    return m_shell_forces;
//  }
  
/////////////////////
//  void ElasticShell::addForce( ElasticShellForce* force )
//  {
//    assert(force != NULL);
//    
//    m_shell_forces.push_back(force);
//  }
  
  void ElasticRodModel::setRadii(Scalar ra, Scalar rb)
  {
/////////////////////
//    for (FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit) 
//    {
//      m_thicknesses[*fit] = thickness;
//      Scalar area = getArea(*fit, false);
//      m_volumes[*fit] = thickness * area;
//    }
    
  }
  
  void ElasticRodModel::setDensity(Scalar density) 
  {
    m_density = density;
  }
  
  void ElasticRodModel::setEdgeUndeformed( const EdgeProperty<Scalar>& undef )/////////////////////
  {
    m_undef_xi = undef;
  }
  
  void ElasticRodModel::setEdgeXis( const EdgeProperty<Scalar>& positions )/////////////////////
  {
    m_xi = positions;
  }
  
  void ElasticRodModel::setEdgeVelocities(const EdgeProperty<Scalar>& velocities)/////////////////////
  {
    m_xi_vel = velocities;
  }
  
  Scalar ElasticRodModel::getThickness(const VertexHandle& vh) const 
  {
/////////////////////
//    Scalar totalA = 0.0;
//    Scalar w;
//    Scalar total = 0.0;
//    for (VertexFaceIterator vfit = m_obj->vf_iter(vh); vfit; ++vfit){
//      w = getArea(*vfit);
//      totalA += w;
//      total += w*m_thicknesses[*vfit];
//    }
//    
//    assert ( totalA > 0.);
//    assert ( total > 0.);
//    return total / totalA;
  }
  
  Scalar ElasticRodModel::getMaxThickness() const 
  {
/////////////////////
//    Scalar maxVal = -1000000;
//    for( FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit ){
//      if (m_thicknesses[*fit] > maxVal ) maxVal = m_thicknesses[*fit];
//    }
//    return maxVal;
  }
  
  Scalar ElasticRodModel::getMinThickness() const 
  {
/////////////////////
//    Scalar minVal = 1000000;
//    for( FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit ){
//      if (m_thicknesses[*fit] < minVal ) minVal = m_thicknesses[*fit];
//    }
//    return minVal;
  }
  
  void ElasticRodModel::getThickness(VertexProperty<Scalar> & vThickness) const
  {
/////////////////////
//    for ( VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit){
//      vThickness[*vit] = getThickness(*vit);
//    }
  }
  
  Scalar ElasticRodModel::getArea(const FaceHandle& f, bool current) const  
  {
/////////////////////
//    FaceVertexIterator fvit = m_obj->fv_iter(f);
//    VertexHandle v0_hnd = *fvit; ++fvit; assert(fvit);
//    VertexHandle v1_hnd = *fvit; ++fvit; assert(fvit);
//    VertexHandle v2_hnd = *fvit; ++fvit; assert(!fvit);
//    
//    //compute triangle areas
//    if(current) 
//    {
//      Vec3d pos0 = m_obj->getVertexPosition(v0_hnd);
//      Vec3d pos1 = m_obj->getVertexPosition(v1_hnd);
//      Vec3d pos2 = m_obj->getVertexPosition(v2_hnd);
//      
//      Vec3d v0 = pos1 - pos0;
//      Vec3d v1 = pos2 - pos0;
//      Vec3d triVec = v0.cross(v1);
//      return 0.5*triVec.norm();
//    }
//    else 
//    {
//      Vec3d pos0 = m_obj->getVertexUndeformedPosition(v0_hnd);
//      Vec3d pos1 = m_obj->getVertexUndeformedPosition(v1_hnd);
//      Vec3d pos2 = m_obj->getVertexUndeformedPosition(v2_hnd);
//      
//      Vec3d v0 = pos1 - pos0;
//      Vec3d v1 = pos2 - pos0;
//      Vec3d triVec = v0.cross(v1);
//      return 0.5*triVec.norm();
//    }
  }
  
  void ElasticRodModel::computeMasses()
  {
/////////////////////
//    //Compute vertex masses in a lumped mass way.
//    
//    m_vertex_masses.assign(0);
//    m_edge_masses.assign(0);
//    
//    Scalar area = 0;
//    
//    //Iterate over all triangles active in this shell and accumulate vertex masses
//    for(FaceIterator f_iter = m_obj->faces_begin(); f_iter != m_obj->faces_end(); ++f_iter) {
//      FaceHandle& f_hnd = *f_iter;
//      if(m_active_faces[f_hnd]) {
//        
//        //get the three vertices
//        FaceVertexIterator fvit = m_obj->fv_iter(f_hnd);
//        VertexHandle v0_hnd = *fvit; ++fvit; assert(fvit);
//        VertexHandle v1_hnd = *fvit; ++fvit; assert(fvit);
//        VertexHandle v2_hnd = *fvit; ++fvit; assert(!fvit);
//        
//        //compute triangle areas
//        Vec3d v0 = getVertexPosition(v1_hnd) - getVertexPosition(v0_hnd);
//        Vec3d v1 = getVertexPosition(v2_hnd) - getVertexPosition(v0_hnd);
//        Vec3d triVec = v0.cross(v1);
//        Scalar area = 0.5*sqrt(triVec.dot(triVec)) / 3.0;
//        Scalar contribution = m_thicknesses[f_hnd] * m_density * area;
//        
//        //accumulate mass to the vertices
//        m_vertex_masses[v0_hnd] += contribution;
//        m_vertex_masses[v1_hnd] += contribution;
//        m_vertex_masses[v2_hnd] += contribution;
//        
//        //also accumulate mass to the edges (this mass computation is probably not consistent with what we want)
//        FaceEdgeIterator feit = m_obj->fe_iter(f_hnd);
//        EdgeHandle e0_hnd = *feit; ++feit; assert(feit);
//        EdgeHandle e1_hnd = *feit; ++feit; assert(feit);
//        EdgeHandle e2_hnd = *feit; ++feit; //assert(feit);
//        
//        m_edge_masses[e0_hnd] += contribution;
//        m_edge_masses[e1_hnd] += contribution;
//        m_edge_masses[e2_hnd] += contribution;
//        
//        //store the current volumes
//        m_volumes[f_hnd] = 3*area*m_thicknesses[f_hnd];
//      }
//    }
//    
//    m_obj->updateVertexMasses();
  }
  
/////////////////////
//  bool ElasticShell::isVertexActive( const VertexHandle& v ) const
//  {
//    //determine if the vertex is on any active face
//    VertexFaceIterator vf = m_obj->vf_iter(v);
//    for(;vf; ++vf) {
//      if(isFaceActive(*vf)) {
//        return true;
//      }
//    }
//    
//    return false;
//  }
//  
//  bool ElasticShell::isEdgeActive( const EdgeHandle& e) const {
//    //if any adjacent face is active, we say this edge is active.
//    EdgeFaceIterator ef = m_obj->ef_iter(e);
//    for(;ef;++ef) {
//      if(isFaceActive(*ef)) {
//        return true;
//      }
//    }
//    
//    return false;
//  }
  
  const Scalar& ElasticRodModel::getDof( const DofHandle& hnd ) const
  {
    //they're all edge Dofs for a rod
    assert(hnd.getType() == DofHandle::EDGE_DOF);
    
    //return reference to the appropriate position in the vector
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    return const_cast<Scalar&>(m_xi[eh]);/////////////////////
  }
  
  void ElasticRodModel::setDof( const DofHandle& hnd, const Scalar& dof )
  {
    //they're all vertex Dofs for a rod
    assert(hnd.getType() == DofHandle::EDGE_DOF);
    
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    m_xi[eh] = dof;/////////////////////
  }
  
  const Scalar& ElasticRodModel::getVel( const DofHandle& hnd ) const
  {
    assert(hnd.getType() == DofHandle::EDGE_DOF);
    
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    return const_cast<Scalar&>(m_xi_vel[eh]);/////////////////////
  }
  
  void ElasticRodModel::setVel( const DofHandle& hnd, const Scalar& vel )
  {
    assert(hnd.getType() == DofHandle::EDGE_DOF);
    
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    m_xi_vel[eh] = vel;/////////////////////
  }
  
  const Scalar& ElasticRodModel::getMass( const DofHandle& hnd ) const
  {
    assert(hnd.getType() == DofHandle::EDGE_DOF);
    
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    return m_edge_masses[eh];/////////////////////
  }
  
  void ElasticRodModel::startStep(Scalar time, Scalar timestep)
  {
    std::cout << "Starting startStep\n";

    //update the damping "reference configuration" for computing viscous forces.
    m_damping_undef_xi = m_xi;/////////////////////
    
    //tell the forces to update anything they need to update
/////////////////////
//    const std::vector<ElasticShellForce*>& forces = getForces();
//    for(unsigned int i = 0; i < forces.size(); ++i)
//      forces[i]->update();

    std::cout << "Done startStep\n";
  }
    
  void ElasticRodModel::endStep(Scalar time, Scalar timestep) 
  {
    std::cout << "Starting endStep.\n";
    bool do_relabel = false;
    
    std::cout << "Vertex count: " << m_obj->nv() << std::endl;
    
    //Adjust thicknesses based on area changes
    updateThickness();/////////////////////
    
    //Update masses based on new areas/thicknesses
    computeMasses();
    
    std::cout << "Completed endStep\n";
  }
    
} //namespace BASim
