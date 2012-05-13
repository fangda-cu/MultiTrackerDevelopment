#include "BASim/src/Physics/DeformableObjects/Rods/ElasticRodModel.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/Math.hh"

namespace BASim 
{  
  ElasticRodModel::ElasticRodModel(DeformableObject* object, const std::vector<EdgeHandle> & rodedges, Scalar timestep) : 
  PhysicalModel(*object), 
//  m_active_edges(rodedges),
  m_edge_stencils(),
  m_joint_stencils(),
  m_obj(object), 
//  m_active_faces(shellFaces), /////////////////////
  m_theta(object), /////////////////////
  m_theta_vel(object),/////////////////////
  m_undef_theta(object),
  m_damping_undef_theta(object),
  m_vertex_masses(object),
  m_edge_masses(object),
  m_radii(object),
  m_volumes(object),
  m_density(1),
  m_properties_edge(object),
  m_properties_edge_tangent(object),
  m_properties_edge_length(object),
  m_properties_reference_director1(object),
  m_properties_reference_director2(object),
  m_properties_material_director1(object),
  m_properties_material_director2(object),
  m_properties_voronoi_length(object),
  m_properties_reference_twist(object),
  m_properties_curvature_binormal(object)
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
  
  void ElasticRodModel::setEdgeThetas( const EdgeProperty<Scalar>& positions )/////////////////////
  {
    m_theta = positions;
  }
  
  void ElasticRodModel::setEdgeThetaVelocities(const EdgeProperty<Scalar>& velocities)/////////////////////
  {
    m_theta_vel = velocities;
  }

  void ElasticRodModel::setEdgeUndeformedThetas( const EdgeProperty<Scalar>& undef )/////////////////////
  {
    m_undef_theta = undef;
  }  

//  Scalar ElasticRodModel::getThickness(const VertexHandle& vh) const 
//  {
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
//  }
//  
//  Scalar ElasticRodModel::getMaxThickness() const 
//  {
/////////////////////
//    Scalar maxVal = -1000000;
//    for( FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit ){
//      if (m_thicknesses[*fit] > maxVal ) maxVal = m_thicknesses[*fit];
//    }
//    return maxVal;
//  }
//  
//  Scalar ElasticRodModel::getMinThickness() const 
//  {
/////////////////////
//    Scalar minVal = 1000000;
//    for( FaceIterator fit = m_obj->faces_begin(); fit != m_obj->faces_end(); ++fit ){
//      if (m_thicknesses[*fit] < minVal ) minVal = m_thicknesses[*fit];
//    }
//    return minVal;
//  }
//  
//  void ElasticRodModel::getThickness(VertexProperty<Scalar> & vThickness) const
//  {
/////////////////////
//    for ( VertexIterator vit = m_obj->vertices_begin(); vit != m_obj->vertices_end(); ++vit){
//      vThickness[*vit] = getThickness(*vit);
//    }
//  }
//  
//  Scalar ElasticRodModel::getArea(const FaceHandle& f, bool current) const  
//  {
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
//  }
  
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
  bool ElasticRodModel::isVertexActive( const VertexHandle& v ) const
  {
//    //determine if the vertex is on any active face
//    VertexFaceIterator vf = m_obj->vf_iter(v);
//    for(;vf; ++vf) {
//      if(isFaceActive(*vf)) {
//        return true;
//      }
//    }
//    
//    return false;
  }
  
  bool ElasticRodModel::isEdgeActive( const EdgeHandle& e) const 
  {
//    //if any adjacent face is active, we say this edge is active.
//    EdgeFaceIterator ef = m_obj->ef_iter(e);
//    for(;ef;++ef) {
//      if(isFaceActive(*ef)) {
//        return true;
//      }
//    }
//    
//    return false;
  }
  
  const Scalar& ElasticRodModel::getDof( const DofHandle& hnd ) const
  {
    //they're all edge Dofs for a rod
    assert(hnd.getType() == DofHandle::EDGE_DOF);
    
    //return reference to the appropriate position in the vector
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    return const_cast<Scalar&>(m_theta[eh]);/////////////////////
  }
  
  void ElasticRodModel::setDof( const DofHandle& hnd, const Scalar& dof )
  {
    //they're all vertex Dofs for a rod
    assert(hnd.getType() == DofHandle::EDGE_DOF);
    
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    m_theta[eh] = dof;/////////////////////
  }
  
  const Scalar& ElasticRodModel::getVel( const DofHandle& hnd ) const
  {
    assert(hnd.getType() == DofHandle::EDGE_DOF);
    
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    return const_cast<Scalar&>(m_theta_vel[eh]);/////////////////////
  }
  
  void ElasticRodModel::setVel( const DofHandle& hnd, const Scalar& vel )
  {
    assert(hnd.getType() == DofHandle::EDGE_DOF);
    
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    m_theta_vel[eh] = vel;/////////////////////
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

    startIteration(time, timestep);
    
    //update the damping "reference configuration" for computing viscous forces.
    m_damping_undef_theta = m_theta;/////////////////////
    
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
//    bool do_relabel = false;/////////////////////
    
    endIteration(time, timestep);
    
    std::cout << "Vertex count: " << m_obj->nv() << std::endl;
    
    //Adjust thicknesses based on area changes
    updateRadii();/////////////////////
    
    //Update masses based on new areas/thicknesses
    computeMasses();
    
    std::cout << "Completed endStep\n";
  }

  void ElasticRodModel::startIteration(Scalar time, Scalar timestep)
  {
    
  }
  
  void ElasticRodModel::endIteration(Scalar time, Scalar timestep)
  {
//    // copy position dofs and edge theta dofs into ElasticRod's local copy
//    
//    
//    // ElasticRod endIteration operation: upcate the derived properties
//    m_elastic_rod.updateProperties();
  }
  
  ////////////////////////////////////////
  void ElasticRodModel::updateRadii()
  {
    
  }

  void ElasticRodModel::upateProperties()
  {
    // This code is adapted from BASim::ElasticRod::updateProperties(). The order of the computation is preserved. No optimization applied.
    DeformableObject & obj = getDefoObj();
    
    // compute edges
    for (size_t i = 0; i < m_edge_stencils.size(); i++)
    {
      EdgeStencil & s = m_edge_stencils[i];
      m_properties_edge[s.e] = obj.getVertexPosition(s.v2) - obj.getVertexPosition(s.v1);
    }      
    
    // compoute reference directors (Time Parallel only; code from BASim::ElasticRod::computeTimeParallel())
    for (size_t i = 0; i < m_edge_stencils.size(); i++)
    {
      EdgeStencil & s = m_edge_stencils[i];
      Vec3d t = getEdge(s.e).normalized();
      Vec3d u = parallel_transport(getReferenceDirector1(s.e), getEdgeTangent(s.e), t);
      u = (u - u.dot(t) * t).normalized();
      m_properties_reference_director1[s.e] = u;
      m_properties_reference_director2[s.e] = t.cross(u);
    }
      
    // compute edge tangents
    for (size_t i = 0; i < m_edge_stencils.size(); i++)
    {
      EdgeStencil & s = m_edge_stencils[i];
      m_properties_edge_tangent[s.e] = getEdge(s.e).normalized();
    }
    
    // compute reference twists (code from BASim::ElasticRod::computeReferenceTwist())
    for (size_t i = 0; i < m_joint_stencils.size(); i++)
    {
      JointStencil & s = m_joint_stencils[i];
      
      Scalar referenceTwist = getReferenceTwist(s.v2);// std::cout << "previous referenceTwist = " << referenceTwist << '\n';
      const Vec3d & u0 = getReferenceDirector1(s.e1);// std::cout << "u0 = " << u0 << '\n';
      const Vec3d & u1 = getReferenceDirector1(s.e2);// std::cout << "u1 = " << u1 << '\n';
      const Vec3d & tangent0 = getEdgeTangent(s.e1);// std::cout << "t0 = " << t0 << '\n';
      const Vec3d & tangent1 = getEdgeTangent(s.e2);// std::cout << "t1 = " << t1 << '\n';
      
      // transport reference frame to next edge
      Vec3d ut = parallel_transport(u0, tangent0, tangent1);// std::cout << "ut = " << ut << '\n';
      
      // rotate by current value of reference twist
      rotateAxisAngle(ut, tangent1, referenceTwist);// std::cout << "ut = " << ut << '\n';
      
      // compute increment to reference twist to align reference frames
      referenceTwist += signedAngle(ut, u1, tangent1);// std::cout << "referenceTwist = " << referenceTwist << '\n';
      
      m_properties_reference_twist[s.v2] = referenceTwist;
    }    
    
    // compute curvature binormals
    for (size_t i = 0; i < m_joint_stencils.size(); i++)
    {
      JointStencil & s = m_joint_stencils[i];
      Vec3d & t0 = getEdgeTangent(s.e1);
      Vec3d & t1 = getEdgeTangent(s.e2);
      assert(approxEq(t0.norm(), 1.0));
      assert(approxEq(t1.norm(), 1.0));
      m_properties_curvature_binormal[s.v2] = 2.0 * t0.cross(t1) / (1.0 + t0.dot(t1));
    }    

    // compute edge lengths    
    for (size_t i = 0; i < m_edge_stencils.size(); i++)
    {
      EdgeStencil & s = m_edge_stencils[i];
      m_properties_edge_length[s.e] = m_properties_edge[s.e].norm();
    }
    
    // compute voronoi lengths (this code is different than BASim::ElasticRod::computeVoronoiLength() due to topology)
    for (VertexIterator i = getDefoObj().vertices_begin(); i != getDefoObj().vertices_end(); ++i) 
    {
      m_properties_voronoi_length[*i] = 0;
    }
    
    for (size_t i = 0; i < m_edge_stencils.size(); i++)
    {
      EdgeStencil & s = m_edge_stencils[i];
      m_properties_voronoi_length[s.v1] += getEdgeLength(s.e) * 0.5;
      m_properties_voronoi_length[s.v2] += getEdgeLength(s.e) * 0.5;
    }    

    // compute material directors (code from BASim::ElasticRod::computeMaterialDirectors())
    for (size_t i = 0; i < m_edge_stencils.size(); i++)
    {
      EdgeStencil & s = m_edge_stencils[i];
      Scalar ca = cos(getEdgeTheta(s.e));
      Scalar sa = sin(getEdgeTheta(s.e));
      const Vec3d& u = getReferenceDirector1(s.e);
      const Vec3d& v = getReferenceDirector2(s.e);
      m_properties_material_director1[s.e] = ca * u + sa * v;
      m_properties_material_director1[s.e] = -sa * u + ca * v;
    }
    
//    updateForceProperties();

  }

} //namespace BASim
