#include "BASim/src/Physics/DeformableObjects/Rods/ElasticRodModel.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Math/Math.hh"
#include "BASim/src/Physics/DeformableObjects/Rods/RodModelStretchingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Rods/RodModelBendingForce.hh"
#include "BASim/src/Physics/DeformableObjects/Rods/RodModelTwistingForce.hh"

namespace BASim 
{  
  ElasticRodModel::ElasticRodModel(DeformableObject * object, const std::vector<EdgeHandle> & rodedges, Scalar timestep) : 
  PhysicalModel(*object), 
//  m_active_edges(rodedges),
  m_edge_stencils(),
  m_joint_stencils(),
  m_edge_active(object),
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
  m_forces(),
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
    // parse the edges in the context of the mesh and generate a list of edge and joint stencils
    // edge stencils
    m_edge_active.assign(0);    
    for (size_t i = 0; i < rodedges.size(); i++)
    {
      m_edge_active[rodedges[i]] = 1;

      // each rod edge forms an edge stencil
      EdgeStencil s;
      s.e = rodedges[i];
      EdgeVertexIterator evit = object->ev_iter(s.e);
      s.v1 = *evit; ++evit;
      s.v2 = *evit; ++evit;
      assert(!evit);
      
      m_edge_stencils.push_back(s);      
    }
    
    // joint stencils
    for (VertexIterator i = object->vertices_begin(); i != object->vertices_end(); ++i)
    {
      // detect all the rod edges that are incident to this vertex
      std::vector<EdgeHandle> active_incident_edges;
      for (VertexEdgeIterator veit = object->ve_iter(*i); veit; ++veit)
        if (m_edge_active[*veit])
          active_incident_edges.push_back(*veit);
      
      if (active_incident_edges.size() >= 2)
      {
        // 2 or more incident edges are rod edges: need to generate one or more joint stencils here
        // 2 edges: generate one joint stencil only
        // >2 edges: generate a ring of joint stencils
        for (size_t j = 0; j < (active_incident_edges.size() == 2 ? 1 : active_incident_edges.size()); j++)
        {
          JointStencil s;
          s.v2 = *i;
          s.e1 = active_incident_edges[j];
          s.e2 = active_incident_edges[(j + 1) % active_incident_edges.size()];

          EdgeVertexIterator evit;
          VertexHandle v1, v2;
          
          evit = object->ev_iter(s.e1);
          v1 = *evit; ++evit;
          v2 = *evit; ++evit;
          assert(!evit);
          s.v1 = (s.v2 == v1 ? v2 : v1);
          assert((s.v1 == v1 && s.v2 == v2) || (s.v1 == v2 && s.v2 == v1));
          
          evit = object->ev_iter(s.e2);
          v1 = *evit; ++evit;
          v2 = *evit; ++evit;
          assert(!evit);
          s.v3 = (s.v2 == v1 ? v2 : v1);
          assert((s.v3 == v1 && s.v2 == v2) || (s.v3 == v2 && s.v2 == v1));
          
          m_joint_stencils.push_back(s);
        }
      }
    }
  }
  
  void ElasticRodModel::setup(
             Scalar youngs, 
             Scalar youngs_damping, 
             Scalar shear, 
             Scalar shear_damping, 
             Scalar timestep, 
             const EdgeProperty<Vec3d> * undeformed_reference_director1)
  {    
    DeformableObject & obj = getDefoObj();

    // swap in the undeformed configuration as current configuration, because rod force initialization code assumes this
    VertexProperty<Vec3d> current_position_copy(obj.getVertexPositions());
    obj.setVertexPositions(obj.getVertexUndeformedPositions());
    EdgeProperty<Scalar> current_theta_copy(m_theta);
    m_theta = m_undef_theta;
    
    // adopt the reference directors specified by user, or generate them if not specified
    if (undeformed_reference_director1)
    {
      m_properties_reference_director1 = *undeformed_reference_director1; // director2 will be computed in updateProperties()
      for (size_t i = 0; i < m_edge_stencils.size(); i++)
      {
        EdgeStencil & s = m_edge_stencils[i];
        Vec3d t = (obj.getVertexPosition(s.v2) - obj.getVertexPosition(s.v1)).normalized();
        if (approxEq(m_properties_reference_director1[s.e].norm(), 0.0, 1e-6) || 
            !approxEq(m_properties_reference_director1[s.e].dot(t), 0.0, 1e-6)) // reject user specified if it's too short or not orthogonal to the edge tangent
          findOrthogonal(m_properties_reference_director1[s.e], t);
      }
    } else
    {
      for (size_t i = 0; i < m_edge_stencils.size(); i++)
      {
        EdgeStencil & s = m_edge_stencils[i];
        Vec3d t = (obj.getVertexPosition(s.v2) - obj.getVertexPosition(s.v1)).normalized();
        findOrthogonal(m_properties_reference_director1[s.e], t);
      }
    }

    // compute the derived properties, based on undeformed configuration and user sepcified undeformed reference directors
    updateProperties();
    
    // set up the three internal forces
    m_stretching_force = new RodModelStretchingForce(*this, m_edge_stencils, youngs, youngs_damping, timestep);
    m_bending_force = new RodModelBendingForce(*this, m_joint_stencils, youngs, youngs_damping, timestep);
    m_twisting_force = new RodModelTwistingForce(*this, m_joint_stencils, shear, shear_damping, timestep);
    
    addForce(m_stretching_force);
    addForce(m_bending_force);
    addForce(m_twisting_force);
    
    // restore the current configuration
    obj.setVertexPositions(current_position_copy);
    m_theta = current_theta_copy;

    // compute the derived properties again, this time based on current (i.e. initial) configuration, and the reference directors will be updated too
    updateProperties();
  }
  
  ElasticRodModel::~ElasticRodModel() 
  {
    delete m_stretching_force;
    delete m_bending_force;
    delete m_twisting_force;
  }
  
  void ElasticRodModel::computeForces(VecXd & force)
  {
/////////////////////
    VecXd curr_force(force.size());
    for (size_t i = 0; i < m_forces.size(); i++) 
    {
      curr_force.setZero();
      m_forces[i]->globalForce(curr_force);
      force += curr_force;
    }    
  }
  
  void ElasticRodModel::computeJacobian(Scalar scale, MatrixBase & J)
  {
/////////////////////
    for (size_t i = 0; i < m_forces.size(); i++)
      m_forces[i]->globalJacobian(scale, J);
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
    for (EdgeIterator i = getDefoObj().edges_begin(); i != getDefoObj().edges_end(); ++i)
    {
      m_radii[*i] = Vec2d(ra, rb);
      m_volumes[*i] = M_PI * ra * rb * getEdgeLength(*i);
    }
  }
  
  void ElasticRodModel::setDensity(Scalar density) 
  {
    m_density = density;
  }
  
  void ElasticRodModel::setEdgeThetas(const EdgeProperty<Scalar>& thetas)/////////////////////
  {
    m_theta = thetas;
  }
  
  void ElasticRodModel::setEdgeThetaVelocities(const EdgeProperty<Scalar>& vels)/////////////////////
  {
    m_theta_vel = vels;
  }

  void ElasticRodModel::setEdgeUndeformedThetas(const EdgeProperty<Scalar>& undef)/////////////////////
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
    //Compute vertex masses in a lumped mass way.
    
    m_vertex_masses.assign(0);
    m_edge_masses.assign(0);
    
    for (size_t i = 0; i < m_edge_stencils.size(); i++)
    {
      EdgeStencil & s = m_edge_stencils[i];
      
      Vec2d r = getRadii(s.e);
      Scalar mass = m_density * M_PI * r(0) * r(1) * getEdgeLength(s.e);
      Scalar inertia = 0.25 * mass * (square(r(0)) + square(r(1)));
      
      m_edge_masses[s.e] = inertia;
      m_vertex_masses[s.v1] += mass * 0.5;
      m_vertex_masses[s.v2] += mass * 0.5;
    }
    
    getDefoObj().updateVertexMasses();
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

    // viscous update
    for (size_t i = 0; i < m_forces.size(); i++)
      m_forces[i]->updateViscousReferenceStrain();
    
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
    
    std::cout << "Vertex count: " << getDefoObj().nv() << std::endl;
    
    // adjust radii based on edge length changes
    updateRadii();/////////////////////
    
    // update masses based on new edge length/radii
    computeMasses();
    
    // update stiffness, since radii have changed
    for (size_t i = 0; i < m_forces.size(); i++)
      m_forces[i]->updateStiffness();
    
    std::cout << "Completed endStep\n";
  }

  void ElasticRodModel::startIteration(Scalar time, Scalar timestep)
  {

  }
  
  void ElasticRodModel::endIteration(Scalar time, Scalar timestep)
  {
    // ElasticRod endIteration operation: upcate the derived properties
    updateProperties();
  }
  
  ////////////////////////////////////////
  void ElasticRodModel::updateRadii()
  {
    for (size_t i = 0; i < m_edge_stencils.size(); i++)
    {
      EdgeStencil & s = m_edge_stencils[i];
      Vec2d old_radii = m_radii[s.e];
      m_radii[s.e] = old_radii * sqrt(m_volumes[s.e] / (getEdgeLength(s.e) * M_PI * old_radii(0) * old_radii(1)));
    }
  }

  void ElasticRodModel::updateProperties()
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
//    for (VertexIterator i = getDefoObj().vertices_begin(); i != getDefoObj().vertices_end(); ++i) 
//    {
//      m_properties_voronoi_length[*i] = 0;
//    }
    m_properties_voronoi_length.assign(0);
    
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
    
    for (size_t i = 0; i < m_forces.size(); i++)
      m_forces[i]->updateProperties();
  }

} //namespace BASim
