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
  m_edge_stencils(),
  m_joint_stencils(),
  m_triedge_stencils(),
  m_edge_active(object),
  m_theta(object),
  m_theta_vel(object),
  m_undef_theta(object),
  m_damping_undef_theta(object),
  m_vertex_masses(object),
  m_edge_masses(object),
  m_radii(object),
  m_volumes(object),
  m_density(1),
  m_forces(),
  m_stretching_force(NULL),
  m_bending_force(NULL),
  m_twisting_force(NULL),
  m_properties_edge(object),
  m_properties_edge_tangent(object),
  m_properties_edge_length(object),
  m_properties_reference_director1(object),
  m_properties_reference_director2(object),
  m_properties_material_director1(object),
  m_properties_material_director2(object),
  m_undeformed_reference_director1(NULL),
  m_undeformed_positions(NULL)
  {
    // parse the edges in the context of the mesh and generate a list of edge and joint stencils
    // edge stencils
    m_edge_active.assign(0);  
    for (size_t i = 0; i < rodedges.size(); i++)
    {
      m_edge_active[rodedges[i]] = 1;
    }
    
    buildStencils();

    
  }
  
  void ElasticRodModel::setUndeformedReferenceDirector1(const EdgeProperty<Vec3d> & undeformed_reference_director1)
  {
    m_undeformed_reference_director1 = new EdgeProperty<Vec3d>(undeformed_reference_director1);
  }
  
  void ElasticRodModel::setUndeformedPositions(const VertexProperty<Vec3d> & undeformed_positions)
  {
    m_undeformed_positions = new VertexProperty<Vec3d>(undeformed_positions);
  }
  
  void ElasticRodModel::setup(Scalar youngs, Scalar youngs_damping, Scalar shear, Scalar shear_damping, Scalar timestep)
  {    
    
    assignStencilDofs();

    DeformableObject & obj = getDefoObj();

    // swap in the undeformed configuration as current configuration, because rod force initialization code assumes this
    VertexProperty<Vec3d> current_position_copy(obj.getVertexPositions());
    if (m_undeformed_positions)
      obj.setVertexPositions(*m_undeformed_positions);
    else
      obj.setVertexPositions(obj.getVertexUndeformedPositions());
    EdgeProperty<Scalar> current_theta_copy(m_theta);
    m_theta = m_undef_theta;
    
    // adopt the reference directors specified by user, or generate them if not specified
    if (m_undeformed_reference_director1)
    {
      m_properties_reference_director1 = *m_undeformed_reference_director1; // director2 will be computed in updateProperties()
      for (size_t i = 0; i < m_edge_stencils.size(); i++)
      {
        EdgeStencil & s = m_edge_stencils[i];
        Vec3d t = (obj.getVertexPosition(s.v2) - obj.getVertexPosition(s.v1)).normalized();
        if (approxEq(m_properties_reference_director1[s.e].norm(), 0.0, 1e-6) || 
            !approxEq(m_properties_reference_director1[s.e].dot(t), 0.0, 1e-6)) // reject user specified if it's too short or not orthogonal to the edge tangent
          findOrthogonal(m_properties_reference_director1[s.e], t);
      }
      delete m_undeformed_reference_director1;
    } else
    {
      for (size_t i = 0; i < m_edge_stencils.size(); i++)
      {
        EdgeStencil & s = m_edge_stencils[i];
        Vec3d t = (obj.getVertexPosition(s.v2) - obj.getVertexPosition(s.v1)).normalized();
        findOrthogonal(m_properties_reference_director1[s.e], t);
      }
    }

    // compute the derived properties, based on undeformed configuration and user specified undeformed reference directors
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
  
  void ElasticRodModel::computeConservativeForcesEnergy(VecXd & force, Scalar& energy)
  {
    VecXd curr_force(force.size());
    for (size_t i = 0; i < m_forces.size(); i++) 
    {
      curr_force.setZero();
      m_forces[i]->globalForce(curr_force);
      force += curr_force;
      energy += m_forces[i]->globalEnergy();
    }    
  }

  void ElasticRodModel::computeForces(VecXd & force)
  {
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
    for (size_t i = 0; i < m_forces.size(); i++)
      m_forces[i]->globalJacobian(scale, J);
  }
  
  void ElasticRodModel::setRadii(Scalar ra, Scalar rb)
  {
    for (EdgeIterator i = getDefoObj().edges_begin(); i != getDefoObj().edges_end(); ++i)
    {
      m_radii[*i] = Vec2d(ra, rb);
      m_volumes[*i] = M_PI * ra * rb * getEdgeLength(*i);
    }
  }

  void ElasticRodModel::setRadii(const EdgeProperty<Vec2d>& radii)
  {
    for (EdgeIterator i = getDefoObj().edges_begin(); i != getDefoObj().edges_end(); ++i)
    {
      if(isEdgeActive(*i)) {
        Vec2d rad = radii[*i];
        m_radii[*i] = rad;
        m_volumes[*i] = M_PI * rad[0] * rad[1] * getEdgeLength(*i);
      }
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

  void ElasticRodModel::computeMasses()
  {
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
    return const_cast<Scalar&>(m_theta[eh]);
  }
  
  void ElasticRodModel::setDof( const DofHandle& hnd, const Scalar& dof )
  {
    //they're all edge Dofs for a rod
    assert(hnd.getType() == DofHandle::EDGE_DOF);
    
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    m_theta[eh] = dof;
  }
  
  const Scalar& ElasticRodModel::getVel( const DofHandle& hnd ) const
  {
    assert(hnd.getType() == DofHandle::EDGE_DOF);
    
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    return const_cast<Scalar&>(m_theta_vel[eh]);
  }
  
  void ElasticRodModel::setVel( const DofHandle& hnd, const Scalar& vel )
  {
    assert(hnd.getType() == DofHandle::EDGE_DOF);
    
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    m_theta_vel[eh] = vel;
  }
  
  const Scalar& ElasticRodModel::getMass( const DofHandle& hnd ) const
  {
    assert(hnd.getType() == DofHandle::EDGE_DOF);
    
    const EdgeHandle& eh = static_cast<const EdgeHandle&>(hnd.getHandle());
    return m_edge_masses[eh];
  }
  
  void ElasticRodModel::startStep(Scalar time, Scalar timestep)
  {
    // viscous update
    for (size_t i = 0; i < m_forces.size(); i++)
      m_forces[i]->updateViscousReferenceStrain();
    
    //update the damping "reference configuration" for computing viscous forces.
    m_damping_undef_theta = m_theta;
    
  }
    
  void ElasticRodModel::endStep(Scalar time, Scalar timestep) 
  {
    // adjust radii based on edge length changes
    updateRadii();
    
    // update masses based on new edge length/radii
    computeMasses();
    
    // update stiffness, since radii have changed
    for (size_t i = 0; i < m_forces.size(); i++)
      m_forces[i]->updateStiffness();
  }

  void ElasticRodModel::startIteration(Scalar time, Scalar timestep)
  {

  }
  
  void ElasticRodModel::endIteration(Scalar time, Scalar timestep)
  {
    // ElasticRod endIteration operation: update the derived properties
    updateProperties();
  }
  
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
    
    // compute edges (from update vertex positions)
    for (size_t i = 0; i < m_edge_stencils.size(); i++)
    {
      EdgeStencil & s = m_edge_stencils[i];
      m_properties_edge[s.e] = obj.getVertexPosition(s.v2) - obj.getVertexPosition(s.v1);
    }      
    
    // compute reference directors (Time Parallel only; code from BASim::ElasticRod::computeTimeParallel())
    for (size_t i = 0; i < m_edge_stencils.size(); i++)
    {
      EdgeStencil & s = m_edge_stencils[i];
      Vec3d t = getEdge(s.e).normalized();
      
      // parallel transport the reference director from the old edge tangent, to the new one, t
      Vec3d u = parallel_transport(getReferenceDirector1(s.e), getEdgeTangent(s.e), t);
      
      // remove the component along the tangent
      u = (u - u.dot(t) * t).normalized();

      //assign the new reference directors
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
      
      Scalar & referenceTwist = s.referenceTwist;
      const Vec3d & u0 = getReferenceDirector1(s.e1);
      const Vec3d & u1 = getReferenceDirector1(s.e2);
      const Vec3d   tangent0 = getEdgeTangent(s.e1) * (s.e1flip ? -1 : 1);
      const Vec3d   tangent1 = getEdgeTangent(s.e2) * (s.e2flip ? -1 : 1);
      
      // transport reference frame to next edge
      Vec3d ut = parallel_transport(u0, tangent0, tangent1);// std::cout << "ut = " << ut << '\n';
      
      // rotate by current value of reference twist
      rotateAxisAngle(ut, tangent1, referenceTwist);// std::cout << "ut = " << ut << '\n';
      
      // compute increment to reference twist to align reference frames
      referenceTwist += signedAngle(ut, u1, tangent1);// std::cout << "referenceTwist = " << referenceTwist << '\n';
    }
    
    // compute curvature binormals
    for (size_t i = 0; i < m_joint_stencils.size(); i++)
    {
      JointStencil & s = m_joint_stencils[i];
      Vec3d t0 = getEdgeTangent(s.e1) * (s.e1flip ? -1 : 1);
      Vec3d t1 = getEdgeTangent(s.e2) * (s.e2flip ? -1 : 1);
      assert(approxEq(t0.norm(), 1.0));
      assert(approxEq(t1.norm(), 1.0));
      s.curvatureBinormal = 2.0 * t0.cross(t1) / (1.0 + t0.dot(t1));
    }    

    // compute edge lengths    
    for (size_t i = 0; i < m_edge_stencils.size(); i++)
    {
      EdgeStencil & s = m_edge_stencils[i];
      m_properties_edge_length[s.e] = m_properties_edge[s.e].norm();
    }
    
    // compute voronoi lengths (this code is different than BASim::ElasticRod::computeVoronoiLength() due to topology)
    for (size_t i = 0; i < m_joint_stencils.size(); i++)
    {
      JointStencil & s = m_joint_stencils[i];
      s.voronoiLength = 0.5 * (getEdgeLength(s.e1) + getEdgeLength(s.e2));
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
      m_properties_material_director2[s.e] = -sa * u + ca * v;
    }
    
    // propagate newly computed stencil data to internal force stencils
    for (size_t i = 0; i < m_edge_stencils.size(); i++)
    {
      EdgeStencil & s = m_edge_stencils[i];
      if (m_stretching_force && m_stretching_force->stencils().size() == m_edge_stencils.size())
        m_stretching_force->stencils()[s.id].copyData(s);
    }
    
    for (size_t i = 0; i < m_joint_stencils.size(); i++)
    {
      JointStencil & s = m_joint_stencils[i];
      if (m_bending_force && m_bending_force->stencils().size() == m_joint_stencils.size())
        m_bending_force->stencils()[s.id].copyData(s);
      if (m_twisting_force && m_twisting_force->stencils().size() == m_joint_stencils.size())
        m_twisting_force->stencils()[s.id].copyData(s);
    }
    
    // forces' updateProperties()
    for (size_t i = 0; i < m_forces.size(); i++)
      m_forces[i]->updateProperties();
  }

  void ElasticRodModel::getScriptedDofs(IntArray & dofIndices, std::vector<Scalar> & dofValues, Scalar time) const
  {
    for (std::vector<std::pair<EdgeHandle, Scalar> >::const_iterator i = m_edge_constraints.begin(); i != m_edge_constraints.end(); i++)
    {
      dofIndices.push_back(getEdgeDofBase(i->first));
      dofValues.push_back(i->second);
    }
    for (std::vector<std::pair<EdgeHandle, Vec3d> >::const_iterator i = m_edge_vel_constraints.begin(); i != m_edge_vel_constraints.end(); i++)
    {
      dofIndices.push_back(getEdgeDofBase(i->first));
      dofValues.push_back(i->second(0) + i->second(1) * (time - i->second(2)));
    }  
  }
  
  // scripting on edge dofs
  void ElasticRodModel::constrainEdge(const EdgeHandle & e, Scalar t)
  {
    m_edge_constraints.push_back(std::pair<EdgeHandle, Scalar>(e, t));
  }
  
  void ElasticRodModel::constrainEdgeVel(const EdgeHandle & e, Scalar init_value, Scalar velocity, Scalar start_time)
  {
    m_edge_vel_constraints.push_back(std::pair<EdgeHandle, Vec3d>(e, Vec3d(init_value, velocity, start_time)));
  }
  
  void ElasticRodModel::releaseEdge(const EdgeHandle & e)
  {
    // assume no duplicates
    for (std::vector<std::pair<EdgeHandle, Scalar> >::iterator i = m_edge_constraints.begin(); i != m_edge_constraints.end(); i++)
      if (i->first == e)
      {
        m_edge_constraints.erase(i);
        return;
      }
    for (std::vector<std::pair<EdgeHandle, Vec3d> >::iterator i = m_edge_vel_constraints.begin(); i != m_edge_vel_constraints.end(); i++)
      if (i->first == e)
      {
        m_edge_vel_constraints.erase(i);
        return;
      }
  }
  
  bool ElasticRodModel::isConstrained(const EdgeHandle & e) const
  {
    for (std::vector<std::pair<EdgeHandle, Scalar> >::const_iterator i = m_edge_constraints.begin(); i != m_edge_constraints.end(); i++)
      if (i->first == e)
        return true;
    for (std::vector<std::pair<EdgeHandle, Vec3d> >::const_iterator i = m_edge_vel_constraints.begin(); i != m_edge_vel_constraints.end(); i++)
      if (i->first == e)
        return true;
    return false;
  }

  void ElasticRodModel::remesh( Scalar min_length, Scalar max_length )
  {
    DeformableObject& obj = getDefoObj();

    //look at each edge.
    //if it's too long, split it
    //if it's too short, collapse together with shorter neighbour
    
    std::vector<EdgeHandle> edges_to_split;
    for(EdgeIterator eit = obj.edges_begin(); eit != obj.edges_end(); ++eit) {
      EdgeHandle eh = *eit;
      
      if(!m_edge_active[eh]) //ignore non-rod edges
        continue;

      Scalar len = getEdgeLength(eh);
      if(len > max_length || len < min_length) {
        edges_to_split.push_back(eh);
      }
    }
    

    for(unsigned int i = 0; i < edges_to_split.size(); ++i) {
      EdgeHandle eh = edges_to_split[i];
      
      //the edge has been deleted
      if(!obj.edgeExists(eh)) 
        continue;
      assert(eh.isValid());

      //compute edge length *from scratch* in case anything has changed due to remeshing
      VertexHandle vh0 = obj.fromVertex(eh); VertexHandle vh1 = obj.toVertex(eh);
      Scalar len = (obj.getVertexPosition(vh0) - obj.getVertexPosition(vh1)).norm();
      
      if(len > max_length) {
        subdivideEdge(eh);
      }
      else if(len < min_length) {
        collapseEdge(eh);
      }
    }

    //construct the stencil set from scratch
    buildStencils();

    //rebuild the dof indices as a result
    obj.computeDofIndexing();

    //next we need to set the dof indices in all the stencils...
    assignStencilDofs();

    //then we need to update the undeformed configurations used by the various forces:
    
    // swap in the undeformed configuration as current configuration, because rod force initialization code assumes this
    VertexProperty<Vec3d> current_position_copy(obj.getVertexPositions());
    obj.setVertexPositions(obj.getVertexUndeformedPositions()); //(remeshing here assumes that rod undeformed positions are same as object ones)
    EdgeProperty<Scalar> current_theta_copy(m_theta);
    m_theta = m_undef_theta;
    
    //compute rod properties based on undeformed state.
    updateProperties();
    
    //adjust rod forces for possibly modified undeformed state!
    if(m_bending_force)
      m_bending_force->updateAfterRemesh();
    if(m_stretching_force)
      m_stretching_force->updateAfterRemesh();
    if(m_twisting_force)
      m_twisting_force->updateAfterRemesh();

    //reload current configuration
    obj.setVertexPositions(current_position_copy);
    m_theta = current_theta_copy;


    //now recompute the various derived quantities
    updateProperties();
    
    //and recompute masses based on new allocations of radius/volume
    computeMasses();

  }

  void ElasticRodModel::subdivideEdge(EdgeHandle& eh) {
    assert(m_edge_active[eh]);
    
    DeformableObject& obj = getDefoObj();

    VertexHandle from_vh = obj.fromVertex(eh);
    VertexHandle to_vh = obj.toVertex(eh);
    VertexHandle new_vh = obj.addVertex();

    EdgeHandle new_eh0 = obj.addEdge(from_vh, new_vh);
    EdgeHandle new_eh1 = obj.addEdge(new_vh, to_vh);
    
    // For now, treat positions by averaging, and 
    // thetas by copying...
    
    //First vertex variables:
    /////////

    // set up new positions
    Vec3d from_pos = obj.getVertexPosition(from_vh);
    Vec3d to_pos = obj.getVertexPosition(to_vh);
    obj.setVertexPosition(new_vh, 0.5*(from_pos + to_pos));
    
    // undeformed positions
    from_pos = obj.getVertexUndeformedPosition(from_vh);
    to_pos = obj.getVertexUndeformedPosition(to_vh);
    obj.setVertexUndeformedPosition(new_vh, 0.5*(from_pos + to_pos));

    // vertex velocities
    Vec3d from_vel = obj.getVertexVelocity(from_vh);
    Vec3d to_vel = obj.getVertexVelocity(to_vh);
    obj.setVertexVelocity(new_vh, 0.5*(from_vel + to_vel));

    // undeformed damping positions
    from_pos = obj.getVertexDampingUndeformedPosition(from_vh);
    to_pos = obj.getVertexDampingUndeformedPosition(to_vh);
    obj.setVertexDampingUndeformedPosition(new_vh, 0.5*(from_pos + to_pos));

    //Then rotation/edge variables:
    /////////

    Scalar theta = getEdgeTheta(eh);
    Scalar theta_vel = getEdgeThetaVelocity(eh);
    Scalar theta_undef = getEdgeUndeformedTheta(eh);
    Scalar theta_undef_damp = getEdgeDampingUndeformedTheta(eh);

    m_theta[new_eh0] = theta;
    m_theta[new_eh1] = theta;

    m_theta_vel[new_eh0] = theta_vel;
    m_theta_vel[new_eh1] = theta_vel;

    m_undef_theta[new_eh0] = theta_undef;
    m_undef_theta[new_eh1] = theta_undef;

    m_damping_undef_theta[new_eh0] = theta_undef_damp;
    m_damping_undef_theta[new_eh1] = theta_undef_damp;

    m_volumes[new_eh0] = 0.5 * getVolume(eh);
    m_volumes[new_eh1] = 0.5 * getVolume(eh);
    
    m_radii[new_eh0] = getRadii(eh);
    m_radii[new_eh1] = getRadii(eh);
    
    //copy over the reference directors
    m_properties_reference_director1[new_eh0] = m_properties_reference_director1[eh];
    m_properties_reference_director1[new_eh1] = m_properties_reference_director1[eh];

    // eliminate the old edge
    obj.deleteEdge(eh, false);

    // activate the new edges
    m_edge_active[new_eh0] = 1;
    m_edge_active[new_eh1] = 1;
    
  }

  void ElasticRodModel::collapseEdge(EdgeHandle& eh) {
    assert(m_edge_active[eh]);

    DeformableObject& obj = getDefoObj();
    assert(obj.edgeExists(eh));
    
    // Get the relevant vertices
    VertexHandle from_vh = obj.fromVertex(eh);
    VertexHandle to_vh = obj.toVertex(eh);

    // Get both adjacent edges, and pick the shorter one to work with
    EdgeHandle from_eh, to_eh;
    VertexEdgeIterator veit = obj.ve_iter(from_vh);
    while(veit && *veit == eh)
      ++veit;
    from_eh = *veit; 

    veit = obj.ve_iter(to_vh);
    while(veit && *veit == eh)
      ++veit;
    to_eh = *veit;

    VertexHandle centre_vh, nbr_vh0, nbr_vh1; //the stencil at/around the vertex to be removed

    // Check if both nbr edges are valid, and if so, pick the shorter one.
    bool from_valid = obj.edgeExists(from_eh) && obj.vertexIncidentEdges(from_vh) == 2;
    bool to_valid = obj.edgeExists(to_eh) && obj.vertexIncidentEdges(to_vh) == 2;
    
    EdgeHandle partner_eh;
    bool choseLeft;
    if(from_valid && to_valid) { //pick the shorter one
      Scalar from_len = getEdgeLengthFromVertices(from_eh);
      Scalar to_len = getEdgeLengthFromVertices(to_eh);
      if(from_len < to_len) {
        choseLeft = true;
        partner_eh = from_eh;
        centre_vh = from_vh;
        nbr_vh0 = to_vh;
      }
      else {
        choseLeft = false;
        partner_eh = to_eh;
        centre_vh = to_vh;
        nbr_vh0 = from_vh;
      }
    }
    else if(from_valid) {
      choseLeft = true;
      partner_eh = from_eh;
      centre_vh = from_vh;
      nbr_vh0 = to_vh;
    }
    else if(to_valid) {
      choseLeft = false;
      partner_eh = to_eh;
      centre_vh = to_vh;
      nbr_vh0 = from_vh;
    }
   
    // Get the nbr vertex in the partner edge
    EdgeVertexIterator evit = obj.ev_iter(partner_eh);
    while(evit && *evit == centre_vh)
      ++evit;
    nbr_vh1 = *evit;
    assert(nbr_vh1.isValid());
    
    //create the new edge
    EdgeHandle new_edge = !choseLeft ? obj.addEdge(nbr_vh0, nbr_vh1) : obj.addEdge(nbr_vh1, nbr_vh0);
    
    // Simply delete the middle vertex (for now...), and perform averaging for the edge variables
    
    // TODO: Work out what new velocities of adjacent vertices should be by preserving linear momentum
    
    // This momentum must be distributed to the neighbours
    Vec3d centre_mom = getMass(centre_vh) * obj.getVertexVelocity(centre_vh);
    Vec3d old_mom_nbr0 = getMass(nbr_vh0) * obj.getVertexVelocity(nbr_vh0);
    Vec3d old_mom_nbr1 = getMass(nbr_vh1) * obj.getVertexVelocity(nbr_vh1);

    Scalar alpha = 0.5; //TODO: Replace with some better mass/length proportional value

    Scalar new_edge_mass = (m_volumes[eh] + m_volumes[partner_eh]) * m_density;
    Scalar old_edge_mass0 = m_volumes[eh] * m_density;
    Scalar old_edge_mass1 = m_volumes[partner_eh] * m_density;

    // Compute new masses by: old vertex mass + newseg - oldseg;
    Scalar new_mass_nbr0 = getMass(nbr_vh0) + new_edge_mass / 2 - old_edge_mass0 / 2; 
    Scalar new_mass_nbr1 = getMass(nbr_vh1) + new_edge_mass / 2 - old_edge_mass1 / 2;

    // Distribute momentum
    Vec3d new_mom_nbr0 = old_mom_nbr0 + alpha * centre_mom;
    Vec3d new_mom_nbr1 = old_mom_nbr1 + (1-alpha) * centre_mom;

    // Assign appropriate new velocities
    obj.setVertexVelocity(nbr_vh0, new_mom_nbr0 / new_mass_nbr0);
    obj.setVertexVelocity(nbr_vh1, new_mom_nbr1 / new_mass_nbr1);
    
    // Incrementally update the masses while we're at it...
    obj.setVertexMass(nbr_vh0, new_mass_nbr0);
    obj.setVertexMass(nbr_vh1, new_mass_nbr1);

    //
    Scalar theta = 0.5*(m_theta[eh] + m_theta[partner_eh]);
    Scalar theta_vel = 0.5*(m_theta_vel[eh] + m_theta_vel[partner_eh]);
    Scalar theta_undef = 0.5*(m_undef_theta[eh] + m_undef_theta[partner_eh]);
    Scalar theta_undef_damp = 0.5*(m_damping_undef_theta[eh] + m_damping_undef_theta[partner_eh]);
    
    m_theta[new_edge] = theta;
    m_theta_vel[new_edge] = theta_vel;
    m_undef_theta[new_edge] = theta_undef;
    m_damping_undef_theta[new_edge] = theta_undef_damp;

    m_volumes[new_edge] = m_volumes[eh] + m_volumes[partner_eh];

    //Find a new set of radii that preserves the volume...
    Scalar edgeLen = (obj.getVertexPosition(nbr_vh0) - obj.getVertexPosition(nbr_vh1)).norm();
    Scalar avgRad = sqrt(m_volumes[new_edge] / M_PI / edgeLen);
    
    //...and the ratio between radii on each axis.
    Vec2d rad0 = m_radii[eh]; 
    Vec2d rad1 = m_radii[partner_eh];
    Vec2d radHalf = 0.5*(rad0+rad1); 
    Scalar ratio = radHalf[0] / radHalf[1];
    m_radii[new_edge] = Vec2d(avgRad * sqrt(ratio), avgRad/sqrt(ratio));
    //m_radii[new_edge] = Vec2d(avgRad, avgRad);
    

    //what to do with reference directors? average them?
    m_properties_reference_director1[new_edge] = m_properties_reference_director1[eh];
    
    // eliminate the old connectivity data
    bool success = obj.deleteEdge(eh, false);
    assert(success);
    success = obj.deleteEdge(partner_eh, false);
    assert(success);
    success = obj.deleteVertex(centre_vh);
    assert(success);

    // activate the new edges
    m_edge_active[new_edge] = 1;

  }

  void ElasticRodModel::buildStencils( )
  {
    m_edge_stencils.clear();
    m_joint_stencils.clear();
    m_triedge_stencils.clear();
    const DeformableObject& object = getDefoObj();
    
    std::vector<EdgeHandle> rodedges;

    int count = 0;
    for(EdgeIterator eit = object.edges_begin(); eit != object.edges_end(); ++eit) {
      EdgeHandle eh = *eit;
      
      if(!m_edge_active[eh]) 
        continue;
      
      rodedges.push_back(eh);
      
      // each rod edge forms an edge stencil
      EdgeStencil s;
      s.id = m_edge_stencils.size();

      s.e = eh;
      EdgeVertexIterator evit = object.ev_iter(s.e);
      s.v1 = *evit; ++evit;
      s.v2 = *evit; ++evit;
      assert(!evit);

      m_edge_stencils.push_back(s);      
    }

    // joint stencils
    for (VertexIterator i = object.vertices_begin(); i != object.vertices_end(); ++i)
    {
      // detect all the rod edges that are incident to this vertex
      std::vector<EdgeHandle> active_incident_edges;
      for (VertexEdgeIterator veit = object.ve_iter(*i); veit; ++veit)
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
          s.id = m_joint_stencils.size();

          s.v2 = *i;
          s.e1 = active_incident_edges[j];
          s.e2 = active_incident_edges[(j + 1) % active_incident_edges.size()];

          EdgeVertexIterator evit;
          VertexHandle v1, v2;

          evit = object.ev_iter(s.e1);
          v1 = *evit; ++evit;
          v2 = *evit; ++evit;
          assert(!evit);
          s.v1 = (s.v2 == v1 ? v2 : v1);
          assert((s.v1 == v1 && s.v2 == v2) || (s.v1 == v2 && s.v2 == v1));

          evit = object.ev_iter(s.e2);
          v1 = *evit; ++evit;
          v2 = *evit; ++evit;
          assert(!evit);
          s.v3 = (s.v2 == v1 ? v2 : v1);
          assert((s.v3 == v1 && s.v2 == v2) || (s.v3 == v2 && s.v2 == v1));

          // detect if the edges' intrinsic orientation (defined by the underlying mesh data structure) matches the stencil's direction
          s.e1flip = (s.v1 == object.toVertex(s.e1));
          s.e2flip = (s.v3 == object.fromVertex(s.e2));

          if (s.e1flip && s.e2flip)
          {
            // both edge misoriented, we should just flip the stencil itself
            EdgeHandle e = s.e1;
            s.e1 = s.e2;
            s.e2 = e;
            VertexHandle v = s.v1;
            s.v1 = s.v3;
            s.v3 = v;
            s.e1flip = false;
            s.e2flip = false;
          }

          s.curvatureBinormal.setZero();
          s.referenceTwist = 0;
          s.voronoiLength = 0;

          m_joint_stencils.push_back(s);
        }
      }
    }

    for (size_t i = 0; i < rodedges.size(); i++)
    {

      // each rod edge forms an edge stencil
      ThreeEdgeStencil s;
      s.id = m_triedge_stencils.size();

      s.e1 = rodedges[i];

      EdgeVertexIterator evit = object.ev_iter(s.e1);
      s.v1 = *evit; ++evit;
      s.v2 = *evit; ++evit;
      assert(!evit);

      VertexEdgeIterator veit = object.ve_iter(s.v1);
      if(*veit == s.e1) ++veit;
      s.e0 = *veit;

      VertexEdgeIterator veit2 = object.ve_iter(s.v2);
      if(*veit2 == s.e1) ++veit2;
      s.e2 = *veit2;

      if (s.e0.isValid() && m_edge_active[s.e0] && s.e2.isValid() && m_edge_active[s.e2] )
      {
        EdgeVertexIterator evit2 = object.ev_iter(s.e0);
        if(*evit2 == s.v1) ++evit2;
        s.v0 = *evit2;

        EdgeVertexIterator evit3 = object.ev_iter(s.e2);
        if(*evit3 == s.v2) ++evit3;
        s.v3 = *evit3;

        m_triedge_stencils.push_back(s);      
      }
    }
  }

  void ElasticRodModel::assignStencilDofs() {
    DeformableObject & obj = getDefoObj();

    // collect the dof indices info in the stencils. this can't be done in constructor because dof indexing hasn't been computed then
    for (size_t i = 0; i < m_edge_stencils.size(); i++)
    {
      EdgeStencil & s = m_edge_stencils[i];

      s.dofindices.resize(6);
      int dofbase = getDefoObj().getPositionDofBase(s.v1);
      s.dofindices[0] = dofbase;
      s.dofindices[1] = dofbase + 1;
      s.dofindices[2] = dofbase + 2;
      dofbase = getDefoObj().getPositionDofBase(s.v2);
      s.dofindices[3] = dofbase;
      s.dofindices[4] = dofbase + 1;
      s.dofindices[5] = dofbase + 2;
    }

    for (size_t i = 0; i < m_joint_stencils.size(); i++)
    {
      JointStencil & s = m_joint_stencils[i];

      s.dofindices.resize(11);
      int dofbase = getDefoObj().getPositionDofBase(s.v1);
      s.dofindices[0] = dofbase;
      s.dofindices[1] = dofbase + 1;
      s.dofindices[2] = dofbase + 2;
      dofbase = getDefoObj().getPositionDofBase(s.v2);
      s.dofindices[4] = dofbase;
      s.dofindices[5] = dofbase + 1;
      s.dofindices[6] = dofbase + 2;
      dofbase = getDefoObj().getPositionDofBase(s.v3);
      s.dofindices[8] = dofbase;
      s.dofindices[9] = dofbase + 1;
      s.dofindices[10] = dofbase + 2;

      s.dofindices[3] = getEdgeDofBase(s.e1);
      s.dofindices[7] = getEdgeDofBase(s.e2);      
    }

    for (size_t i = 0; i < m_triedge_stencils.size(); i++)
    {
      ThreeEdgeStencil & s = m_triedge_stencils[i];

      s.dofindices.resize(12);
      int dofbase;
      if(s.e0.isValid()) {
        dofbase= getDefoObj().getPositionDofBase(s.v0);
        s.dofindices[0] = dofbase;
        s.dofindices[1] = dofbase + 1;
        s.dofindices[2] = dofbase + 2;
      }

      dofbase = getDefoObj().getPositionDofBase(s.v1);
      s.dofindices[3] = dofbase;
      s.dofindices[4] = dofbase + 1;
      s.dofindices[5] = dofbase + 2;
      dofbase = getDefoObj().getPositionDofBase(s.v2);
      s.dofindices[6] = dofbase;
      s.dofindices[7] = dofbase + 1;
      s.dofindices[8] = dofbase + 2;

      if(s.e2.isValid()) {
        dofbase = getDefoObj().getPositionDofBase(s.v3);
        s.dofindices[9] = dofbase;
        s.dofindices[10] = dofbase + 1;
        s.dofindices[11] = dofbase + 2;
      }
    }
  }

  void ElasticRodModel::setInflow(VertexHandle& vh, Vec3d velocity, Scalar radius) {
     m_inflow_velocity = velocity;
     m_inflow_radius = radius;
     m_inflow_vertex = vh;
     m_inflow_active = true;
     
     //store the position
     Vec3d pos = getDefoObj().getVertexPosition(vh);
     m_inflow_position = pos;

     //assign a constraint velocity
     getDefoObj().constrainVertex(vh, new FixedVelocityConstraint(pos, velocity, 0));
  }

  void ElasticRodModel::extendMesh(Scalar current_time) {
     //create a new vertex
     DeformableObject& defo = getDefoObj();
     VertexHandle newVert = defo.addVertex();
     EdgeHandle newEdge = defo.addEdge(newVert, m_inflow_vertex);
     
     //change the inflow vertex
     m_inflow_vertex = newVert;

     //set vertex position, edge length, etc.
     defo.setVertexPosition(newVert, m_inflow_position);
     defo.setVertexVelocity(newVert, m_inflow_velocity);
     setRadii(newEdge, Vec2d(m_inflow_radius,m_inflow_radius));

     getDefoObj().constrainVertex(newVert, new FixedVelocityConstraint(m_inflow_position, m_inflow_velocity, current_time));
     
     //TODO: Update all properties for the new guys!
     
  }

} //namespace BASim
