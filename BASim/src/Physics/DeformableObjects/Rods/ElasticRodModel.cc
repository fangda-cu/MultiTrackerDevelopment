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

      // each rod edge forms an edge stencil
      EdgeStencil s;
      s.id = m_edge_stencils.size();
      
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
          s.id = m_joint_stencils.size();
          
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
          
          // detect if the edges' intrinsic orientation (defined by the underlying mesh data structure) matches the stencil's direction
          s.e1flip = (s.v1 == object->toVertex(s.e1));
          s.e2flip = (s.v3 == object->fromVertex(s.e2));
          
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

      EdgeVertexIterator evit = object->ev_iter(s.e1);
      s.v1 = *evit; ++evit;
      s.v2 = *evit; ++evit;
      assert(!evit);

      VertexEdgeIterator veit = object->ve_iter(s.v1);
      if(*veit == s.e1) ++veit;
      s.e0 = *veit;

      VertexEdgeIterator veit2 = object->ve_iter(s.v2);
      if(*veit2 == s.e1) ++veit2;
      s.e2 = *veit2;

      EdgeVertexIterator evit2 = object->ev_iter(s.e0);
      if(*evit2 == s.v1) ++evit2;
      s.v0 = *evit2;

      EdgeVertexIterator evit3 = object->ev_iter(s.e2);
      if(*evit3 == s.v2) ++evit3;
      s.v3 = *evit3;

      m_triedge_stencils.push_back(s);      
    }

    
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

    // collect the dof indices info in the stencils. this can't be done in constructor because dof indexing hasn't been computed then
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
    
    // compute edges
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
      if (m_stretching_force)
        m_stretching_force->stencils()[s.id].copyData(s);
    }
    
    for (size_t i = 0; i < m_joint_stencils.size(); i++)
    {
      JointStencil & s = m_joint_stencils[i];
      if (m_bending_force)
        m_bending_force->stencils()[s.id].copyData(s);
      if (m_twisting_force)
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

} //namespace BASim
