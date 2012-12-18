//
//  PositionDofsModel.cc
//  BASim
//
//  Created by Fang Da (fang@cs.columbia.edu) on 5/11/12.
//  Copyright (c) 2012 Columbia. All rights reserved.
//

#include "BASim/src/Physics/DeformableObjects/PositionDofsModel.hh"
#include "BASim/src/Physics/DeformableObjects/DefoObjForce.hh"

namespace BASim
{

  PositionDofsModel::~PositionDofsModel() {
    for(unsigned int i = 0; i < m_position_forces.size(); ++i) {
      delete m_position_forces[i];
    }
  }

  void PositionDofsModel::startStep(Scalar time, Scalar timestep)
  {
    m_damping_undeformed_positions = m_positions;
  }

  void PositionDofsModel::endStep(Scalar time, Scalar timestep)
  {
    
  }
  
//  void PositionDofsModel::accumulateMasses(const VertexProperty<Scalar> &masses)
//  { 
//    for (VertexIterator i = getDefoObj().vertices_begin(); i != getDefoObj().vertices_end(); ++i) 
//      m_vertex_masses[*i] += masses[*i]; 
//  }

  void PositionDofsModel::constrainVertex(const VertexHandle & v, const Vec3d& pos)
  {
    m_constrained_vertices.push_back(v);
    PositionConstraint* c = new FixedPositionConstraint(pos);
    m_constraint_positions.push_back(c);
  }
  
  void PositionDofsModel::constrainVertex(const VertexHandle & v, PositionConstraint * c)
  {
    m_constrained_vertices.push_back(v);
    m_constraint_positions.push_back(c);
  }
  
  void PositionDofsModel::releaseVertex(const VertexHandle & v)
  {    
    bool deletedVertex = true;
    while(deletedVertex) {
      deletedVertex = false;
      //can only have one constraint or things get broken anyways(right?), so no need to search for multiple
      int index = -1;
      for(unsigned int i = 0; i < m_constraint_positions.size(); ++i) {
        if(m_constrained_vertices[i] == v) {
          index = i;
          deletedVertex = true;
          break;
        }
      }
      
      //remove the constraint
      if(index != -1) {
        delete m_constraint_positions[index];
        m_constraint_positions.erase(m_constraint_positions.begin()+index);
        m_constrained_vertices.erase(m_constrained_vertices.begin()+index);
      }
    }
  }

  void PositionDofsModel::releaseAllVertices()
  {
    assert(m_constrained_vertices.size() == m_constraint_positions.size());
    for (size_t i = 0; i < m_constraint_positions.size(); i++)
      delete m_constraint_positions[i];
    m_constrained_vertices.clear();
    m_constraint_positions.clear();
  }
  
  bool PositionDofsModel::isConstrained(const VertexHandle & v) const 
  {
    for(unsigned int i = 0; i < m_constrained_vertices.size(); ++i)
      if(m_constrained_vertices[i] == v)
        return true;
    return false;
  }

  void PositionDofsModel::getScriptedDofs(IntArray & dofIndices, std::vector<Scalar> & dofValues, Scalar time) const
  {
    for(unsigned int i = 0; i < m_constrained_vertices.size(); ++i) 
    {  
      int dofBase = getVertexDofBase(m_constrained_vertices[i]);
      Vec3d pos = m_constraint_positions[i]->operator()(time);
      if(m_constraint_positions[i]->xEnabled) {
        dofIndices.push_back(dofBase); 
        dofValues.push_back(pos[0]);
      }
      if(m_constraint_positions[i]->yEnabled) {
        dofIndices.push_back(dofBase+1); 
        dofValues.push_back(pos[1]);
      }
      if(m_constraint_positions[i]->zEnabled) {
        dofIndices.push_back(dofBase+2); 
        dofValues.push_back(pos[2]);
      }
    }
  }

  void PositionDofsModel::computeForces(VecXd& force) {
    const std::vector<DefoObjForce*>& forces = m_position_forces;
    std::vector<DefoObjForce*>::const_iterator fIt;

    VecXd curr_force(force.size());
    for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
      curr_force.setZero();
      (*fIt)->globalForce(curr_force);

      force += curr_force;
    }
  }

  void PositionDofsModel::computeJacobian(Scalar scale, MatrixBase& J) {
    
    const std::vector<DefoObjForce*>& forces = m_position_forces;
    std::vector<DefoObjForce*>::const_iterator fIt;

    for (fIt = forces.begin(); fIt != forces.end(); ++fIt)
      (*fIt)->globalJacobian(scale, J);
    
  }
}
