
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"

namespace BASim {

DeformableObject::DeformableObject() :
  m_dt(1), m_models(0) 
{

}

DeformableObject::~DeformableObject()
{

}


void DeformableObject::setTimeStep(Scalar dt) {
  m_dt = dt;
}

Scalar DeformableObject::getTimeStep() {
  return m_dt;
}
void DeformableObject::setTime(Scalar time) {
  m_time = time;
}

Scalar DeformableObject::getTime() {
  return m_dt;
}


void DeformableObject::computeForces(VecXd& force) {
  
  std::vector<PhysicalModel*>::iterator model_it;
  VecXd curr_force(force.size());
  for(model_it = m_models.begin(); model_it != m_models.end(); ++model_it) {
    curr_force.setZero();
    (*model_it)->computeForces(curr_force);
    force += curr_force;
  }

}


void DeformableObject::computeJacobian(Scalar scale, MatrixBase& J) {

  std::vector<PhysicalModel*>::iterator model_it;
  for(model_it = m_models.begin(); model_it != m_models.end(); ++model_it) {
    (*model_it)->computeJacobian(scale, J);
  }

}

void DeformableObject::addModel( PhysicalModel* model )
{
  m_models.push_back(model);
}

const Scalar& DeformableObject::getDof(int i) const {
  int model = m_dofModels[i];
  DofHandle hnd = m_dofHandles[i];
  return m_models[model]->getDof(hnd);
}

void DeformableObject::setDof(int i, const Scalar& dof) {
  int model = m_dofModels[i];
  DofHandle hnd = m_dofHandles[i];
  m_models[model]->setDof(hnd, dof);
}

const Scalar& DeformableObject::getVel(int i) const {
  int model = m_dofModels[i];
  DofHandle hnd = m_dofHandles[i];
  return m_models[model]->getVel(hnd);
}

void DeformableObject::setVel(int i, const Scalar& vel) {
  int model = m_dofModels[i];
  DofHandle hnd = m_dofHandles[i];
  m_models[model]->setVel(hnd, vel);
}

const Scalar& DeformableObject::getMass(int i) const {
  int model = m_dofModels[i];
  DofHandle hnd = m_dofHandles[i];
  return m_models[model]->getMass(hnd);
}

void DeformableObject::computeDofIndexing()
{
  //Number each degree of freedom that exists in the active models

  //Central DOF numbering scheme:
  //for each simplex type
    //for each such simplex 
      //for each model (that is active for this simplex)
        //for each DOF the model requests
          //index++;

  m_dofHandles.clear();
  m_dofModels.clear();

  //Vertex DOF's
  int dofIndex = 0;
  VertexIterator vert_it;
  for(vert_it = vertices_begin(); vert_it != vertices_end(); ++vert_it) {
    for(unsigned int m = 0; m < m_models.size(); ++m) {
      if(m_models[m]->isVertexActive(*vert_it)) {
        m_models[m]->setVertexDofBase(*vert_it, dofIndex);
        for(int d = 0; d < m_models[m]->numVertexDofs(); ++d) {
          DofHandle h(dofIndex);
          h.setNum(d);
          h.setType(DofHandle::VERTEX_DOF);
          h.setHandle(*vert_it);
          m_dofHandles.push_back(h);
          m_dofModels.push_back(m);
          dofIndex++;
        }
      }
    }
  }

  //Edge DOF's
  EdgeIterator edge_it;
  for(edge_it = edges_begin(); edge_it != edges_end(); ++edge_it) {
    for(unsigned int m = 0; m < m_models.size(); ++m) {
      if(m_models[m]->isEdgeActive(*edge_it)) {
        m_models[m]->setEdgeDofBase(*edge_it, dofIndex);
        for(int d = 0; d < m_models[m]->numEdgeDofs(); ++d) {
          DofHandle h(dofIndex);
          h.setNum(d);
          h.setType(DofHandle::EDGE_DOF);
          h.setHandle(*edge_it);
          m_dofHandles.push_back(h);
          m_dofModels.push_back(m);
          dofIndex++;
        }
      }
    }
  }

  //Face DOF's
  FaceIterator face_it;
  for(face_it = faces_begin(); face_it != faces_end(); ++face_it) {
    for(unsigned int m = 0; m < m_models.size(); ++m) {
      if(m_models[m]->isFaceActive(*face_it)) {
        m_models[m]->setFaceDofBase(*face_it, dofIndex);
        for(int d = 0; d < m_models[m]->numFaceDofs(); ++d) {
          DofHandle h(dofIndex);
          h.setNum(d);
          h.setType(DofHandle::FACE_DOF);
          h.setHandle(*face_it);
          m_dofHandles.push_back(h);
          m_dofModels.push_back(m);
          dofIndex++;
        }
      }
    }
  }

  //Tet DOF's
  TetIterator tet_it;
  for(tet_it = tets_begin(); tet_it != tets_end(); ++tet_it) {
    for(unsigned int m = 0; m < m_models.size(); ++m) {
      if(m_models[m]->isTetActive(*tet_it)) {
        m_models[m]->setTetDofBase(*tet_it, dofIndex);
        for(int d = 0; d < m_models[m]->numTetDofs(); ++d) {
          DofHandle h(dofIndex);
          h.setNum(d);
          h.setType(DofHandle::TET_DOF);
          h.setHandle(*tet_it);
          m_dofHandles.push_back(h);
          m_dofModels.push_back(m);
          dofIndex++;
        }
      }
    }
  }
  m_ndof = dofIndex;
  
}

void DeformableObject::getScriptedDofs( IntArray& dofIndices, std::vector<Scalar>& dofValues, Scalar time ) const
{
  for(unsigned int i = 0; i < m_models.size(); ++i)
    m_models[i]->getScriptedDofs(dofIndices, dofValues, time);
}



} //namespace BASim