
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Physics/DeformableObjects/PositionDofsModel.hh"

namespace BASim {

DeformableObject::DeformableObject() :
  m_dt(1), m_models(0), m_posdofsmodel(NULL)
{
  m_posdofsmodel = new PositionDofsModel(this);
  addModel(m_posdofsmodel);
}

DeformableObject::~DeformableObject()
{
  if (m_posdofsmodel) delete m_posdofsmodel;
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

// All DOFS at once
const VertexProperty<Vec3d>& DeformableObject::getVertexPositions() const                   { return m_posdofsmodel->getPositions(); }
const VertexProperty<Vec3d>& DeformableObject::getVertexVelocities() const                  { return m_posdofsmodel->getVelocities(); }
const VertexProperty<Vec3d>& DeformableObject::getVertexUndeformedPositions() const         { return m_posdofsmodel->getUndeformedPositions(); }
const VertexProperty<Vec3d>& DeformableObject::getVertexDampingUndeformedPositions() const  { return m_posdofsmodel->getDampingUndeformedPositions(); }

void DeformableObject::setVertexPositions                 (const VertexProperty<Vec3d>& pos) { m_posdofsmodel->setPositions(pos); }
void DeformableObject::setVertexVelocities                (const VertexProperty<Vec3d>& vel) { m_posdofsmodel->setVelocities(vel); }
void DeformableObject::setVertexUndeformedPositions       (const VertexProperty<Vec3d>& pos) { m_posdofsmodel->setUndeformedPositions(pos); }
void DeformableObject::setVertexDampingUndeformedPositions(const VertexProperty<Vec3d>& pos) { m_posdofsmodel->setDampingUndeformedPositions(pos); }

//Individual DOFs
Vec3d DeformableObject::getVertexPosition                 (const VertexHandle& v) const { return m_posdofsmodel->getPosition(v); }
Vec3d DeformableObject::getVertexVelocity                 (const VertexHandle& v) const { return m_posdofsmodel->getVelocity(v); }
Scalar DeformableObject::getVertexMass                    (const VertexHandle& v) const { return m_posdofsmodel->getMass(v); }
Vec3d DeformableObject::getVertexUndeformedPosition       (const VertexHandle& v) const { return m_posdofsmodel->getUndeformedPosition(v); }
Vec3d DeformableObject::getVertexDampingUndeformedPosition(const VertexHandle& v) const { return m_posdofsmodel->getDampingUndeformedPosition(v); }

void DeformableObject::setVertexPosition                  (const VertexHandle& v, const Vec3d& pos) { m_posdofsmodel->setPosition(v, pos); }
void DeformableObject::setVertexVelocity                  (const VertexHandle& v, const Vec3d& vel) { m_posdofsmodel->setVelocity(v, vel); }
void DeformableObject::setVertexMass                      (const VertexHandle& v, Scalar m)         { m_posdofsmodel->setMass(v, m); }
void DeformableObject::setVertexUndeformedPosition        (const VertexHandle& v, const Vec3d& pos) { m_posdofsmodel->setUndeformedPosition(v, pos); }
void DeformableObject::setVertexDampingUndeformedPosition (const VertexHandle& v, const Vec3d& pos) { m_posdofsmodel->setDampingUndeformedPosition(v, pos); }
  
void DeformableObject::clearVertexMasses() { m_posdofsmodel->clearMasses(); }
void DeformableObject::accumulateVertexMasses(const VertexProperty<Scalar>& masses) { m_posdofsmodel->accumulateMasses(masses); }
void DeformableObject::accumulateVertexMass(const VertexHandle& v, Scalar mass) { m_posdofsmodel->accumulateMass(v, mass); }

void DeformableObject::updateVertexMasses()
{
  clearVertexMasses();
  for (size_t i = 0; i < m_models.size(); i++)
    accumulateVertexMasses(m_models[i]->getVertexMasses());
}
  
int DeformableObject::getPositionDofBase(const VertexHandle& vh) const { return m_posdofsmodel->getVertexDofBase(vh); }


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

void DeformableObject::startStep() 
{ 
  for(unsigned int i = 0; i < m_models.size(); ++i) 
    m_models[i]->startStep(m_time, m_dt); 
}

void DeformableObject::endStep() 
{ 
  for(unsigned int i = 0; i < m_models.size(); ++i) 
    m_models[i]->endStep(m_time, m_dt); 
}

void DeformableObject::startIteration() 
{ 
  for(unsigned int i = 0; i < m_models.size(); ++i) 
    m_models[i]->startIteration(m_time, m_dt); 
}

void DeformableObject::endIteration() 
{ 
  for(unsigned int i = 0; i < m_models.size(); ++i) 
    m_models[i]->endIteration(m_time, m_dt); 
}


} //namespace BASim