#include "BASim/src/Collisions/ElTopo/broadphasegrid.hh"
#include "BASim/src/Physics/DeformableObjects/Solids/ElasticSolid.hh"
#include "BASim/src/Physics/DeformableObjects/Solids/ElasticSolidForce.hh"
#include "BASim/src/Physics/DeformableObjects/DeformableObject.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"
#include "BASim/src/Math/Math.hh"
#include "BASim/src/Core/TopologicalObject/TopObjUtil.hh"

#include <algorithm>
#include <numeric>

namespace BASim {

ElasticSolid::ElasticSolid(DeformableObject* object, const TetProperty<char>& solidTets, Scalar timestep) : 
  PhysicalModel(*object), m_obj(object), 
    m_active_tets(solidTets), 
    m_vertex_masses(object),
    m_density(1)
{
}

ElasticSolid::~ElasticSolid() {
}

void ElasticSolid::computeForces( VecXd& force )
{
  const std::vector<ElasticSolidForce*>& forces = getForces();
  std::vector<ElasticSolidForce*>::const_iterator fIt;
  
  VecXd curr_force(force.size());
  for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
    curr_force.setZero();
    (*fIt)->globalForce(curr_force);
    
    force += curr_force;
  }

}

void ElasticSolid::computeJacobian( Scalar scale, MatrixBase& J )
{
  const std::vector<ElasticSolidForce*>& forces = getForces();
  std::vector<ElasticSolidForce*>::const_iterator fIt;

  for (fIt = forces.begin(); fIt != forces.end(); ++fIt)
    (*fIt)->globalJacobian(scale, J);
}

const std::vector<ElasticSolidForce*>& ElasticSolid::getForces() const
{
  return m_solid_forces;
}

void ElasticSolid::addForce( ElasticSolidForce* force )
{
  assert(force != NULL);

  m_solid_forces.push_back(force);
}


void ElasticSolid::setDensity(Scalar density) {
  m_density = density;
}

void ElasticSolid::computeMasses()
{
  //Compute vertex masses in a lumped mass way.

  m_vertex_masses.assign(0);
  
  Scalar area = 0;

  //Iterate over all tets active in this solid and accumulate vertex masses
  for(TetIterator t_iter = m_obj->tets_begin(); t_iter != m_obj->tets_end(); ++t_iter) {
    TetHandle& t_hnd = *t_iter;
    if(m_active_tets[t_hnd]) {

      //get the four vertices
      TetVertexIterator tvit = m_obj->tv_iter(t_hnd);
      VertexHandle v0_hnd = *tvit; ++tvit; assert(tvit);
      VertexHandle v1_hnd = *tvit; ++tvit; assert(tvit);
      VertexHandle v2_hnd = *tvit; ++tvit; assert(tvit);
      VertexHandle v3_hnd = *tvit; ++tvit; assert(!tvit);

      //compute tets volumes
      Vec3d v0 = getVertexPosition(v1_hnd) - getVertexPosition(v0_hnd);
      Vec3d v1 = getVertexPosition(v2_hnd) - getVertexPosition(v0_hnd);
      Vec3d v2 = getVertexPosition(v3_hnd) - getVertexPosition(v0_hnd);
      
      //Volume per the formula V = | (a-d).((b-d)x(c-d)) | / 6

      Scalar volVec = fabs(v0.dot(v1.cross(v2))) / 6;
      Scalar contribution = volVec * m_density;
      
      //accumulate mass to the vertices
      m_vertex_masses[v0_hnd] += contribution;
      m_vertex_masses[v1_hnd] += contribution;
      m_vertex_masses[v2_hnd] += contribution;
      m_vertex_masses[v3_hnd] += contribution;

    }
  }
 
  m_obj->updateVertexMasses();
}

bool ElasticSolid::isVertexActive( const VertexHandle& v ) const
{
  //determine if the vertex is on any active face
  VertexFaceIterator vf = m_obj->vf_iter(v);
  for(;vf; ++vf) {
    if(isFaceActive(*vf)) {
      return true;
    }
  }
  
  return false;
}

bool ElasticSolid::isEdgeActive( const EdgeHandle& e) const {
  //if any adjacent face is active, we say this edge is active.
  EdgeFaceIterator ef = m_obj->ef_iter(e);
  for(;ef;++ef) {
    if(isFaceActive(*ef)) {
      return true;
    }
  }
  
  return false;
}

bool ElasticSolid::isFaceActive( const FaceHandle& f) const {
  //if any adjacent face is active, we say this edge is active.
  FaceTetIterator ft = m_obj->ft_iter(f);
  for(;ft;++ft) {
    if(isTetActive(*ft)) {
      return true;
    }
  }

  return false;
}

const Scalar& ElasticSolid::getDof( const DofHandle& hnd ) const
{
  //no DoFs, all handled by positionDofsModel
  assert(false);
  static Scalar dummy=0;
  return dummy;
}

void ElasticSolid::setDof( const DofHandle& hnd, const Scalar& dof )
{
  //no DoFs, all handled by positionDofsModel
 
}

const Scalar& ElasticSolid::getVel( const DofHandle& hnd ) const
{
  //no DoFs, all handled by positionDofsModel
  static Scalar dummy = 0;
  return dummy;
}

void ElasticSolid::setVel( const DofHandle& hnd, const Scalar& vel )
{
  //no DoFs, all handled by positionDofsModel
}

const Scalar& ElasticSolid::getMass( const DofHandle& hnd ) const
{
  //no DoFs, all handled by positionDofsModel
  static Scalar dummy = 0;
  return dummy;
}

void ElasticSolid::getScriptedDofs( IntArray& dofIndices, std::vector<Scalar>& dofValues, Scalar time ) const
{
  // position dof scripting is moved to PositionDofsModel.
}

void ElasticSolid::startStep(Scalar time, Scalar timestep)
{
  std::cout << "Starting startStep\n";

  //tell the forces to update anything they need to update
  const std::vector<ElasticSolidForce*>& forces = getForces();
  for(unsigned int i = 0; i < forces.size(); ++i) {
    forces[i]->update();
  }
  std::cout << "Done startStep\n";
}

void ElasticSolid::endStep(Scalar time, Scalar timestep) {

  std::cout << "Starting endStep.\n";
  bool do_relabel = false;

  std::cout << "Completed endStep\n";

}


} //namespace BASim
