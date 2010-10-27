/**
 * \file RodCollisionTimeStepper.hh
 *
 * \author acoull@wetafx.co.nz
 * \date 12/01/2009
 */

#ifndef RODCOLLISIONTIMESTEPPER_HH
#define RODCOLLISIONTIMESTEPPER_HH

#include <tr1/unordered_map>

#include "../../Collisions/CollisionMeshData.hh"
#include "../../Core/ObjectControllerBase.hh"
#include "RodTimeStepper.hh"
#include "RodPenaltyForce.hh"

namespace BASim {

class RodPenaltyForce;
class CollisionMeshData;

typedef std::tr1::unordered_map<size_t, CollisionMeshData *> CollisionMeshDataHashMap;
typedef std::tr1::unordered_map<size_t, CollisionMeshData *>::iterator CollisionMeshDataHashMapIterator;

/** Class to time step a rod, detecting and responding to collisions. */
class RodCollisionTimeStepper : public ObjectControllerBase
{
public:
  RodCollisionTimeStepper(RodTimeStepper* rodTimeStepper, ElasticRod* rod);

  ~RodCollisionTimeStepper();

  ObjectControllerBase* getTimeStepper()
  {
    return m_rodTimeStepper;
  }

  void setVertexPositionPenalty(int vertex_id, Vec3d& target_position, double stiffness, short type = 0 );
  RodVertexConstraint *setVertexPositionPenalty2(int vertex_id, Vec3d& target_position, double stiffness, short type);


  // id vertex_id = -1, delete all
  void clearVertexPositionPenalty(int vertex_id = -1);

  void setTimeStep(Scalar dt)
  {
    if (is_multiple_stepper) {
      m_dt = dt; return;
    } 
  
    m_rodTimeStepper->setTimeStep(dt);
  }
  
  void setBoundaryCondition(RodBoundaryCondition* bc)
  {
    m_rodTimeStepper->setBoundaryCondition(bc);
  }

  RodBoundaryCondition* getBoundaryCondition()
  {
    return m_rodTimeStepper->getBoundaryCondition();
  }
  
  void setCollisionMeshesMap(CollisionMeshDataHashMap* collisionMeshes)
  {
      m_collisionMeshes = collisionMeshes;
  }

  void shouldDoCollisions(bool doCollisionsFlag)
  {
    m_collisionsEnabled = doCollisionsFlag;
  }

  void initialiseCollisionsAndApplyObjectCollisionForces()
  {
    if (!m_enabled)
      return;

    m_rod->setCollisionStartPositions();
    
    if ( !m_collisionsEnabled || m_collisionMeshes == NULL || m_collisionMeshes->size()==0 )
    {
        return;
    }
    
    getProximities(*m_collisionMeshes);    
  }

  void collisionsBegin()
  {
    if (!m_enabled)
      return; 

    if (is_multiple_stepper) {
      m_rod->collisionsBegin(m_dt);
      return;
    }
   
    m_rod->collisionsBegin(m_rodTimeStepper->getTimeStep());
  }

  void updateEndPositions()
  {
    if (!m_enabled)
      return; 

    m_rod->updateEndPositions(m_rodTimeStepper->getTimeStep());    
  }

  void respondToObjectCollisions()
  {
    if (!m_enabled)
      return;

	if ( !m_collisionsEnabled || m_collisionMeshes == NULL || m_collisionMeshes->size()==0 )
    {
        return;
    }

    if (is_multiple_stepper) {
      respondObjectCollisions(*m_collisionMeshes, m_dt);
      return;
    } 

    respondObjectCollisions(*m_collisionMeshes, m_rodTimeStepper->getTimeStep());
  }

  void tidyUpCollisionStructuresForNextStep()
  {
     if (!m_enabled)
      return;

     if (is_multiple_stepper) {
       m_rod->collisionsEnd(m_dt);
       return;
     } 

     m_rod->collisionsEnd(m_rodTimeStepper->getTimeStep());

  }

  bool execute();
  
  ElasticRod* getRod() { return m_rod; }


  static void getProximities(vector<ElasticRod*> &rods);
  static void respondRodCollisions(vector<ElasticRod*> &rods, Scalar dt, int maxIterations, Scalar COR);
  
  static void getClumpingPairs(vector<ElasticRod*> &rods);

  void setClumping(bool flag, Scalar coeff = 0.0);

  bool impulse_enabled;

    // the following two functions are temp, just to test how to
  // change the constraints in one sim
//  void addVertexPositionContstraints( VertexConstraintMapIter &i_iter )
//  {
//      m_vertexContraints.push_back( i_iter );
//  }

  void updateVertexPositionConstraints();

protected:
  void getProximities(CollisionMeshDataHashMap &collisionMeshes);
  void respondObjectCollisions(CollisionMeshDataHashMap &collisionMeshes, Real dt);

  bool m_collisionsEnabled;
  CollisionMeshDataHashMap* m_collisionMeshes;
  RodPenaltyForce* m_rodPenaltyForce;
  RodTimeStepper* m_rodTimeStepper;
  ElasticRod* m_rod;

  bool is_multiple_stepper; // no penalty method, no m_rodTimeStepper
  double m_dt; // only used in multiple stepping

  std::vector< VertexConstraintMapIter > m_vertexContraints;

};


} // namespace BASim

#endif // RODCOLLISIONTIMESTEPPER_HH
