/**
 * \file RodCollisionTimeStepper.hh
 *
 * \author acoull@wetafx.co.nz
 * \date 12/01/2009
 */

#ifndef RODCOLLISIONTIMESTEPPER_HH
#define RODCOLLISIONTIMESTEPPER_HH

#include <tr1/unordered_map>

#include <BASim/src/Collisions/CollisionMeshData.hh>
//#include "RodTimeStepper.hh"

namespace BASim {

class RodPenaltyForce;
class CollisionMeshData;

typedef std::tr1::unordered_map<size_t, CollisionMeshData *> CollisionMeshDataHashMap;
typedef std::tr1::unordered_map<size_t, CollisionMeshData *>::iterator CollisionMeshDataHashMapIterator;

/** Class to time step a rod, detecting and responding to collisions. */
class RodCollisionTimeStepper : public ObjectControllerBase
{
public:
  RodCollisionTimeStepper(ObjectControllerBase* rodTimeStepper, ElasticRod* rod);

  ~RodCollisionTimeStepper();

  ObjectControllerBase* getTimeStepper()
  {
    return m_rodTimeStepper;
  }
  
  void execute()
  {
    m_rodTimeStepper->execute();
  }

  void execute(CollisionMeshDataHashMap &collisionMeshes, Scalar dt)
  {
    m_rod->setCollisionStartPositions();
    getProximities(collisionMeshes);
    m_rodTimeStepper->execute();
    m_rod->collisionsBegin(dt);
    
 //   if (m_fullSelfCollisionsEnabled)
   //   Rod::respondRodCollisions(_rods, collisionDT, m_selfCollisionsIterations,
   //                             m_collisionsCoefficientOfRestitution);
    respondObjectCollisions(collisionMeshes, dt);
    m_rod->collisionsEnd(dt);
  }

  void getProximities(CollisionMeshDataHashMap &collisionMeshes);
  void respondObjectCollisions(CollisionMeshDataHashMap &collisionMeshes, Real dt);


protected:
  RodPenaltyForce* m_rodPenaltyForce;
  ObjectControllerBase* m_rodTimeStepper;
  ElasticRod* m_rod;
};


} // namespace BASim

#endif // RODCOLLISIONTIMESTEPPER_HH
