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
class RodTimeStepper;

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

  void setTimeStep(Scalar dt)
  {
    m_dt = dt;   
  }
  
  void setCollisionMeshesMap(CollisionMeshDataHashMap* collisionMeshes)
  {
      m_collisionMeshes = collisionMeshes;
  }

  void doCollisions(bool doCollisionsFlag)
  {
    m_collisionsEnabled = doCollisionsFlag;
  }

  void execute()
  {
    if (!m_enabled)
      return;

    // Sanity check, m_dt MUST EQUAL m_rodTimeStepper->getDt();
    // get it from m_rodTimeStepper rather than using the one here.
   // if ( !m_collisionsEnabled || m_collisionMeshes == NULL || m_collisionMeshes->size()==0 )
    {
        m_rodTimeStepper->execute();
    }
    /*
    m_rod->setCollisionStartPositions();
    getProximities(*m_collisionMeshes);
    m_rodTimeStepper->execute();
    m_rod->collisionsBegin(m_dt);
    
    //   if (m_fullSelfCollisionsEnabled)
    //   Rod::respondRodCollisions(_rods, collisionDT, m_selfCollisionsIterations,
    //                             m_collisionsCoefficientOfRestitution);
    respondObjectCollisions(*m_collisionMeshes, m_dt);
    m_rod->collisionsEnd(m_dt);*/
  }

protected:
  void getProximities(CollisionMeshDataHashMap &collisionMeshes);
  void respondObjectCollisions(CollisionMeshDataHashMap &collisionMeshes, Real dt);

  bool m_collisionsEnabled;
  CollisionMeshDataHashMap* m_collisionMeshes;
  RodPenaltyForce* m_rodPenaltyForce;
  ObjectControllerBase* m_rodTimeStepper;
  ElasticRod* m_rod;
  Scalar m_dt;

  // Tracks which vertices are fixed for this rod.
 // vector<bool>& m_vertexFixed;
};


} // namespace BASim

#endif // RODCOLLISIONTIMESTEPPER_HH
