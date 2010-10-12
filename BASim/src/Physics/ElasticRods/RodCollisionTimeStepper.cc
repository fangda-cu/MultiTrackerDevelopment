
#include "RodCollisionTimeStepper.hh"
#include "RodTimeStepper.hh"
#include "RodPenaltyForce.hh"

namespace BASim {

RodCollisionTimeStepper::RodCollisionTimeStepper(RodTimeStepper* rodTimeStepper, ElasticRod* rod)
  : m_collisionsEnabled(false), m_rodPenaltyForce(NULL), m_rodTimeStepper(rodTimeStepper), m_rod(rod), is_multiple_stepper(false)
{
  if (m_rodTimeStepper) {
    m_rodPenaltyForce = new RodPenaltyForce();
    dynamic_cast<RodTimeStepper*>(m_rodTimeStepper)->addExternalForce(m_rodPenaltyForce);  
    m_rod->setPenaltyForce(m_rodPenaltyForce);
  } else {
    // multiple rod time stepper is being used
    
    is_multiple_stepper = true;
    m_rodPenaltyForce = NULL;
    
  }
  
  m_collisionMeshes = NULL;
  impulse_enabled = true;
  
}

//void RodCollisionTimeStepper::updateVertexPositionConstraints()
//{
//    //this is dummy test, meaningless at all
//    for( size_t s = 0; s < m_vertexContraints.size(); s++ )
//    {
//        Vec3d direction = m_vertexContraints[ s ]->second.m_target;
//        double norm = direction.norm();
//        if( norm > 0 )
//            direction = direction * ( 1. / norm );
//
//        //make the target position step forward a little bit...
//        m_vertexContraints[ s ]->second.m_target =
//                m_vertexContraints[ s ]->second.m_target + 0.1 * direction;
//    }
//
//}

void RodCollisionTimeStepper::setVertexPositionPenalty(int vertex_id, Vec3d& target_position, double stiffness, short type)
{
  if (m_rodPenaltyForce) {
    if (vertex_id < m_rod->nv() && vertex_id >= 0)
    {
      m_vertexContraints.push_back(m_rodPenaltyForce->setVertexPositionPenalty( m_rod,
                vertex_id, target_position, stiffness, type));
               
    } else {
      std::cout << "invalid vertex id\n";
    }
  } else {
    std::cout << "add penatly force first\n";
  }
}

RodVertexConstraint *RodCollisionTimeStepper::setVertexPositionPenalty2(int vertex_id, Vec3d& target_position, double stiffness, short type)
{
  if (m_rodPenaltyForce) {
    if (vertex_id < m_rod->nv() && vertex_id >= 0)
    {
      m_vertexContraints.push_back(m_rodPenaltyForce->setVertexPositionPenalty( m_rod,
                vertex_id, target_position, stiffness, type));
      return &(m_vertexContraints.back()->second);

    } else {
      std::cout << "invalid vertex id\n";
    }
  } else {
    std::cout << "add penatly force first\n";
  }
}

// id vertex_id = -1, delete all
void RodCollisionTimeStepper::clearVertexPositionPenalty(int vertex_id) 
{
  if (m_rodPenaltyForce) {
    if (vertex_id == -1 || (vertex_id< m_rod->nv() && vertex_id >= 0)) {
      m_rodPenaltyForce->clearVertexPositionPenalty(vertex_id);
    } else {
      std::cout << "invalid vertex id\n";
    }
  } else {
    std::cout << "add penatly force first\n";
  }
}

  
void RodCollisionTimeStepper::setClumping(bool flag, Scalar coeff) {
	if (m_rodPenaltyForce) {
		m_rodPenaltyForce->setClumping(flag, coeff);
	}
}

bool RodCollisionTimeStepper::execute()
  {
    if (!m_enabled)
      return false;
      /*
    if (m_rod) {
			std::cout << "\n\n\n\n\n\nBEFORE dynamic execute vertex \n" << m_rod->quasistatic() << "\n";
		  for(int i=0; i < (int)m_rod->nv(); i++) {
			  std::cout << m_rod->getVertex(i) << " ";
			  if (i < (int)m_rod->nv() - 1) std::cout << m_rod->getTheta(i) << "       ";
			}
			std::cout << "\n";
			std::cout << "execute vel\n";
		  for(int i=0; i < (int)m_rod->nv(); i++) {
			  std::cout << m_rod->getVelocity(i) << " ";
			  if (i < (int)m_rod->nv() - 1) std::cout << m_rod->getThetaDot(i) << "       ";
			}
		  for(int i=0; i < (int)m_rod->nv() * 4 - 1; i++) {
			  std::cout << m_rod->getMass(i) << " ";
			}
			std::cout << "\n";
			std::cout << "\n\n";
    } else {
			std::cout << "m_rod null\n\n";
    
    }
      */

  //std::cout << "COLLISION execute\n";
    return m_rodTimeStepper->execute();
/*
		std::cout << "after dynamic execute vertex\n";
	  for(int i=0; i < (int)m_rod->nv(); i++) {
		  std::cout << m_rod->getVertex(i) << " ";
		  if (i < (int)m_rod->nv() - 1) std::cout << m_rod->getTheta(i) << "       ";
		}
		std::cout << "\n";
		std::cout << "execute vel\n";
	  for(int i=0; i < (int)m_rod->nv(); i++) {
		  std::cout << m_rod->getVelocity(i) << " ";
		  if (i < (int)m_rod->nv() - 1) std::cout << m_rod->getThetaDot(i) << "       ";
		}
		std::cout << "\n";
		*/

  }


RodCollisionTimeStepper::~RodCollisionTimeStepper()
{
  if (m_rodTimeStepper != NULL)
  {
    delete m_rodTimeStepper;
  }

  // Don't delete penalty forces here any more, they are removed in RodTimeStepper  
}

/* static */ void RodCollisionTimeStepper::getProximities(vector<ElasticRod*> &rods)
{
    // Get list of edge-edge proximities
    //
    UniformGrid grid;
    Collisions collisions;
    grid.getProximities(rods, collisions);

    for (CollisionsIterator cItr=collisions.begin(); cItr!=collisions.end(); ++cItr)
    {
        // Add penalty forces for each proximity
        // It must be added to BOTH rods, since a force computation for a rod
        // can only be for that particular rod, a little wasteful but not many options
        // without changing the framework
        //
        ElasticRod *rod1 = dynamic_cast<ElasticRod *>(cItr->getFirstObject());
        ElasticRod *rod2 = dynamic_cast<ElasticRod *>(cItr->getSecondObject());
        rod1->getPenaltyForce()->addRodPenaltyForce(cItr->getFirstPrimitiveIndex(0),
                                          rod2, cItr->getSecondPrimitiveIndex(0));
        rod2->getPenaltyForce()->addRodPenaltyForce(cItr->getSecondPrimitiveIndex(0),
                                          rod1, cItr->getFirstPrimitiveIndex(0));
    }
}


void RodCollisionTimeStepper::getClumpingPairs(vector<ElasticRod*> &rods)
{

  const Scalar max_clump_size = 5.0;
	
  for (int i=0; i<(int) rods.size(); i++)   
  {
    ElasticRod *rod1 = rods[i];
    Vec3d x1 = rod1->getVertex(0);
  	
	for (int j=i+1; j<(int) rods.size(); j++) 		
    {
	  ElasticRod *rod2 = rods[j];
	  Vec3d x2 = rod2->getVertex(0);

	  Scalar dist = (x1 - x2).norm();
			
	  if (dist < max_clump_size) 
      {
	    for (int v1=0; v1<rod1->nv(); v1++) 
        {
		  if (!rod1->vertFixed(v1)) 
          {
		    for (int v2=0; v2<rod2->nv(); v2++) 
            {
			  if (!rod2->vertFixed(v2)) 
              {
			    rod1->getPenaltyForce()->addRodClumpingForce(v1, rod2, v2);
				rod2->getPenaltyForce()->addRodClumpingForce(v2, rod1, v1);
			  }
			}
		  }
		}
	  }	  
	}
  }      
}


void RodCollisionTimeStepper::getProximities(CollisionMeshDataHashMap &collisionMeshes)
{
  if (is_multiple_stepper) return;

    if (m_rod == NULL)
    {
      cerr << "Rod is null in getProximities!\n";
      return;
    }
    //else
    //  cerr << "Doing collision proximity check\n";

    ElasticRod *rod = m_rod;

    Collisions colls;
    std::vector<uint> cands;
    
    Scalar radius_scale = rod->getRadiusScale();

		// Jungseock - This might be necessary for the future.
		m_rodPenaltyForce->clearPenaltyForces();

    // Check for proximities against every collision mesh
    //
    for (CollisionMeshDataHashMapIterator cmItr =collisionMeshes.begin();
                                          cmItr!=collisionMeshes.end(); ++cmItr)
    {
        CollisionMeshData *cmData = cmItr->second;

        // Check every vertex of the rod
        //
        for (int i=0; i<rod->nv(); ++i)
        {
            // Ignore fixed vertices (you can't apply a force to them anyway)
            //
            if (!rod->vertFixed(i))
            {
                // Build a bounding box around the vertex, extruding it by
                // the collision mesh thickness plus the radius, which is our proximity
                //
                Vec3d vMin = *(Vec3d *)&rod->getVertex(i);
                Vec3d vMax = *(Vec3d *)&rod->getVertex(i);
                Vec3d thickness((cmData->getThickness() + radius_scale * rod->radius()),
                                (cmData->getThickness() + radius_scale * rod->radius()),
                                (cmData->getThickness() + radius_scale * rod->radius()));
                vMin += -thickness;
                vMax += thickness;

                // We're re-using the uint vector (to save memory allocs), so clear out
                // old stuff
                //
                cands.clear();

                // Rasterizes our bounding box (vMin, vMax) to the grid, getting a list
                // of cells that contain it, then returns all the other elements
                // (in this case collision mesh triangles) that are also contained in
                // each of those cells, also makes sure the bounding boxes overlap
                //
                cmData->_grid.findOverlappingElements(vMin, vMax, cands);
                for (std::vector<uint>::iterator uiItr=cands.begin(); uiItr!=cands.end(); ++uiItr)
                {
                    // Do a full check that the vertex and triangle are actually within
                    // proximity
                    //
                    CandidateCollision cand(rod, i, cmData, *uiItr, VERTEX_TRIANGLE);
                    if (cand.getProximity(colls))
                    {
                        //cerr << "Found candidate Collision!\n";
                        // If they are, add a penalty force for this vertex-triangle pair
                        //
                        m_rodPenaltyForce->addRodPenaltyForce(cand.getFirstPrimitive(),
                                               cmData, cand.getSecondPrimitive());
                    }
                    //else
                   //     cerr << "FAILED to find candidate Collision!\n";
                }
            }
        }
    }
}


void RodCollisionTimeStepper::respondObjectCollisions(CollisionMeshDataHashMap &collisionMeshes, Real dt)
{
    // Mostly the same as the getProximities function above, except where noted
    //

//  return;
  
  if (!m_collisionsEnabled || !impulse_enabled) return;

    ElasticRod *rod = m_rod;
    Scalar radius_scale = rod->getRadiusScale();
    
    std::vector<uint> cands;
    for (CollisionMeshDataHashMapIterator cmItr=collisionMeshes.begin(); cmItr!=collisionMeshes.end(); ++cmItr)
    {
        CollisionMeshData *cmData = cmItr->second;

        CandidateCollisionSet candCollisions;
        for (int i=0; i<rod->nv(); ++i)
        {
            if (!rod->vertFixed(i))
            {
                // We're finding continuous time collisions, so build bounds (min and max)
                // around both the start-of-timestep and the end-of-timestep positions
                // Inflate them by some small epsilon amount
                //
                Vec3d vMin, vMax;
                minmax(rod->getStartPositions()[i],
                       rod->getEndPositions()[i],
                       vMin, vMax);
                Vec3d thickness((cmData->getThickness() + radius_scale * rod->radius()),
                                (cmData->getThickness() + radius_scale * rod->radius()),
                                (cmData->getThickness() + radius_scale * rod->radius()));
                vMin += -thickness;
                vMax += thickness;
                
                cands.clear();
                cmData->_grid.findOverlappingElements(vMin, vMax, cands);
                for (std::vector<uint>::iterator uiItr=cands.begin(); uiItr!=cands.end(); ++uiItr)
                    candCollisions.push_back(CandidateCollision(rod, i, cmData, *uiItr, VERTEX_TRIANGLE));
            }
        }

        // Full collisions means we check the rod segments against the grid as well
        //
        if (cmData->getFullCollisions())
        {
            for (int i=0; i<rod->nv()-1; ++i)
            {
                if ( !rod->vertFixed(i) && !rod->vertFixed((i+1)) )
                {
                    // Start and end positions for both rod vertices, so four points in all
                    //
                    Vec3d vMin, vMax;
                    minmax(rod->getStartPositions()[i],
                           rod->getEndPositions()[i],
                           rod->getStartPositions()[(i+1)],
                           rod->getEndPositions()[(i+1)],
                           vMin, vMax);
                    Vec3d thickness((cmData->getThickness() + radius_scale * rod->radius()),
                                    (cmData->getThickness() + radius_scale * rod->radius()),
                                    (cmData->getThickness() + radius_scale * rod->radius()));
                    vMin += -thickness;
                    vMax += thickness;

                    cands.clear();
                    cmData->_grid.findOverlappingElements(vMin, vMax, cands);
                    for (std::vector<uint>::iterator uiItr=cands.begin(); uiItr!=cands.end(); ++uiItr)
                    {
                        // We are looking for edge-edge collisions (vertex-triangle are covered)
                        // So check the rod segment against each of the tri's 3 edges
                        //
                        for (uint j=0; j<3; ++j)
                            candCollisions.push_back(CandidateCollision(rod, i, cmData, cmData->_triangleEdgeIndices[*uiItr][j], EDGE_EDGE));
                    }
                }
            }
            
            // This is how David wrote it for Rodney
            /*for (int i=0; i<rod->ne(); ++i)
            {
                if (!rod->vertFixed(i) && !rod->vertFixed((i+1)%rod->nv()))
                {
                    Can I just remove the % nv in here?
                    
                    // Start and end positions for both rod vertices, so four points in all
                    //
                    Vec3d vMin, vMax;
                    minmax(rod->getStartPositions()[i],
                           rod->getEndPositions()[i],
                           rod->getStartPositions()[(i+1)%rod->nv()],
                           rod->getEndPositions()[(i+1)%rod->nv()],
                           vMin, vMax);
                    Vec3d thickness((cmData->getThickness() + rod->radius()),
                                    (cmData->getThickness() + rod->radius()),
                                    (cmData->getThickness() + rod->radius()));
                    vMin += -thickness;
                    vMax += thickness;

                    cands.clear();
                    cmData->_grid.findOverlappingElements(vMin, vMax, cands);
                    for (std::vector<uint>::iterator uiItr=cands.begin(); uiItr!=cands.end(); ++uiItr)
                    {
                        // We are looking for edge-edge collisions (vertex-triangle are covered)
                        // So check the rod segment against each of the tri's 3 edges
                        //
                        for (uint j=0; j<3; ++j)
                            candCollisions.push_back(CandidateCollision(rod, i, cmData, cmData->_triangleEdgeIndices[*uiItr][j], EDGE_EDGE));
                    }
                }
            }*/
        }
        
        Collisions colls;
        for (CandidateCollisionSetIterator ccsItr=candCollisions.begin();
                                           ccsItr!=candCollisions.end(); ++ccsItr)
        {
            // We have the candidates, find those which are actual collisions
            //
            ccsItr->getContinuousTime(dt, colls);
        }

//        std::cout << "Collision time stepper : detected collisions in respondObjectCollisions : " << colls.size() << "\n";
        // For all collisions, apply an inelastic impulse
        //
        for (CollisionsIterator cItr=colls.begin(); cItr!=colls.end(); ++cItr) {
            cItr->applyImpulse(0.01);
        }
    }
}

void RodCollisionTimeStepper::respondRodCollisions(vector<ElasticRod*> &rods, Real dt, int maxIterations, Real COR)
{
    int iterations = 0;

    // Below is the code that applies iterative impulses//

    // Only do anything if full (continuous time) collisions are turned on and
    // the maximum number of iterations is positive
    //
    Collisions collisions;
    if ( maxIterations > 0)
    {
        // Loop until there are no more collisions are the maximum number of iterations
        // is reached
        //
        do
        {
            // TODO: Re-use the grid somehow?
            //
            UniformGrid grid;
            collisions.clear();
            grid.getContinuousTimeCollisions(rods, dt, collisions);

            for (CollisionsIterator cItr=collisions.begin(); cItr!=collisions.end(); ++cItr)
                cItr->applyImpulse(COR);

            ++iterations;
        } while (collisions.size() && iterations < maxIterations);
    
        if (collisions.size() && iterations == maxIterations)
        {
            // TODO: Fail-safe call would go here
            //
            std::cout << "WARNING: Maximum number of iterations reached. There may be unresolved collisions." << std::endl;
        }

        if (collisions.size())
        {
            // If we actually modified any velocities, update the candidate end positions
            // for any collision detection done after this
            //
            for (vector<ElasticRod*>::iterator rItr=rods.begin(); rItr!=rods.end(); ++rItr)
                (*rItr)->updateEndPositions(dt);
        }
    }
}


} // namespace BASim
