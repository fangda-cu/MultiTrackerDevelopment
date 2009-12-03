
#include "RodCollisionTimeStepper.hh"
#include "RodPenaltyForce.hh"

namespace BASim {

RodCollisionTimeStepper::RodCollisionTimeStepper(ObjectControllerBase* rodTimeStepper, ElasticRod* rod)
    : m_rodTimeStepper(rodTimeStepper), m_rod(rod)
{
  m_rodPenaltyForce = new RodPenaltyForce();
  dynamic_cast<RodTimeStepper*>(m_rodTimeStepper)->addExternalForce(m_rodPenaltyForce);  
}

RodCollisionTimeStepper::~RodCollisionTimeStepper()
{
  if (m_rodTimeStepper != NULL)
    delete m_rodTimeStepper;
  if (m_rodPenaltyForce != NULL)
    delete m_rodPenaltyForce;
}

void RodCollisionTimeStepper::getProximities(CollisionMeshDataHashMap &collisionMeshes)
{
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
                Vec3d thickness((cmData->getThickness() + rod->radius()),
                                (cmData->getThickness() + rod->radius()),
                                (cmData->getThickness() + rod->radius()));
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
                       // cerr << "Found candidate Collision!\n";
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
    
    ElasticRod *rod = m_rod;

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
                Vec3d thickness((cmData->getThickness() + rod->radius()),
                                (cmData->getThickness() + rod->radius()),
                                (cmData->getThickness() + rod->radius()));
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
      //  if (cmData->getFullCollisions())
        {
            for (int i=0; i<rod->ne(); ++i)
            {
                if (!rod->vertFixed(i) && !rod->vertFixed((i+1)%rod->nv()))
                {
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
            }
        }

        Collisions colls;
        for (CandidateCollisionSetIterator ccsItr=candCollisions.begin();
                                           ccsItr!=candCollisions.end(); ++ccsItr)
        {
            // We have the candidates, find those which are actual collisions
            //
            ccsItr->getContinuousTime(dt, colls);
        }

        // For all collisions, apply an inelastic impulse
        //
        for (CollisionsIterator cItr=colls.begin(); cItr!=colls.end(); ++cItr)
            cItr->applyImpulse(0.01);
    }
}


} // namespace BASim
