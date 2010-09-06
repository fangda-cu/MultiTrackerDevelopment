/*
 *
 * See http://twiki.wetafx.co.nz/Weta/VolumetricHairCollisions for more documentation
 *
 */

#ifndef __VOLUMETRIC_COLLISIONS_GPU_HH__
#define __VOLUMETRIC_COLLISIONS_GPU_HH__

#include <vector>
#include <iostream>
#include <cmath>

#include "objects/Rod.hh"
#include "collisions/CollisionMeshData.hh"
#include "collisions/Collision.hh"
#include "collisions/VolumetricCollisions.hh"

#include <ocustorage/grid3d.h>

using namespace std;

class VolumetricCollisionsGPU : public VolumetricCollisions
{

    typedef ocu::Grid3DHostF Grid3DHostF;
    typedef ocu::Grid3DDeviceF Grid3DDeviceF;
    typedef ocu::Grid3DDeviceB Grid3DDeviceB;

public:
    VolumetricCollisionsGPU();
    VolumetricCollisionsGPU(Rods &rods);
    ~VolumetricCollisionsGPU();

 
    void respondVolumetricCollisions(Rods &rods, Real targetEdgeDensity, 
                                     Real volumetricRadius, Real gridDx, Vec3<Real> separationCondition,
                                     CollisionMeshDataHashMap &collisionMeshes);

    void draw(bool displayGrid, Real displayGridVelocitiesMultiplier, Real maxDisplayDensity,
              bool displayPsiN, bool displayPsiD);

private:
    
    void reinitialiseRods(Rods &rods);
    void deleteRods();


    bool _initialised;
    uint _numEdges;
    uint _numVertices;
    uint _numOriginalVertices;

    float _radius;

    // CPU data
    float* _hPos; 
    float* _hVel;

   
    // GPU data
    float* _dPos;
    float* _dVel;
    float* _dSortedPos;
    float* _dSortedVel;
    
    // grid data for sorting method
    uint* _dGridHash;
    uint* _dGridIndex;
    uint _gridSortBits;

    bool printdbg;

    
};

#endif //__VOLUMETRIC_COLLISIONS_GPU_HH__
