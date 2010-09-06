/* 
 *
 * See http://twiki.wetafx.co.nz/Weta/VolumetricHairCollisions for more documentation
 *
 */


#ifndef __VOLUMETRIC_COLLISIONS_HH__
#define __VOLUMETRIC_COLLISIONS_HH__

//#include "../Beaker.hh"

#ifdef WETA
#include <weta/Wfigaro/Core/EigenIncludes.hh>
#include <weta/Wfigaro/Physics/ElasticRods/ElasticRod.hh>
#include <weta/Wfigaro/Collisions/CollisionMeshData.hh>
#include <weta/Wfigaro/Collisions/Collision.hh>
#else
#include <BASim/src/Core/EigenIncludes.hh>
#include <BASim/src/Physics/ElasticRod.hh>
#include <BASim/src/Collisions/CollisionMeshData.hh>
#include <BASim/src/Collisions/Collision.hh>
#endif

#include <vector>
#include <iostream>
#include <cmath>

#include "../WmFigRodGroup.hh"

#include <Geometry/SEGMENTED_CURVE.h>
#include <Grids/GRID_3D.h>
#include <Incompressible_Flows/INCOMPRESSIBLE_UNIFORM.h>
#include <Matrices_And_Vectors/VECTOR.h>

using namespace std;

class VolumetricCollisions
{

public:
    VolumetricCollisions() {}
    virtual ~VolumetricCollisions() {}


    
    virtual void respondVolumetricCollisions( RodDataMap &rodDataMap, Real targetEdgeDensity, 
				     Real volumetricRadius, Real gridDx, Vec3d separationCondition,  double flip,
                     double slip, CollisionMeshDataHashMap &collisionMeshes) = 0;

    virtual void draw( bool displayGrid, Real displayGridVelocitiesMultiplier, 
		               Real maxDisplayDensity,
		               bool displayPsiN, bool displayPsiD ) = 0;

private:

    virtual void reinitialiseRods(  RodDataMap &rodDataMap ) = 0;
    virtual void deleteRods() = 0;

};

#endif //__VOLUMETRIC_COLLISIONS_HH__
