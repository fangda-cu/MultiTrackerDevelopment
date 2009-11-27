// CollisionObject.hh
//
// Anything that wants to collide and be collided must implement
// this interface. All collision detection and response routines
// take object pointers of this base type
//

#ifndef COLLISIONOBJECT_HH
#define COLLISIONOBJECT_HH

#include <BASim/Math>
#include "BridsonArray.hh"
//#include "BridsonVec.hh"

namespace BASim {

typedef std::vector<Vec3d> Positions;
typedef std::vector<Vec3d> Velocities;
typedef std::vector<unsigned int> Indices;

// This should be unsigned int, doubt it matters but should change and check results.
typedef std::vector<Vec3i> Vec3Indices;
typedef double Real;
    
class CollisionObject
{
public:

    virtual ~CollisionObject() {};
    virtual Positions&  getPositions()=0;
    virtual Velocities& getVelocities()=0;

    virtual Indices& getEdgeIndices()=0;
    virtual Indices& getTriangleIndices()=0;

    virtual Real getMassInverse(uint vertex)=0;

    virtual Real getFrictionCoefficient()=0;

    virtual Real getThickness()=0;

};

}

#endif

