// Collision.hh
//
// This class allows the user to respond to collisions between
// two primitives
//

#ifndef COLLISION_H
#define COLLISION_H

#include "CollisionObject.hh"

namespace BASim {
    
// All possible collision types, more may be added.
// is currently not implemented
//
enum CollisionType { VERTEX_TRIANGLE=0, EDGE_EDGE };


class Collision
{
public:
    Collision();
    Collision(uint index1, uint index2, uint index3, uint index4,
              Real dist, Vec3d &n, CollisionType type);
    ~Collision();

    CollisionObject *getFirstObject() const
    { return _object1; }

    CollisionObject *getSecondObject() const
    { return _object2; }

    uint getFirstPrimitiveIndex(int index) const;
    uint getSecondPrimitiveIndex(int index) const;

    void setFirstPrimitiveIndex(int index, uint val)
    { _primitiveIndex1[index] = val; }

    void setSecondPrimitiveIndex(int index, uint val)
    { _primitiveIndex2[index] = val; }

    void setFirstObject(CollisionObject *object1)
    { _object1 = object1; }
    
    void setSecondObject(CollisionObject *object2)
    { _object2 = object2; }

    void setBarycentricCoordinates(Real u, Real v, Real w=0.0)
    { _u = u; _v = v; _w = w; }

    Real getBarycentricCoordinateU()
    { return _u; }

    Real getBarycentricCoordinateV()
    { return _v; }

    Real getBarycentricCoordinateW()
    { return _w; }

    void setDistance(Real dist)
    { _distance = dist; }

    Real getDistance()
    { return _distance; }

    void setNormal(Vec3d &n)
    { _normal = n; }

    void setType(CollisionType type)
    { _type = type; }

    // Currently the only response supported, applies an impulse
    // with the specified coefficient of restitution (COR)
    //
    void applyImpulse(Real COR);

protected:

    CollisionObject *_object1;
    CollisionObject *_object2;

    uint _primitiveIndex1[2];
    uint _primitiveIndex2[3];

    Real _u;
    Real _v;
    Real _w;

    Real _distance;

    Vec3d _normal;

    CollisionType _type;

private:

};

typedef std::vector<Collision> Collisions;
typedef std::vector<Collision>::iterator CollisionsIterator;

}

#endif

