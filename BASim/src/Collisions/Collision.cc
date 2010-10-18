// Collision.cc
//

#ifndef _COLLISION_CC_
#define _COLLISION_CC_

#include "Collision.hh"

//#include "objects/Rod.hh"

namespace BASim {

Collision::Collision()
 : _object1(0), _object2(0), _u(0.0), _v(0.0), _w(0.0)
{
}

Collision::Collision(uint index1, uint index2, uint index3, uint index4,
                     Real dist, Vec3d &n, CollisionType type)
 : _object1(0), _object2(0), _distance(dist), _normal(n), _type(type)
{
    switch (_type)
    {
        case VERTEX_TRIANGLE:
        {
            _primitiveIndex1[0] = index1;
            _primitiveIndex2[0] = index2;
            _primitiveIndex2[1] = index3;
            _primitiveIndex2[2] = index4;

            break;
        }
        case EDGE_EDGE:
        {
            _primitiveIndex1[0] = index1;
            _primitiveIndex1[1] = index2;
            _primitiveIndex2[0] = index3;
            _primitiveIndex2[1] = index4;

            break;
        }
    }
}

Collision::~Collision()
{
}

uint Collision::getFirstPrimitiveIndex(int index) const
{
    if (_type == VERTEX_TRIANGLE)
        assert(index >= 0 && index < 1);
    else if (_type == EDGE_EDGE)
        assert(index >= 0 && index < 2);

    return (_primitiveIndex1[index]);
}

uint Collision::getSecondPrimitiveIndex(int index) const
{
    if (_type == VERTEX_TRIANGLE)
        assert(index >= 0 && index < 3);
    else if (_type == EDGE_EDGE)
        assert(index >= 0 && index < 2);

    return (_primitiveIndex2[index]);
}

void Collision::applyImpulse(Real COR)
{
    // First get the position and velocity vectors from the objects colliding
    //
    Positions&  x1 = _object1->getPositions();
    Positions&  x2 = _object2->getPositions();
    Velocities& v1 = _object1->getVelocities();
    Velocities& v2 = _object2->getVelocities();

    // Friction here is the maximum of the two objects' frictions, there are different
    // ways to do this so you may want to play around with it. I.e., take the average
    //
    Real friction = std::max(_object1->getFrictionCoefficient(), _object2->getFrictionCoefficient());

    switch (_type)
    {
        case VERTEX_TRIANGLE:
        {
            // Consider the object as stationary by using the relative velocity
            // If you ever need vertex-triangle response between two simulated
            // objects this code will have to be modified
            //
            assert(_primitiveIndex2[0] < v2.size());
            assert(_primitiveIndex2[1] < v2.size());
            assert(_primitiveIndex2[2] < v2.size());
            assert(_primitiveIndex1[0] < v1.size() && _primitiveIndex1[0] < x1.size());

            Vec3d& vVertex = v1[_primitiveIndex1[0]];
            Vec3d& vObj1   = v2[_primitiveIndex2[0]];
            Vec3d& vObj2   = v2[_primitiveIndex2[1]];
            Vec3d& vObj3   = v2[_primitiveIndex2[2]];

            Vec3d vTriangle = (_u * vObj1 + _v * vObj2 + _w * vObj3);

            Vec3d vRel = (vVertex - vTriangle);
            Real normalRelVelReal = vRel.dot(_normal);
            Vec3d normalRelVel = normalRelVelReal * _normal;

            // Negative relative velocity means they are approaching
            //
            if (normalRelVelReal < 0.0)
            {
                // Vertex gets the whole kick since triangle is assumed to be fixed
                //
                vVertex -= (1.0 + COR) * normalRelVel;
  
                // Calcuate the new relative velocity and make sure our kick worked
                //
                Real normalRelVelNewReal = (vVertex - vTriangle).dot(_normal);
                assert(normalRelVelNewReal > -1e-6);
            
                if (friction > 0.0)
                {
                    // Tangent (sliding) velocity is the non-normal component of
                    // the original relative velocity
                    //
                    Vec3d vTan  = vRel - normalRelVel;
                    Real tangentVel = vTan.norm();

                    // Only apply friction if there is sliding motion
                    //
                    if (tangentVel > 1e-6)
                    {
                        // Normalize tangent vector and find friction using Coulomb's Law
                        //
                        vTan /= tangentVel;
                        Real Ic = std::min(friction * std::fabs(normalRelVelNewReal - normalRelVelReal),
                                           tangentVel);
                       
                        vVertex -= Ic * vTan;
                    }
                }
            }
            break;
        }
        case EDGE_EDGE:
        {
            Vec3d& v00 = v1[_primitiveIndex1[0]];
            Vec3d& v01 = v1[_primitiveIndex1[1]];
            Vec3d& v10 = v2[_primitiveIndex2[0]];
            Vec3d& v11 = v2[_primitiveIndex2[1]];

            Vec3d vRel = ((1.0 - _u) * v00 + _u * v01) -
                         ((1.0 - _v) * v10 + _v * v11);

            Vec3d vRelN = vRel.dot(_normal) * _normal;

            Real vRelReal = vRel.dot(_normal);
            if (vRelReal < 0.0)
            {
                Real m00 = _object1->getMassInverse(_primitiveIndex1[0]);
                Real m01 = _object1->getMassInverse(_primitiveIndex1[1]);
                Real m10 = _object2->getMassInverse(_primitiveIndex2[0]);
                Real m11 = _object2->getMassInverse(_primitiveIndex2[1]);

                // Weighted impulse, see [Bridson et al 2002] and Bridson's
                // SIGGRAPH slides on collisions
                //
                Real J = ((1.0 + COR) * vRelReal) / ((1.0 - _u) * (1.0 - _u) * m00 +
                                                     (      _u) * (      _u) * m01 +
                                                     (1.0 - _v) * (1.0 - _v) * m10 +
                                                     (      _v) * (      _v) * m11);

                v00 -= (1.0 - _u) * J * m00 * _normal;
                v01 -= (      _u) * J * m01 * _normal;
                v10 += (1.0 - _v) * J * m10 * _normal;
                v11 += (      _v) * J * m11 * _normal;

                // Get new relative velocity in the normal direction
                //
                Vec3d vRelNew = ((1.0 - _u) * v00 + _u * v01) -
                                ((1.0 - _v) * v10 + _v * v11);
                Real vRelNnew = vRelNew.dot(_normal);

                if (friction > 0.0)
                {
                    // Tangent (sliding) direction is component of relative velocity
                    // that isn't in the normal direction (vRel = vRelN + vTan)
                    //
                    Vec3d vTan  = vRel - vRelN;
                    Real tangentVel = vTan.norm();
                    if (tangentVel > 1e-6)
                    {
                        // Normalize tangent direction so we can apply impulse
                        //
                        vTan /= tangentVel;

                        // Coulomb's Law
                        //
                        J  = std::min(friction * std::fabs(vRelNnew - vRelReal), tangentVel);
                        J /= ((1.0 - _u) * (1.0 - _u) * m00 +
                              (      _u) * (      _u) * m01 +
                              (1.0 - _v) * (1.0 - _v) * m10 +
                              (      _v) * (      _v) * m11);

                        v00 -= (1.0 - _u) * J * m00 * _normal;
                        v01 -= (      _u) * J * m01 * _normal;
                        v10 += (1.0 - _v) * J * m10 * _normal;
                        v11 += (      _v) * J * m11 * _normal;
                    }
                }
            }

            break;
        }
        default:
        {
            assert(0);
        }
    }
}

}

#endif

