// CandidateCollision.hh
//
// Two primitives might be colliding, this class does the low-level
// tests to see if they actually are
//

#ifndef _CANDIDATE_COLLISION_H_
#define _CANDIDATE_COLLISION_H_

#include <tr1/unordered_map>
#include "Collision.hh"
//#include <ext/hash_set>

//class ElasticRod;

namespace BASim {

class CandidateCollision
{
public:
    CandidateCollision(CollisionObject *object1, uint primitive1,
                       CollisionObject *object2, uint primitive2, CollisionType type);

    ~CandidateCollision();

    CollisionObject *getFirstObject() const
    { return _object1; }

    CollisionObject *getSecondObject() const
    { return _object2; }

    uint getFirstPrimitive() const
    { return _primitive1; }
    
    uint getSecondPrimitive() const
    { return _primitive2; }

    CollisionType getCollisionType() const
    { return _type; }

    // Checks if two primitives intersect during the interval [0,dt]
    //
    bool getContinuousTime(Real dt, Collisions &collisions) const;

    // Checks if two primitives are within proximity, where proximity
    // is defined as the sum of their two object's thicknesses
    //
    bool getProximity(Collisions &collisions) const;

	bool operator<(const CandidateCollision &cc) const;
	bool operator==(const CandidateCollision &cc) const;
    
  

protected:

    CollisionObject *_object1;
    CollisionObject *_object2;

    uint _primitive1;
    uint _primitive2;

    CollisionType _type;

protected:

    // Once you know what type of collision you are looking at
    // (types are defined in Collision.hh), do these low-level tests
    //

    bool getContinuousTimeVertexTriangle(Vec3d &x0,
                                         Vec3d &x1,
                                         Vec3d &x2,
                                         Vec3d &x3,
                                         Vec3d &v0,
                                         Vec3d &v1,
                                         Vec3d &v2,
                                         Vec3d &v3,
                                         Real dt, Collision &collision) const;

    bool getProximityVertexTriangle(Vec3d &v,
                                    Vec3d &t1, Vec3d &t2, Vec3d &t3,
                                    Real h, Collision &collision) const;


    Real getClosestPointsVertexTriangle(const Vec3d& v0, const Vec3d& v1,
                                        const Vec3d& v2, const Vec3d& v3,
                                        Real &t1, Real &t2, Real &t3) const;

    bool getContinuousTimeEdgeEdge(Vec3d &x00,
                                   Vec3d &x01,
                                   Vec3d &x10,
                                   Vec3d &x11,
                                   Vec3d &v00,
                                   Vec3d &v01,
                                   Vec3d &v10,
                                   Vec3d &v11,
                                   Real dt, Collision &collision) const;

    bool getProximityEdgeEdge(Vec3d &e11, Vec3d &e12,
                              Vec3d &e21, Vec3d &e22,
                              Real h, Collision &collision) const;


    Real getClosestPointsEdgeEdge(const Vec3d& e11, const Vec3d& e12,
                                  const Vec3d& e21, const Vec3d& e22,
                                  Real &s, Real &t) const;

    void getCoplanarityTimes(const Vec3d &x0, const Vec3d &x1,
                             const Vec3d &x2, const Vec3d &x3,
                             const Vec3d &v0, const Vec3d &v1,
                             const Vec3d &v2, const Vec3d &v3,
                             std::vector<Real> &times, std::vector<Real> &errors) const;

    Real signed_volume(const Vec3d &x0, const Vec3d &x1,
                       const Vec3d &x2, const Vec3d &x3) const;

private:
    CandidateCollision();

};

class CCHash
{
public:
    size_t operator()(const CandidateCollision &cc) const
    {
        return (cc.getFirstPrimitive() * cc.getSecondPrimitive());
    }
};

//typedef __gnu_cxx::hash_set<CandidateCollision, CCHash> CandidateCollisionSet;
//typedef __gnu_cxx::hash_set<CandidateCollision, CCHash>::iterator CandidateCollisionSetIterator;
typedef std::vector<CandidateCollision> CandidateCollisionSet;
typedef std::vector<CandidateCollision>::iterator CandidateCollisionSetIterator;

}

#endif

