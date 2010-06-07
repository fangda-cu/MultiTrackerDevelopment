// CandidateCollision.cc
//

#include "CandidateCollision.hh"
#include "../Math/Math.hh"

#include <cmath>

namespace BASim {

#define COLLISION_EPSILON 1e-6

using namespace bridson;

CandidateCollision::CandidateCollision(CollisionObject *object1, uint primitive1,
                                       CollisionObject *object2, uint primitive2, CollisionType type)
 : _object1(object1), _object2(object2), _primitive1(primitive1), _primitive2(primitive2), _type(type)
{
}

CandidateCollision::~CandidateCollision()
{
}

bool CandidateCollision::operator<(const CandidateCollision &cc) const
{
    if (_type != cc._type)
        return (_type < cc._type);
    else if (_primitive1 != cc._primitive1)
        return (_primitive1 < cc._primitive1);
    else if (_primitive2 != cc._primitive2)
        return (_primitive2 < cc._primitive2);
    else
        return false;
}

bool CandidateCollision::operator==(const CandidateCollision &cc) const
{
    return ((_type       == cc._type)       &&
            (_primitive1 == cc._primitive1) &&
            (_primitive2 == cc._primitive2) &&
            (_object1    == cc._object1) &&
            (_object2    == cc._object2));
}

bool CandidateCollision::getContinuousTime(Real dt, Collisions &collisions) const
{
    Positions&  x1 = _object1->getPositions();
    Positions&  x2 = _object2->getPositions();
    Velocities& v1 = _object1->getVelocities();
    Velocities& v2 = _object2->getVelocities();

    switch (_type)
    {
        case VERTEX_TRIANGLE:
        {
            Indices& indices2 = _object2->getTriangleIndices();

            Collision collision;
            if (getContinuousTimeVertexTriangle(x1[_primitive1],
                                                x2[indices2[(3 * _primitive2)    ]],
                                                x2[indices2[(3 * _primitive2) + 1]],
                                                x2[indices2[(3 * _primitive2) + 2]],
                                                v1[_primitive1],
                                                v2[indices2[(3 * _primitive2)    ]],
                                                v2[indices2[(3 * _primitive2) + 1]],
                                                v2[indices2[(3 * _primitive2) + 2]],
                                                dt, collision))
            {
                // Fill out the rest of the collision info
                //
                collision.setFirstPrimitiveIndex(0, _primitive1);
                collision.setSecondPrimitiveIndex(0, indices2[(3 * _primitive2)    ]);
                collision.setSecondPrimitiveIndex(1, indices2[(3 * _primitive2) + 1]);
                collision.setSecondPrimitiveIndex(2, indices2[(3 * _primitive2) + 2]);
                collision.setFirstObject(_object1);
                collision.setSecondObject(_object2);
                collisions.push_back(collision);

                return true;
            }

            return false;
        }
        case VERTEX_EDGE:
        {
            assert(0);
            break;
        }
        case EDGE_EDGE:
        {
            Indices& indices1 = _object1->getEdgeIndices();
            Indices& indices2 = _object2->getEdgeIndices();

            Collision collision;
            if (getContinuousTimeEdgeEdge(x1[indices1[(2 * _primitive1)    ]],
                                          x1[indices1[(2 * _primitive1) + 1]],
                                          x2[indices2[(2 * _primitive2)    ]],
                                          x2[indices2[(2 * _primitive2) + 1]],
                                          v1[indices1[(2 * _primitive1)    ]],
                                          v1[indices1[(2 * _primitive1) + 1]],
                                          v2[indices2[(2 * _primitive2)    ]],
                                          v2[indices2[(2 * _primitive2) + 1]],
                                          dt, collision))
            {
//                if (collision.getBarycentricCoordinateU() > 0.0 &&
//                    collision.getBarycentricCoordinateU() < 1.0)
                {
                    collision.setFirstPrimitiveIndex( 0, indices1[(2 * _primitive1)    ]);
                    collision.setFirstPrimitiveIndex( 1, indices1[(2 * _primitive1) + 1]);
                    collision.setSecondPrimitiveIndex(0, indices2[(2 * _primitive2)    ]);
                    collision.setSecondPrimitiveIndex(1, indices2[(2 * _primitive2) + 1]);
                    collision.setFirstObject(_object1);
                    collision.setSecondObject(_object2);
                    collisions.push_back(collision);

                    return true;
                }
            }

            return false;
        }
        default:
        {
            assert(0);
            break;
        }
    }

    return false;
}

bool CandidateCollision::getProximity(Collisions &collisions) const
{
    Positions&  x1 = _object1->getPositions();
    Positions&  x2 = _object2->getPositions();
    Velocities& v1 = _object1->getVelocities();
    Velocities& v2 = _object2->getVelocities();

    Real h = _object1->getThickness() + _object2->getThickness();

    switch (_type)
    {
        case VERTEX_TRIANGLE:
        {
            Indices& indices2 = _object2->getTriangleIndices();

            Collision collision;
            if (getProximityVertexTriangle(x1[_primitive1],
                                           x2[indices2[(3 * _primitive2)    ]],
                                           x2[indices2[(3 * _primitive2) + 1]],
                                           x2[indices2[(3 * _primitive2) + 2]],
                                           h, collision))
            {
                collision.setFirstPrimitiveIndex(0, _primitive1);
                collision.setSecondPrimitiveIndex(0, indices2[(3 * _primitive2)    ]);
                collision.setSecondPrimitiveIndex(1, indices2[(3 * _primitive2) + 1]);
                collision.setSecondPrimitiveIndex(2, indices2[(3 * _primitive2) + 2]);
                collision.setFirstObject(_object1);
                collision.setSecondObject(_object2);
                collisions.push_back(collision);

                return true;
            }

            return false;
        }
        case VERTEX_EDGE:
        {
            assert(0);
            break;
        }
        case EDGE_EDGE:
        {
            Indices& indices1 = _object1->getEdgeIndices();
            Indices& indices2 = _object2->getEdgeIndices();

            Collision collision;
            if (getProximityEdgeEdge(x1[indices1[(2 * _primitive1)    ]],
                                     x1[indices1[(2 * _primitive1) + 1]],
                                     x2[indices2[(2 * _primitive2)    ]],
                                     x2[indices2[(2 * _primitive2) + 1]],
                                     h, collision))
            {
                collision.setFirstPrimitiveIndex( 0, indices1[(2 * _primitive1)    ]);
                collision.setFirstPrimitiveIndex( 1, indices1[(2 * _primitive1) + 1]);
                collision.setSecondPrimitiveIndex(0, indices2[(2 * _primitive2)    ]);
                collision.setSecondPrimitiveIndex(1, indices2[(2 * _primitive2) + 1]);
                collision.setFirstObject(_object1);
                collision.setSecondObject(_object2);
                collisions.push_back(collision);

                return true;
            }

            return false;
        }
        default:
        {
            assert(0);
            break;
        }
    }

    return false;
}

bool CandidateCollision::getContinuousTimeVertexTriangle(Vec3d &x0,
                                                         Vec3d &x1,
                                                         Vec3d &x2,
                                                         Vec3d &x3,
                                                         Vec3d &v0,
                                                         Vec3d &v1,
                                                         Vec3d &v2,
                                                         Vec3d &v3,
                                                         Real dt, Collision &collision) const
{
    // Get the times when the vertex is coplanar to the triangle, a necessary condition
    // (but not sufficient) for a collision
    //
    std::vector<Real> times, errors;
    getCoplanarityTimes(x0, x1, x2, x3, x0+v0*dt, x1+v1*dt, x2+v2*dt, x3+v3*dt, times, errors);

    // Check each time in order
    //
    for (uint a=0; a<times.size(); ++a)
    {
        // Times are given between [0,1], so scale by timestep
        //
        Real t = times[a] * dt;

        Vec3d xt0 = x0 + t * v0;
        Vec3d xt1 = x1 + t * v1;
        Vec3d xt2 = x2 + t * v2;
        Vec3d xt3 = x3 + t * v3;

        Vec3d normal;
        Real s1, s2, s3;
        Real distance = std::sqrt(getClosestPointsVertexTriangle(xt0, xt1, xt2, xt3, s1, s2, s3));

        // Was this a real collision or a false positive?
        //
        if (distance < COLLISION_EPSILON)
        {
//            normal = xt0 - (s1 * xt1 + s2 * xt2 + s3 * xt3);

            // now figure out a decent normal
            //
//           if (distance < (1e-2 * COLLISION_EPSILON))
            {
                // if we don't trust the normal...
                // first try the triangle normal at collision time
                //
                normal = (xt2 - xt1).cross(xt3 - xt1);

                Real m = normal.norm();
                if (m > COLLISION_EPSILON)
                {
                    normal /= m;
                }
                else
                {
                    // if that didn't work, try triangle normal at start
                    //
                    normal = (x2 - x1).cross(x3 - x1);
                    m = normal.norm();
                    if (m > COLLISION_EPSILON)
                    {
                        normal /= m;
                    }
                    else
                    {
                        // if that didn't work, try vector between points at the start
                        //
                        normal = x0 - (s1 * x1 + s2 * x2 + s3 * x3);
                        m = normal.norm();
                        if (m > COLLISION_EPSILON)
                        {
                            normal /= m;
                        }
                        else
                        {
                            // if that didn't work, boy are we in trouble; just get any non-parallel vector
                            //
                            Vec3d dx = xt2 - xt1;
                            if (dx[0] != 0 || dx[1] != 0)
                            {
                                normal = Vec3d(dx[1], -dx[0], 0);
                                normal.normalize();
                            }
                            else
                            {
                                dx = xt3 - xt1;
                                if (dx[0] != 0 || dx[1] != 0)
                                {
                                    normal = Vec3d(dx[1], -dx[0], 0);
                                    normal.normalize();
                                }
                                else
                                {
                                    normal = Vec3d(0, 1, 0); // the last resort
                                }
                            }
                        }
                    }
                }
            }

            collision.setBarycentricCoordinates(s1, s2, s3);
            collision.setNormal(normal);
            collision.setDistance(distance);
            collision.setType(VERTEX_TRIANGLE);

            return true;
        }
    }

    return false;
}

bool CandidateCollision::getProximityVertexTriangle(Vec3d &v,
                                                    Vec3d &t1, Vec3d &t2, Vec3d &t3,
                                                    Real h, Collision &collision) const
{
    Real a1, a2, a3;
    Real distance = getClosestPointsVertexTriangle(v, t1, t2, t3, a1, a2, a3);

    if (distance < h * h && distance > COLLISION_EPSILON)
    {
//        Vec3d normal = (v - (a1 * t1 + a2 * t2 + a3 * t3));
        Vec3d normal = (t2 - t1).cross(t3 - t1);
        normal.normalize();

        collision.setNormal(normal);
        collision.setBarycentricCoordinates(a1, a2, a3);
        collision.setDistance(std::sqrt(distance));
        collision.setType(VERTEX_TRIANGLE);

        return true;
    }

    return false;
}

bool CandidateCollision::getContinuousTimeEdgeEdge(Vec3d &x00,
                                                   Vec3d &x01,
                                                   Vec3d &x10,
                                                   Vec3d &x11,
                                                   Vec3d &v00,
                                                   Vec3d &v01,
                                                   Vec3d &v10,
                                                   Vec3d &v11,
                                                   Real dt, Collision &collision) const
{
    std::vector<Real> times, errors;
    getCoplanarityTimes(x00, x01, x10, x11, x00+v00*dt, x01+v01*dt, x10+v10*dt, x11+v11*dt, times, errors);

    for (uint a=0; a<times.size(); ++a)
    {
        Real t = times[a] * dt;
       
        Vec3d xt00 = x00 + t * v00;
        Vec3d xt01 = x01 + t * v01;
        Vec3d xt10 = x10 + t * v10;
        Vec3d xt11 = x11 + t * v11;

        Vec3d normal;
        Real s1, s2;
        Real distance = std::sqrt(getClosestPointsEdgeEdge(xt00, xt01, xt10, xt11, s1, s2));

        if (distance < COLLISION_EPSILON)
        {
//            normal = ((1.0 - s1) * xt00 + s1 * xt01) -
//                     ((1.0 - s2) * xt10 + s2 * xt11);
//
//            // now figure out a decent normal
//            //
//            if (distance < 1e-10)
            {
                normal = (xt11 - xt10).cross(xt01 - xt00);

                Real m = normal.norm();
                if (m > COLLISION_EPSILON)
                {
                    normal /= m;
                }
                else
                {
                    normal = ((1.0 - s1) * x00 + s1 * x01) -
                             ((1.0 - s2) * x10 + s2 * x11);

                    m = normal.norm();
                    if (m > COLLISION_EPSILON)
                    {
                        normal /= m;
                    }
                    else
                    {
                        // if that didn't work, boy are we in trouble; just get any non-parallel vector
                        //
                        Vec3d dx = xt01 - xt00;
                        if (dx[0] != 0 || dx[1] != 0)
                        {
                            normal = Vec3d(dx[1], -dx[0], 0);
                            normal.normalize();
                        }
                        else
                        {
                            dx = xt11 - xt10;
                            if (dx[0] != 0 || dx[1] != 0)
                            {
                                normal = Vec3d(dx[1], -dx[0], 0);
                                normal.normalize();
                            }
                            else
                            {
                                normal = Vec3d(0, 1, 0); // the last resort
                            }
                        }
                    }
                }
            }
//            else
//                normalize(normal);


            if (!(normal.dot(( (1.0 - s1) * v00 + s1 * v01) - ((1.0 - s2) * v10 + s2 * v11)) < 0.0))
                normal *= -1;

            collision.setBarycentricCoordinates(s1, s2);
            collision.setNormal(normal);
            collision.setDistance(distance);
            collision.setType(EDGE_EDGE);

            return true;
        }
    }

    return false;
}

bool CandidateCollision::getProximityEdgeEdge(Vec3d &e11, Vec3d &e12,
                                              Vec3d &e21, Vec3d &e22,
                                              Real h, Collision &collision) const
{
    Real alpha, beta;
    Real distance = getClosestPointsEdgeEdge(e11, e12, e21, e22, alpha, beta);

    if (distance < h * h && distance > COLLISION_EPSILON)
    {
        Vec3d normal = ((1.0 - alpha) * e11 + alpha * e12) -
                            ((1.0 -  beta) * e21 +  beta * e22);
        normal.normalize();

        collision.setNormal(normal);
        collision.setBarycentricCoordinates(alpha, beta);
        collision.setDistance(std::sqrt(distance));
        collision.setType(EDGE_EDGE);

        return true;
    }

    return false;
}

Real CandidateCollision::getClosestPointsVertexTriangle(const Vec3d& v0, const Vec3d& v1,
                                                        const Vec3d& v2, const Vec3d& v3,
                                                        Real &t1, Real &t2, Real &t3) const

{
    const double *p = v0.data();
    const double *a = v1.data();
    const double *b = v2.data();
    const double *c = v3.data();

    Real ab[3], ac[3], ap[3], bp[3];

    ab[0] = b[0] - a[0];
    ab[1] = b[1] - a[1];
    ab[2] = b[2] - a[2];

    ac[0] = c[0] - a[0];
    ac[1] = c[1] - a[1];
    ac[2] = c[2] - a[2];

    ap[0] = p[0] - a[0];
    ap[1] = p[1] - a[1];
    ap[2] = p[2] - a[2];

    Real d1 = ab[0]*ap[0] + ab[1]*ap[1] + ab[2]*ap[2];
    Real d2 = ac[0]*ap[0] + ac[1]*ap[1] + ac[2]*ap[2];

    if ((d1 <= 0.0f) && (d2 <= 0.0f))
    {
        t1 = 1.0f;
        t2 = 0.0f;
        t3 = 0.0f;

        return ((p[0]-a[0])*(p[0]-a[0]) + (p[1]-a[1])*(p[1]-a[1]) + (p[2]-a[2])*(p[2]-a[2]));
    }

    bp[0] = p[0] - b[0];
    bp[1] = p[1] - b[1];
    bp[2] = p[2] - b[2];

    Real d3 = ab[0]*bp[0] + ab[1]*bp[1] + ab[2]*bp[2];
    Real d4 = ac[0]*bp[0] + ac[1]*bp[1] + ac[2]*bp[2];

    if ((d3 >= 0.0f) && (d4 <= d3))
    {
        t1 = 0.0f;
        t2 = 1.0f;
        t3 = 0.0f;

        return ((p[0]-b[0])*(p[0]-b[0]) + (p[1]-b[1])*(p[1]-b[1]) + (p[2]-b[2])*(p[2]-b[2]));
    }

    Real vc = d1*d4 - d3*d2;

    if ((vc <= 0.0f) && (d1 >= 0.0f) && (d3 <= 0.0f))
    {
        Real v = d1 / (d1 - d3);

        t1 = 1-v;
        t2 = v;
        t3 = 0;

        Real vec[3];
        vec[0] = p[0] - (a[0]+v*ab[0]);
        vec[1] = p[1] - (a[1]+v*ab[1]);
        vec[2] = p[2] - (a[2]+v*ab[2]);

        return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    }

    Real cp[3];
    cp[0] = p[0] - c[0];
    cp[1] = p[1] - c[1];
    cp[2] = p[2] - c[2];

    Real d5 = ab[0]*cp[0] + ab[1]*cp[1] + ab[2]*cp[2];
    Real d6 = ac[0]*cp[0] + ac[1]*cp[1] + ac[2]*cp[2];

    if ((d6 >= 0.0f) && (d5 <= d6))
    {
        t1 = 0;
        t2 = 0;
        t3 = 1;

        return ((p[0]-c[0])*(p[0]-c[0]) + (p[1]-c[1])*(p[1]-c[1]) + (p[2]-c[2])*(p[2]-c[2]));
    }

    Real vb = d5*d2 - d1*d6;

    if ((vb <= 0.0f) && (d2 >= 0.0f) && (d6 <= 0.0f))
    {
        Real w = d2 / (d2 - d6);

        t1 = 1-w;
        t2 = 0;
        t3 = w;

        Real vec[3];
        vec[0] = p[0] - (a[0]+w*ac[0]);
        vec[1] = p[1] - (a[1]+w*ac[1]);
        vec[2] = p[2] - (a[2]+w*ac[2]);

        return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    }

    Real va = d3*d6 - d5*d4;

    if ((va <= 0.0f) && ((d4-d3) >= 0.0f) && ((d5-d6) >= 0.0f))
    {
        Real w = (d4 - d3) / ((d4 - d3) + (d5 - d6));

        t1 = 0;
        t2 = 1-w;
        t3 = w;

        Real vec[3];
        vec[0] = p[0] - (b[0]+w*(c[0]-b[0]));
        vec[1] = p[1] - (b[1]+w*(c[1]-b[1]));
        vec[2] = p[2] - (b[2]+w*(c[2]-b[2]));

        return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
    }

    Real denom = 1.0f / (va + vb + vc);
    Real v = vb * denom;
    Real w = vc * denom;
    Real u = 1.0 - v - w;

    t1 = u;
    t2 = v;
    t3 = w;

    Real vec[3];
    vec[0] = p[0] - (u*a[0] + v*b[0] + w*c[0]);
    vec[1] = p[1] - (u*a[1] + v*b[1] + w*c[1]);
    vec[2] = p[2] - (u*a[2] + v*b[2] + w*c[2]);

    return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
}

Real CandidateCollision::getClosestPointsEdgeEdge(const Vec3d& e11, const Vec3d& e12,
                                                  const Vec3d& e21, const Vec3d& e22,
                                                  Real &s, Real &t) const
{
    const double *p1 = e11.data();
    const double *q1 = e12.data();
    const double *p2 = e21.data();
    const double *q2 = e22.data();

	Real d1[3], d2[3], r[3], a, e, f;
	Real c1[3], c2[3];

	d1[0] = q1[0] - p1[0];
	d1[1] = q1[1] - p1[1];
	d1[2] = q1[2] - p1[2];

	d2[0] = q2[0] - p2[0];
	d2[1] = q2[1] - p2[1];
	d2[2] = q2[2] - p2[2];

	r[0] = p1[0] - p2[0];
	r[1] = p1[1] - p2[1];
	r[2] = p1[2] - p2[2];

	a = d1[0]*d1[0] + d1[1]*d1[1] + d1[2]*d1[2];
	e = d2[0]*d2[0] + d2[1]*d2[1] + d2[2]*d2[2];
	f = d2[0]*r[0] + d2[1]*r[1] + d2[2]*r[2];

	// check if either or both segments degenerate into points
	//
	if ((a <= COLLISION_EPSILON) && (e <= COLLISION_EPSILON))
	{
		s = t = 0.0f;
		c1[0] = p1[0]; c1[1] = p1[1]; c1[2] = p1[2];
		c2[0] = p2[0]; c2[1] = p2[1]; c2[2] = p2[2];

		return ((c1[0]-c2[0])*(c1[0]-c2[0]) + (c1[1]-c2[1])*(c1[1]-c2[1]) + (c1[2]-c2[2])*(c1[2]-c2[2]));
	}

	if (a <= COLLISION_EPSILON)
	{
		// first segment degenerates into a point
		//
		s = 0.0f;
		t = f / e;
		if (t<0.0f) t = 0.0f;
		if (t>1.0f) t = 1.0f;
	}
	else
	{
		Real c = d1[0]*r[0] + d1[1]*r[1] + d1[2]*r[2];

		if (e <= COLLISION_EPSILON)
		{
			// second segment degenerates into a point
			//
			t = 0.0f;
			s = -c / a;
			if (s<0.0f) s = 0.0f;
			if (s>1.0f) s = 1.0f;
		}
		else
		{
			// nondegenerate case
			//
			Real b = d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2];
			Real denom = a*e - b*b;

			if (denom != 0.0f)
			{
				s = (b*f - c*e) / denom;
				if (s<0.0f) s = 0.0f;
				if (s>1.0f) s = 1.0f;
			}
			else
				s = 0.0f;

			Real tnom = b*s + f;
			if (tnom < 0.0f)
			{
				t = 0.0f;
				s = -c / a;
				if (s<0.0f) s = 0.0f;
				if (s>1.0f) s = 1.0f;
			}
			else if (tnom > e)
			{
				t = 1.0f;
				s = (b - c) / a;
				if (s<0.0f) s = 0.0f;
				if (s>1.0f) s = 1.0f;
			}
			else
				t = tnom / e;
		}
	}

	c1[0] = p1[0] + d1[0] * s;
	c1[1] = p1[1] + d1[1] * s;
	c1[2] = p1[2] + d1[2] * s;

	c2[0] = p2[0] + d2[0] * t;
	c2[1] = p2[1] + d2[1] * t;
	c2[2] = p2[2] + d2[2] * t;

	return ((c1[0]-c2[0])*(c1[0]-c2[0]) + (c1[1]-c2[1])*(c1[1]-c2[1]) + (c1[2]-c2[2])*(c1[2]-c2[2]));

}

void CandidateCollision::getCoplanarityTimes(const Vec3d &x0, const Vec3d &x1,
                                             const Vec3d &x2, const Vec3d &x3,
                                             const Vec3d &xnew0, const Vec3d &xnew1,
                                             const Vec3d &xnew2, const Vec3d &xnew3,
                                             std::vector<Real> &times, std::vector<Real> &errors) const
{
    const Real tol = 1e-8;
    times.clear();
    errors.clear();

    // cubic coefficients, A*t^3+B*t^2+C*t+D (for t in [0,1])
    Vec3d x03=x0-x3, x13=x1-x3, x23=x2-x3;
    Vec3d v03=(xnew0-xnew3)-x03, v13=(xnew1-xnew3)-x13, v23=(xnew2-xnew3)-x23;

    double A=triple(v03,v13,v23),
           B=triple(x03,v13,v23)+triple(v03,x13,v23)+triple(v03,v13,x23),
           C=triple(x03,x13,v23)+triple(x03,v13,x23)+triple(v03,x13,x23),
           D=triple(x03,x13,x23);

//	Real x21[3], v21[3], x31[3], v31[3], x41[3], v41[3];
//	for (int i=0; i<3; ++i)
//	{
//		x21[i] = x1[i] - x0[i];
//		x31[i] = x2[i] - x0[i];
//		x41[i] = x3[i] - x0[i];
//		v21[i] = v1[i] - v0[i];
//		v31[i] = v2[i] - v0[i];
//		v41[i] = v3[i] - v0[i];
//	}

	// The polynomial coefficients
	//
//	Real A = -v21[2]*v31[1]*v41[0] + v21[1]*v31[2]*v41[0] + v21[2]*v31[0]*v41[1] -
//			  v21[0]*v31[2]*v41[1] - v21[1]*v31[0]*v41[2] + v21[0]*v31[1]*v41[2];
//
//	Real B = -v31[2]*v41[1]*x21[0] + v31[1]*v41[2]*x21[0] + v31[2]*v41[0]*x21[1] -
//			  v31[0]*v41[2]*x21[1] - v31[1]*v41[0]*x21[2] + v31[0]*v41[1]*x21[2] +
//			  v21[2]*v41[1]*x31[0] - v21[1]*v41[2]*x31[0] - v21[2]*v41[0]*x31[1] +
//			  v21[0]*v41[2]*x31[1] + v21[1]*v41[0]*x31[2] - v21[0]*v41[1]*x31[2] -
//			  v21[2]*v31[1]*x41[0] + v21[1]*v31[2]*x41[0] + v21[2]*v31[0]*x41[1] -
//			  v21[0]*v31[2]*x41[1] - v21[1]*v31[0]*x41[2] + v21[0]*v31[1]*x41[2];
//
//	Real C = -v41[2]*x21[1]*x31[0] + v41[1]*x21[2]*x31[0] + v41[2]*x21[0]*x31[1] -
//			  v41[0]*x21[2]*x31[1] - v41[1]*x21[0]*x31[2] + v41[0]*x21[1]*x31[2] +
//			  v31[2]*x21[1]*x41[0] - v31[1]*x21[2]*x41[0] - v31[2]*x21[0]*x41[1] +
//			  v31[0]*x21[2]*x41[1] + v21[2]*x31[0]*x41[1] - v21[0]*x31[2]*x41[1] +
//			  v31[1]*x21[0]*x41[2] - v31[0]*x21[1]*x41[2] - v21[1]*x31[0]*x41[2] -
//			  v21[2]*x31[1]*x41[0] + v21[1]*x31[2]*x41[0] + v21[0]*x31[1]*x41[2];
//
//	Real D = -x21[2]*x31[1]*x41[0] + x21[1]*x31[2]*x41[0] + x21[2]*x31[0]*x41[1] -
//			  x21[0]*x31[2]*x41[1] - x21[1]*x31[0]*x41[2] + x21[0]*x31[1]*x41[2];

/*
	int leadCoeff = 3;
	Real coeffs[4];
	if (a == 0.0)
	{
		leadCoeff = 2;
		if (b == 0.0)
		{
			leadCoeff = 1;
			if (c == 0.0)
			{
				// Degenerate polynomial
				//
				return;
			}
			else
			{
                times.push_back(-d / c);
                return;
//				coeffs[0] = c;
//				coeffs[1] = d;
			}
		}
		else
		{
			coeffs[0] = b;
			coeffs[1] = c;
			coeffs[2] = d;
		}
	}
	else
	{
		coeffs[0] = a;
		coeffs[1] = b;
		coeffs[2] = c;
		coeffs[3] = d;
	}

    RootFinder rf;
    Real realRoots[3];
    Real imagRoots[3];
 
    int count = rf.rpoly(coeffs, leadCoeff, realRoots, imagRoots);
 
	for (int i=0; i<count; ++i)
	{
		if (!imagRoots[i])
            times.push_back(realRoots[i]);
	}

    std::sort(times.begin(), times.end());
*/

    const double convergence_tol = tol*(std::fabs(A)+std::fabs(B)+std::fabs(C)+std::fabs(D));

    // find intervals to check, or just solve it if it reduces to a quadratic =============================
    std::vector<double> interval_times;
    double discriminant=B*B-3*A*C; // of derivative of cubic, 3*A*t^2+2*B*t+C, divided by 4 for convenience
    if(discriminant<=0){ // monotone cubic: only one root in [0,1] possible
        // so we just 
        interval_times.push_back(0);
        interval_times.push_back(1);
    }
    else
    { // positive discriminant, B!=0
        if(A==0)
        { // the cubic is just a quadratic, B*t^2+C*t+D ========================================
            discriminant=C*C-4*B*D; // of the quadratic
            if(discriminant<=0)
            {
                double t=-C/(2*B);
                if(t>=-tol && t<=1+tol)
                {
                    t=clamp(t, 0., 1.);
                    Real val = std::fabs(signed_volume((1-t)*x0+t*xnew0,
                                                       (1-t)*x1+t*xnew1,
                                                       (1-t)*x2+t*xnew2,
                                                       (1-t)*x3+t*xnew3));
                    if (val < convergence_tol)
                    {
                        times.push_back(t);
                    }
                }
            }
            else
            { // two separate real roots
                double t0, t1;
                if(C>0) t0=(-C-std::sqrt(discriminant))/(2*B);
                else    t0=(-C+std::sqrt(discriminant))/(2*B);
                t1=D/(B*t0);
                if(t1<t0) std::swap(t0,t1);
                if(t0>=-tol && t0<=1+tol)
                {
                    times.push_back(clamp(t0, 0., 1.));
                }
                if(t1>=-tol && t1<=1+tol)
                {
                    addUnique(times, clamp(t1, 0., 1.));
                }
            }

            for (size_t i=0; i<times.size(); ++i)
            {
                Real ti = times[i];
                Real val = std::fabs(signed_volume((1-ti)*x0+ti*xnew0,
                            (1-ti)*x1+ti*xnew1,
                            (1-ti)*x2+ti*xnew2,
                            (1-ti)*x3+ti*xnew3));
                errors.push_back(val);
            }

            return;
        }
        else
        { // cubic is not monotone: divide up [0,1] accordingly =====================================
            double t0, t1;
            if(B>0)
                t0=(-B-std::sqrt(discriminant))/(3*A);
            else
                t0=(-B+std::sqrt(discriminant))/(3*A);
            t1=C/(3*A*t0);
            if(t1<t0)
                std::swap(t0,t1);
            interval_times.push_back(0);
            if(t0>0 && t0<1)
                interval_times.push_back(t0);
            if(t1>0 && t1<1)
                interval_times.push_back(t1);

            interval_times.push_back(1);
        }
    }

    // look for roots in indicated intervals ==============================================================
    // evaluate coplanarity more accurately at each endpoint of the intervals
    std::vector<double> interval_values(interval_times.size());
    for(unsigned int i=0; i<interval_times.size(); ++i){
        double t=interval_times[i];
        interval_values[i]=signed_volume((1-t)*x0+t*xnew0, (1-t)*x1+t*xnew1, (1-t)*x2+t*xnew2, (1-t)*x3+t*xnew3);
    }
    // first look for interval endpoints that are close enough to zero, without a sign change
    for(unsigned int i=0; i<interval_times.size(); ++i){
        if(interval_values[i]==0)
        {
            times.push_back(interval_times[i]);
        }
        else if(std::fabs(interval_values[i])<convergence_tol)
        {
            if((i==0 || (interval_values[i-1]>=0 && interval_values[i]>=0) ||
                        (interval_values[i-1]<=0 && interval_values[i]<=0)) &&
               (i==interval_times.size()-1 || (interval_values[i+1]>=0 && interval_values[i]>=0) ||
                                              (interval_values[i+1]<=0 && interval_values[i]<=0)))
            {
                times.push_back(interval_times[i]);
            }
        }
    }
    // and then search in intervals with a sign change
    for(unsigned int i=1; i<interval_times.size(); ++i)
    {
        double tlo=interval_times[i-1], thi=interval_times[i], tmid;
        double vlo=interval_values[i-1], vhi=interval_values[i], vmid;
        if((vlo<0 && vhi>0) || (vlo>0 && vhi<0)){
            // start off with secant approximation (in case the cubic is actually linear)
            double alpha=vhi/(vhi-vlo);
            tmid=alpha*tlo+(1-alpha)*thi;
            for(int iteration=0; iteration<50; ++iteration){
                vmid=signed_volume((1-tmid)*x0+tmid*xnew0, (1-tmid)*x1+tmid*xnew1,
                        (1-tmid)*x2+tmid*xnew2, (1-tmid)*x3+tmid*xnew3);
                if(std::fabs(vmid)<1e-2*convergence_tol) break;
                if((vlo<0 && vmid>0) || (vlo>0 && vmid<0)){ // if sign change between lo and mid
                    thi=tmid;
                    vhi=vmid;
                }else{ // otherwise sign change between hi and mid
                    tlo=tmid;
                    vlo=vmid;
                }
                if(iteration%2) alpha=0.5; // sometimes go with bisection to guarantee we make progress
                else alpha=vhi/(vhi-vlo); // other times go with secant to hopefully get there fast
                tmid=alpha*tlo+(1-alpha)*thi;
            }
            times.push_back(tmid);
        }
    }
    std::sort(times.begin(), times.end());

    for (size_t i=0; i<times.size(); ++i)
    {
        Real ti = times[i];
        Real val = std::fabs(signed_volume((1-ti)*x0+ti*xnew0,
                                           (1-ti)*x1+ti*xnew1,
                                           (1-ti)*x2+ti*xnew2,
                                           (1-ti)*x3+ti*xnew3));
        errors.push_back(val);
    }
}

Real CandidateCollision::signed_volume(const Vec3d &x0, const Vec3d &x1,
                                         const Vec3d &x2, const Vec3d &x3) const
{
    // Equivalent to triple(x1-x0, x2-x0, x3-x0), six times the signed volume of the tetrahedron.
    // But, for robustness, we want the result (up to sign) to be independent of the ordering.
    // And want it as accurate as possible...
    // But all that stuff is hard, so let's just use the common assumption that all coordinates are >0,
    // and do something reasonably accurate in fp.

    // This formula does almost four times too much multiplication, but if the coordinates are non-negative
    // it suffers in a minimal way from cancellation error.
    return ( x0[0]*(x1[1]*x3[2]+x3[1]*x2[2]+x2[1]*x1[2])
            +x1[0]*(x2[1]*x3[2]+x3[1]*x0[2]+x0[1]*x2[2])
            +x2[0]*(x3[1]*x1[2]+x1[1]*x0[2]+x0[1]*x3[2])
            +x3[0]*(x1[1]*x2[2]+x2[1]*x0[2]+x0[1]*x1[2]) )

        - ( x0[0]*(x2[1]*x3[2]+x3[1]*x1[2]+x1[1]*x2[2])
                +x1[0]*(x3[1]*x2[2]+x2[1]*x0[2]+x0[1]*x3[2])
                +x2[0]*(x1[1]*x3[2]+x3[1]*x0[2]+x0[1]*x1[2])
                +x3[0]*(x2[1]*x1[2]+x1[1]*x0[2]+x0[1]*x2[2]) );
}

}

