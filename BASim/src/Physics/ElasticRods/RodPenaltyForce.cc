// RodPenaltyForce.cc
//

//#include "BASim/src/Collisions/CollisionMeshData.hh"
#include "RodPenaltyForce.hh"
#ifdef WETA
#include "../../Math/Math.hh"
#else
#include "BASim/src/Math/Math.hh"
#endif

namespace BASim
{

#define COLLISION_EPSILON 1e-6

RodPenaltyForce::RodPenaltyForce( double penaltyThicknessFraction ) : m_penaltyThicknessFraction(penaltyThicknessFraction)
{
}

RodPenaltyForce::~RodPenaltyForce()
{
}

void RodPenaltyForce::computeForceDX(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J)
{
    MatXd localJ(3, 3);
    IntArray indices(3);
    int nv = rod.nv();
    Scalar r = rod.radius();

    for (int i = 0; i < (int) vertex_face_collisions.size(); i++)
    {
        int vertex = vidx[i];
        Vec3d v0 = rod.getVertex(vertex);

	double penaltyThickness =  vertex_face_collisions[i]->h / 10.0;

        Scalar distance = (v0 - vertex_face_collisions[i]->cp).dot(vertex_face_collisions[i]->m_normal)
                - (vertex_face_collisions[i]->r0 + vertex_face_collisions[i]->r1 + penaltyThickness);

	if (distance > 0) continue;

        localJ.setZero();
        localJacobian(localJ, vertex_face_collisions[i]->k, vertex_face_collisions[i]->m_normal);

        for (int i = 0; i < 3; ++i)
        {
            indices[i] = baseindex + rod.vertIdx(vertex, i);
        }
        localJ *= scale;
        J.add(indices, indices, localJ);
    }

}

void RodPenaltyForce::localJacobian(MatXd& J, const Scalar stiffness, const Vec3d& normal)
{
    Mat3d M = -stiffness * outerProd(normal, normal);

    for (int j = 0; j < 3; ++j)
    {
        for (int k = 0; k < 3; ++k)
        {
            J(j, k) += M(j, k);
        }
    }
}

void RodPenaltyForce::clearPenaltyForces()
{
    // std::cerr << "Clearing penalty forces" << std::endl;
    vidx.clear();
    vertex_face_collisions.clear();
}

void RodPenaltyForce::computeForce(const ElasticRod& rod, VecXd& F)
{
  //  VecXd beforeF = F;
  //  std::cout << "Forces (BEFORE) = \n " << F << std::endl;

    //  ElasticRod& rod = const_cast<ElasticRod&>(const_rod);

    for (int i = 0; i < (int) vertex_face_collisions.size(); i++)
    {
        int vertex = vidx[i];
        Vec3d v0 = rod.getVertex(vertex);

	double penaltyThickness =  vertex_face_collisions[i]->h / 10.0;

        Scalar distance = (v0 - vertex_face_collisions[i]->cp).dot(vertex_face_collisions[i]->m_normal)
                - (vertex_face_collisions[i]->r0 + vertex_face_collisions[i]->r1 + penaltyThickness);

	if (distance > 0) continue;

        Vec3d force = -vertex_face_collisions[i]->k * distance * vertex_face_collisions[i]->m_normal;

        for (int i = 0; i < 3; ++i)
        {
            F[rod.vertIdx(vertex, i)] += force[i];
        }

     //   std::cout << "Collision, applying force to vertex " << vertex << std::endl;
     //   std::cout << "Distance = " << distance << std::endl;
      //  std::cout << "Normal = " << vertex_face_collisions[i]->m_normal << std::endl;
     //   std::cout << "Force = " << force << std::endl;
     //   std::cout << "Stiffness = " << vertex_face_collisions[i]->k << std::endl;

    }
  //  std::cerr << "Penalty forces added = \n " << F - beforeF << std::endl;
}

void RodPenaltyForce::addRodPenaltyForce(int vertex, VertexFaceProximityCollision* vfpcol)
{
    //	std::cout << "pntl added " << vertex << " " << cllsn.n << " " << cllsn.cp << "\n";
    vidx.push_back(vertex);
    vertex_face_collisions.push_back(vfpcol);
}

// Helper routines for getting closest distance between vertex-triangle pairs and
// edge-edge pairs
//
/*
 Scalar RodPenaltyForce::getClosestPointsVertexTriangle(const Vec3d& v0, const Vec3d& v1,
 const Vec3d& v2, const Vec3d& v3,
 Scalar &t1, Scalar &t2, Scalar &t3) const

 {
 const double *p = v0.data();
 const double *a = v1.data();
 const double *b = v2.data();
 const double *c = v3.data();

 Scalar ab[3], ac[3], ap[3], bp[3];

 ab[0] = b[0] - a[0];
 ab[1] = b[1] - a[1];
 ab[2] = b[2] - a[2];

 ac[0] = c[0] - a[0];
 ac[1] = c[1] - a[1];
 ac[2] = c[2] - a[2];

 ap[0] = p[0] - a[0];
 ap[1] = p[1] - a[1];
 ap[2] = p[2] - a[2];

 Scalar d1 = ab[0]*ap[0] + ab[1]*ap[1] + ab[2]*ap[2];
 Scalar d2 = ac[0]*ap[0] + ac[1]*ap[1] + ac[2]*ap[2];

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

 Scalar d3 = ab[0]*bp[0] + ab[1]*bp[1] + ab[2]*bp[2];
 Scalar d4 = ac[0]*bp[0] + ac[1]*bp[1] + ac[2]*bp[2];

 if ((d3 >= 0.0f) && (d4 <= d3))
 {
 t1 = 0.0f;
 t2 = 1.0f;
 t3 = 0.0f;

 return ((p[0]-b[0])*(p[0]-b[0]) + (p[1]-b[1])*(p[1]-b[1]) + (p[2]-b[2])*(p[2]-b[2]));
 }

 Scalar vc = d1*d4 - d3*d2;

 if ((vc <= 0.0f) && (d1 >= 0.0f) && (d3 <= 0.0f))
 {
 Scalar v = d1 / (d1 - d3);

 t1 = 1-v;
 t2 = v;
 t3 = 0;

 Scalar vec[3];
 vec[0] = p[0] - (a[0]+v*ab[0]);
 vec[1] = p[1] - (a[1]+v*ab[1]);
 vec[2] = p[2] - (a[2]+v*ab[2]);

 return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
 }

 Scalar cp[3];
 cp[0] = p[0] - c[0];
 cp[1] = p[1] - c[1];
 cp[2] = p[2] - c[2];

 Scalar d5 = ab[0]*cp[0] + ab[1]*cp[1] + ab[2]*cp[2];
 Scalar d6 = ac[0]*cp[0] + ac[1]*cp[1] + ac[2]*cp[2];

 if ((d6 >= 0.0f) && (d5 <= d6))
 {
 t1 = 0;
 t2 = 0;
 t3 = 1;

 return ((p[0]-c[0])*(p[0]-c[0]) + (p[1]-c[1])*(p[1]-c[1]) + (p[2]-c[2])*(p[2]-c[2]));
 }

 Scalar vb = d5*d2 - d1*d6;

 if ((vb <= 0.0f) && (d2 >= 0.0f) && (d6 <= 0.0f))
 {
 Scalar w = d2 / (d2 - d6);

 t1 = 1-w;
 t2 = 0;
 t3 = w;

 Scalar vec[3];
 vec[0] = p[0] - (a[0]+w*ac[0]);
 vec[1] = p[1] - (a[1]+w*ac[1]);
 vec[2] = p[2] - (a[2]+w*ac[2]);

 return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
 }

 Scalar va = d3*d6 - d5*d4;

 if ((va <= 0.0f) && ((d4-d3) >= 0.0f) && ((d5-d6) >= 0.0f))
 {
 Scalar w = (d4 - d3) / ((d4 - d3) + (d5 - d6));

 t1 = 0;
 t2 = 1-w;
 t3 = w;

 Scalar vec[3];
 vec[0] = p[0] - (b[0]+w*(c[0]-b[0]));
 vec[1] = p[1] - (b[1]+w*(c[1]-b[1]));
 vec[2] = p[2] - (b[2]+w*(c[2]-b[2]));

 return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
 }

 Scalar denom = 1.0f / (va + vb + vc);
 Scalar v = vb * denom;
 Scalar w = vc * denom;
 Scalar u = 1.0 - v - w;

 t1 = u;
 t2 = v;
 t3 = w;

 Scalar vec[3];
 vec[0] = p[0] - (u*a[0] + v*b[0] + w*c[0]);
 vec[1] = p[1] - (u*a[1] + v*b[1] + w*c[1]);
 vec[2] = p[2] - (u*a[2] + v*b[2] + w*c[2]);

 return (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);
 } */

}
