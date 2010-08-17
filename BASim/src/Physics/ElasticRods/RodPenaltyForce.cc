// RodPenaltyForce.cc
//

//#include "BASim/src/Collisions/CollisionMeshData.hh"
#include "RodPenaltyForce.hh"

namespace BASim {

#define COLLISION_EPSILON 1e-6



RodPenaltyForce::RodPenaltyForce()
{
}

RodPenaltyForce::~RodPenaltyForce()
{
}
/*
  void RodPenaltyForce::computeEnergy(Scalar& e)
  {
  for (VertexObjectMapIterator voItr=_vertexObjects.begin(); voItr!=_vertexObjects.end(); ++voItr)
  {
  CollisionMeshData *cmData = voItr->second.first;

  int vertex   = voItr->first;
  int triangle = voItr->second.second;

  // Get distance between vertex and triangle
  //
  Scalar t1, t2, t3;
  Scalar distance = std::sqrt(getClosestPointsVertexTriangle(_rod.vertex(vertex),
  cmData->prevPositions[cmData->triangleIndices[(3 * triangle)    ]],
  cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 1]],
  cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 2]],
  t1, t2, t3));

  if (distance < (cmData->getThickness() + _rod.Radius()))
  {
  Scalar stiffness = cmData->getSeparationStrength();
  Scalar separation = (cmData->getThickness() + _rod.Radius()) - distance;
  e += 0.5 * stiffness * separation * separation;
  }
  }

  for (EdgeRodMapIterator erItr=_edgeRods.begin(); erItr!=_edgeRods.end(); ++erItr)
  {
  Rod *otherRod = erItr->second.first;

  int edge1 = erItr->first;
  int edge2 = erItr->second.second;

  // Get distance between two edges
  //
  Scalar s1, s2;
  Scalar distance = std::sqrt(getClosestPointsEdgeEdge(_rod.vertex(edge1),
  _rod.vertex((edge1+1)%_rod.nv()),
  otherRod->vertex(edge2),
  otherRod->vertex((edge2+1)%otherRod->nv()),
  s1, s2));

  if (distance < (_rod.Radius() + otherRod->Radius()))
  {
  Scalar stiffness = _rod.getSeparationStrength();
  Scalar separation = (_rod.Radius() + otherRod->Radius()) - distance;

  e += 0.5 * stiffness * separation * separation;
  }
  }
  }
*/


void RodPenaltyForce::computeForceDX(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J)
{
  MatXd localJ(3, 3);
  IntArray indices(3);
  int nv = rod.nv();
  Scalar r = rod.radius();
  
  for (VertexObjectMapIterator voItr=_vertexObjects.begin(); voItr!=_vertexObjects.end(); ++voItr)
  {
    CollisionMeshData *cmData = voItr->second.first;

    int vertex   = voItr->first;
    int triangle = voItr->second.second;

//  	Scalar thickness = cmData->getThickness() + rod.radius();

    Vec3d n = (cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 1]] -
               cmData->prevPositions[cmData->triangleIndices[(3 * triangle)    ]]).cross(
               cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 2]] - 
               cmData->prevPositions[cmData->triangleIndices[(3 * triangle)    ]]);

    n.normalize();
    
    Scalar stiffness = cmData->getSeparationStrength();

    localJ.setZero();
    localJacobian(localJ, stiffness, n);
  
	  for (int i = 0; i < 3; ++i) {
	    indices[i] = baseindex + rod.vertIdx(vertex,i);
	  }
	  localJ *= scale;
	  J.add(indices, indices, localJ);
	    
  }
}



void RodPenaltyForce::localJacobian(MatXd& J, const Scalar stiffness, const Vec3d& normal)
{
	Mat3d M = -stiffness * outerProd(normal, normal);

	for (int j = 0; j < 3; ++j) {
		for (int k = 0; k < 3; ++k) {
			J(j,k) += M(j,k);
		}
	}

}


void RodPenaltyForce::clearPenaltyForces() {

	_vertexObjects.clear();


}


void RodPenaltyForce::computeForce(const ElasticRod& const_rod, VecXd& F)
{
  VecXd beforeF = F;

//  cerr << "Forces (BEFORE) = \n " << F << endl;

  ElasticRod& rod = const_cast<ElasticRod&>(const_rod);
  // Record pairs whose distance is greatern than the influence of the
  // force we will delete these later
  //
  std::vector<std::pair<int, std::pair<CollisionMeshData *, int> > > toDelete;
  for (VertexObjectMapIterator voItr=_vertexObjects.begin(); voItr!=_vertexObjects.end(); ++voItr)
  {
    CollisionMeshData *cmData = voItr->second.first;

    int vertex   = voItr->first;
    int triangle = voItr->second.second;

//    std::cout << "pf " << vertex << "\n";
    
    // Get distance between vertex and triangle
    //
    Scalar t1, t2, t3;
    Scalar distance = std::sqrt(getClosestPointsVertexTriangle(rod.getVertex(vertex),
                                                               cmData->prevPositions[cmData->triangleIndices[(3 * triangle)    ]],
                                                               cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 1]],
                                                               cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 2]],
                                                               t1, t2, t3));
    /*Scalar distance = std::sqrt(getClosestPointsVertexTriangle(rod.getVertex(vertex),
                                                               cmData->prevPositions[cmData->triangleIndices[(3 * triangle)    ]],
                                                               cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 1]],
                                                               cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 2]],
                                                               t1, t2, t3));*/

    //Scalar distance = (v0 - vertex_face_collisions[i].cp).dot(vertex_face_collisions[i].n) - (vertex_face_collisions[i].r0 + vertex_face_collisions[i].r1 + vertex_face_collisions[i].h);

    //Vec3d force = -vertex_face_collisions[i].k * distance * vertex_face_collisions[i].n;

    
    if (1 || distance < (cmData->getThickness() + rod.radius()))
    {
    	Scalar thickness = cmData->getThickness() + rod.radius();

      Vec3d n = (cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 1]] -
                 cmData->prevPositions[cmData->triangleIndices[(3 * triangle)    ]]).cross(
                 cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 2]] - 
                 cmData->prevPositions[cmData->triangleIndices[(3 * triangle)    ]]);

      Vec3d normal = rod.getVertex(vertex) -
        (t1 * cmData->prevPositions[cmData->triangleIndices[(3 * triangle)    ]] +
         t2 * cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 1]] +
         t3 * cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 2]]);
				
	//					if (m_damping > 0) {
	//						mag += -m_damping * (v.dot(normal));
	//					}
				
//						Vec3d f = mag * normal;

//						for (int j = 0; j < 3; ++j) {
//							force(j) += f(j);
//						}					

      // Vertex is inside object or the distance is too small to trust the normal
      //
//      if (n.dot(normal) < 0.0 || distance < 1e-6)
      //if (distance < 1e-6)   continue;

//      normal.normalize();

//      Vec3d relVel = rod.getVelocity(vertex) -
//        (t1 * cmData->velocities[cmData->triangleIndices[(3 * triangle)    ]] +
//         t2 * cmData->velocities[cmData->triangleIndices[(3 * triangle) + 1]] +
//         t3 * cmData->velocities[cmData->triangleIndices[(3 * triangle) + 2]]);

      // TODO: Should the COR automatically be set so that the separating velocity
      // is exactly zero when the distance is equal to thickness + radius?
      //
      Scalar e = 1.0;
//      if (!(relVel.dot(normal) < 0.0))
//        e = cmData->getCoefficientOfRestitution();

      Scalar stiffness = cmData->getSeparationStrength();

			// surface normal			
			n.normalize();
			
//      Vec3d force = e * stiffness * ((cmData->getThickness() + rod.radius()) - distance) * normal;
			Vec3d force = -e * stiffness * (n.dot(normal) - thickness) * n;

//      std::cout << vertex << " VERT " << force << "\n";

      // TODO: It is very hard to achieve static friction with a friction force,
      // in fact, I don't know of a way to do it (maybe somebody has figured it out)
      // As such, keep the friction coeff low, otherwise you may get weird behavior
      // Try using a velocity based damping for friction, it may work better
      //
      /*
      double frc = rod.getFrictionCoefficient();
      Scalar friction = std::max(frc, cmData->getFrictionCoefficient());
      if (friction > 0.0)
      {
        Vec3d tangent = relVel - (relVel.dot(normal)) * normal;
        Scalar l = tangent.norm();
        if (l > 1e-6)
        {
          tangent /= l;

          force -= force.norm() * friction * tangent;
        }
      }
      */

      for (int i=0; i<3; ++i)
        F[rod.vertIdx(vertex, i)] += force[i];

      //cerr << "Collision, applying force to vertex " << vertex << ": " << force << endl;
      
    }
    else
    {
    //  cerr << "Scheduling force on vertex " << vertex << " for deletion\n";
      toDelete.push_back(std::make_pair(vertex, std::make_pair(cmData, triangle)));
    }
  }

  /*cerr << "about to delete forces, num forces = " << _vertexObjects.size() << endl;
  cerr << "toDelete.size() = " << toDelete.size() << endl;
  */
  // Kill all pairs flagged for deletion
  //
  
  std::vector<std::pair<int, std::pair<CollisionMeshData *, int> > >::iterator itr=toDelete.begin();
  for( ; itr!=toDelete.end(); ++itr)
  {
    std::pair<VertexObjectMapIterator, VertexObjectMapIterator> p=_vertexObjects.equal_range(itr->first);
    for (VertexObjectMapIterator i=p.first; i!=p.second; ++i)
    {
      if (i->second.first == itr->second.first && i->second.second == itr->second.second)
      {
        _vertexObjects.erase(i);
        break;
      }
    }
  }

  return;
  
  /*cerr << "done deleting: " << endl;
  cerr << "num forces = " << _vertexObjects.size() << endl;
  cerr << "toDelete.size() = " << toDelete.size() << endl;
  */
  
  

  // Now compute all edge-edge forces
  //
  std::vector<std::pair<int, std::pair<ElasticRod *, int> > > rodsToDelete;
  for (EdgeRodMapIterator erItr=_edgeRods.begin(); erItr!=_edgeRods.end(); ++erItr)
  {
    ElasticRod *otherRod = erItr->second.first;

    int edge1 = erItr->first;
    int edge2 = erItr->second.second;

    // Get distance between two edges
    //
    Scalar s1, s2;
    Scalar distance = std::sqrt(getClosestPointsEdgeEdge(rod.getVertex(edge1),
                                                         rod.getVertex((edge1+1)%rod.nv()),
                                                         otherRod->getVertex(edge2),
                                                         otherRod->getVertex((edge2+1)%otherRod->nv()),
                                                         s1, s2));

    if (distance < (rod.radius() + otherRod->radius()))
    {
      if (distance < 1e-6)
        continue;

      Vec3d normal = ((1.0 - s1) * rod.getVertex(edge1) + s1 * rod.getVertex((edge1+1)%rod.nv())) -
        ((1.0 - s2) * otherRod->getVertex(edge2) + s2 * otherRod->getVertex((edge1+1)%otherRod->nv()));
      normal.normalize();

      Vec3d relVel = ((1.0 - s1) * rod.getVelocity(edge1) + s1 * rod.getVelocity((edge1+1)%rod.nv())) -
        ((1.0 - s2) * otherRod->getVelocity(edge2) + s2 * otherRod->getVelocity((edge1+1)%otherRod->nv()));

      Scalar e = 1.0;
      if (!(relVel.dot(normal) < 0.0))
        e = const_cast<ElasticRod&>(rod).getCoefficientOfRestitution();

      Scalar stiffness = rod.getSeparationStrength();
      Vec3d force = e * stiffness * ((rod.radius() + otherRod->radius()) - distance) * normal;

      Scalar friction = std::max(rod.getFrictionCoefficient(), otherRod->getFrictionCoefficient());
      if (friction > 0.0)
      {
        Vec3d tangent = relVel - relVel.dot(normal) * normal;
        Scalar l = tangent.norm();
        if (l > 1e-6)
        {
          tangent /= l;

          // TODO: Friction should be capped by how much relative motion
          // there; this is hard to do for a frictional force
          //
          force -= force.norm() * friction * tangent;
        }
      }

      for (int i=0; i<3; ++i)
      {
        F[rod.vertIdx(edge1, i)] += (1.0 - s1) * force[i];
        F[rod.vertIdx((edge1+1)%rod.nv(), i)] += (      s1) * force[i];
      }
    }
    else
      rodsToDelete.push_back(std::make_pair(edge1, std::make_pair(otherRod, edge2)));
  }

  for (std::vector<std::pair<int, std::pair<ElasticRod *, int> > >::iterator itr=rodsToDelete.begin();
       itr!=rodsToDelete.end(); ++itr)
  {
    std::pair<EdgeRodMapIterator, EdgeRodMapIterator> p=_edgeRods.equal_range(itr->first);
    for (EdgeRodMapIterator i=p.first; i!=p.second; ++i)
    {
      if (i->second.first == itr->second.first && i->second.second == itr->second.second)
      {
        _edgeRods.erase(i);
        break;
      }
    }
  }


// Apply clumping force

  if (clumping_enbld) {
	std::vector<std::pair<int, std::pair<ElasticRod *, int> > > rodsToDelete;
	for (VertexRodMapIterator vrItr=_clumpingVerts.begin(); vrItr!=_clumpingVerts.end(); ++vrItr)
	{
	  ElasticRod *otherRod = vrItr->second.first;

	  int v1 = vrItr->first;
	  int v2 = vrItr->second.second;
			
		if (rod.vertFixed(v1) || otherRod->vertFixed(v2)) {
	       rodsToDelete.push_back(std::make_pair(v1, std::make_pair(otherRod, v2)));
		   continue;
		}
			
		Vec3d d = rod.getVertex(v1) - otherRod->getVertex(v2);
		Scalar distance = d.norm();
			
		if (true || distance > (rod.radius() + otherRod->radius()) + 1e-4) {
           Vec3d force = - clumping_coeff * d / (distance * distance);
				
		   F.segment(rod.vertIdx(v1, 0), 3) += force;			
		
		} else {
		   rodsToDelete.push_back(std::make_pair(v1, std::make_pair(otherRod, v2)));
		   continue;
		}
    }

	for (std::vector<std::pair<int, std::pair<ElasticRod *, int> > >::iterator itr=rodsToDelete.begin();
	     itr!=rodsToDelete.end(); ++itr)
	{
	  std::pair<VertexRodMapIterator, VertexRodMapIterator> p=_clumpingVerts.equal_range(itr->first);
	  for (VertexRodMapIterator i=p.first; i!=p.second; ++i)
	  {
	    if (i->second.first == itr->second.first && i->second.second == itr->second.second)
	    {
	      _clumpingVerts.erase(i);
	      break;
	    }
	  }
	}
  }

  // cerr << "Forces (AFTER) = \n " << F - beforeF << endl;
}

void RodPenaltyForce::addRodPenaltyForce(int vertex, CollisionMeshData *cmData, int triangle)
{
  // Adding a vertex-triangle penalty force, where the triangle is part of a collision
  // object mesh. First make sure it's not already in our list
  //
  std::pair<VertexObjectMapIterator, VertexObjectMapIterator> p = _vertexObjects.equal_range(vertex);

  for (VertexObjectMapIterator i=p.first; i!=p.second; ++i)
  {
    if (i->second.first == cmData && i->second.second == triangle)
      return;
  }

  _vertexObjects.insert(VertexObjectMap::value_type(vertex, std::make_pair(cmData, triangle)));
}

void RodPenaltyForce::addRodPenaltyForce(int edge, ElasticRod *rod, int otherEdge)
{
  // Adding a penalty force between two edges of two rods. Make sure its
  // not already on the list
  //
  std::pair<EdgeRodMapIterator, EdgeRodMapIterator> p = _edgeRods.equal_range(edge);

  for (EdgeRodMapIterator i=p.first; i!=p.second; ++i)
  {
    if (i->second.first == rod && i->second.second == otherEdge)
      return;
  }

  _edgeRods.insert(EdgeRodMap::value_type(edge, std::make_pair(rod, otherEdge)));
}


// Add a force to make clumps, for possible pairs of vertices 
void RodPenaltyForce::addRodClumpingForce(int vertex, ElasticRod *rod, int otherVertex)
{
  if(clumping_enbld) {
    std::pair<VertexRodMapIterator, VertexRodMapIterator> p = _clumpingVerts.equal_range(vertex);

	for (VertexRodMapIterator i=p.first; i!=p.second; ++i)
	{
	  if (i->second.first == rod && i->second.second == otherVertex)
        return;
	}

	_clumpingVerts.insert(VertexRodMap::value_type(vertex, std::make_pair(rod, otherVertex)));
  }
}

// Helper routines for getting closest distance between vertex-triangle pairs and
// edge-edge pairs
//
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
}

Scalar RodPenaltyForce::getClosestPointsEdgeEdge(const Vec3d& e11, const Vec3d& e12,
                                                 const Vec3d& e21, const Vec3d& e22,
                                                 Scalar &s, Scalar &t) const
{
  const double *p1 = e11.data();
  const double *q1 = e12.data();
  const double *p2 = e21.data();
  const double *q2 = e22.data();

  Scalar d1[3], d2[3], r[3], a, e, f;
  Scalar c1[3], c2[3];

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
    Scalar c = d1[0]*r[0] + d1[1]*r[1] + d1[2]*r[2];

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
      Scalar b = d1[0]*d2[0] + d1[1]*d2[1] + d1[2]*d2[2];
      Scalar denom = a*e - b*b;

      if (denom != 0.0f)
      {
        s = (b*f - c*e) / denom;
        if (s<0.0f) s = 0.0f;
        if (s>1.0f) s = 1.0f;
      }
      else
        s = 0.0f;

      Scalar tnom = b*s + f;
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

}
