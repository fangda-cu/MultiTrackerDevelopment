// RodPenaltyForce.cc
//

//#include "BASim/src/Collisions/CollisionMeshData.hh"

//#include <vector>

#include "RodPenaltyForce.hh"

namespace BASim {

  
#define COLLISION_EPSILON 1e-6

RodPenaltyForce::RodPenaltyForce()
{
}

RodPenaltyForce::~RodPenaltyForce()
{
}


VertexConstraintMapIter RodPenaltyForce::setVertexPositionPenalty(const ElasticRod* rod,
        int vertex_id, Vec3d& target_position, double stiffness, short type)
{
    Vec3d direction = rod->getVertex( vertex_id ) - target_position ;
    RodVertexConstraint newConstraint(target_position, stiffness, direction.norm(), type );
    
    return m_vertex_position_penalties.insert(VertexConstraintMap::value_type(
        std::make_pair(vertex_id, newConstraint) ));
}

// id vertex_id = -1, delete all
void RodPenaltyForce::clearVertexPositionPenalty(int vertex_id)
{
  if (vertex_id == -1) {
    m_vertex_position_penalties.clear();
  } else {
    std::pair<VertexConstraintMapIter, VertexConstraintMapIter> p = m_vertex_position_penalties.equal_range(vertex_id);

    for (VertexConstraintMapIter i=p.first; i!=p.second; ++i)
    {
      if (i->first == vertex_id) {
        m_vertex_position_penalties.erase(i);
      }
    }
  }
}


void RodPenaltyForce::computeForceDX(int baseindex, const ElasticRod& rod, Scalar scale, MatrixBase& J)
{
    
  // vertex - target position penalty
  for (VertexConstraintMapIter it=m_vertex_position_penalties.begin(); it!=m_vertex_position_penalties.end(); ++it)
  {
    int vertex    = it->first;
    Scalar stiff  = it->second.m_stiff;

    for (int i = 0; i < 3; ++i) {
      J.add(baseindex+rod.vertIdx(vertex, i), baseindex+rod.vertIdx(vertex, i), -stiff * scale);
    }
  }
  
    /*
  {
    IntArray indices(3);
    
    for(int i=0; i<(int)m_bulk_springs.size(); i++) {
      RodGroupManager::RodRodSpring spr = m_bulk_springs[i];
      
      Vec3d p(0,0,0);
      
      for (int j=spr.start_v; j<spr.start_v + spr.nv; j++) {
        p += rod.getVertex(j);
      }
    
      double m = (double) spr.nv;
      p /= m;
      
      Vec3d vecA = p - spr.opposite_end;
      
      Mat3d M = 2.0 * outerProd(vecA, vecA);
      double d = vecA.dot(vecA) - spr.rest_length * spr.rest_length;
      
      // d * identity matrix
      for (int j = 0; j < 3; ++j) {
        M(j,j) += d;
      }
      
      M *= - spr.stiff / (m*m);
      M *= scale;
      
      for (int j=spr.start_v; j<spr.start_v + spr.nv; j++) {
        for (int i = 0; i < 3; ++i) {
          indices[i] = baseindex + rod.vertIdx(j,i);
        }
        J.add(indices, indices, M);
      }
    }
  }*/
  
  
  MatXd localJ(3, 3);
  IntArray indices(3);
  
  uint id = 0;
  
  for (VertexObjectMapIterator voItr=_vertexObjects.begin(); voItr!=_vertexObjects.end(); ++voItr, ++id)
  {
    CollisionMeshData *cmData = voItr->second.first;

    int vertex   = voItr->first;
    int triangle = voItr->second.second;
    
    Vec3d n = surface_normals[id];

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
  surface_normals.clear();
}


void RodPenaltyForce::computeForce(const ElasticRod& const_rod, VecXd& F)
{
  ElasticRod& rod = const_cast<ElasticRod&>(const_rod);

  //find distinct vertex ids attached to constraints
  std::vector< int > vertexIds;
  for (VertexConstraintMapIter it=m_vertex_position_penalties.begin(),
          itEnd = m_vertex_position_penalties.end(); it!=itEnd; ++it)
  {
      if( std::find( vertexIds.begin(), vertexIds.end(), it->first )
              == vertexIds.end() )
      {
          vertexIds.push_back( it->first );
      }
  }

  for( std::vector< int >::iterator it = vertexIds.begin(),
          endIt = vertexIds.end(); it != endIt; ++it )
  {

      int vertex_id = (*it);
      std::pair<VertexConstraintMapIter, VertexConstraintMapIter> p
              = m_vertex_position_penalties.equal_range( vertex_id );

      Vec3d vertex = rod.getVertex( vertex_id );

      for (VertexConstraintMapIter i=p.first; i!=p.second; ++i)
      {
          if (i->first == vertex_id) 
          {
              Vec3d target = i->second.m_target;
              double stiff = i->second.m_stiff;
              double distance = i->second.m_restDistance;

              Vec3d n = vertex - target;

//              if( distance < 0 )
//              {
//                  i->second.m_restDistance = n.norm();
//
//              }

              if( i->second.m_type == RodVertexConstraint::kRest )
              {
                  double norm = n.norm();
                  Vec3d normalized = n * ( 1. / norm );
                  Vec3d newTarget = target + distance * normalized;
                  n = vertex - newTarget;
              }
              
              
              Vec3d force = - stiff * n;

              // add force together
              for (int i=0; i<3; ++i)
              {
                  F[rod.vertIdx(vertex_id, i)] += force[i];
              }
          }
      }
  }
  /*
  //std::vector<
  // vertex - target position penalty
  for (VertexPositionMapIterator it=m_vertex_position_penalties.begin(); it!=m_vertex_position_penalties.end(); ++it)
  {
    int vertex    = it->first;
    Vec3d& target = it->second.first;
    Scalar stiff  = it->second.second;
    
    Vec3d n = rod.getVertex(vertex) - target;
    
    Vec3d force = - stiff * n;
    
//    std::cout << n << "\n" << rod.getVertex(vertex) << "\n" << target <<"\n";

    for (int i=0; i<3; ++i) {
      F[rod.vertIdx(vertex, i)] += force[i];
    }
  }
  */
  
  /*
  for(int i=0; i<(int)m_bulk_springs.size(); i++) {
    RodGroupManager::RodRodSpring spr = m_bulk_springs[i];
    
    Vec3d p(0,0,0);
    
    for (int j=spr.start_v; j<spr.start_v + spr.nv; j++) {
      p += rod.getVertex(j);
    }
  
    double m = (double) spr.nv;
    p /= m;
    
    Vec3d vecA = p - spr.opposite_end;
    Vec3d force = - spr.stiff / m * (vecA.dot(vecA) - spr.rest_length * spr.rest_length) * vecA;
    
    for (int j=spr.start_v; j<spr.start_v + spr.nv; j++) {
      for (int k=0; k<3; ++k) {
        F[rod.vertIdx(j, k)] += force[k];
      }
    }
  }*/
  
  // Record pairs whose distance is greatern than the influence of the
  // force we will delete these later
  // => Don't delete (even though the distance is bigger than the influence) cause we need to maintain consistency of equations during Newton solver
  //
  //std::vector<std::pair<int, std::pair<CollisionMeshData *, int> > > toDelete;
  uint id = 0;
  for (VertexObjectMapIterator voItr=_vertexObjects.begin(); voItr!=_vertexObjects.end(); ++voItr, ++id)
  {
    CollisionMeshData *cmData = voItr->second.first;

    int vertex   = voItr->first;
    int triangle = voItr->second.second;
    
//    if (1 || distance < (cmData->getThickness() + rod.radius()))
    {
      Scalar thickness = cmData->getThickness() + rod.radius() * rod.getRadiusScale();
      Scalar stiffness = cmData->getSeparationStrength();

      Vec3d n = surface_normals[id];

      Scalar nnormal = (rod.getVertex(vertex) - cmData->prevPositions[cmData->triangleIndices[(3 * triangle)    ]]).dot(n); 
      Vec3d force = - stiffness * (nnormal - thickness) * n;

      for (int i=0; i<3; ++i)
        F[rod.vertIdx(vertex, i)] += force[i];

    }
//    else
//    {
      // we don't delete
      //toDelete.push_back(std::make_pair(vertex, std::make_pair(cmData, triangle)));
//    }
  }

  /*cerr << "about to delete forces, num forces = " << _vertexObjects.size() << endl;
  cerr << "toDelete.size() = " << toDelete.size() << endl;
  */
  // Kill all pairs flagged for deletion
  //
  
/*  std::vector<std::pair<int, std::pair<CollisionMeshData *, int> > >::iterator itr=toDelete.begin();
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
  }*/

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
  
  Vec3d n = (cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 1]] -
      cmData->prevPositions[cmData->triangleIndices[(3 * triangle)    ]]).cross(
      cmData->prevPositions[cmData->triangleIndices[(3 * triangle) + 2]] - 
      cmData->prevPositions[cmData->triangleIndices[(3 * triangle)    ]]);
  n.normalize();
  
  surface_normals.push_back(n);

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
}*/

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



// BULKY SPRING STUFFS

/*
RodGroupManager::RodGroupManager( std::vector<ElasticRod*> & rods, int ngs, double gap, double stf ) {

  m_rods = rods;
  nr = (int) m_rods.size();
  ng = ngs;

  layer_gap = gap;
  stiffness = stf;

  Setup();

}

int RodGroupManager::ClusterRods( int K ) {
  std::vector<Vec3d> group_mean;
  const int max_it = 20;
  m_rod_groups.clear();
  m_group_nrs.clear();

  if (K > nr) K = nr;

  group_mean.clear();
  m_rod_groups.resize(nr);
  m_group_nrs.resize(K, 1);
  m_rod_locs.resize(nr);

  for(int i=0; i<nr; i++) {
    Vec3d v = (m_rods[i]->getVertex(0) + m_rods[i]->getVertex(1)) * 0.5;
    m_rod_locs[i] = v;
  }
    
  for(int i=0; i<K; i++) {
    m_rod_groups[i] = i;
    group_mean.push_back( m_rod_locs[i] );
  }

  if (K == nr) return K;

  int it = 0;
  while (it < max_it) {
    it++;
    
//    for(int i=0; i<nr; i++) {
//      std::cout << m_rod_groups[i] << " ";
//    }  
//    cout << "\n";
    
    for(int i=0; i<nr; i++) {
      int best_group = -1;
      double min_d = 0;
        
      for( int j=0; j<K; j++) {
        double d = (m_rod_locs[i] - group_mean[j]).norm();
        if (best_group == -1 || d < min_d) {
          best_group = j;
          min_d = d;
        }
      }
    
      m_rod_groups[i] = best_group;
    }
  
    for (int i=0; i<K; i++) {
      group_mean[i] = Vec3d(0,0,0);
      m_group_nrs[i] = 0;
    }
  
  
    for(int i=0; i<nr; i++) {
      Vec3d v = m_rod_locs[i];
      group_mean[ m_rod_groups[i] ] += v;
      m_group_nrs[ m_rod_groups[i] ]++;
    }

    for (int i=0; i<K; i++) {
      group_mean[i] /= (double) m_group_nrs[i];
    }
  
  }
  
  return K;

}

void RodGroupManager::ActivateSpring(double dt) {
  for(int i=0; i<m_rods.size(); i++) {
    m_rods[i]->getPenaltyForce()->clearBulkSprings();
  }
  
  int nacts=0;
  
  for(int i=0; i<m_springs.size(); i++) {
    RodRodSpring spr = m_springs[i];
    
    if (spr.nv == 0) continue;
    
    // compute average point 
    Vec3d p1(0,0,0);
    Vec3d p2(0,0,0);
    Vec3d v1(0,0,0);
    Vec3d v2(0,0,0);
    
    for (int j=spr.start_v; j<spr.start_v + spr.nv; j++) {
      p1 += spr.rod1->getVertex(j);
      p2 += spr.rod2->getVertex(j);
      v1 += spr.rod1->getVelocity(j);
      v2 += spr.rod2->getVelocity(j);
    }
  
    p1 /= (double) spr.nv;
    p2 /= (double) spr.nv;
    v1 /= (double) spr.nv;
    v2 /= (double) spr.nv;
    
    double length_t0 = (p1-p2).norm();
    double length_t1 = (p1+v1*dt - p2+v2*dt).norm();
  
//    std::cout << spr.start_v << " " << spr.nv << " " << length_t0 << " " << spr.rest_length <<"\n";
    
    // if the spring is in compressive state
    if (length_t0 < spr.rest_length) {
  //    if (length_t0 >= length_t1) { // do we need this?
//        spr.opposite_end = p2 + v2 * dt; 
      spr.opposite_end = p2;
      spr.rod1->getPenaltyForce()->activateBulkSpring(spr);
        
//        spr.opposite_end = p1 + v1 * dt; 
//        spr.rod2->getPenaltyForce()->activateBulkSpring(spr);
        
      nacts++;
    }
    //}
    
  }
  
  std::cout << nacts << " springs are activated\n";
}


void RodGroupManager::SetupSpring() {
  double r1 = layer_gap;
  int n_spr_per_rod = 3; // from one rod
    
  for(int i=0; i<ng; i++) {
    std::vector<int> group;
  
    int n = 0;
  
    for(int j=0; j<nr; j++) {
      if (m_rod_groups[j] == i) {
        group.push_back(j);
        n++;
      }
    }
  
    if (n == 1) continue;
  
    std::vector<int> rod_nspr(n, 0);
    std::vector< std::pair<int, int> > spr_list;
  
    for(int a=0; a<n; a++) {
      std::vector<RodRodSpring> springs;
    
      for(int b=0; b<n; b++) {
        if (a == b) continue;
    
        RodRodSpring spr;
          
        spr.g = i;
        spr.id1 = a;
        spr.id2 = b;
        spr.rod1 = m_rods[group[a]];
        spr.rod2 = m_rods[group[b]];
        spr.start_v = 0;
        spr.nv = 2;
        spr.rest_length = (m_rod_locs[group[a]] - m_rod_locs[group[b]]).norm();
        spr.stiff = stiffness;
        springs.push_back(spr);
      }
    
      sort(springs.begin(), springs.end(), RodGroupManager::RodRodSpringCompare());
    
      int jend = n_spr_per_rod;
      for(int j=0; j<jend && j<springs.size(); j++) {
        pair <int, int> rp = make_pair (springs[j].id1, springs[j].id2);
        pair <int, int> rp2 = make_pair (springs[j].id2, springs[j].id1);
        if (std::find(spr_list.begin(), spr_list.end(), rp)==spr_list.end() && std::find(spr_list.begin(), spr_list.end(), rp2)==spr_list.end()) {
          rod_nspr[springs[j].id1]++;
          rod_nspr[springs[j].id2]++;
          spr_list.push_back(make_pair (springs[j].id1, springs[j].id2));
          m_springs.push_back(springs[j]);
        } else {
          jend++;
        }
      }
    }
  }
}


void RodGroupManager::PrecomputeSpring() {

  const uint nv_seg = 3;
  int ns = (int)m_springs.size();
      
  for(int i=0; i<ns; i++) {
    RodRodSpring spr = m_springs[i];
  
    int end_v = min(spr.rod1->nv() - 1, spr.rod2->nv() - 1);
    
    if (end_v < 2) {
      spr.nv = 0;
      m_springs[i] = spr;
      
    } else {
      
  //    int nv = min(spr.rod1->nv() - 2, spr.rod2->nv() - 2);
    
      spr.start_v = 2;
      spr.nv = min(nv_seg, end_v - spr.start_v + 1);
    
      Vec3d v1 (0,0,0);
      Vec3d v2 (0,0,0);
    
      for (int j=spr.start_v; j<spr.start_v + spr.nv; j++) {
        v1 += spr.rod1->getVertex(j);
        v2 += spr.rod2->getVertex(j);
      }
    
      v1 /= (double) spr.nv;
      v2 /= (double) spr.nv;
    
      spr.rest_length = (v1-v2).norm();
      
      m_springs[i] = spr;
      
      spr.start_v += spr.nv;
      
      while (spr.start_v < end_v) {
        spr.nv = min(nv_seg, end_v - spr.start_v + 1);
        
        Vec3d v1 (0,0,0);
        Vec3d v2 (0,0,0);
    
        for (int j=spr.start_v; j<spr.start_v + spr.nv; j++) {
          v1 += spr.rod1->getVertex(j);
          v2 += spr.rod2->getVertex(j);
        }
    
        v1 /= (double) spr.nv;
        v2 /= (double) spr.nv;
    
        spr.rest_length = (v1-v2).norm();
      
        m_springs.push_back(spr);
        
        spr.start_v += spr.nv;
      }
    }
  }
}

void RodGroupManager::Setup() {
  ns = 0;
  m_springs.clear();

  m_rod_groups.clear();
  ng = ClusterRods( ng );

  SetupSpring();

  PrecomputeSpring();

}

RodGroupSpringForce::RodGroupSpringForce(std::vector<ElasticRod*>& rods)
  : ns(0), nr(0), groupA(), groupB(), restL(), is_activated(), stiffness(30.0), m_rods(rods), m_rod_locs(), m_base_id(), m_group_id(), m_group_rods()
{
  nr = (int) m_rods.size();
  
  int base = 0;
  for( int i = 0; i < nr; ++i ) {
    m_base_id.push_back(base);
    base += m_rods[i]->ndof();
  }
  
  for( int i = 0; i < nr; ++i ) {
    Vec3d v = m_rods[i]->getVertex ( m_rods[i]->nv() * 3 / 4 );
//    Vec3d v = (m_rods[i]->getVertex ( 0 ) + m_rods[i]->getVertex ( 1 )) * 0.5;
    m_rod_locs.push_back( v );
  }
  
  // Add rod-rod springs
  
  // clustering
  std::vector<int> idx;
  for( int i = 0; i < nr; ++i ) idx.push_back(i);
  m_group_id.resize(nr, 0);
  
  ng = binaryClustering(idx, 0);
  std::cout << ng << " groups \n";
  
  // Add group-group springs
  if (ng > 1) {
    double r = 10.0;
    int nvmax = 4;  // maximum number of vertices in each segment (in each rod) 
    int nsmax = 2;  // maximum number of segments (in each rod) 
    
    std::vector <Vec3d> avg;
    
    for(int i=0; i<ng; i++) {
      Vec3d v(0,0,0);
      
      for( int j = 0; j < m_group_rods[i].size(); ++j ) {        v += m_rod_locs[m_group_rods[i][j]];      }
      v /= (double)m_group_rods[i].size();
      
      avg.push_back(v);
    }
    
    for(int i=0; i<ng; i++) {
      for(int j=i+1; j<ng; j++) {
        if ((avg[i] - avg[j]).norm() < r) {
          // add spring
          RodVertexVec gA;
          RodVertexVec gB;
          
//          Vec3d va(0,0,0);
//          Vec3d vb(0,0,0);
          
          int nv = nvmax;
          
          for(int k=0; k<m_group_rods[i].size(); k++) {
            for(int l=1; l<m_rods[m_group_rods[i][k]]->nv() && l<1+nv; l++) {
              gA.push_back(make_pair(m_group_rods[i][k], l));
//              va += m_rods[m_group_rods[i][k]]->getVertex(l);
            }
          }
      
          for(int k=0; k<m_group_rods[j].size(); k++) {
            for(int l=1; l<m_rods[m_group_rods[j][k]]->nv() && l<1+nv; l++) {
              gB.push_back(make_pair(m_group_rods[j][k], l));
//              vb += m_rods[m_group_rods[j][k]]->getVertex(l);
            }
          }
      
          if (gA.size() < 1 || gB.size() < 1) continue;
          
//          va /= (double)gA.size();
//          vb /= (double)gB.size();
          
          Vec3d va = computeCenterPosition(gA);
          Vec3d vb = computeCenterPosition(gB);
          
          double rest_length = (va - vb).norm();
          
          // Add spring into global vector
          
          groupA.push_back(gA);
          groupB.push_back(gB);
          restL.push_back(rest_length);
          
        }
      }
    }
    
  }
  
  ns = restL.size();
  is_activated.resize(ns, false);
  
  std::cout << ns << " springs\n";
  

  
}

int RodGroupSpringForce::binaryClustering(std::vector<int>& idx, int gid) {
  // if all rods (their average positions) are inside of a cirle (sphere) with radius r, they will consist the one same group.
  // Otherwise, divide.
  
  double r = 3.0;
  
  int n = (int) idx.size();
  
  if (n < 1) return gid;
  
  double dmax = 0;
  if (n > 1) {
    
    Vec3d v(0,0,0);
    for( int i=0; i<n; i++) {    v += m_rod_locs[idx[i]];  }
    v /= (double) n;
    
    for( int i=0; i<n; i++) {   
      if ( (m_rod_locs[idx[i]] - v).norm() > dmax ) dmax = (m_rod_locs[idx[i]] - v).norm();
    }
    
  }
    
  if (n==1 || dmax < r) {
    // the same group
    for( int i=0; i<n; i++) {   
      m_group_id[idx[i]] = gid;
      m_rods[idx[i]]->draw_cl = gid+1;
    }
    
    m_group_rods.push_back(idx);
    
    // done
    return gid+1;
    
  } else {
    std::vector<int> ida;
    std::vector<int> idb;
    Vec3d va = m_rod_locs[idx[0]];
    Vec3d vb = m_rod_locs[idx[n-1]];
    
    int it = 30;
    while(it-- > 0) {
      Vec3d nva(0,0,0);
      Vec3d nvb(0,0,0);
      ida.clear();
      idb.clear();
      
      for(int i=0; i<n; i++) {
        double da = (m_rod_locs[idx[i]] - va).norm();
        double db = (m_rod_locs[idx[i]] - vb).norm();
        
        if (da < db) {
          ida.push_back(idx[i]);
          nva += m_rod_locs[idx[i]];
        } else {
          idb.push_back(idx[i]);
          nvb += m_rod_locs[idx[i]];
        }
      }
      
      va = nva / (double)(ida.size());
      vb = nvb / (double)(idb.size());
    }
    
    int gida = binaryClustering(ida, gid);
    return binaryClustering(idb, gida);
    
    
  }
}

RodGroupSpringForce::~RodGroupSpringForce() 
{
}

Vec3d& RodGroupSpringForce::computeCenterPosition( RodVertexVec& rv ) {
  Vec3d v(0,0,0);
  
  for(int i=0; i<(int)rv.size(); i++) {
    v += m_rods[rv[i].first]->getVertex(rv[i].second);
  }
  
  v /= (double) rv.size();
  
  return v;
  
}

void RodGroupSpringForce::checkActivatingCondition() {
  int na = 0;
      
  for(int i=0; i<ns; i++) {
    RodVertexVec& gA = groupA[i];
    RodVertexVec& gB = groupB[i];
    
    Vec3d va = computeCenterPosition(gA);
    Vec3d vb = computeCenterPosition(gB);
          
    double len = (va - vb).norm();
  
    if (len < restL[i]) {
      is_activated[i] = true; na++;
    } else {
      is_activated[i] = false;
    }
  }
  
  std::cout << na << " springs active\n";
}

void RodGroupSpringForce::computeForce( VecXd& force ) {
  
  for(int i=0; i<ns; i++) {
    if (is_activated[i]) {
      RodVertexVec& gA = groupA[i];
      RodVertexVec& gB = groupB[i];
      
      double nA = (double) gA.size();
      
      Vec3d va = computeCenterPosition(gA);
      Vec3d vb = computeCenterPosition(gB);
            
      double len = (va - vb).norm();
    
      Vec3d fa = -stiffness * (len - restL[i]) / (nA * len) * (va - vb);
      
      for (RodVertexVecIterator rvIt = gA.begin(); rvIt!=gA.end(); ++rvIt)
      {
        int rid = (*rvIt).first;
        int vid = (*rvIt).second;
        
        for( int coord = 0; coord < 3; ++coord ) {
          force(m_base_id[rid] + m_rods[rid]->vertIdx(vid,coord)) += fa(coord);
        }
      }
    }
  }
  
}

void RodGroupSpringForce::computeForceDX( Scalar scale, MatrixBase& J ){
  
//  J.finalizeNonzeros();
  
  for(int i=0; i<ns; i++) {
    if (1 || is_activated[i]) {
      RodVertexVec& gA = groupA[i];
      RodVertexVec& gB = groupB[i];
      
      double nA = (double) gA.size();
      
      Vec3d va = computeCenterPosition(gA);
      Vec3d vb = computeCenterPosition(gB);
            
      Vec3d vs = (va - vb);
      double len = (va - vb).norm();
    
      Mat3d jac = - stiffness / (nA*nA) * ( vs * vs.transpose() * restL[i] / (len*len*len) + (len - restL[i]) / len * Mat3d::Identity());
      jac *= -scale;
      
      std::vector<std::vector<int>> indicesList;
      for (RodVertexVecIterator rvIt = gA.begin(); rvIt!=gA.end(); ++rvIt)
      {
        int rid = (*rvIt).first;
        int vid = (*rvIt).second;
        
        std::vector<int> indicesA;
        for( int coord = 0; coord < 3; ++coord ) 
        {
          indicesA.push_back(m_base_id[rid] + m_rods[rid]->vertIdx(vid,coord));
        }
        indicesList.push_back(indicesA);
      }
    
      for(int j=0; j<indicesList.size(); j++) {
        for(int k=0; k<indicesList.size(); k++) {
          IndexArray indicesA(3);
          IndexArray indicesB(3);
          
          for (int l=0; l<indicesList[j].size(); l++) {
            std::cout << indicesList[j][l] << " ";
            indicesA(l) = indicesList[j][l];
          }
          for (int l=0; l<indicesList[k].size(); l++) {
            std::cout << indicesList[k][l] << " ";
            indicesB(l) = indicesList[k][l];
          }
          std::cout << "\n";
          
          //J.add(indicesList[j], indicesList[k], jac);
          J.add(indicesA, indicesB, jac);
        }
      }
      
    }
  }
  
}

*/


