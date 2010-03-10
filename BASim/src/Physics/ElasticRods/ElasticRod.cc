/**
 * \file ElasticRod.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 08/29/2009
 */

#include "ElasticRod.hh"
#include "../../Math/Math.hh"
#include "RodForce.hh"

namespace BASim {

using namespace Util;

ElasticRod::ElasticRod(int numVertices, bool closed)
{
  // create vertices
  for (int i = 0; i < numVertices; ++i) addVertex();

  // create edges
  vertex_iter vit = vertices_begin();
  vertex_handle vh = *vit;
  for (++vit; vit != vertices_end(); ++vit) {
    addEdge(vh, *vit);
    vh = *vit;
  }

  add_property(m_forces, "forces");
  add_property(m_quasistatic, "quasistatic", true);
  add_property(m_refFrameType, "reference-frame", TimeParallel);
  add_property(m_density, "density", 1.0);
  add_property(m_fixed, "fixed_dofs");
  add_property(m_fixedVerts, "fixed_verts");
  add_property(m_fixedEdges, "fixed_edges");
  add_property(m_YoungsModulus, "Young's modulus", 1.0);
  add_property(m_ShearModulus, "shear modulus", 1.0);
  add_property(m_viscosity, "dynamic viscosity", 0.0);
  add_property(m_dt, "rod's time step size", 0.1);

  add_property(m_vertexPositions, "vertex_positions", Vec3d(0,0,0));
  add_property(m_vertexVelocities, "vertex velocities", Vec3d(0,0,0));
  add_property(m_voronoiLengths, "voronoi legths", 0.0);
  add_property(m_vertexMasses, "vertex masses", 0.0);
  add_property(m_vertexFixed, "is vertex fixed", false);
  add_property(m_phi, "phi", 0.0);
  add_property(m_curvatureBinormal, "curvature binormal", Vec3d(0,0,0));
  add_property(m_vertIdx, "vertex index", 0);

  add_property(m_theta, "theta", 0.0);
  add_property(m_thetaDot, "theta dot", 0.0);
  add_property(m_edgeRadius, "edge radius");
  add_property(m_edgeInertias, "edge inertia", 0.0);
  add_property(m_referenceDirectors, "reference directors", Util::pair<Vec3d,Vec3d>(Vec3d(1,0,0), Vec3d(0,1,0)));
  add_property(m_materialDirectors, "material directors", Util::pair<Vec3d,Vec3d>(Vec3d(1,0,0), Vec3d(0,1,0)));
  add_property(m_edges, "edges", Vec3d(0,0,0));
  add_property(m_tangents, "tangents", Vec3d(0,0,0));
  add_property(m_edgeLengths, "edge lengths", 0.0);
  add_property(m_edgeFixed, "is edge fixed", false);
  add_property(m_edgeIdx, "edge index", 0);

  m_friction = 0;
  m_cor = 0;
  m_separationStrength = 1;

  setupDofIndices();

  // For collisions
  m_currentPositions.resize(numVertices);
  m_previousPositions.resize(numVertices);
  m_currentVelocities.resize(numVertices);
  m_edgeIndices.resize(2 * ne());
  for (int i=0; i<ne(); ++i)
  {
    m_edgeIndices[2*i  ] = i;
    m_edgeIndices[2*i+1] = (i+1)%numVertices;
  }
}

ElasticRod::~ElasticRod()
{
  RodForces& forces = getForces();
  RodForces::iterator it;
  for (it = forces.begin(); it != forces.end(); ++it) {
    delete (*it);
  }
}

void ElasticRod::setup()
{
  computeEdges();
  computeTangents();
  computeCurvatureBinormals();
  computeEdgeLengths();
  computeVoronoiLengths();
  findOrthogonal(const_cast<Vec3d&>(getReferenceDirector1(0)), getTangent(0));
  computeSpaceParallel();
  computeMaterialDirectors();
  computeVertexMasses();
  if (!quasistatic()) computeEdgeInertias();
  property(m_ndof) = 3 * nv() + ne();
}

void ElasticRod::addForce(RodForce* force)
{
  assert(force != NULL);

  getForces().push_back(force);
}

void ElasticRod::computeForces(VecXd& force)
{
  RodForces& forces = getForces();
  RodForces::iterator fIt;
  for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
    (*fIt)->globalForce(force);
  }
}

void ElasticRod::computeJacobian(MatrixBase& J)
{
  RodForces& forces = getForces();
  RodForces::iterator fIt;
  for (fIt = forces.begin(); fIt != forces.end(); ++fIt)
    (*fIt)->globalJacobian(J);
}

void ElasticRod::computeEdges()
{
  for (int j = 0; j < ne(); ++j)
    setEdge(j, getVertex(j + 1) - getVertex(j));
}

void ElasticRod::computeTangents()
{
  for (int j = 0; j < ne(); ++j)
    setTangent(j, getEdge(j).normalized());
}

void ElasticRod::computeEdgeLengths()
{
  for (int j = 0; j < ne(); ++j)
    setEdgeLength(j, getEdge(j).norm());
}

void ElasticRod::computeVoronoiLengths()
{
  setVoronoiLength(0, 0.5 * getEdgeLength(0));

  for (int i = 1; i < nv() - 1; ++i) {
    setVoronoiLength(i, 0.5 * (getEdgeLength(i - 1) + getEdgeLength(i)));
  }

  setVoronoiLength(nv() - 1, 0.5 * getEdgeLength(ne() - 1));
}

void ElasticRod::computeReferenceDirectors()
{
  if (refFrameType() == SpaceParallel) computeSpaceParallel();
  else computeTimeParallel();
}

void ElasticRod::computeSpaceParallel()
{
  // transport first edge in time
  edge_iter eit = edges_begin(), end = edges_end();
  edge_handle eh = *eit;
  Vec3d t0 = getEdge(eh).normalized();
  Vec3d u = parallelTransport(getReferenceDirector1(eh), getTangent(eh), t0);
  u = (u - u.dot(t0) * t0).normalized();
  setReferenceDirector1(eh, u);
  setReferenceDirector2(eh, t0.cross(u));

  // transport along centerline (Bishop frame)
  for (++eit; eit != end; ++eit) {
    eh = *eit;
    Vec3d t1 = getEdge(eh).normalized();
    u = parallelTransport(u, t0, t1);
    u = (u - u.dot(t1) * t1).normalized();
    setReferenceDirector1(eh, u);
    setReferenceDirector2(eh, t1.cross(u));
    t0 = t1;
  }
}

void ElasticRod::computeTimeParallel()
{
  edge_iter eit, end = edges_end();
  for (eit = edges_begin(); eit != end; ++eit) {
    edge_handle& eh = *eit;
    Vec3d t = getEdge(eh).normalized();
    Vec3d u = parallelTransport(getReferenceDirector1(eh), getTangent(eh), t);
    u = (u - u.dot(t) * t).normalized();
    setReferenceDirector1(eh, u);
    setReferenceDirector2(eh, t.cross(u));
  }
}

void ElasticRod::computeCurvatureBinormals()
{
  for (int i = 1; i < nv() - 1; ++i) {
    Vec3d& kb = property(m_curvatureBinormal)[i];
    computeCurvatureBinormal(kb, getTangent(i - 1), getTangent(i));
  }
}

void ElasticRod::computePhi()
{
  if (refFrameType() == SpaceParallel) return;

  for (int i = 1; i < nv() - 1; ++i) {
    const Vec3d& u0 = getReferenceDirector1(i - 1);
    const Vec3d& u1 = getReferenceDirector1(i);
    const Vec3d& tangent = getTangent(i);
    Scalar& phi = property(m_phi)[i];

    // transport reference frame to next edge
    Vec3d ut = parallelTransport(u0, getTangent(i - 1), tangent);

    // rotate by current value of phi
    rotateAxisAngle(ut, tangent, phi);

    // compute increment to phi to align reference frames
    phi += signedAngle(ut, u1, tangent);
  }
}

void ElasticRod::computeMaterialDirectors()
{
  for (int j = 0; j < ne(); ++j) {
    Scalar c = cos(getTheta(j));
    Scalar s = sin(getTheta(j));
    const Vec3d& u = getReferenceDirector1(j);
    const Vec3d& v = getReferenceDirector2(j);
    setMaterial1(j,  c * u + s * v);
    setMaterial2(j, -s * u + c * v);
  }
}

void ElasticRod::computeVertexMasses()
{
  for (int i = 0; i < nv(); ++i) {
    Scalar mass = 0;
    if (i > 0) {
      mass +=
        computeMass(density(), radiusA(i - 1), radiusB(i - 1),
                    0.5 * getEdgeLength(i - 1));
    }
    if (i < nv() - 1) {
      mass +=
        computeMass(density(), radiusA(i), radiusB(i),
                    0.5 * getEdgeLength(i));
    }
    setVertexMass(i, mass);
  }
}

void ElasticRod::computeEdgeInertias()
{
  for (int j = 0; j < ne(); ++j) {
    Scalar a = radiusA(j);
    Scalar b = radiusB(j);
    Scalar mass = computeMass(density(), a, b, getEdgeLength(j));
    setEdgeInertia(j, 0.25 * mass * (square(a) + square(b)));
  }
}

Scalar ElasticRod::computeMass(Scalar density, Scalar a, Scalar b, Scalar h)
{
  return (density * M_PI * a * b * h);
}

void ElasticRod::setRadius(const Scalar& r)
{
  for (int j = 0; j < ne(); ++j) {
    setRadius(j, r, r);
  }
}

void ElasticRod::setRadius(const Scalar& a, const Scalar& b)
{
  for (int j = 0; j < ne(); ++j) {
    setRadius(j, a, b);
  }
}

void ElasticRod::setRadius(int j, const Scalar& r)
{
  setRadius(j, r, r);
}

void ElasticRod::setRadius(const edge_handle& eh, const Scalar& r)
{
  setRadius(eh.idx(), r, r);
}

void ElasticRod::setRadius(int j, const Scalar& a, const Scalar& b)
{
  Util::pair<Scalar, Scalar>& radius = property(m_edgeRadius)[j];
  radius.first = a;
  radius.second = b;
}

const Scalar& ElasticRod::radius() const
{
  return property(m_edgeRadius)[0].first;
}

const Scalar& ElasticRod::radiusA(int j) const
{
  return property(m_edgeRadius)[j].first;
}

const Scalar& ElasticRod::radiusA(const edge_handle& eh) const
{
  return property(m_edgeRadius)[eh].first;
}

const Scalar& ElasticRod::radiusB(int j) const
{
  return property(m_edgeRadius)[j].second;
}

const Scalar& ElasticRod::radiusB(const edge_handle& eh) const
{
  return property(m_edgeRadius)[eh].second;
}

void ElasticRod::setMaterial1(int j, const Vec3d& m1)
{
  property(m_materialDirectors)[j].first = m1;
}

void ElasticRod::setMaterial2(int j, const Vec3d& m2)
{
  property(m_materialDirectors)[j].second = m2;
}

void ElasticRod::updateProperties()
{
  computeEdges();
  computeReferenceDirectors();
  computeTangents();
  computePhi();
  computeCurvatureBinormals();
  computeEdgeLengths();
  computeVoronoiLengths();
  computeMaterialDirectors();

  RodForces& forces = getForces();
  RodForces::iterator fIt;
  for (fIt = forces.begin(); fIt != forces.end(); ++fIt)
    (*fIt)->updateProperties();

#ifdef DEBUG
  verifyProperties();
#endif // DEBUG
}

void ElasticRod::updateReferenceProperties()
{
  RodForces forces = getForces();
  RodForces::iterator it;
  for (it = forces.begin(); it != forces.end(); ++it) {
    (*it)->updateUndeformedStrain();
    (*it)->updateStiffness();
  }
}

void ElasticRod::verifyProperties()
{
  for (int j = 0; j < ne(); ++j) {
    assert(approxEq((getVertex(j + 1) - getVertex(j)).eval(), getEdge(j)));
    assert(approxEq(getTangent(j).norm(), 1.0));
    assert(approxEq((getEdgeLength(j) * getTangent(j)).eval(), getEdge(j)));

    assert(approxEq(getReferenceDirector1(j).norm(), 1.0));
    assert(approxEq(getReferenceDirector1(j).dot(getEdge(j)), 0.0));

    assert(approxEq(getReferenceDirector2(j).norm(), 1.0));
    assert(approxEq(getReferenceDirector2(j).dot(getEdge(j)), 0.0));

    assert(approxEq(getMaterial1(j).norm(), 1.0));
    assert(approxEq(getMaterial1(j).dot(getEdge(j)), 0.0));

    assert(approxEq(getMaterial2(j).norm(), 1.0));
    assert(approxEq(getMaterial2(j).dot(getEdge(j)), 0.0));
  }

  for (int i = 1; i < nv() - 1; ++i) {
    assert(approxEq(2.0 * getVoronoiLength(i),
                    getEdgeLength(i-1)+getEdgeLength(i)));
    assert(approxEq(getCurvatureBinormal(i).dot(getEdge(i - 1)), 0.0));
    assert(approxEq(getCurvatureBinormal(i).dot(getEdge(i)), 0.0));
    assert(approxEq(getCurvatureBinormal(i).norm(),
                    2.0 * tan(angle(getEdge(i - 1), getEdge(i)) / 2.0)));
  }

  RodForces& forces = getForces();
  RodForces::iterator it;
  for (it = forces.begin(); it != forces.end(); ++it) {
    (*it)->verifyProperties();
  }
}

void ElasticRod::viscousUpdate()
{
  RodForces forces = getForces();
  RodForces::iterator it;
  for (it = forces.begin(); it != forces.end(); ++it) {
    if ((*it)->viscous()) (*it)->updateUndeformedStrain();
  }
}

void ElasticRod::setupDofIndices()
{
  vertex_iter vit, vend = vertices_end();
  int i = 0;
  for (vit = vertices_begin(); vit != vend; ++vit, ++i) {
    property(m_vertIdx)[*vit] = 4 * i;
    for (int j = 0; j < 3; ++j) {
      DofHandle handle(4 * i + j);
      handle.setType(DofHandle::VERTEX_DOF);
      handle.setHandle(*vit);
      handle.setNum(j);
      property(m_map).addMapping(handle, 4 * i + j);
    }
  }

  edge_iter eit, eend = edges_end();
  i = 0;
  for (eit = edges_begin(); eit != eend; ++eit, ++i) {
    property(m_edgeIdx)[*eit] = 4 * i + 3;
    DofHandle handle(4 * i + 3);
    handle.setType(DofHandle::EDGE_DOF);
    handle.setHandle(*eit);
    handle.setNum(0);
    property(m_map).addMapping(handle, 4 * i + 3);
  }
}

const Scalar& ElasticRod::getVertexDof(const vertex_handle& vh, int num) const
{
  return const_cast<Vec3d&>(property(m_vertexPositions)[vh])[num];
}

void
ElasticRod::setVertexDof(const vertex_handle& vh, int num, const Scalar& dof)
{
  property(m_vertexPositions)[vh][num] = dof;
}

const Scalar& ElasticRod::getEdgeDof(const edge_handle& eh, int num) const
{
  return property(m_theta)[eh];
}

void
ElasticRod::setEdgeDof(const edge_handle& eh, int num, const Scalar& dof)
{
  property(m_theta)[eh] = dof;
}

const Scalar& ElasticRod::getVertexVel(const vertex_handle& vh, int num) const
{
  return const_cast<Vec3d&>(property(m_vertexVelocities)[vh])[num];
}

void
ElasticRod::setVertexVel(const vertex_handle& vh, int num, const Scalar& vel)
{
  property(m_vertexVelocities)[vh][num] = vel;
}

const Scalar& ElasticRod::getEdgeVel(const edge_handle& eh, int num) const
{
  return property(m_thetaDot)[eh];
}

void
ElasticRod::setEdgeVel(const edge_handle& eh, int num, const Scalar& vel)
{
  property(m_thetaDot)[eh] = vel;
}

const Scalar& ElasticRod::getVertexMass(const vertex_handle& vh, int num) const
{
  return property(m_vertexMasses)[vh];
}

void
ElasticRod::setVertexMass(const vertex_handle& vh, int num, const Scalar& mass)
{
  property(m_vertexMasses)[vh]= mass;
}

const Scalar& ElasticRod::getEdgeMass(const edge_handle& eh, int num) const
{
  return property(m_edgeInertias)[eh];
}

void
ElasticRod::setEdgeMass(const edge_handle& eh, int num, const Scalar& mass)
{
  property(m_edgeInertias)[eh] = mass;
}

void ElasticRod::setTimeStep(Scalar dt)
{
  property(m_dt) = dt;
  RodForces& forces = getForces();
  RodForces::iterator fIt;
  for (fIt = forces.begin(); fIt != forces.end(); ++fIt) {
    if ((*fIt)->viscous()) (*fIt)->updateStiffness();
  }
}

} // namespace BASim
