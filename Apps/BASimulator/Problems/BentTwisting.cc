/**
 * \file BentTwisting.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 09/09/2009
 */

#include "BentTwisting.hh"

BentTwisting::BentTwisting()
  : Problem("Bent Twisting", "A bent anisotropic rod is twisted")
  , rod(NULL)
  , stepper(NULL)
  , m_maxTwist(10)
  , m_twistRate(1)
  , m_currentTwist(0)
{
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();

  AddOption("shape-radius", "radius of circular arc defining shape of rod",
            20.0);
  AddOption("twist-rate", "rate of twisting of end of rod", m_twistRate);
  AddOption("twist-max", "maximum amount of twist", m_maxTwist);

  // default to no gravity
  GetVecOpt("gravity") = Vec3d::Zero();
}

void BentTwisting::Setup()
{
  loadDynamicsProps();

  RodOptions opts;
  getRodOptions(opts);

  m_twistRate = GetScalarOpt("twist-rate");
  m_maxTwist = GetScalarOpt("twist-max");

  Scalar radius = GetScalarOpt("shape-radius");
  std::vector<Vec3d> vertices, undeformed;
  for (int i = 0; i < opts.numVertices; ++i) {
    vertices.push_back(Vec3d(radius * cos(i * M_PI / (opts.numVertices - 1)),
                             radius * sin(i * M_PI / (opts.numVertices - 1)),
                             0));
  }

  rod = setupRod(opts, vertices, vertices);

  rod->fixVert(0);
  rod->fixVert(1);
  rod->fixVert(rod->nv() - 2);
  rod->fixVert(rod->nv() - 1);
  rod->fixEdge(0);
  rod->fixEdge(rod->ne() - 1);

  stepper = getRodTimeStepper(*rod);

  m_world->addObject(rod);
  m_world->addController(stepper);
}

void BentTwisting::AtEachTimestep()
{
  if (m_maxTwist > 0 && m_currentTwist >= m_maxTwist) return;

  Scalar twistIncrement = getDt() * m_twistRate;
  int edge = rod->ne() - 1;
  rod->setTheta(edge, rod->getTheta(edge) - twistIncrement);
  rod->updateProperties();
  m_currentTwist += twistIncrement;
}
