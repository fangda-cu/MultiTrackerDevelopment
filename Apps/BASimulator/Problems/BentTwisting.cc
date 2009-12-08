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

BentTwisting::~BentTwisting()
{
  if (rod != NULL) delete rod;
  if (stepper != NULL) delete stepper;
}

void BentTwisting::Setup()
{
  /*BandMatrix mat(100, 100, 7, 7);
  int num = 1;
  for (int i = 0; i < mat.rows(); ++i) {
    int lower = std::max(i - 7, 0);
    int upper = std::min(i + 7 + 1, mat.cols());
    for (int j = lower; j < upper; ++j) {
      mat(i, j) = num;
      ++num;
    }
    }*/

  /*BandMatrix A(2, 2, 1, 1);
  A(0, 0) = 4;
  A(0, 1) = 1;
  A(1, 0) = 1;
  A(1, 1) = 3;
  A.print();

  VecXd b(2);
  b(0) = 1;
  b(1) = 2;
  std::cout << b << std::endl;

  VecXd x(2);
  x(0) = 2;
  x(1) = 1;

  ConjugateGradient cg(A);
  cg.solve(x, b);

  std::cout << x << std::endl;

  exit(0);*/
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
