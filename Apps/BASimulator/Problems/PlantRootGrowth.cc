/**
 * \file PlantRootGrowth.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 05/13/2010
 */

#include "PlantRootGrowth.hh"
#include <cstdlib>

using namespace std;

PlantRootGrowth::PlantRootGrowth()
  : Problem("Plant root growth", "Simulations of a growing plant root")
  , rod(NULL)
  , stepper(NULL)
  , gelForce(NULL)
  , rootLength(3.0)
  , gelStiffness(100.0)
{
  // add default options associated with dynamics, rods, and time stepping
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();

  // default to no gravity
  GetVecOpt("gravity") = Vec3d::Zero();

  // default to solving statics instead of dynamics
  GetStringOpt("integrator") = "statics";

  // change default tolerance of nonlinear solver
  GetScalarOpt("atol") = 1e-5;

  // default to circular cross-section
  GetScalarOpt("major-radius") = 0.1;
  GetScalarOpt("minor-radius") = 0.1;

  // default to no velocity-dependent forces
  GetScalarOpt("viscosity") = 0.0;
  GetScalarOpt("mass-damping") = 0.0;

  // add problem-specific options
  AddOption("root-length", "length of the root", rootLength);
  AddOption("gel-stiffness", "stiffness of gel surrounding root", gelStiffness);
}

PlantRootGrowth::~PlantRootGrowth()
{
  if (rod != NULL) delete rod;
  if (stepper != NULL) delete stepper;
}

void PlantRootGrowth::Setup()
{
  // loads dynamics options that were parsed from the command line
  loadDynamicsProps();

  // loads options associated to the elastic rod from the command line.
  // RodOptions structure defined in BASim/src/Physics/ElasticRods/RodUtils.hh
  RodOptions opts;
  getRodOptions(opts);

  // load problem-specific options
  rootLength = GetScalarOpt("root-length");
  gelStiffness = GetScalarOpt("gel-stiffness");

  // set up undeformed configuration of elastic rod
  std::vector<Vec3d> undeformed;
  for (int i = 0; i < opts.numVertices; ++i) {
    undeformed.push_back( Vec3d( 0,
                                 i * rootLength / (opts.numVertices - 1),
                                 0 ) );
  }

  // set up initial configuration of elastic rod
  std::vector<Vec3d> initial = undeformed;
  // add some amount of noise to the initial configuration to break
  // symmetry
  for (int i = 0; i < opts.numVertices; ++i) {
    Scalar x = 0.01 * (rand() - RAND_MAX / 2) / RAND_MAX;
    Scalar z = 0.01 * (rand() - RAND_MAX / 2) / RAND_MAX;
    initial[i] += Vec3d( x, 0, z );
  }

  // create the rod and the time stepper
  rod = setupRod(opts, initial, undeformed);
  stepper = getRodTimeStepper(*rod);

  // add external force that simulates the resistance of the gel
  // surrounding the root
  gelForce = new RodGelForce(gelStiffness);
  stepper->addExternalForce(gelForce);

  // set up the boundary conditions for the rod
  RodBoundaryCondition* boundary = stepper->getBoundaryCondition();
  int nv = rod->nv();
  int ne = rod->ne();
  boundary->setDesiredVertexPosition(0, Vec3d(0, 0, 0));
  //boundary->setDesiredVertexPosition(1, Vec3d(0, rootLength / ne, 0));
  //boundary->setDesiredEdgeAngle(0, rod->getTheta(0));
  Vec3d last( 0, (1.0 - 0.5 / ne) * rootLength, 0 );
  Vec3d prev( 0, (1.0 - 1.5 / ne) * rootLength, 0 );
  boundary->setDesiredVertexPosition(nv - 2, prev);
  boundary->setDesiredVertexPosition(nv - 1, last);
  boundary->setDesiredEdgeAngle(ne - 1, rod->getTheta(ne - 1));

  // add the rod and the time stepper to the world so it gets
  // simulated
  m_world->addObject(rod);
  m_world->addController(stepper);
}

void PlantRootGrowth::AtEachTimestep()
{
  // create a new rod based on old rod

  // create a new stepper for the new rod

  // add gel force to new rod

  // set up boundary conditions for new rod

  // add new rod and stepper to the world to get simulated

  // remove old rod and stepper from the world and delete the memory
}
