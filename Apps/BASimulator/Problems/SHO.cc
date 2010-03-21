/**
 * \file SHO.cc
 *
 * \author miklos@cs.columbia.edu
 * \date 09/09/2009
 */

#include "SHO.hh"

SHO::SHO()
  : Problem("SHO", "test problem for the rods")
  , rod(NULL)
  , stepper(NULL)
{
  addDynamicsProps();

  AddOption("integrator", "type of integrator to use for the rod",
            "implicit");
}

void SHO::Setup()
{
  loadDynamicsProps();

  RodOptions opts;
  opts.numVertices = 3;
  opts.YoungsModulus = 1000;
  opts.ShearModulus = 375;
  opts.radiusA = 1;
  opts.radiusB = 1;
  std::vector<Vec3d> vertices, undeformed;
  for (int i = 0; i < opts.numVertices; ++i) {
    vertices.push_back(Vec3d(0, -10 * i, 0));
  }

  rod = setupRod(opts, vertices, vertices);
  stepper = new RodTimeStepper(*rod);
  stepper->getBoundaryCondition()->setDesiredVertexPosition(0, rod->getVertex(0));
  stepper->getBoundaryCondition()->setDesiredVertexPosition(1, rod->getVertex(1));
  stepper->getBoundaryCondition()->setDesiredEdgeAngle(0, rod->getTheta(0));
  //stepper->addExternalForce(new RodMassDamping(10.0));
  if (getGravity().norm() > 0)
    stepper->addExternalForce(new RodGravity(getGravity()));
  std::string& method = GetStringOpt("integrator");
  if (method == "implicit")
    stepper->setDiffEqSolver(RodTimeStepper::IMPL_EULER);
  else if (method == "symplectic")
    stepper->setDiffEqSolver(RodTimeStepper::SYMPL_EULER);
  stepper->setTimeStep(0.001);
  m_world->addObject(rod);
  m_world->addController(stepper);
}

void SHO::AtEachTimestep()
{
  //VecXd F(rod->ndof());
  //F.setZero();
  //stepper->evaluatePDot(F);
  //std::cout << F << std::endl;
  /*
  std::vector<RodForce*>& forces = rod->getForces();
  std::cout << "internal forces" << std::endl;
  for (size_t i = 0; i < forces.size(); ++i) {
    F.setZero();
    RodForce& force = *forces[i];
    force.globalForce(F);
    std::cout << F << std::endl;
  }
  std::vector<RodExternalForce*>& external = stepper->getExternalForces();
  std::cout << "external forces" << std::endl;
  for (size_t i = 0; i < external.size(); ++i) {
    F.setZero();
    external[i]->computeForce(*rod,F);
    std::cout << F << std::endl;
  }
  */
/*
  Scalar dt = 0.1;
  if (rod->getTheta(0) < 10) {
    rod->setTheta(0, rod->getTheta(0) + dt);
    rod->updateProperties();
    }*/
}
