#include "ImplicitEulerTestsOne.hh"

ImplicitEulerTestsOne::ImplicitEulerTestsOne()
: Problem("Implicit Euler Tests", "Tests the stability of implicit euler by pulling a rod with a fixed edge.")
//, m_br_stpr(NULL)
{
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();

  // Values for the default options
  GetStringOpt("integrator") = "implicit";
  GetVecOpt("gravity") = Vec3d(0,0,0);
  GetScalarOpt("dt") = 0.001;
  // Maximum number of implicit solver itations
  GetIntOpt("iterations") = 1000; 
  GetIntOpt("nv") = 23;
  GetScalarOpt("mass-damping") = 0;
  GetScalarOpt("viscosity") = 0;
  GetScalarOpt("radius-scale") = 10.0;

  GetBoolOpt("quasistatic") = false;
  
  //AddOption("radius-scale", "Amount to scale major/minor radii by for rendering and collision detection/response.", 10.0);
  //AddOption("straight-undeformed", "True generates a straight undeformed config, false generates a helix.", true);
}

ImplicitEulerTestsOne::~ImplicitEulerTestsOne()
{
}

void ImplicitEulerTestsOne::Setup()
{
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  //opts.radiusScale = GetScalarOpt("radius-scale");
  
  // 12 inches == 60.96 cm
  double L = 60.96/3.0;  
  double dh = L/((double)opts.numVertices-1);
  
  Vec3d baseposition(0.0,0.0,0.0);

  std::vector<Vec3d> deformed;
  deformed.push_back(baseposition);
  for( int k = 0; k < opts.numVertices-1; ++k )
  {
    Vec3d vert = deformed.back();
    if( k%2 == 0 ) vert += Vec3d(1.0/sqrt(2.0),1.0/sqrt(2.0),0.0);
    else vert += Vec3d(1.0/sqrt(2.0),-1.0/sqrt(2.0),0.0);
    deformed.push_back(vert);
  }

  std::vector<Vec3d> undeformed;
  undeformed.push_back(baseposition);
  for( int k = 0; k < opts.numVertices-1; ++k )
  {
    Vec3d vert = undeformed.back();
    vert += Vec3d(1.0,0.0,0.0);
    undeformed.push_back(vert);
  }
  
  
  ElasticRod* newrod = setupRod(opts,deformed,undeformed);

  newrod->getBoundaryCondition()->setDesiredVertexPosition(0,newrod->getVertex(0));
  newrod->getBoundaryCondition()->setDesiredVertexPosition(1,newrod->getVertex(1));
  newrod->getBoundaryCondition()->setDesiredVertexPosition(9,newrod->getVertex(9));
  newrod->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
  
  m_world->addObject(newrod);
  m_rods.push_back(newrod);

  RodTimeStepper* newstepper = getRodTimeStepper(*newrod);
  m_world->addController(newstepper);
  
  std::cout << "Solver Options: " << std::endl;
  std::cout << "   integrator: " << GetStringOpt("integrator") << std::endl;
  std::cout << "   dt:         " << GetScalarOpt("dt") << std::endl;
  std::cout << "   iterations: " << GetIntOpt("iterations") << std::endl;
  std::cout << "   solver: " << SolverUtils::instance()->getSolverName() << std::endl;
  
  std::cout << "\"World\" Options: " << std::endl;
  std::cout << "   mass-damping: " << GetScalarOpt("mass-damping") << std::endl;
  
  std::cout << "Rod Options: " << std::endl;
  std::cout << "   numVertices:   " << opts.numVertices << std::endl;
  std::cout << "   density:       " << opts.density << std::endl;
  std::cout << "   radiusA:       " << opts.radiusA << std::endl;
  std::cout << "   radiusB:       " << opts.radiusB << std::endl;
  std::cout << "   radiusScale:   " << opts.radiusScale << std::endl;
  std::cout << "   YoungsModulus: " << opts.YoungsModulus << std::endl;
  std::cout << "   ShearModulus:  " << opts.ShearModulus << std::endl;
  std::cout << "   viscosity:     " << opts.viscosity << std::endl;
  std::cout << "   anisotropic:   " << opts.anisotropic << std::endl;
  std::cout << "   elastic:       " << opts.elastic << std::endl;
  std::cout << "   quasistatic:   " << opts.quasistatic << std::endl;
  std::cout << "   inextensible:  " << opts.inextensible << std::endl;
  std::cout << "   reframe:       ";
  switch(opts.refFrame) 
  {
    case ElasticRod::TimeParallel:
      std::cout << "TimeParallel" << std::endl;
      break;
    case ElasticRod::SpaceParallel:
      std::cout << "SpaceParallel" << std::endl;
      break;
    default:
      std::cout << "WARNING NONHANDLED CASE" << std::endl;
      break;
  }  
}


void ImplicitEulerTestsOne::AtEachTimestep()
{
//  // cm / second
//  double pull_rate = 5.0;
//
//  if( fmod((getTime()/10.0),4.0) < 1.0 )
//  {
//    Vec3d dx = Vec3d(pull_rate*GetScalarOpt("dt"),0.0,0.0);
//    for( int i = 0; i < (int) m_rods.size(); ++i ) 
//    {
//      m_rods[i]->setVertex(0,m_rods[i]->getVertex(0)+dx);
//      m_rods[i]->setVertex(1,m_rods[i]->getVertex(1)+dx);
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(0,m_rods[i]->getVertex(0));
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(1,m_rods[i]->getVertex(1));
//      m_rods[i]->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
//    }
//  }
//  else if( fmod((getTime()/10.0),4.0) < 2.0 )
//  {
//    Vec3d dz = -Vec3d(0.0,0.0,pull_rate*GetScalarOpt("dt"));
//    for( int i = 0; i < (int) m_rods.size(); ++i ) 
//    {
//      m_rods[i]->setVertex(0,m_rods[i]->getVertex(0)+dz);
//      m_rods[i]->setVertex(1,m_rods[i]->getVertex(1)+dz);
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(0,m_rods[i]->getVertex(0));
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(1,m_rods[i]->getVertex(1));
//      m_rods[i]->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
//    }
//  }
//  else if( fmod((getTime()/10.0),4.0) < 3.0 )
//  {
//    Vec3d dx = -Vec3d(pull_rate*GetScalarOpt("dt"),0.0,0.0);
//    for( int i = 0; i < (int) m_rods.size(); ++i ) 
//    {
//      m_rods[i]->setVertex(0,m_rods[i]->getVertex(0)+dx);
//      m_rods[i]->setVertex(1,m_rods[i]->getVertex(1)+dx);
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(0,m_rods[i]->getVertex(0));
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(1,m_rods[i]->getVertex(1));
//      m_rods[i]->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
//    }
//  }
//  else
//  {
//    Vec3d dz = Vec3d(0.0,0.0,pull_rate*GetScalarOpt("dt"));
//    for( int i = 0; i < (int) m_rods.size(); ++i ) 
//    {
//      m_rods[i]->setVertex(0,m_rods[i]->getVertex(0)+dz);
//      m_rods[i]->setVertex(1,m_rods[i]->getVertex(1)+dz);
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(0,m_rods[i]->getVertex(0));
//      m_rods[i]->getBoundaryCondition()->setDesiredVertexPosition(1,m_rods[i]->getVertex(1));
//      m_rods[i]->getBoundaryCondition()->setDesiredEdgeAngle(0,0.0);
//    }
//  }
}





