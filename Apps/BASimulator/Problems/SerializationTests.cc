#include "SerializationTests.hh"

SerializationTests::SerializationTests()
: Problem("Serialization Tests", "A collection of tests to ensure serialization functions correctly")
, m_rod(NULL)
, m_stepper(NULL)
{
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();
 
  // Values for the default options
  GetStringOpt("integrator") = "implicit";
  GetVecOpt("gravity") = Vec3d(0,0,0);
  GetScalarOpt("dt") = 0.01;
  GetScalarOpt("mass-damping") = 0;
  GetScalarOpt("viscosity") = 0;
  GetScalarOpt("radius-scale") = 10.0;
  
  GetBoolOpt("quasistatic") = false;

  AddOption("SerializationTest", "Which serialization test to execute", 0);
}

SerializationTests::~SerializationTests()
{
  if( m_rod != NULL )
  {
    delete m_rod;
    m_rod = NULL;
  }
  if( m_stepper != NULL )
  {
    delete m_stepper;
    m_stepper = NULL;
  }
}

void SerializationTests::setupTwoVertexNoMovement()
{
  assert( m_rod == NULL );
  assert( m_stepper == NULL );
  assert( !GetBoolOpt("quasistatic") );

  loadDynamicsProps();

  RodOptions opts;
  getRodOptions(opts);
  opts.numVertices = 2;

  std::vector<Vec3d> undeformed;
  undeformed.push_back(Vec3d(0.0,0.0,0.0));
  undeformed.push_back(Vec3d(2.0,0.0,0.0));
  std::vector<Vec3d> deformed = undeformed;

  m_rod = setupRod(opts,deformed,undeformed);
  assert( m_rod != NULL );
  m_stepper = getRodTimeStepper(*m_rod);
  assert( m_stepper != NULL );

  m_world->addObject(m_rod);
  m_world->addController(m_stepper);
  
  insertBreakpoint( 10.0 );
  insertBreakpoint( 20.0 );
}

// This test has a difference in D0D0D0D0, which is PROBLEMBASEMAGIC
void SerializationTests::setupTwoVertexConstantVelocityNoForcesExcited()
{
  assert( m_rod == NULL );
  assert( m_stepper == NULL );
  assert( !GetBoolOpt("quasistatic") );
  
  loadDynamicsProps();
  
  RodOptions opts;
  getRodOptions(opts);
  opts.numVertices = 2;
  
  std::vector<Vec3d> undeformed;
  undeformed.push_back(Vec3d(0.0,0.0,0.0));
  undeformed.push_back(Vec3d(2.0,0.0,0.0));
  std::vector<Vec3d> deformed = undeformed;
  
  m_rod = setupRod(opts,deformed,undeformed);
  assert( m_rod != NULL );
  m_stepper = getRodTimeStepper(*m_rod);
  assert( m_stepper != NULL );
  
  m_rod->setVelocity(0,Vec3d(0.0,0.0,0.2));
  m_rod->setVelocity(1,Vec3d(0.0,0.0,0.2));
  
  m_world->addObject(m_rod);
  m_world->addController(m_stepper);
  
  insertBreakpoint( 1.0 );
  insertBreakpoint( 2.0 );
}


void SerializationTests::Setup()
{
  switch (GetIntOpt("SerializationTest")) 
  {
    case 0:
    {
      std::cout << "Executing test: Two Vertex No Movement" << std::endl;
      setupTwoVertexNoMovement();
      break;
    }
    case 1:
    {
      std::cout << "Executing test: Two Vertex Constant Velocity No Forces Excited" << std::endl;
      setupTwoVertexConstantVelocityNoForcesExcited();
      break;
    }
    default:
    {
      std::cout << "UNKOWN TEST CASE: " << GetIntOpt("ViscosityTest") << ", exiting." << std::endl;
      exit(1);
      break;
    }
  }
}

void SerializationTests::serialize( std::ofstream& of )
{
  assert( of.is_open() );
  
  // Serialize the parent class
  Problem::serializeProblem(of);
  
  int SERIALIZATIONTESTMAGIC = 0xC0C0C0C0;
  of.write((char*)&SERIALIZATIONTESTMAGIC,sizeof(int));
  of.write((char*)&SERIALIZATIONTESTMAGIC,sizeof(int));
  
  // Serialize the rod
  TopologicalObjectSerializer objserializer;
  objserializer.appendTopologicalObjectToFile(*m_rod,of);

  of.write((char*)&SERIALIZATIONTESTMAGIC,sizeof(int));
  of.write((char*)&SERIALIZATIONTESTMAGIC,sizeof(int));


//  objserializer.saveTopologicalObject( *m_rod, "tempfile.bin" );
//  delete m_rod;
//  m_rod = NULL;
//  objserializer.loadTopologicalObject( &m_rod, "tempfile.bin" );
//  assert( m_rod != NULL );

  // RodTimeStepper* m_stepper; < Just explicitly recreate this
}

void SerializationTests::resumeFromfile( std::ifstream& ifs )
{
  assert( m_rod == NULL );
  assert( m_stepper == NULL );

  Problem::resumeProblem(ifs);
  
  int byteeater;
  ifs.read((char*)&byteeater,sizeof(int));
  ifs.read((char*)&byteeater,sizeof(int));  
  
  // Load the serialized rod
  TopologicalObjectSerializer objserializer;
  objserializer.loadTopologicalObjectFromFile(&m_rod,ifs);
  assert( m_rod != NULL );

  //RodState statebackup;
  //statebackup.copyState( *rod );
  //statebackup.print(*rod);

  ifs.read((char*)&byteeater,sizeof(int));
  ifs.read((char*)&byteeater,sizeof(int));  
  
  m_stepper = getRodTimeStepper(*m_rod);
  assert( m_stepper != NULL );

  m_world->addObject(m_rod);
  m_world->addController(m_stepper);
}








