/**
 * \file HandTest.hh
 *
 * \author 
 * \date 
 */

#include "HandTest.hh"

#include <iostream>
#include <fstream>
#include <string>
using namespace std;

HandController::HandController( TriangleMesh& mesh,double time, double dt, std::vector<Vec3d>& scripted_translation, std::vector<double>& scripted_time )
: ScriptingController(time,dt)
, m_mesh(mesh)
, m_time(time)
, m_dt(dt)
, hand_translation(scripted_translation)
, hand_time(scripted_time)
{
}

bool HandController::execute()
{
	int time_range = -1;
	
	for(int i=0; i<(int)hand_time.size() - 1; i++) {
		if (hand_time[i] <= m_time && m_time <= hand_time[i+1]) {
			time_range = i;
			break;
		}
	}		
	
//	std::cout << m_time << " " << m_dt << "\n";
	
	if (time_range >= 0) {
		Vec3d x0 = hand_translation[time_range];
		Vec3d x1 = hand_translation[time_range + 1];

		Scalar dt = hand_time[time_range + 1] - hand_time[time_range];		
		
		Vec3d v = (x1 - x0) / dt;
		Vec3d dx = v * m_dt;
		
		for( TriangleMesh::vertex_iter itr = m_mesh.vertices_begin(); itr != m_mesh.vertices_end(); ++itr )
		{
	    m_mesh.getVertex(*itr) += dx;
		}			
	}
		  
  return true;
}

void HandController::setTime( double time )
{
  m_time = time;
}

void HandController::setDt( double dt )
{
  m_dt = dt;
}

double HandController::getTime()
{
  return m_time;
}

double HandController::getDt()
{
  return m_dt;
}

HandTest::HandTest()
: Problem("Hand Test", "Hand Test.")
, m_rods()
, m_tri_objs()
, m_controllers()
, m_br_stepper(NULL)
, m_scripting_controllers()
{
	hand_mesh = NULL;
	
  addDynamicsProps();
  addRodOptions();
  addRodTimeStepperOptions();
  
  AddOption("penalty-collision", "using implicit penalty collision", true);
  AddOption("penalty-collision-stiff", "stiffness of implicit penalty force", 2.0);
  AddOption("penalty-collision-extrathick", "extra thickness of implicit penalty force", 0.3);
  
  AddOption("hand-input-filename", "The name of file containing hand data", "assets/handtest/hand.obj");
  AddOption("rod-input-filename", "The name of file containing rod data", "assets/handtest/fur.dat");
  AddOption("rod-input-rodnum", "The number of rods to be used", 1);
  
//  AddOption("vertical-gravity", "vertical-gravity", -981.0);

  AddOption("hand-relative-speed", "hand-relative-speed", 1.0);
  AddOption("rod-input-vertex-density", "One vertex will be used per this number, the others will be skipped", 2);
  
  // Set Defaults for built-in options
  
//  GetVecOpt("gravity") = Vec3d(0, GetScalarOpt("vertical-gravity") ,0); // 981 cm/s^2 = 9.81 m/s^2
  GetVecOpt("gravity") = Vec3d(0, 0, 0);
  GetScalarOpt("dt") = 0.002;
    
}

HandTest::~HandTest()
{
  for( int i = 0; i < (int) m_rods.size(); ++i )
  {
    assert( m_rods[i] != NULL );
    delete m_rods[i];
    m_rods[i] = NULL;
  }
  
  for( int i = 0; i < (int) m_controllers.size(); ++i )
  {
    assert( m_controllers[i] != NULL );
    delete m_controllers[i];
    m_controllers[i] = NULL;
  }
}

void HandTest::Setup()
{
//  GetVecOpt("gravity") = Vec3d(0, GetScalarOpt("vertical-gravity") ,0); // 981 cm/s^2 = 9.81 m/s^2
  
  loadDynamicsProps();
  RodOptions opts;
  
  opts.viscosity = 0.0;
  opts.radiusScale = 5;
  
	std::string name = GetStringOpt ("rod-input-filename");
	std::ifstream myfile;
	
	myfile.open(name.c_str() , std::ifstream::in);

	int nr = 0;
	int max_r = GetIntOpt("rod-input-rodnum");
	
	if (myfile.is_open())
	{
	  while ( (!myfile.eof()) && (nr<max_r || max_r <= 0) )
	  {
	  	int nv;
	  	std::vector<Vec3d> vertices;
	  	
	  	if (!(myfile >> nv)) break;
	  	
	  	int vcount = 0;
	  	
	  	for(int i=0; i<nv; i++) {
				Scalar x, y, z;
				
				myfile >> x >> y >> z ;
				
				if (vcount % GetIntOpt("rod-input-vertex-density") != 0) {
					vcount++;
					continue;
				}
				vcount++;
				Vec3d pos(x,y,z);
				
				vertices.push_back(pos);
				
				std::cout << vertices[i] << "\n";
			}
			
			nv = (int) vertices.size();

			opts.numVertices = nv;

			ElasticRod* newrod = setupRod(opts, vertices, vertices);
			
			newrod->getBoundaryCondition()->setDesiredVertexPosition(0, newrod->getVertex(0));
			newrod->getBoundaryCondition()->setDesiredVertexPosition(1, newrod->getVertex(1));
		  newrod->getBoundaryCondition()->setDesiredEdgeAngle(0, newrod->getTheta(0));
		  
			RodTimeStepper* newstepper = getRodTimeStepper(*newrod);
			m_rods.push_back(newrod);
			m_controllers.push_back(newstepper);

			m_world->addObject(newrod);
//				m_world->addController(newstepper);

			nr++;
	  }
	  myfile.close();
	}

	else {
		cout << "Unable to open file"; 
  	exit(1);
	}
  
  setupHand();
  
  m_br_stepper = new BARodStepper( m_rods, m_tri_objs, m_scripting_controllers, m_controllers, GetScalarOpt("dt") );
  m_world->addController(m_br_stepper);
  
  m_br_stepper->enableIterativeInelasticImpulses();
//  m_br_stepper->setEdgeEdgePenalty(100.0);
	if (GetBoolOpt("penalty-collision")) {
	  m_br_stepper->enableImplicitPenaltyImpulses();
		m_br_stepper->setImplicitPenaltyStiffness(GetScalarOpt("penalty-collision-stiff"));  
		m_br_stepper->setImplicitPenaltyThickness(GetScalarOpt("penalty-collision-extrathick"));
	
	} else {
	  m_br_stepper->disableImplicitPenaltyImpulses();
	}
   
}

void HandTest::setupHand() {
  std::vector<Vec3d> hand_translation;
  std::vector<double> hand_time;
	
	Scalar t = GetScalarOpt("hand-relative-speed") / 24.0;
	
  hand_translation.clear();
  hand_time.clear();
  
  hand_time.push_back(0);
  hand_translation.push_back(Vec3d(29.956095599946128, 7.5180718557333197, -16.580964140147159));

  hand_time.push_back(9.0 * t);
  hand_translation.push_back(Vec3d(46.53932836691812, -1.0532844046782408, -16.580964140147159));

  hand_time.push_back(11.0 * t);
  hand_translation.push_back(Vec3d(46.53932836691812, -1.9677718745971124, -17.382847474214532));

  hand_time.push_back(14.0 * t);
  hand_translation.push_back(Vec3d(46.53932836691812, -2.4085183530238519, -20.436171976519496));

  hand_time.push_back(19.0 * t);
  hand_translation.push_back(Vec3d(46.53932836691812, -2.5785176276134809, -28.146586204550552));

  hand_time.push_back(24.0 * t);
  hand_translation.push_back(Vec3d(46.53932836691812, -2.5860296998671419, -38.059309814220462));

  hand_time.push_back(29.0 * t);
  hand_translation.push_back(Vec3d(46.53932836691812, -2.2016343462671468, -47.42262456498986));

  hand_time.push_back(40.0 * t);
  hand_translation.push_back(Vec3d(46.53932836691812, 9.0020439198977655, -61.312565478462119));

	Vec3d base = hand_translation[3];

	for(int i=0; i<(int)hand_translation.size(); i++) {
		hand_translation[i] -= base;
	}

  ObjParser objparser;
  hand_mesh = new TriangleMesh();
  objparser.loadTriangularMesh( GetStringOpt ("hand-input-filename"), *hand_mesh );
  for( TriangleMesh::vertex_iter itr = hand_mesh->vertices_begin(); itr != hand_mesh->vertices_end(); ++itr )
  {
    hand_mesh->getVertex(*itr) += hand_translation[0] - hand_translation[3];
  }
  m_world->addObject(hand_mesh);  

	m_tri_objs.push_back(hand_mesh);
	
	hand_controller = new HandController ( *hand_mesh, getTime(), GetScalarOpt("dt"), hand_translation, hand_time);
  m_scripting_controllers.push_back(hand_controller);
	
}


void HandTest::AtEachTimestep()
{

  for( int i = 0; i < (int) m_controllers.size(); ++i ) m_controllers[i]->setTimeStep(GetScalarOpt("dt"));
  if( m_br_stepper != NULL ) m_br_stepper->setDt(GetScalarOpt("dt"));
  hand_controller->setDt(GetScalarOpt("dt"));
  hand_controller->setTime(getTime());
  
}
