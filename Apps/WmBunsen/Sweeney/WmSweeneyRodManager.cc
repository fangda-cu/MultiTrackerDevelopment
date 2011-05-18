#include "WmSweeneyRodManager.hh"

using namespace std;
using namespace BASim;

WmSweeneyRodManager::WmSweeneyRodManager()
{
    m_bridsonStepper = NULL;
    m_rods.clear();
    m_rodTimeSteppers.clear();
    m_triangleMeshes.clear();
    m_scriptingControllers.clear();
    m_rodRenderers.clear();
}

WmSweeneyRodManager::~WmSweeneyRodManager()
{
    delete m_bridsonStepper;
}

bool WmSweeneyRodManager::addRod( const std::vector< BASim::Vec3d >& i_vertices, 
                                  const double i_time,
                                  const double i_youngsModulus,
                                  const double i_shearModulus, const double i_viscosity, 
                                  const double i_density, const double i_radiusA, 
                                  const double i_radiusB, 
                                  const BASim::ElasticRod::RefFrameType i_referenceFrame,
                                  const double i_massDamping, 
                                  const BASim::Vec3d i_gravity,
                                  const BASim::RodTimeStepper::Method i_solverType )
{
    cerr << "WmSweeneyRodManager::addRod: About to create rod\n";
    RodOptions rodOptions;
    rodOptions.YoungsModulus = i_youngsModulus; /* megapascal */
    rodOptions.ShearModulus = i_shearModulus;   /* megapascal */
    rodOptions.viscosity = i_viscosity;         /* poise */
    rodOptions.density = i_density;             /* grams per cubic centimeter */
    rodOptions.radiusA = i_radiusA;             /* millimeter */
    rodOptions.radiusB = i_radiusB;             /* millimeter */
    rodOptions.refFrame = i_referenceFrame;
    rodOptions.numVertices = ( int )( i_vertices.size() );
        
    // Use the rod helper function to build the rod
    ElasticRod* rod = setupRod( rodOptions,
                                i_vertices,
                                i_vertices );

    cerr << "WmSweeneyRodManager::addRod: setupRod returned\n";

    // We need a rod renderer to draw the rod in OpenGL
    RodRenderer* rodRenderer = new RodRenderer( *rod );
    
    // Create a timeStepper to simulate the rod forward in time
    RodTimeStepper* stepper = new RodTimeStepper( *rod );
	stepper->setDiffEqSolver( i_solverType );
	    
    // Add a damping force to the rod
    stepper->addExternalForce( new RodMassDamping( i_massDamping ) );
    
    // Lock the first two vertices of the rod
    RodBoundaryCondition* boundary = stepper->getBoundaryCondition();
                
    // Set the velocity to be zero as we're grooming static hair
    boundary->setDesiredVertexPosition( 0, i_time, rod->getVertex( 0 ), BASim::Vec3d( 0.0, 0.0, 0.0 ) );
    boundary->setDesiredVertexPosition( 1, i_time, rod->getVertex( 1 ), BASim::Vec3d( 0.0, 0.0, 0.0 ) );
 
    // If the magnitude of gravity is 0 then don't bother adding the force
    if ( i_gravity.norm() > 0 )
    {
        stepper->addExternalForce( new RodGravity( i_gravity ) );        
    }
            
    // Add a force class that we will use to pass in forces from Maya
    RodMayaForces* rodMayaForces = new RodMayaForces( rod );
    stepper->addExternalForce( rodMayaForces );

    // Reverse hairdo is still experimental and optional...
    if ( 0 )
    {
        cerr << "Doing reverse hairdo!\n";
        rod->doReverseHairdo(stepper);
    }
    
    // Arbitrarily scale the rod up so it can be seen
    rod->setRadiusScale( 10.0 );
    
    // Store all the things we need to control the rod or add it to a BridsonStepper
    m_rods.push_back( rod );
    m_rodTimeSteppers.push_back( stepper );
    m_rodRenderers.push_back( rodRenderer) ;
    
    cerr << "WmSweeneyRodManager::addRod: Created rod\n";
    
    return true;
}

void WmSweeneyRodManager::addCollisionMesh( BASim::TriangleMesh* i_triangleMesh,
                                            WmFigMeshController* i_scriptingController )
{
    m_triangleMeshes.push_back( i_triangleMesh );
    m_scriptingControllers.push_back( i_scriptingController );
}

void WmSweeneyRodManager::initialiseSimulation( const double i_timeStep, const double i_startTime )
{
    // FIXME: pass in timestep from Maya
    m_bridsonStepper = new BridsonStepper( m_rods, m_triangleMeshes, m_scriptingControllers, 
                                           m_rodTimeSteppers, i_timeStep, i_startTime );
                                           
}

void WmSweeneyRodManager::takeStep()
{
    // Check if anything has actually been initialised yet. We may still be being loaded by Maya.
    if ( m_bridsonStepper == NULL )
        return;
    
    // Force self collisions to be off whilst testing
    m_bridsonStepper->skipRodRodCollisions( true );

    // TEST: DO WE NEED TO DO THIS EVERY FRAME
   
  /*  for ( size_t r =0; r< m_rods.size(); ++r )
    {
        ElasticRod* rod = m_rods[ r ];
        
        // Lock the first two vertices of the rod
        RodBoundaryCondition* boundary = m_rodTimeSteppers[ r ]->getBoundaryCondition();
                    
        for ( int v=0; v<rod->nv(); ++v )
        {
            // Probably most of these were not fixed, this may be slow to do.
            // ??? Benchmark this and check what is going on.
            boundary->releaseVertex( v );
        }
        
        // Set the velocity to be zero as we're grooming static hair
        boundary->setDesiredVertexPosition( 0, m_rodTimeSteppers[ r ]->getTime(), rod->getVertex( 0 ), BASim::Vec3d( 0.0, 0.0, 0.0 ) );
        boundary->setDesiredVertexPosition( 1, m_rodTimeSteppers[ r ]->getTime(), rod->getVertex( 1 ), BASim::Vec3d( 0.0, 0.0, 0.0 ) );
    }*/

    // Ensure the rod stays stuck on the head
	//updateAllBoundaryConditions();                                   
    
    m_bridsonStepper->execute();
}

void WmSweeneyRodManager::drawAllRods()
{
    for ( size_t r = 0; r < m_rodRenderers.size(); ++r )
    {
        m_rodRenderers[ r ]->render();
    }
}
