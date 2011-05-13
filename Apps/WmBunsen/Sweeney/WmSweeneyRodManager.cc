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
                                  const double i_youngsModulus,
                                  const double i_shearModulus, const double i_viscosity, 
                                  const double i_density, const double i_radiusA, 
                                  const double i_radiusB, 
                                  const BASim::ElasticRod::RefFrameType i_referenceFrame,
                                  const double i_massDamping, 
                                  const BASim::Vec3d i_gravity,
                                  const BASim::RodTimeStepper::Method i_solverType )
{
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

    // We need a rod renderer to draw the rod in OpenGL
    RodRenderer* rodRenderer = new RodRenderer( *rod );
    
    // Create a timeStepper to simulate the rod forward in time
    RodTimeStepper* stepper = new RodTimeStepper( *rod );
	stepper->setDiffEqSolver( i_solverType );
	    
    // Add a damping force to the rod
    stepper->addExternalForce( new RodMassDamping( i_massDamping ) );
        
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
    
    // Store all the things we need to control the rod or add it to a BridsonStepper
    m_rods.push_back( rod );
    m_rodTimeSteppers.push_back( stepper );
    m_rodRenderers.push_back( rodRenderer) ;
    
    return true;
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
