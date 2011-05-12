#include "WmSweeneyRodManager.hh"

WmSweeneyRodManager::WmSweeneyRodManager()
{
    
}

WmSweeneyRodManager::~WmSweeneyRodManager()
{
    
}

bool wmSweeneyRodManager::addRod( const double i_youngsModulus,
                                  const double i_shearModulus,
                                  const double i_viscosity,
                                  const double i_density,
                                  const double i_radiusA,
                                  const double i_radiusB,
                                  const  i_referenceFrame,
                                  const double i_massDamping,
                                  const double i_gravity,
                                  const RodTimeStepper::Method i_solverType,
                                  const vector< BASim::Vec3d >& i_vertices )
{
    
    RodOptions rodOptions;
    rodOptions.YoungsModulus = 1000.0; /* megapascal */
    rodOptions.ShearModulus = 340.0;   /* megapascal */
    rodOptions.viscosity = 10.0;       /* poise */
    rodOptions.density = 1.3;          /* grams per cubic centimeter */
    rodOptions.radiusA = 0.05;         /* millimeter */
    rodOptions.radiusB = 0.05;         /* millimeter */
    rodOptions.refFrame = BASim::ElasticRod::TimeParallel;
    rodOptions.numVertices = 10;

    double massDamping = 10.0;
    BASim::Vec3d gravity( 0,0 -980.0, 0.0 );
    RodTimeStepper::Method solverType( RodTimeStepper::IMPL_EULER );
    
    bool doReverseHairdo = false;

    vector< BASim::Vec3d > vertices;
    vertices.resize( 10 );

    for ( unsigned int r = 0; r < 5; ++r )
    {
        for ( size_t v = 0; v < 10; ++v )
        {
            vertices[ v ] = BASim::Vec3d( r, 0.0, v );        
            int rodIndex = m_rodGroup.addRod( vertices, rodOptions, massDamping, 
                                              gravity, solverType, false, false );
                                
            // Lock the first edge in place
            m_rodGroup.addKinematicEdge( rodIndex, 0 );
        }
    }
    
    m_rods.clear();
    m_rodTimeSteppers.clear();
    m_triangleMeshes.clear();
    m_scriptingControllers.clear();  
    
    int numRods = m_rodGroup.numberOfRods();
    for ( int r=0; r<numRods; r++ )
    {        
        cerr << "Adding rod " << r << " to world\n";
        m_rods.push_back( m_rodGroup.elasticRod( r ) );
        m_world->addObject( m_rodGroup.elasticRod( r ) );
        
        cerr << "Adding rod time stepper " << r << " to world\n";
        m_rodTimeSteppers.push_back( m_rodGroup.stepper( r ) );
    }
    
    // FIXME: pass in timestep from Maya
    m_bridsonStepper = new BridsonStepper( m_rods, m_triangleMeshes, m_scriptingControllers, 
                                           m_rodTimeSteppers, (1.0/24.0/10.0), startTime);
                                           
    */
}
