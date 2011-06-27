#include "WmSweeneyNode.hh"
#include "WmFigConnectionNode.hh"
#include "../WmBunsenCollisionMeshNode.hh"

#include <maya/MFnMatrixAttribute.h>
#include <maya/MPlugArray.h>
#include <maya/MQuaternion.h>

#include <sys/stat.h>

#include "../../../BASim/src/Physics/ElasticRods/RodBendingForceSym.hh"
#include "../../../BASim/src/Physics/ElasticRods/RodTwistingForceSym.hh"

using namespace BASim;

// Required by Maya to identify the node
/* static */ MTypeId WmSweeneyNode::typeID ( 0x001135, 0xF6  );
/* static */ MString WmSweeneyNode::typeName( "wmSweeneyNode" );

// Input attributes
/* static */ MObject WmSweeneyNode::ia_time;
/* static */ MObject WmSweeneyNode::ia_startTime;

// Hair Property Attributes
/* static */ MObject WmSweeneyNode::ia_length;
/* static */ MObject WmSweeneyNode::ia_edgeLength;
/* static */ MObject WmSweeneyNode::ia_verticesPerRod;
/* static */ MObject WmSweeneyNode::ia_rodRadius;
/* static */ MObject WmSweeneyNode::ia_rodAspectRatio;
/* static */ MObject WmSweeneyNode::ia_rodRotation;
/* static */ MObject WmSweeneyNode::ia_curlRadius;
/* static */ MObject WmSweeneyNode::ia_curlPitch;
/* static */ MObject WmSweeneyNode::ia_curlStart;
/* static */ MObject WmSweeneyNode::ia_rodPitch;
/* static */ MObject WmSweeneyNode::ia_rodDamping;

// Barbershop specific inputs
/*static*/ MObject WmSweeneyNode::ia_strandVertices;
/*static*/ MObject WmSweeneyNode::ia_strandRootFrames;
/*static*/ MObject WmSweeneyNode::ia_verticesPerStrand;

// Output to the Barbershop guide curve deformer
/*static*/ MObject WmSweeneyNode::oa_simulatedNurbs;

// Sync attributes
/* static */ MObject WmSweeneyNode::ca_rodPropertiesSync;

// Collision meshes
/* static */ MObject WmSweeneyNode::ia_collisionMeshes;

//Solver Tolerances
/* static */ MObject WmSweeneyNode::ia_stol;
/* static */ MObject WmSweeneyNode::ia_atol;
/* static */ MObject WmSweeneyNode::ia_rtol;
/* static */ MObject WmSweeneyNode::ia_inftol;
/* static */ MObject WmSweeneyNode::ia_numLineSearchIters;

// Performance Tuning
//GeneralParameters
/* static */ MObject WmSweeneyNode::ia_enablePenaltyResponse;
/* static */ MObject WmSweeneyNode::ia_implicitThickness;
/* static */ MObject WmSweeneyNode::ia_implicitStiffness;
/* static */ MObject WmSweeneyNode::ia_inextensibilityThreshold;

//Failuredetection
/* static */ MObject WmSweeneyNode::ia_maxNumOfSolverIters;
/* static */ MObject WmSweeneyNode::ia_maxNumOfCollisionIters;
/* static */ MObject WmSweeneyNode::ia_enableExplosionDetection;
/* static */ MObject WmSweeneyNode::ia_explosionDampening;
/* static */ MObject WmSweeneyNode::ia_explosionThreshold;
/* static */ MObject WmSweeneyNode::ia_stretchingThreshold;

//FailureResponse
/* static */ MObject WmSweeneyNode::ia_solverFailure;
/* static */ MObject WmSweeneyNode::ia_collisionFailure;
/* static */ MObject WmSweeneyNode::ia_explosionFailure;
/* static */ MObject WmSweeneyNode::ia_stretchingFailure;
/* static */ MObject WmSweeneyNode::ia_maxNumSolverSubsteps;
/* static */ MObject WmSweeneyNode::ia_maxNumCollisionSubsteps;
/* static */ MObject WmSweeneyNode::ia_maxNumExplosionSubsteps;
/* static */ MObject WmSweeneyNode::ia_maxNumStretchingSubsteps;

// Debug drawing
/* static */ MObject WmSweeneyNode::ia_shouldDrawVelocity;

WmSweeneyNode::WmSweeneyNode() : m_rodManager( NULL )
{
}

WmSweeneyNode::~WmSweeneyNode()
{
}

WmSweeneyRodManager* WmSweeneyNode::rodManager()
{
    return m_rodManager;
}

MStatus WmSweeneyNode::compute( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
	MStatus status;

	if ( i_plug == ca_rodPropertiesSync )
	{
	    m_currentTime = i_dataBlock.inputValue( ia_time ).asTime().value();
		m_startTime = i_dataBlock.inputValue( ia_startTime ).asDouble();

		// Hair properties
		m_length = i_dataBlock.inputValue( ia_length ).asDouble();
		m_edgeLength = i_dataBlock.inputValue( ia_edgeLength ).asDouble();
	    m_rodRadius = i_dataBlock.inputValue( ia_rodRadius ).asDouble();
		m_rodAspectRatio = i_dataBlock.inputValue( ia_rodAspectRatio ).asDouble();
		m_rodRotation = i_dataBlock.inputValue( ia_rodRotation ).asDouble();
	    m_curlRadius = i_dataBlock.inputValue( ia_curlRadius ).asDouble();
		m_curlPitch = i_dataBlock.inputValue( ia_curlPitch ).asDouble();
		m_curlStart = i_dataBlock.inputValue( ia_curlStart ).asDouble();
		m_rodPitch = i_dataBlock.inputValue( ia_rodPitch ).asDouble();
		m_rodDamping = i_dataBlock.inputValue( ia_rodDamping ).asBool();

        bool shouldDrawVelocity = i_dataBlock.inputValue( ia_shouldDrawVelocity ).asBool();
        if ( m_rodManager != NULL )
        {
            m_rodManager->setRodsDrawVelocity( shouldDrawVelocity );
        }     

		MObject strandVerticesObj = i_dataBlock.inputValue( ia_strandVertices ).data();
		MFnVectorArrayData strandVerticesArrayData( strandVerticesObj, &status );
		CHECK_MSTATUS( status );

		MVectorArray strandVertices = strandVerticesArrayData.array( &status );
		CHECK_MSTATUS( status );

		MDataHandle strandRootFramesHandle = i_dataBlock.inputValue( ia_strandRootFrames, & status );
		CHECK_MSTATUS( status );
		MObject strandRootFramesObj = strandRootFramesHandle.data();
		MFnVectorArrayData strandRootFramesArrayData( strandRootFramesObj, &status );
		CHECK_MSTATUS( status );

		MVectorArray strandRootFrames = strandRootFramesArrayData.array( &status );
		CHECK_MSTATUS( status );

		//for ( unsigned int i = 0 ; i < strandRootFrames.length() ; i++ )
		//	cout << "ROD FRAME # " << i << " " << strandRootFrames[ i ].x << " " << strandRootFrames[ i ].y << " " << strandRootFrames[ i ].z << endl;

		int verticesPerRod = i_dataBlock.inputValue( ia_verticesPerRod ).asInt();

		int numberOfVerticesPerStrand = i_dataBlock.inputValue( ia_verticesPerStrand ).asInt();

		if ( m_currentTime == m_startTime )
		{
			// We can't use the assignment operator because strandVertices is actually
			// a reference to the MVectorArrayData array and it will go out of scope
			// and we'll be left with a reference to nothing and obviously a crash.
			status = m_strandVertices.copy( strandVertices );
			CHECK_MSTATUS( status );

			status = m_strandRootFrames.copy( strandRootFrames );
			CHECK_MSTATUS( status );

			m_numberOfVerticesPerStrand = numberOfVerticesPerStrand;

			m_verticesPerRod = verticesPerRod;

			initialiseRodFromBarberShopInput( i_dataBlock );
		}
		else
		{
		    if ( m_rodManager != NULL )
		    {
		    	bool update_rod = false;
				// If the rods exist then just update their undeformed configuration but keep running
				// the simulation.

				// Apply kinetic damping
				m_rodManager->setUseKineticDamping( m_rodDamping );

				 for (size_t i = 0; i < m_rodManager->m_rods.size(); ++i)
				 {
				     // TODO: Add code to update undeformed configuration
				     // for just now, recreate the rod
				     // initialiseRodFromBarberShopInput( i_dataBlock );

				     // get total rod length for scaling
				     Scalar curl_len = 0;
				     int curl_resolution = 1;
				     // TODO (sainsley) : create method inside ElasticRod that returns a total length
				     // this can be computed within ElasticRod::computeEdgeLengths and stored locally
				     // so we do not need to recompute it here

				    // Compute new edge length
				    Scalar updated_edge_length = m_length / m_verticesPerRod;
				    std::vector<Scalar> rest_lengths;

					// Adjust edge lengths and compute the resulting curl length
				    for ( ElasticRod::edge_iter eh = m_rodManager->m_rods[i]->edges_begin();
				          eh != m_rodManager->m_rods[i]->edges_end(); ++eh )
				    {
				    	// TODO: remove this check for the fixed edge--it doesn't matter
				    	if ( updated_edge_length != m_rodManager->m_rods[i]->getEdgeLength( *eh ))
				    	{
				    		update_rod = true;
				    	}
						rest_lengths.push_back( updated_edge_length );
						//cout << "WmSweeneyNode::compute::edges: idx = " << eh->idx() << " new edge length = " << updated_edge_length << " current edge length " << m_rodManager->m_rods[i]->getEdgeLength(  *eh  )  << endl;
						if ( eh->idx() >= m_curlStart*( m_verticesPerRod - 1 ) )
						{
							curl_len += m_rodManager->m_rods[i]->getEdgeLength( *eh );
						        curl_resolution++;
						}
					}
				    assert( m_curlStart == 1.0 || curl_len != 0 );

				    // pass remaining lengths to the forces
				    m_rodManager->m_rods[i]->setRestLengths( rest_lengths );

				    // parametric variable for walking the rod lengthwise
				    Scalar t = 0;
				    int j = 0;

					// adjust the rod size
					Scalar radius_a = m_rodRadius;
					Scalar radius_b = radius_a;
					// apply apsect ratio : flip axis if aspect ratio is less than 1 to preserve radius scale
					if ( m_rodAspectRatio > 1.0 )
					{
					    radius_b *= m_rodAspectRatio;
					}
					else
					{
						radius_a *= 1.0/m_rodAspectRatio;
					}



					m_rodManager->m_rods[i]->setRadius( radius_a, radius_b );

					// TODO(sainsley) : strip all of this logic out of ElasticRod and RodBendingForce
					// m_rodManager->m_rods[i]->setBaseRotation( m_rodRotation*M_PI );

					// set initial rotation
					m_rodManager->m_rods[i]->setTheta( 0, m_rodRotation * M_PI );
					Scalar c = cos( m_rodManager->m_rods[i]->getTheta( 0 ) );
					Scalar s = sin( m_rodManager->m_rods[i]->getTheta( 0 ) );
					const Vec3d& u = m_rodManager->m_rods[i]->getReferenceDirector1( 0 );
					const Vec3d& v = m_rodManager->m_rods[i]->getReferenceDirector2( 0 );
					m_rodManager->m_rods[i]->setMaterial1( 0,  c * u + s * v );
					m_rodManager->m_rods[i]->setMaterial2( 0, -s * u + c * v );

					m_rodManager->m_rods[i]->updateStiffness();

					for ( ElasticRod::vertex_iter vh = m_rodManager->m_rods[i]->vertices_begin();
	                          vh != m_rodManager->m_rods[i]->vertices_end(); ++vh )
	                {

						// curl curvature and torsion

	                    assert( m_rodManager->m_rods[i]->m_bendingForce != NULL );

	                    // curl curvature and torsion

	                    Scalar curvature = 0.0;
	                    Scalar torsion = 0.0;

	                    if ( t > 0 )
	                    {
	                    	curvature = m_curlRadius * curl_len/curl_resolution;
	                    	torsion = m_curlPitch * M_PI * curl_len/curl_resolution;
	                    }

	                    m_rodManager->m_rods[i]->m_twistingForce->setUndeformedTwist( *vh, torsion );
	                    m_rodManager->m_rods[i]->m_bendingForce->setKappaBar( *vh, Vec2d( curvature, 0 ) );
	                    //m_rodManager->m_rods[i]->m_bendingForce->setKappaBar(
						//	*vh, Vec2d(  curvature * cos( torsion * t ),
						//		curvature * sin( torsion * t ) ) );

						//cout << "WmSweeneyNode::compute::simulate: idx = " << vh->idx() << " parametric var = " << t << " curvature " <<  m_rodManager->m_rods[i]->m_bendingForce->getKappaBar( *vh ) << " bending stiffness " <<  m_rodManager->m_rods[i]->m_bendingForce->getB( *vh ) << " vertex mass " << m_rodManager->m_rods[i]->getVertexMass( vh->idx() ) << endl;

						if ( vh->idx() >= m_curlStart*(m_verticesPerRod)  && vh != m_rodManager->m_rods[i]->vertices_end() )
						{
							t += m_rodManager->m_rods[i]->getEdgeLength( j++ )/curl_len;
						}
	                }

				 }
				 updateCollisionMeshes( i_dataBlock );
				 m_rodManager->takeStep();


				 double atol = powf(10, -i_dataBlock.inputValue( ia_atol).asDouble());
			 	 double stol = powf(10, -i_dataBlock.inputValue( ia_stol).asDouble());
				 double rtol  = powf(10, -i_dataBlock.inputValue( ia_rtol).asDouble());
				 double inftol  = powf(10, -i_dataBlock.inputValue( ia_inftol).asDouble());
				 int numLineSearchIters = i_dataBlock.inputValue( ia_numLineSearchIters).asInt();

				 double stiffness = i_dataBlock.inputValue( ia_implicitStiffness).asDouble();
				 m_rodManager->updateSolverSettings( atol, stol, rtol, inftol, numLineSearchIters, stiffness );
			}
		}
		i_dataBlock.setClean( i_plug );
	}
	else if ( i_plug == oa_simulatedNurbs )
	{
		compute_oa_simulatedNurbs( i_plug, i_dataBlock );
	}
	else
	{
		return MS::kUnknownParameter;
	}
	return MS::kSuccess;
}

void WmSweeneyNode::initialiseCollisionMeshes( MDataBlock &i_data )
{
    MStatus status;

    MArrayDataHandle inArrayH = i_data.inputArrayValue( ia_collisionMeshes, &status );
    CHECK_MSTATUS( status );
    size_t numMeshesConnected = inArrayH.elementCount();

    for ( unsigned int i=0; i < numMeshesConnected; i++ )
    {
        // Even if we don't use it, grab the data so Maya knows to evaluate the node
        inArrayH.jumpToElement(i);
        MDataHandle collisionMeshH = inArrayH.inputValue( &status );
        CHECK_MSTATUS( status );

        MPlug plug( thisMObject(), ia_collisionMeshes );
        CHECK_MSTATUS( status );
        if ( plug.isArray( &status ) )
        {
            MPlug indxPlug = plug.elementByLogicalIndex( i, &status );
            CHECK_MSTATUS( status );
            if ( indxPlug.isConnected( &status ) )
            {
                MPlugArray inPlugArr;
                indxPlug.connectedTo( inPlugArr, true, false, &status );
                CHECK_MSTATUS( status );

                // Since we asked for the destination there can only be one plug in the array
                MPlug meshPlug = inPlugArr[0];
                MObject collisionMeshNodeObj = meshPlug.node( &status );
                CHECK_MSTATUS( status );
                MFnDependencyNode collisionMeshNodeFn( collisionMeshNodeObj );
                WmBunsenCollisionMeshNode* collisionMeshNode = (WmBunsenCollisionMeshNode*)collisionMeshNodeFn.userNode();

                TriangleMesh* triangleMesh = NULL;
                WmFigMeshController* figMeshController = NULL;

                collisionMeshNode->initialise( NULL, i, &triangleMesh, &figMeshController );

                // Now add the mesh to the rod manager
                m_rodManager->addCollisionMesh( triangleMesh, figMeshController->currentLevelSet(),
                                                figMeshController );
            }
            else
            {
                CHECK_MSTATUS( status );
            }
        }
    }
}

void WmSweeneyNode::updateCollisionMeshes( MDataBlock& i_dataBlock )
{
    MStatus status;

    MArrayDataHandle inArrayH = i_dataBlock.inputArrayValue( ia_collisionMeshes, &status );
    CHECK_MSTATUS( status );
    size_t numMeshesConnected = inArrayH.elementCount();

    for ( unsigned int i=0; i < numMeshesConnected; i++ )
    {
        // All we need to do is ask Maya for the data and it will pull the attr,
        // causing a compute in the collision mesh node which will directly
        // update the collision mesh data in the RodManager.
        inArrayH.jumpToElement(i);
        MDataHandle collisionMeshH = inArrayH.inputValue( &status );
        CHECK_MSTATUS( status );
    }
}

void WmSweeneyNode::initialiseRodFromBarberShopInput( MDataBlock& i_dataBlock )
{
    // We need to be able to run a sim from pre-groomed barbershop strands too. That seems
    // pretty simple as we just take the vertices from the barbershop input and create the
    // rod from that. Then when taking a time step we set the first edge to match the barbershop
    // input. It looks like Sweeney code easily work for grooming or for dynamic sims.

    cerr << "initialiseRodFromBarberShopInput() - About to create rods from Barbershop input\n";

    // Reset the manager and remove all rods before adding more
    delete m_rodManager;

    m_rodManager = new WmSweeneyRodManager();

    cerr << "initialiseRodFromBarberShopInput() - Deleted and created a new WmSweeneyRodManager\n";

    if ( m_strandVertices.length() == 0 )
    {
        cerr << "initialiseRodFromBarberShopInput() - no input strands so can't create any rods";
        return;
    }

    // First, get all the collision mesh data organised
    initialiseCollisionMeshes( i_dataBlock );

    // Create one rod for each barbershop strand. Ignore the strand shape or length but do
    // take its initial direction as a follicle angle
    unsigned int currentVertexIndex = 0;
    unsigned int numberOfStrands = m_strandVertices.length() / m_numberOfVerticesPerStrand;


    cerr << "initialiseRodFromBarberShopInput() - m_strandVertices.length() = " << m_strandVertices.length() << endl;
    cerr << "initialiseRodFromBarberShopInput() - number of barberShop strands = " << numberOfStrands << endl;
    cerr << "initialiseRodFromBarberShopInput() - number of vertices per barberShop strand = " << m_numberOfVerticesPerStrand<< endl;

    vector< BASim::Vec3d > vertices;

    for ( unsigned int inputStrandNumber = 0; inputStrandNumber < numberOfStrands; ++inputStrandNumber )
    {
    	MVector direction = m_strandVertices[ currentVertexIndex + 1 ]
                                - m_strandVertices[ currentVertexIndex ];
        direction.normalize();

        constructRodVertices( vertices, direction, m_strandVertices[ currentVertexIndex ] );

        cout << "initialiseRodFromBarberShopInput() - check for root frames for " << inputStrandNumber << endl;

        if ( m_strandRootFrames.length() != 0 )
        {

        	BASim::Vec3d m1 = Vec3d(  m_strandRootFrames[ 3*inputStrandNumber ].x,
       				m_strandRootFrames[ 3*inputStrandNumber ].y,
       				m_strandRootFrames[ 3*inputStrandNumber ].z  );

        	m1.normalize();

        	m_rodManager->addRod( vertices, m_startTime, m1 );

        }
        else
        {
        	m_rodManager->addRod( vertices, m_startTime );
        }

		cerr << "Creating rod at time " << m_startTime << endl;

		currentVertexIndex += m_numberOfVerticesPerStrand;
	}

    //Set performancetuningparameters to pass through
    PerformanceTuningParameters perfParams;
    perfParams.m_enable_penalty_response=i_dataBlock.inputValue( ia_enablePenaltyResponse).asBool();
    perfParams.m_implicit_thickness=i_dataBlock.inputValue( ia_implicitThickness).asDouble();
    perfParams.m_implicit_stiffness=i_dataBlock.inputValue( ia_implicitStiffness).asDouble();
    perfParams.m_inextensibility_threshold=i_dataBlock.inputValue( ia_inextensibilityThreshold).asInt();
    perfParams.m_solver.m_max_iterations=i_dataBlock.inputValue( ia_maxNumOfSolverIters).asInt();
    perfParams.m_collision.m_max_iterations=i_dataBlock.inputValue( ia_maxNumOfCollisionIters).asInt();
    perfParams.m_enable_explosion_detection=i_dataBlock.inputValue( ia_enableExplosionDetection).asBool();
    perfParams.m_explosion_damping=i_dataBlock.inputValue( ia_explosionDampening).asDouble();
    perfParams.m_explosion_threshold=i_dataBlock.inputValue( ia_explosionThreshold).asDouble();
    perfParams.m_stretching_threshold=i_dataBlock.inputValue( ia_stretchingThreshold).asDouble();
    perfParams.m_solver.m_in_case_of= (BASim::FailureMode::ResponseSeverity) i_dataBlock.inputValue( ia_solverFailure).asInt();
    perfParams.m_collision.m_in_case_of=(BASim::FailureMode::ResponseSeverity) i_dataBlock.inputValue( ia_collisionFailure).asInt();
    perfParams.m_explosion.m_in_case_of=(BASim::FailureMode::ResponseSeverity) i_dataBlock.inputValue( ia_explosionFailure).asInt();
    perfParams.m_stretching.m_in_case_of=(BASim::FailureMode::ResponseSeverity) i_dataBlock.inputValue( ia_stretchingFailure).asInt();
    perfParams.m_solver.m_max_substeps=i_dataBlock.inputValue( ia_maxNumSolverSubsteps).asInt();
    perfParams.m_collision.m_max_substeps=i_dataBlock.inputValue( ia_maxNumCollisionSubsteps).asInt();
    perfParams.m_explosion.m_max_substeps=i_dataBlock.inputValue( ia_maxNumExplosionSubsteps).asInt();
    perfParams.m_stretching.m_max_substeps=i_dataBlock.inputValue( ia_maxNumStretchingSubsteps).asInt();

	double m_atol = powf( 10, -i_dataBlock.inputValue( ia_atol ).asDouble() );
	double m_stol = powf( 10, -i_dataBlock.inputValue( ia_stol ).asDouble() );
	double m_rtol = powf( 10, -i_dataBlock.inputValue( ia_rtol ).asDouble() );
	double m_inftol = powf( 10, -i_dataBlock.inputValue( ia_inftol ).asDouble() );
	int  m_numLineSearchIters=i_dataBlock.inputValue( ia_numLineSearchIters).asInt();

	cerr << "initialiseRodFromBarberShopInput() - About to initialise simulation\n";
	    m_rodManager->initialiseSimulation( 1 / 24.0, m_startTime, perfParams, m_atol, m_stol, m_rtol, m_inftol,
	                                        m_numLineSearchIters );
	cerr << "initialiseRodFromBarberShopInput() - Simulation initialised at time " << m_startTime << endl;
}

void WmSweeneyNode::constructRodVertices( vector< BASim::Vec3d >& o_rodVertices, const MVector& i_direction,
                                  const MVector& i_rootPosition )
{
    // Construct a list of vertices for a rod with its root at i_rootPosition and going in direction
    // i_direction

    o_rodVertices.clear();

    MVector edge = i_direction;
    edge.normalize();

	cerr << "constructRodVertices(): m_length = " << m_length << endl;
	cerr << "constructRodVertices(): m_verticesPerRod = " << m_verticesPerRod << endl;
	cerr << "constructRodVertices(): i_direction = " << i_direction << endl;

	edge *= m_length / m_verticesPerRod;

	cerr << "constructRodVertices(): edgeLength = " << edge.length() << "\n";

	MVector currentVertex( i_rootPosition );

	//Scalar a = 0.0;
	//if ( m_curlRadius != 0.0 )
	//	a = 1.0 / (m_curlRadius * m_length / m_verticesPerRod);
	//	Scalar b = 0.5 * m_length / m_verticesPerRod;

	for ( int v = 0; v < m_verticesPerRod; ++v )
    {
        // Straight rods as twist is controlled by the rod properties

        o_rodVertices.push_back( BASim::Vec3d( currentVertex.x, currentVertex.y, currentVertex.z ) );

        currentVertex += edge;

        //MVector newPoint( a * cos( (double)v ),
        	//		b * (double)v, a * sin( (double)v ) );
        //

        // For testing, force a straight rod
        // MVector newPoint( 0.0, v, 0.0 );

        // The helix is created with the y-axis as the centre, rotate it
        // so that it has i_direction as the centre
        //MQuaternion rotationQ( MVector( 0.0, 1.0, 0.0 ), i_direction );
        //newPoint = newPoint.rotateBy( rotationQ );

        // Now move the point to sit where the Barbershop input strand comes from
        //newPoint += i_rootPosition;

        //o_rodVertices.push_back( BASim::Vec3d( newPoint.x, newPoint.y, newPoint.z ) );
	}

	cerr << "constructRodVertices(): Finished constructing rod vertices\n";
}

void WmSweeneyNode::compute_oa_simulatedNurbs( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    MStatus status;

    cerr << "compute_oa_simulatedNurbs()\n";

    // First pull all the inputs to make sure we're up to date.
    i_dataBlock.inputValue( ca_rodPropertiesSync, &status ).asBool();
    CHECK_MSTATUS( status );

    // The above may have been clean so just make sure we actually read time
    i_dataBlock.inputValue( ia_time, &status ).asTime().value();
    CHECK_MSTATUS( status );

    MArrayDataHandle simulatedNurbsArrayHandle = i_dataBlock.outputArrayValue( oa_simulatedNurbs, &status );
    CHECK_MSTATUS( status );

    MArrayDataBuilder simulatedNurbsArrayDataBuilder( & i_dataBlock, oa_simulatedNurbs,
                                                      (unsigned int)m_rodManager->numberOfRods(), & status);
    CHECK_MSTATUS( status );

    for ( size_t r = 0; r < m_rodManager->numberOfRods(); ++r )
    {
        MPointArray nurbsEditPoints;

        ElasticRod* rod = m_rodManager->rod( r );

        for ( unsigned int v = 0; v < rod->nv(); v++ )
        {
            BASim::Vec3d vertex = rod->getVertex( v );
            MPoint nurbsEditPoint( vertex[ 0 ], vertex[ 1 ], vertex[ 2 ] );
            nurbsEditPoints.append( nurbsEditPoint );
        }

        MFnNurbsCurveData nurbsDataFn;
        MObject nurbsDataObj = nurbsDataFn.create();
        MFnNurbsCurve nurbsFn;
        MObject nurbsObj = nurbsFn.createWithEditPoints( nurbsEditPoints, 1, MFnNurbsCurve::kOpen,
            false /*not 2d*/, false /*not rational*/, true /*uniform params*/, nurbsDataObj, & status );
        CHECK_MSTATUS( status );

        MDataHandle simulatedNurbsHandle = simulatedNurbsArrayDataBuilder.addElement( (unsigned int) r, & status );
        CHECK_MSTATUS( status );

        status = simulatedNurbsHandle.set( nurbsDataObj );
        CHECK_MSTATUS( status );
    }

    simulatedNurbsArrayHandle.set( simulatedNurbsArrayDataBuilder );
    simulatedNurbsArrayHandle.setAllClean();

    i_dataBlock.setClean( i_plug );
}

void WmSweeneyNode::draw( M3dView& i_view, const MDagPath& i_path,
                            M3dView::DisplayStyle i_style,
                            M3dView::DisplayStatus i_status )
{
    MStatus status;

    // Pull on the sync plugs to cause compute() to be called if any
    // of the rod properties or time has changed.
    double d;

    MPlug propertiesSyncPlug( thisMObject(), ca_rodPropertiesSync );
    propertiesSyncPlug.getValue( d );

    i_view.beginGL();
    //glPushAttrib( GL_CURRENT_BIT | GL_POINT_BIT | GL_LINE_BIT );
	glPushAttrib( GL_ALL_ATTRIB_BITS );

	// draw dynamic Hair
    if ( m_rodManager != NULL )
    {
        m_rodManager->drawAllRods();
    }

	glPopAttrib();

	i_view.endGL();
}

MStatus WmSweeneyNode::connectionMade( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc )
{
    MStatus stat;
    MStatus retVal(MS::kUnknownParameter );

    return retVal;
}

MStatus WmSweeneyNode::connectionBroken( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc )
{
    MStatus stat;
    MStatus retVal( MS::kUnknownParameter );

    return retVal;
}

bool WmSweeneyNode::isBounded() const
{
	return false;
}

void* WmSweeneyNode::creator()
{
	return new WmSweeneyNode();
}

/*static */ MStatus WmSweeneyNode::addNumericAttribute( MObject& i_attribute, MString i_longName,
    MString i_shortName, MFnNumericData::Type i_type, double i_defaultValue, bool i_isInput,
    bool i_isArray )
{
    // Creates a numeric attribute with default attributes
    MStatus stat = MS::kSuccess;

    MFnNumericAttribute nAttr;
    i_attribute = nAttr.create( i_longName, i_shortName, i_type, i_defaultValue, &stat);
    if ( !stat )
    {
        cerr << "Failed to create attribute " << i_longName << endl;
        return stat;
    }
    if ( i_isInput )
    {
        nAttr.setWritable( true );
    }
    else
    {
        nAttr.setWritable( false );
    }

    if ( i_isArray )
    {
        nAttr.setArray( true );
    }

    stat = addAttribute( i_attribute );
    if ( !stat ) { stat.perror( "addAttribute " + i_longName ); return stat; }

    return stat;
}

/* static */ MStatus WmSweeneyNode::initialize()
{
    MStatus status;

    addNumericAttribute( ca_rodPropertiesSync, "rodPropertiesSync", "rps", MFnNumericData::kBoolean, false, false );

    {
        MFnUnitAttribute uAttr;
        ia_time = uAttr.create( "time", "t", MTime( 0.0 ), &status );
        if ( !status)
        {
            status.perror("create ia_time attribute");
            return status;
        }
        CHECK_MSTATUS( uAttr.setWritable(true) );
        CHECK_MSTATUS( uAttr.setConnectable(true) );
        CHECK_MSTATUS( uAttr.setStorable(false) );
        status = addAttribute( ia_time );
        if ( !status ) { status.perror( "addAttribute ia_time" ); return status; }
    }
    status = attributeAffects( ia_time, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_time->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_startTime, "startTime", "stt", MFnNumericData::kDouble, 1.0, true );
	status = attributeAffects( ia_startTime, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_startTime->ca_rodPropertiesSync" ); return status; }

	//addNumericAttribute( ia_length, "length", "len", MFnNumericData::kDouble, 10.0, true );
	//status = attributeAffects( ia_length, ca_rodPropertiesSync );
	//if ( !status ) { status.perror( "attributeAffects ia_length->ca_rodPropertiesSync" ); return status; }

    {
        MFnNumericAttribute numericAttr;
        ia_length = numericAttr.create( "length", "len", MFnNumericData::kDouble, 10.0, &status );
        CHECK_MSTATUS( status );
        CHECK_MSTATUS( numericAttr.setReadable( true ) );
        CHECK_MSTATUS( numericAttr.setWritable( true ) );
        CHECK_MSTATUS( numericAttr.setMin( 1.0 ) );
        CHECK_MSTATUS( numericAttr.setMax( 100.0 ) );
        status = addAttribute( ia_length );
        CHECK_MSTATUS( status );

        status = attributeAffects( ia_length, ca_rodPropertiesSync );
        if ( !status ) { status.perror( "attributeAffects ia_length->ca_rodPropertiesSync" ); return status; }
    }

        // TODO : remove this? i don't think we are using it
	addNumericAttribute( ia_edgeLength, "edgeLength", "ele", MFnNumericData::kDouble, 1.0, true );
	status = attributeAffects( ia_edgeLength, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_edgeLength->ca_rodPropertiesSync" ); return status; }

	//addNumericAttribute( ia_rodRadius, "rodRadius", "ror", MFnNumericData::kDouble, 0.0, true );
    {
        MFnNumericAttribute numericAttr;
        ia_rodRadius = numericAttr.create( "rodRadius", "ror", MFnNumericData::kDouble, 0.005, &status );
        CHECK_MSTATUS( status );
        CHECK_MSTATUS( numericAttr.setReadable( true ) );
        CHECK_MSTATUS( numericAttr.setWritable( true ) );
        CHECK_MSTATUS( numericAttr.setMin( 0.001 ) );
        CHECK_MSTATUS( numericAttr.setMax( 0.1 ) );
        status = addAttribute( ia_rodRadius );
        CHECK_MSTATUS( status );

        status = attributeAffects( ia_rodRadius, ca_rodPropertiesSync );
        if ( !status ) { status.perror( "attributeAffects ia_rodRadius->ca_rodPropertiesSync" ); return status; }
    }

    {
        MFnNumericAttribute numericAttr;
        ia_rodAspectRatio = numericAttr.create( "rodAspectRatio", "roar", MFnNumericData::kDouble, 1.0, &status );
        CHECK_MSTATUS( status );
        CHECK_MSTATUS( numericAttr.setReadable( true ) );
        CHECK_MSTATUS( numericAttr.setWritable( true ) );
        CHECK_MSTATUS( numericAttr.setMin( 0.1 ) );
        CHECK_MSTATUS( numericAttr.setMax( 10.0 ) );
        status = addAttribute( ia_rodAspectRatio );
        CHECK_MSTATUS( status );

        status = attributeAffects( ia_rodAspectRatio, ca_rodPropertiesSync );
        if ( !status ) { status.perror( "attributeAffects ia_rodAspectRatio->ca_rodPropertiesSync" ); return status; }
    }

    {
        MFnNumericAttribute numericAttr;
        ia_rodRotation = numericAttr.create( "rodRotation", "rorot", MFnNumericData::kDouble, 0.0, &status );
        CHECK_MSTATUS( status );
        CHECK_MSTATUS( numericAttr.setReadable( true ) );
        CHECK_MSTATUS( numericAttr.setWritable( true ) );
        CHECK_MSTATUS( numericAttr.setMin( -1.0 ) );
        CHECK_MSTATUS( numericAttr.setMax( 1.0 ) );
        status = addAttribute( ia_rodRotation );
        CHECK_MSTATUS( status );

        status = attributeAffects( ia_rodRotation, ca_rodPropertiesSync );
        if ( !status ) { status.perror( "attributeAffects ia_rodRotation->ca_rodPropertiesSync" ); return status; }
    }

    {
        MFnNumericAttribute numericAttr;
        ia_curlRadius = numericAttr.create( "curlTightness", "crlrad", MFnNumericData::kDouble, 0.0, &status );
        CHECK_MSTATUS( status );
        CHECK_MSTATUS( numericAttr.setReadable( true ) );
        CHECK_MSTATUS( numericAttr.setWritable( true ) );
        CHECK_MSTATUS( numericAttr.setMin( -2.0 ) );
        CHECK_MSTATUS( numericAttr.setMax( 2.0 ) );
        status = addAttribute( ia_curlRadius );
        CHECK_MSTATUS( status );

        status = attributeAffects( ia_curlRadius, ca_rodPropertiesSync );
        if ( !status ) { status.perror( "attributeAffects ia_curlRadius->ca_rodPropertiesSync" ); return status; }
    }

    {
    	MFnNumericAttribute numericAttr;
        ia_curlPitch = numericAttr.create( "curlSpacing", "crlptch", MFnNumericData::kDouble, 0.0, &status );
        CHECK_MSTATUS( status );
        CHECK_MSTATUS( numericAttr.setReadable( true ) );
        CHECK_MSTATUS( numericAttr.setWritable( true ) );
        CHECK_MSTATUS( numericAttr.setMin( 0.0 ) );
        CHECK_MSTATUS( numericAttr.setMax( 5.0 ) );
        status = addAttribute( ia_curlPitch );
        CHECK_MSTATUS( status );

        status = attributeAffects( ia_curlPitch, ca_rodPropertiesSync );
        if ( !status ) { status.perror( "attributeAffects ia_curlPitch->ca_rodPropertiesSync" ); return status; }
    }

    {
    	MFnNumericAttribute numericAttr;
        ia_curlStart = numericAttr.create( "curlStart", "crlstrt", MFnNumericData::kDouble, 0.0, &status );
        CHECK_MSTATUS( status );
        CHECK_MSTATUS( numericAttr.setReadable( true ) );
        CHECK_MSTATUS( numericAttr.setWritable( true ) );
        CHECK_MSTATUS( numericAttr.setMin( 0.0 ) );
        CHECK_MSTATUS( numericAttr.setMax( 1.0 ) );
        status = addAttribute( ia_curlStart );
        CHECK_MSTATUS( status );

        status = attributeAffects( ia_curlStart, ca_rodPropertiesSync );
        if ( !status ) { status.perror( "attributeAffects ia_curlStart->ca_rodPropertiesSync" ); return status; }
    }

    addNumericAttribute( ia_rodDamping, "rodDamping", "roddamp", MFnNumericData::kBoolean, true, true );
    	status = attributeAffects( ia_rodDamping, ca_rodPropertiesSync );
    if ( !status ) { status.perror( "attributeAffects ia_rodDamping->ca_rodPropertiesSync" ); return status; }


    addNumericAttribute( ia_rodPitch, "rodPitch", "rop", MFnNumericData::kDouble, 0.5, true );
    status = attributeAffects( ia_rodPitch, ca_rodPropertiesSync );
    if ( !status ) { status.perror( "attributeAffects ia_rodPitch->ca_rodPropertiesSync" ); return status; }

    {
    	MFnTypedAttribute tAttr;
    	ia_strandVertices = tAttr.create( "strandVertices", "sve",
                                          MFnData::kVectorArray, & status );
    	CHECK_MSTATUS( status );
    	CHECK_MSTATUS( tAttr.setReadable( false ) );
    	CHECK_MSTATUS( tAttr.setWritable( true ) );
    	CHECK_MSTATUS( tAttr.setStorable( false ) );
    	status = addAttribute( ia_strandVertices );
    	CHECK_MSTATUS( status );
    }
    status = attributeAffects( ia_strandVertices, ca_rodPropertiesSync );
    if ( !status ) { status.perror( "attributeAffects ia_strandVertices->ca_rodPropertiesSync" ); return status; }

	{
		MFnTypedAttribute tAttr;
		ia_strandRootFrames = tAttr.create( "strandRootFrames", "srf",
											MFnData::kVectorArray, &status );
		CHECK_MSTATUS( status );
		CHECK_MSTATUS( tAttr.setReadable( false ) );
		CHECK_MSTATUS( tAttr.setWritable( true ) );
		CHECK_MSTATUS( tAttr.setStorable( false ) );
		status = addAttribute( ia_strandRootFrames );
		CHECK_MSTATUS( status );
	}
	status = attributeAffects( ia_strandRootFrames, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_strandRootFrames->ca_rodPropertiesSync" ); return status; }

	{
		MFnTypedAttribute tAttr;
		oa_simulatedNurbs = tAttr.create( "simulatedNurbs", "sns",
			MFnData::kNurbsCurve, & status );
		CHECK_MSTATUS( status );
	    CHECK_MSTATUS( tAttr.setArray( true ) );
	    CHECK_MSTATUS( tAttr.setReadable( true ) );
	    CHECK_MSTATUS( tAttr.setWritable( false ) );
	    CHECK_MSTATUS( tAttr.setConnectable( true ) );
	    CHECK_MSTATUS( tAttr.setUsesArrayDataBuilder( true ) );
	    status = addAttribute( oa_simulatedNurbs );
	    CHECK_MSTATUS( status );
	}
	status = attributeAffects( ia_time, oa_simulatedNurbs );
	if ( !status ) { status.perror( "attributeAffects ia_time->oa_simulatedNurbs" ); return status; }
	status = attributeAffects( ca_rodPropertiesSync, oa_simulatedNurbs );
	if ( !status ) { status.perror( "attributeAffects ca_rodPropertiesSync->oa_simulatedNurbs" ); return status; }

	{
		MFnNumericAttribute nAttr;
		ia_collisionMeshes = nAttr.create( "collisionMeshes", "com", MFnNumericData::kBoolean, false, &status );
		if (!status)
		{
			status.perror( "create ia_collisionMeshes attribute" );
		    return status;
		}
		nAttr.setWritable( true );
		nAttr.setReadable( false );
		nAttr.setConnectable( true );
		nAttr.setDisconnectBehavior( MFnAttribute::kDelete );
		nAttr.setArray( true );
		status = addAttribute( ia_collisionMeshes );
		if ( !status ) { status.perror( "addAttribute ia_collisionMeshes" ); return status; }
	}
	status = attributeAffects( ia_collisionMeshes, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_collisionMeshes->ca_rodPropertiesSync" ); return status; }

	//Solver settings
	addNumericAttribute( ia_stol, "stol", "stl", MFnNumericData::kDouble, 99, true );
		status = attributeAffects( ia_stol, ca_rodPropertiesSync );
	    if ( !status ) { status.perror( "attributeAffects ia_stol->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_atol, "atol", "atl", MFnNumericData::kDouble, 8, true );
		status = attributeAffects( ia_atol, ca_rodPropertiesSync );
	    if ( !status ) { status.perror( "attributeAffects ia_atol->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_rtol, "rtol", "rtl", MFnNumericData::kDouble, 99, true );
		status = attributeAffects( ia_rtol, ca_rodPropertiesSync );
	    if ( !status ) { status.perror( "attributeAffects ia_rtol->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_inftol, "inftol", "itl", MFnNumericData::kDouble, 8, true );
	status = attributeAffects( ia_inftol, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_inftol->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_numLineSearchIters, "numLineSearchIters", "nlsi", MFnNumericData::kInt, 2, true );
	status = attributeAffects( ia_numLineSearchIters, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_numLineSearchIters->ca_rodPropertiesSync" ); return status; }

	//General parameters
	addNumericAttribute( ia_enablePenaltyResponse, "enablePenaltyResponse", "epr", MFnNumericData::kBoolean, true, true );
	status = attributeAffects( ia_enablePenaltyResponse, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_enablePenaltyResponse->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_implicitThickness, "implicitThickness", "imt", MFnNumericData::kDouble, 0.10, true );
	status = attributeAffects( ia_implicitThickness, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_implicitThickness->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_implicitStiffness, "implicitStiffness", "ims", MFnNumericData::kDouble, 1.0, true );
	status = attributeAffects( ia_implicitStiffness, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_implicitStiffness->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_inextensibilityThreshold, "inextensibilityThreshold", "ixf", MFnNumericData::kInt, 0, true );
	status = attributeAffects( ia_inextensibilityThreshold, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_inextensibilityThreshold->ca_rodPropertiesSync" ); return status; }

    // Debug drawing
    addNumericAttribute( ia_shouldDrawVelocity, "shouldDrawVelocity", "sdv", MFnNumericData::kBoolean, false, true );
	status = attributeAffects( ia_shouldDrawVelocity, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_shouldDrawVelocity->ca_rodPropertiesSync" ); return status; }

	//Failure  Detection
	addNumericAttribute( ia_maxNumOfSolverIters, "maxNumOfSolverIters", "mnsi", MFnNumericData::kInt, 250, true );
	status = attributeAffects( ia_maxNumOfSolverIters, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_maxNumOfSolverIters->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_maxNumOfCollisionIters, "maxNumOfCollisionIters", "mnci", MFnNumericData::kInt, 0, true );
	status = attributeAffects( ia_maxNumOfCollisionIters, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_maxNumOfCollisionIters->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_enableExplosionDetection, "enableExplosionDetection", "eex", MFnNumericData::kBoolean, true, true );
	status = attributeAffects( ia_enableExplosionDetection, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_enableExplosionDetection->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_explosionDampening, "explosionDampening", "exd", MFnNumericData::kDouble, 100.0, true );
	status = attributeAffects( ia_explosionDampening, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_iexplosionDampening->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_explosionThreshold, "explosionThreshold", "ext", MFnNumericData::kDouble, 0.5, true );
	status = attributeAffects( ia_explosionThreshold, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_explosionThreshold->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_stretchingThreshold, "stretchingThreshold", "sxt", MFnNumericData::kDouble, 2.0, true );
	status = attributeAffects( ia_stretchingThreshold, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_stretchingThreshold->ca_rodPropertiesSync" ); return status; }

	{
		MFnEnumAttribute enumAttrFn;
	    ia_solverFailure = enumAttrFn.create( "ifSolverStillFails", "svf", (short) FailureMode::IgnoreError, & status );
	    CHECK_MSTATUS( status );
	    enumAttrFn.addField( "Ignore error",   (short) FailureMode::IgnoreError );
	    enumAttrFn.addField( "Kill the rod",  (short) FailureMode::KillTheRod );
	    enumAttrFn.addField( "Halt simulation",  (short) FailureMode::HaltSimulation );
	    enumAttrFn.setKeyable( false );
	    enumAttrFn.setStorable( true );
	    enumAttrFn.setWritable( true );
	    enumAttrFn.setReadable( true );
	    status = addAttribute( ia_solverFailure );
	    CHECK_MSTATUS( status );
	}
	status = attributeAffects( ia_solverFailure, ca_rodPropertiesSync );
    if (!status) { status.perror( "attributeAffects ia_solverFailure->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_maxNumSolverSubsteps, "maxNumSolverSubsteps", "mnss", MFnNumericData::kInt, 0, true );
	status = attributeAffects( ia_maxNumSolverSubsteps, ca_rodPropertiesSync );
    if ( !status ) { status.perror( "attributeAffects ia_maxNumSolverSubsteps->ca_rodPropertiesSync" ); return status; }

	{
		MFnEnumAttribute enumAttrFn;
		ia_collisionFailure = enumAttrFn.create( "ifCollisionStillFails", "clf", (short) FailureMode::IgnoreError, & status );
		CHECK_MSTATUS( status );
		enumAttrFn.addField( "Ignore error",   (short) FailureMode::IgnoreError );
		enumAttrFn.addField( "Kill the rod",  (short) FailureMode::KillTheRod );
		enumAttrFn.addField( "Halt simulation",  (short) FailureMode::HaltSimulation );
		enumAttrFn.setKeyable( false );
		enumAttrFn.setStorable( true );
		enumAttrFn.setWritable( true );
		enumAttrFn.setReadable( true );
		status = addAttribute( ia_collisionFailure );
		CHECK_MSTATUS( status );
	}
	status = attributeAffects( ia_collisionFailure, ca_rodPropertiesSync );
	if (!status) { status.perror( "attributeAffects ia_collisionFailure->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_maxNumCollisionSubsteps, "maxNumCollisionSubsteps", "mncs", MFnNumericData::kInt, 0, true );
	status = attributeAffects( ia_maxNumCollisionSubsteps, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_maxNumCollisionSubsteps->ca_rodPropertiesSync" ); return status; }

	{
		MFnEnumAttribute enumAttrFn;
		ia_explosionFailure = enumAttrFn.create( "ifExplosionStillFails", "exf", (short) FailureMode::IgnoreError, & status );
	    CHECK_MSTATUS( status );
	    enumAttrFn.addField( "Ignore error",   (short) FailureMode::IgnoreError );
	    enumAttrFn.addField( "Kill the rod",  (short) FailureMode::KillTheRod );
	    enumAttrFn.addField( "Halt simulation",  (short) FailureMode::HaltSimulation );
	    enumAttrFn.setKeyable( false );
	    enumAttrFn.setStorable( true );
	    enumAttrFn.setWritable( true );
	    enumAttrFn.setReadable( true );
	    status = addAttribute( ia_explosionFailure );
	    CHECK_MSTATUS( status );
	}
	status = attributeAffects( ia_explosionFailure, ca_rodPropertiesSync );
	if (!status) { status.perror( "attributeAffects ia_explosionFailure->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_maxNumExplosionSubsteps, "maxNumExplosionSubsteps", "mnes", MFnNumericData::kInt, 7, true );
	status = attributeAffects( ia_maxNumExplosionSubsteps, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_maxNumExplosionSubsteps->ca_rodPropertiesSync" ); return status; }

	{
		MFnEnumAttribute enumAttrFn;
		ia_stretchingFailure = enumAttrFn.create( "ifStretchingStillFails", "sxf", (short) FailureMode::KillTheRod, & status );
		CHECK_MSTATUS( status );
	    enumAttrFn.addField( "Ignore error",   (short) FailureMode::IgnoreError );
	    enumAttrFn.addField( "Kill the rod",  (short) FailureMode::KillTheRod );
	    enumAttrFn.addField( "Halt simulation",  (short) FailureMode::HaltSimulation );
	    enumAttrFn.setKeyable( false );
	    enumAttrFn.setStorable( true );
	    enumAttrFn.setWritable( true );
	    enumAttrFn.setReadable( true );
	    status = addAttribute( ia_stretchingFailure );
	    CHECK_MSTATUS( status );
	}
	status = attributeAffects( ia_stretchingFailure, ca_rodPropertiesSync );
	if (!status) { status.perror( "attributeAffects ia_stretchingFailure->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_maxNumStretchingSubsteps, "maxNumStretchingSubsteps", "mnts", MFnNumericData::kInt, 0, true );
	status = attributeAffects( ia_maxNumStretchingSubsteps, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_maxNumStretchingSubsteps->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_verticesPerStrand, "verticesPerStrand", "vps", MFnNumericData::kInt, 12, true );
	status = attributeAffects( ia_verticesPerStrand, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects verticesPerStrand->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_verticesPerRod, "verticesPerRod", "cpr", MFnNumericData::kInt, 10, true );
	status = attributeAffects( ia_verticesPerRod, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_verticesPerRod->ca_rodPropertiesSync" ); return status; }

	return MS::kSuccess;
}
