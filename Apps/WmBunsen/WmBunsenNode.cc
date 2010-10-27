#include "WmBunsenNode.hh"
#include "WmBunsenCollisionMeshNode.hh"
#include "constraints/WmFigConstraintNode.hh"
#include "WmFigRodComponentList.hh"
#include <maya/MFnMessageAttribute.h>

MTypeId WmBunsenNode::typeID( 0x001135, 0x18 ); 
MString WmBunsenNode::typeName( "wmFigaroNode" );
MObject WmBunsenNode::ca_syncAttrs;
MObject WmBunsenNode::ia_time;
MObject WmBunsenNode::ia_fps;
MObject WmBunsenNode::ia_substeps;
MObject WmBunsenNode::ia_maxDt;
MObject WmBunsenNode::ia_maxIter;
MObject WmBunsenNode::ia_startTime;
MObject WmBunsenNode::ia_rodsNodes;
MObject WmBunsenNode::ia_gravity;
MObject WmBunsenNode::ia_numberOfThreads;
MObject WmBunsenNode::ia_solver;
MObject WmBunsenNode::ia_enabled;
MObject WmBunsenNode::ia_collisionMeshes;
MObject WmBunsenNode::ia_collisionsEnabled;
MObject WmBunsenNode::ia_plasticDeformations;
MObject WmBunsenNode::ia_solverType;

MObject WmBunsenNode::ia_msgConstraintNodes;

// Clumping
/* static */ MObject WmBunsenNode::ia_isClumpingEnabled;
/* static */ MObject WmBunsenNode::ia_clumpingCoefficient;

MObject WmBunsenNode::ia_selfCollisionPenaltyForces;
MObject WmBunsenNode::ia_fullSelfCollisions;
MObject WmBunsenNode::ia_fullSelfCollisionCOR;
MObject WmBunsenNode::ia_fullSelfCollisionIterations;
MObject WmBunsenNode::ia_timingsFile;
MObject WmBunsenNode::ia_timingEnabled;
MObject WmBunsenNode::oa_simStepTaken;

// Volumetric Collisions
/* static */ MObject WmBunsenNode::ia_volumetricCollisions;
/* static */ MObject WmBunsenNode::ia_gridDX;
/* static */ MObject WmBunsenNode::ia_targetEdgeDensity;
/* static */ MObject WmBunsenNode::ia_volumetricRadius;
/* static */ MObject WmBunsenNode::ia_flip;
/* static */ MObject WmBunsenNode::ia_slip;
/* static */ MObject WmBunsenNode::ia_separationCondition;
/* static */ MObject WmBunsenNode::ia_displayGrid;
/* static */ MObject WmBunsenNode::ia_displayAirBoundary;
/* static */ MObject WmBunsenNode::ia_displayCollisionBoundary;
/* static */ MObject WmBunsenNode::ia_displayGridVelocitiesMultiplier;
/* static */ MObject WmBunsenNode::ia_maxDisplayDensity;

// Drawing
/* static */ MObject WmBunsenNode::ia_drawSubSteppedVertices;

// XML output
/* static */ MObject WmBunsenNode::ia_writeToXMLFile;
/* static */ MObject WmBunsenNode::ia_XMLFilePath;

WmBunsenNode::WmBunsenNode() : m_initialised( false ), m_enabled( true ), m_beaker( NULL ),
    m_solverType( RodTimeStepper::NONE ), m_writeXMLData( false ), m_xmlFilePath( "" )
{
    m_beaker = new Beaker();
}

WmBunsenNode::~WmBunsenNode()
{
    if ( m_beaker != NULL )
        delete m_beaker;
}

void WmBunsenNode::pullOnAllRodNodes( MDataBlock& i_dataBlock )
{
    MStatus stat;
    
    // Pull all on all input rod nodes, causing them to update the rod data owned by beaker
    // that they each have pointers to.
    MArrayDataHandle inArrayH = i_dataBlock.inputArrayValue( ia_rodsNodes, &stat );
    CHECK_MSTATUS(stat);
    size_t numRodsConnected = inArrayH.elementCount();
  
    for ( unsigned int r=0; r < numRodsConnected; r++ ) 
    {
        inArrayH.jumpToElement( r );
        inArrayH.inputValue( &stat );
        CHECK_MSTATUS( stat );
        
        // and thats it! The rod node will get the signal that it needs to update its output
        // and will directly change the data in m_beaker. It's dumb to pass it along Maya connections
        // to here then to beaker. So we cut out the middle man.
    }
}  

void WmBunsenNode::addRodsToWorld( MDataBlock& i_dataBlock )
{
    MStatus stat;
    
    // Run through each attached rod node and create the associated rod data structure inside 
    // Beaker. This will get called after all the rods have been deleted as Maya has
    // just been moved to start time.
  
    MPlug rodPlugArray( thisMObject(), ia_rodsNodes );
    CHECK_MSTATUS( stat );
    size_t numRodsConnected = rodPlugArray.numConnectedElements( &stat );
    CHECK_MSTATUS( stat );

    for ( unsigned int r=0; r < numRodsConnected; r++ )
    {
        if ( rodPlugArray.isArray( &stat ) )
        {
            MPlug rodPlug = rodPlugArray.elementByLogicalIndex( r, &stat );
            CHECK_MSTATUS( stat );
            if ( rodPlug.isConnected( &stat ) )
            {
                MPlugArray inPlugArr;
                rodPlug.connectedTo( inPlugArr, true, false, &stat );
                CHECK_MSTATUS( stat );
                
                // Since we asked for the destination there can only be one plug in the array
                MPlug rodNodePlug = inPlugArr[0];
                MObject rodNodeObj = rodNodePlug.node( &stat );
                CHECK_MSTATUS( stat );
                MFnDependencyNode rodNodeFn( rodNodeObj );
                WmFigRodNode* wmFigRodNode = ( WmFigRodNode* )rodNodeFn.userNode();
    
                // Now the rod node has been initialised the undeformed positions for the rods
                // it owns. Since we are resetting the sim we need to actually now create the
                // rods and add them to the world.
                
                m_beaker->addRodsToWorld( r, wmFigRodNode->rodGroup() );
            }
            else
                CHECK_MSTATUS( stat );
        }
    }

    if ( m_writeXMLData )
    {        
        std::string filename =  MFileIO::currentFile().asChar();
        m_beaker->startXMLLogging( m_xmlFilePath, filename );
        m_writeXMLData = true;
    }
}


void WmBunsenNode::updateAllCollisionMeshes( MDataBlock &data )
{
    //
    // Kinematic Objects
    // Each Kinematic Object in m_derManager holds a pointer to the mesh data in the collisionMeshNode.
    // This means that we don't need to update the mesh data in Kinematic objects, they will always 
    // have the up to date info. 
    // However - when Maya loads a scene it connects attributes in an arbitrary fashion and so
    // it is possible that the Kinematic object is created before the collision Mesh node had
    // been created. In that situation we can't pass that pointer when this node is connected to 
    // collision mesh nodes ( in connectionMade() ). So rather than having code to connect the ones
    // we can connect there in connectionMade() and also code here we just do it here each time 
    // step. It's a super low overhead test so it is ok to do it once per frame.
    // What we do - for each connected kinematic object check if it has been initialised
    // and if it hasn't and we now have valid data from the collision mesh node then initialise
    // the kinematic object.
    MStatus stat;
    
    MArrayDataHandle inArrayH = data.inputArrayValue( ia_collisionMeshes, &stat );
    CHECK_MSTATUS(stat);
    size_t numMeshesConnected = inArrayH.elementCount();
  
    for ( unsigned int i=0; i < numMeshesConnected; i++ ) 
    {

        // Even if we don't use it, grab the data so Maya knows to evaluate the node    
        inArrayH.jumpToElement(i);
        MDataHandle collisionMeshH = inArrayH.inputValue( &stat);
        CHECK_MSTATUS(stat);

        MObject collisionMeshDataObj = collisionMeshH.data();

        // First check if the kinematic object has been initialised, if it has then don't bother
        // continuing as we would be updating the pointer with the same pointer.
        if ( m_beaker->collisionMeshInitialised( i ) )
            continue;

        MPlug plug( thisMObject(), ia_collisionMeshes );
        CHECK_MSTATUS( stat );
        if ( plug.isArray( &stat ) )
        {
            MPlug indxPlug = plug.elementByLogicalIndex( i, &stat );
            CHECK_MSTATUS( stat );
            if ( indxPlug.isConnected( &stat ) ) 
            {
                MPlugArray inPlugArr;
                indxPlug.connectedTo( inPlugArr, true, false, &stat );
                CHECK_MSTATUS( stat );
                
                // Since we asked for the destination there can only be one plug in the array
                MPlug meshPlug = inPlugArr[0];
                MObject collisionMeshNodeObj = meshPlug.node( &stat );
                CHECK_MSTATUS( stat );
                MFnDependencyNode collisionMeshNodeFn( collisionMeshNodeObj );
                WmBunsenCollisionMeshNode* collisionMeshNode = (WmBunsenCollisionMeshNode*)collisionMeshNodeFn.userNode();
                BASim::CollisionMeshData* collisionMeshData = collisionMeshNode->collisionMeshData();

                // If the connection was made on file load then it is possible the collision mesh
                // node had not yet created the data from the mesh. If so then skip it and we'll 
                // get it when we are finished loading and time moves.
                if (collisionMeshData != NULL)
                    m_beaker->initialiseCollisionMesh(collisionMeshData, i);
            }
            else
                CHECK_MSTATUS( stat );            
        }
    }
}

void WmBunsenNode::updateAllRodNodes( MDataBlock &i_dataBlock )
{
    MStatus stat;
    
    // Pull all on all input rod nodes, causing them to update the rod data owned by beaker
    // that they each have pointers to.
    MArrayDataHandle inArrayH = i_dataBlock.inputArrayValue( ia_rodsNodes, &stat );
    CHECK_MSTATUS(stat);
    size_t numRodsConnected = inArrayH.elementCount();
  
    for ( unsigned int r=0; r < numRodsConnected; r++ ) 
    {
        inArrayH.jumpToElement( r );
        inArrayH.inputValue( &stat );
        CHECK_MSTATUS( stat );

        // and thats it! The rod node will get the signal that it needs to update its output
        // and will directly change the data in m_beaker. It's dumb to pass it along Maya connections
        // to here then to beaker. So we cut out the middle man.
    }
}

void WmBunsenNode::addAllConstraints( MDataBlock &i_dataBlock )
{
	MStatus stat;

	MArrayDataHandle msgConstraintNodesH = i_dataBlock.inputArrayValue( ia_msgConstraintNodes, &stat );
	CHECK_MSTATUS(stat);
	size_t numConstraints = msgConstraintNodesH.elementCount();

	//MGlobal::displayInfo( MString("Num constraints: ") + numConstraints );

	MPlug msgConstraintNodesPlug( thisMObject(), ia_msgConstraintNodes );
	//MGlobal::displayInfo( MString("Num constraints2: ") + msgConstraintNodesPlug.numElements() );

	unsigned int i;
	for( i=0; i < numConstraints; i++ )
	{
		MPlug msgConstraintNodePlug = msgConstraintNodesPlug.elementByPhysicalIndex( i );
		if( !msgConstraintNodePlug.isNull() )
		{
			//MGlobal::displayInfo( MString("found element: ") + i );

			MPlugArray srcPlugs;
			msgConstraintNodePlug.connectedTo( srcPlugs, true, false );
			if( srcPlugs.length() )
			{
				MObject figConstraintNode = srcPlugs[0].node();
				MFnDependencyNode figConstraintFn( figConstraintNode );
				//MGlobal::displayInfo( "Constraint node: " + figConstraintFn.name() );

		        MPlug enablePlug = figConstraintFn.findPlug( "enable" );
		        if( !enablePlug.asBool() )
		        	continue;

				WmFigConstraintNode &constraint = *static_cast<WmFigConstraintNode*>( figConstraintFn.userNode() );
				constraint.rodVertexConstraints.clear();

		        MPlug stiffnessPlug = figConstraintFn.findPlug( "stiffness" );
			    double stiffness = stiffnessPlug.asInt();

		        MPlug targetWorldPositionPlug = figConstraintFn.findPlug( "targetWorldPosition" );
			    MObject targetWorldPositionData = targetWorldPositionPlug.asMObject();
			    MFnNumericData numericDataFn( targetWorldPositionData );
		        MPoint targetWorldPosition;
			    numericDataFn.getData3Double( targetWorldPosition.x, targetWorldPosition.y, targetWorldPosition.z );

//		        MPlug rodIdPlug = figConstraintFn.findPlug( "rodId" );
//			    int rodId = rodIdPlug.asInt();
//
//			    MPlug vertexIdPlug = figConstraintFn.findPlug( "vertexId" );
//			    int vertexId = vertexIdPlug.asInt();

			    MPlug rodVerticesPlug = figConstraintFn.findPlug( "rodVertices" );
			    MString rodVerticesTxt( rodVerticesPlug.asString() );

			    MPlug constraintTypePlug = figConstraintFn.findPlug( "constraintType" );
				int constraintType = constraintTypePlug.asInt();

				MObject figRodNode;
				MPlug figRodNodeMsgPlug = figConstraintFn.findPlug( "figRodNodeMsg" );
				figRodNodeMsgPlug.connectedTo( srcPlugs, true, false );
				if( srcPlugs.length() )
				{
					figRodNode = srcPlugs[0].node();
				}
				MFnDependencyNode figRodNodeFn( figRodNode );

#if 0
				MGlobal::displayInfo( "Setting vertex position penalty. Constraint node: " + figConstraintFn.name() + " rodId: " + rodId + " vertexId: " + vertexId + " constraintType: " + constraintType + \
									  " worldPosition: " + worldPosition.x + ", " + worldPosition.y + ", " + worldPosition.z + \
									  " figRodNode: " + figRodNodeFn.name() );
#endif

#if 0
				WmFigRodNode *rodNode = static_cast< WmFigRodNode*>( figRodNodeFn.userNode());
			    WmFigRodGroup*  rodGroup = rodNode->rodGroup();
		        if( rodId < rodGroup->numberOfRods() )
		        {
		            if( vertexId < 0 )
		            	vertexId = rodGroup->elasticRod( rodId )->nv() -1 ;

		            if( vertexId < rodGroup->elasticRod( rodId )->nv() )
		            {
		                //Vec3d target_position = Vec3d( i_pos.x, i_pos.y, i_pos.z );
		                Vec3d target_position = Vec3d( targetWorldPosition.x, targetWorldPosition.y, targetWorldPosition.z );

		                BASim::RodVertexConstraint *rodVertexConstraint = rodGroup->collisionStepper( rodId )->setVertexPositionPenalty2( vertexId, target_position, stiffness, constraintType );
		                constraint.rodVertexConstraints.push_back( rodVertexConstraint );
		            }
		        }
#endif

		        WmFigRodComponentList rodComponentList;
		        rodComponentList.unserialise( rodVerticesTxt );

				WmFigRodNode *rodNode = static_cast< WmFigRodNode*>( figRodNodeFn.userNode());
			    WmFigRodGroup*  rodGroup = rodNode->rodGroup();
			    unsigned int iRod, rodId, iVert, vertexId, nVerts;
			    for( iRod=0; iRod < rodGroup->numberOfRods(); iRod++ ) {
			    	rodId = iRod;

			    	nVerts = rodGroup->elasticRod( rodId )->nv();
			    	for( iVert=0; iVert < nVerts; iVert++ ) {
			    		vertexId = iVert;

			    		if( rodComponentList.containsRodVertex( rodId, vertexId ) ) {
			                Vec3d target_position = Vec3d( targetWorldPosition.x, targetWorldPosition.y, targetWorldPosition.z );

			                BASim::RodVertexConstraint *rodVertexConstraint = rodGroup->collisionStepper( rodId )->setVertexPositionPenalty2( vertexId, target_position, stiffness, constraintType );
			                constraint.rodVertexConstraints.push_back( rodVertexConstraint );
			    		}
			    	}
			    }
			}
		}
	}
}

void WmBunsenNode::updateAllConstraints( MDataBlock &i_dataBlock )
{
	MPlug msgConstraintNodesPlug( thisMObject(), ia_msgConstraintNodes );
	const unsigned int numConstraints = msgConstraintNodesPlug.numElements();

	unsigned int i;
	for( i=0; i < numConstraints; i++ )
	{
		MPlug msgConstraintNodePlug = msgConstraintNodesPlug.elementByPhysicalIndex( i );
		if( !msgConstraintNodePlug.isNull() )
		{
			//MGlobal::displayInfo( MString("found element: ") + i );

			MPlugArray srcPlugs;
			msgConstraintNodePlug.connectedTo( srcPlugs, true, false );
			if( srcPlugs.length() )
			{
				MObject figConstraintNode = srcPlugs[0].node();
				MFnDependencyNode figConstraintFn( figConstraintNode );
				//MGlobal::displayInfo( "Updating Constraint node: " + figConstraintFn.name() );

		        MPlug enablePlug = figConstraintFn.findPlug( "enable" );
		        if( !enablePlug.asBool() )
		        	continue;

				WmFigConstraintNode &constraint = *static_cast<WmFigConstraintNode*>( figConstraintFn.userNode() );

		        MPlug stiffnessPlug = figConstraintFn.findPlug( "stiffness" );
			    double stiffness = stiffnessPlug.asInt();

		        MPlug worldPositionPlug = figConstraintFn.findPlug( "worldPosition" );
			    MObject worldPositionData = worldPositionPlug.asMObject();
			    MFnNumericData numericDataFn( worldPositionData );
		        MPoint worldPosition;
			    numericDataFn.getData3Double( worldPosition.x, worldPosition.y, worldPosition.z );

#if 0
			    MGlobal::displayInfo( "Updating vertex position penalty. Constraint node: " + figConstraintFn.name() + " stiffness: " + stiffness + \
									  " worldPosition: " + worldPosition.x + ", " + worldPosition.y + ", " + worldPosition.z  );
#endif

			    BASim::RodVertexConstraint *rodVertexConstraint;
			    std::list<BASim::RodVertexConstraint *>::iterator it;
			    for( it=constraint.rodVertexConstraints.begin(); it != constraint.rodVertexConstraints.end(); it++ )
			    {
			    	rodVertexConstraint = *it;
			    	rodVertexConstraint->m_stiff = stiffness;
			    	rodVertexConstraint->m_target = Vec3d( worldPosition.x, worldPosition.y, worldPosition.z );
			    }
			}
		}
	}
}

//void RodCollisionTimeStepper::updateVertexPositionConstraints()
//{
//    //this is dummy test, meaningless at all
//    for( size_t s = 0; s < m_vertexContraints.size(); s++ )
//    {
//        Vec3d direction = m_vertexContraints[ s ]->second.m_target;
//        double norm = direction.norm();
//        if( norm > 0 )
//            direction = direction * ( 1. / norm );
//
//        //make the target position step forward a little bit...
//        m_vertexContraints[ s ]->second.m_target =
//                m_vertexContraints[ s ]->second.m_target + 0.1 * direction;
//    }
//
//}



MStatus WmBunsenNode::compute( const MPlug& i_plug, MDataBlock& i_dataBlock ) 
{
    MStatus stat;

    //cerr << "WmBunsenNode::compute() with i_plug = " << i_plug.name() << endl;
    
    if ( i_plug == ca_syncAttrs )
    {
        m_enabled = i_dataBlock.inputValue( ia_enabled, &stat ).asBool();
        CHECK_MSTATUS( stat );
        
        m_previousTime = m_currentTime;
        m_currentTime = i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
        CHECK_MSTATUS( stat );

        m_startTime = i_dataBlock.inputValue( ia_startTime, &stat ).asDouble();
        CHECK_MSTATUS( stat );

        m_framedt = 1.0 / i_dataBlock.inputValue(ia_fps, &stat).asDouble();
        CHECK_MSTATUS( stat );

        int substeps = i_dataBlock.inputValue(ia_substeps, &stat).asInt();
        CHECK_MSTATUS( stat );

        //m_beaker->setDt(i_dataBlock.inputValue(ia_maxDt, &stat ).asDouble());
        //CHECK_MSTATUS( stat );	

        const double3 &gravity = i_dataBlock.inputValue( ia_gravity, &stat ).asDouble3();
        CHECK_MSTATUS( stat );
        m_beaker->setGravity( BASim::Vec3d( gravity[0], gravity[1], gravity[2] ) );

    	//m_beaker->setMaxIter(i_dataBlock.inputValue(ia_maxIter, &stat).asInt());
        //CHECK_MSTATUS( stat );
    
        int numberOfThreads = i_dataBlock.inputValue( ia_numberOfThreads, &stat ).asInt();
        CHECK_MSTATUS( stat );

        bool collisionsEnabled = i_dataBlock.inputValue( ia_collisionsEnabled, &stat ).asBool();
        CHECK_MSTATUS( stat );
        bool selfCollisionPenaltyForces = i_dataBlock.inputValue( ia_selfCollisionPenaltyForces, &stat ).asBool();
        CHECK_MSTATUS( stat );
        bool fullSelfCollisions = i_dataBlock.inputValue( ia_fullSelfCollisions, &stat ).asBool();
        CHECK_MSTATUS( stat );
        int selfCollisionsIters = i_dataBlock.inputValue( ia_fullSelfCollisionIterations, &stat ).asInt();
        CHECK_MSTATUS( stat );
        double selfCollisionsCOR = i_dataBlock.inputValue( ia_fullSelfCollisionCOR, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        
        // Clumping
        m_beaker->clumpingEnabled( i_dataBlock.inputValue( ia_isClumpingEnabled, &stat ).asBool() );
        CHECK_MSTATUS( stat );
        m_beaker->clumpingCoefficient( i_dataBlock.inputValue( ia_clumpingCoefficient, &stat ).asDouble() );
        CHECK_MSTATUS( stat );

        
        ////////////////////////////////////////////////////////////////////////////////////////////////
        //
        // Volumetric Collision attributes
        //
        ////////////////////////////////////////////////////////////////////////////////////////////////
    
        getAllVolumetricCollisionAttributes( i_dataBlock );

        int solver = i_dataBlock.inputValue( ia_solver, &stat ).asInt();

        bool plasticDeformations = i_dataBlock.inputValue( ia_plasticDeformations, &stat ).asBool();
        CHECK_MSTATUS( stat );
        m_beaker->setPlasticDeformations( plasticDeformations );

        m_beaker->shouldDrawSubsteppedVertices( i_dataBlock.inputValue( ia_drawSubSteppedVertices, &stat ).asBool() ) ;

        MString xmlPath = i_dataBlock.inputValue( ia_XMLFilePath, &stat ).asString();
        CHECK_MSTATUS( stat );

        m_xmlFilePath = xmlPath.asChar();

        bool writeXMLData = i_dataBlock.inputValue( ia_writeToXMLFile, &stat ).asBool();
        CHECK_MSTATUS( stat );
        
        if ( writeXMLData == false && m_writeXMLData == true )
        {
            m_beaker->writeXMLFileToDisk();
        }

        m_writeXMLData = writeXMLData;

        if ( m_enabled )
        {
            if ( m_currentTime == m_startTime )
            {
                m_beaker->resetEverything();

                updateAllRodNodes( i_dataBlock );
                updateAllCollisionMeshes( i_dataBlock );
                addRodsToWorld( i_dataBlock );

                //MGlobal::displayInfo( "COMPUTE AT START TIME" );
                addAllConstraints( i_dataBlock );
            }
            else
            {
                updateAllRodNodes( i_dataBlock );
                updateAllCollisionMeshes( i_dataBlock );
                updateAllConstraints( i_dataBlock );
            }
            
        }
        
        m_solverType = ( RodTimeStepper::Method) ( i_dataBlock.inputValue( ia_solverType, &stat ).asInt() );
        CHECK_MSTATUS( stat );
        
        if ( m_enabled && ( ( m_solverType == RodTimeStepper::STATICS ) || ( m_currentTime > m_previousTime && m_currentTime > m_startTime ) ) ) 
        {
            m_beaker->takeTimeStep( numberOfThreads, m_framedt, substeps, collisionsEnabled, 
                                    selfCollisionPenaltyForces, fullSelfCollisions,
                                    selfCollisionsIters, selfCollisionsCOR ); 
        }
        
        bool timingEnabled = i_dataBlock.inputValue( ia_timingEnabled, &stat ).asBool();
        CHECK_MSTATUS( stat );
        MString timingsFile = i_dataBlock.inputValue( ia_timingsFile, &stat ).asString();
        CHECK_MSTATUS( stat );
        m_beaker->setTimingsFile( timingsFile.asChar() );
        m_beaker->setTimingEnabled( timingEnabled );
    
        MDataHandle outputData = i_dataBlock.outputValue ( ca_syncAttrs, &stat );
        if ( !stat )
        {
            stat.perror("WmBunsenNode::compute get ca_syncAttrs");
            return stat;
        }
	
        // We don't even need to put anything in the output handle as nothing uses it.
        // Just tell Maya it's clean so it doesn't repeatedly evaluate it.
	
        stat = i_dataBlock.setClean( i_plug );
        if ( !stat )
        {
            stat.perror("WmBunsenNode::compute setClean");
            return stat;
        }
    }
    else if ( i_plug == oa_simStepTaken )
    {
        // Get time so Maya knows we care about it.
        i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
        CHECK_MSTATUS(stat);
        
        // Get ca_StepTime so that we know it has moved the sim forward
        i_dataBlock.inputValue( ca_syncAttrs, &stat ).asTime().value();
        CHECK_MSTATUS(stat);
        
        //////////////////////////////////////////////////////////////////
        //
        // We don't actually output any data here. For speed we have the
        // connection that is really a message attribute. It tells the
        // node on the other end that things have changed. It will then
        // get a pointer to this class and grab the data it needs.
        // This could be done with compound attributes or an MPxData class
        // but that is a pain in the ass for such a simple task and we never
        // want to save the data so we don't need to do it.
        //
        //////////////////////////////////////////////////////////////////
        
        stat = i_dataBlock.setClean( i_plug );
        if ( !stat ) 
        {
            stat.perror( "WmBunsenNode::compute oa_simStepTaken setClean()" );
            return stat;
        }
    }
    else
    {
	return MS::kUnknownParameter;
    }

    return MS::kSuccess;
}

void WmBunsenNode::draw( M3dView& i_view, const MDagPath& i_path,
                         M3dView::DisplayStyle i_style,
                         M3dView::DisplayStatus i_status )
{ 

    MStatus stat;
    MObject thisNode = thisMObject();
    
    MPlug syncPlug( thisNode, ca_syncAttrs );
    double d; 
    stat = syncPlug.getValue( d );
    if ( !stat )
    {
        stat.perror( "WmBunsenNode::draw getting ca_syncAttrs" );
        return;
    }
    
    i_view.beginGL(); 
    glPushAttrib( GL_CURRENT_BIT | GL_POINT_BIT | GL_LINE_BIT | GL_ENABLE_BIT |  GL_LIGHTING_BIT );
    
    // Only draw beaker if it's not start time as what it draws only makes sense
    // when the sim is running.
    if ( m_currentTime != m_startTime )
        m_beaker->draw();
    
    // What did this line do? it was here from the devkit example. Is it to with point colouring
    //view.setDrawColor ( WmBunsenNode );
    
    glPopAttrib();
    i_view.endGL();
}

void WmBunsenNode::getAllVolumetricCollisionAttributes( MDataBlock& i_dataBlock )
{
    MStatus stat;

    double flip = i_dataBlock.inputValue( ia_flip, &stat).asDouble();
    CHECK_MSTATUS(stat);
    m_beaker->setFlip( flip );

    double slip = i_dataBlock.inputValue( ia_slip, &stat).asDouble();
    CHECK_MSTATUS(stat);
    m_beaker->setSlip( slip );

    bool doVolumetricCollisions = i_dataBlock.inputValue( ia_volumetricCollisions, &stat).asBool();
    CHECK_MSTATUS(stat);
    m_beaker->setDoVolumetricCollisions( doVolumetricCollisions );

    double targetEdgeDensity = i_dataBlock.inputValue( ia_targetEdgeDensity, &stat).asDouble();
    CHECK_MSTATUS(stat);
    m_beaker->setTargetEdgeDensity( targetEdgeDensity );

    double volumetricRadius = i_dataBlock.inputValue( ia_volumetricRadius, &stat).asDouble();
    CHECK_MSTATUS(stat);
    m_beaker->setVolumetricRadius( volumetricRadius );

    double gridDX = i_dataBlock.inputValue( ia_gridDX, &stat).asDouble();
    CHECK_MSTATUS( stat );
    m_beaker->setGridDX( gridDX );

    const double3 &separationCondition = i_dataBlock.inputValue( ia_separationCondition, &stat).asDouble3();    
    CHECK_MSTATUS(stat);
    m_beaker->setSeparationCondition( separationCondition[ 0 ], separationCondition[ 1 ],
                                      separationCondition[ 2 ] );
    
    bool displayGrid = i_dataBlock.inputValue( ia_displayGrid, &stat).asBool();
    CHECK_MSTATUS(stat);
    m_beaker->setDisplayGrid( displayGrid );

    double displayGridVelocitiesMultiplier = i_dataBlock.inputValue(ia_displayGridVelocitiesMultiplier, &stat).asDouble();
    CHECK_MSTATUS(stat);
    m_beaker->setDisplayGridVelocitiesMultiplier( displayGridVelocitiesMultiplier );

    double maxDisplayDensity = i_dataBlock.inputValue(ia_maxDisplayDensity, &stat).asDouble();
    CHECK_MSTATUS(stat);
    m_beaker->setMaxDisplayDensity( maxDisplayDensity );

    bool displayCollisionBoundary = i_dataBlock.inputValue(ia_displayCollisionBoundary, &stat).asBool();
    CHECK_MSTATUS(stat);
    m_beaker->setDisplayCollisionBoundary( displayCollisionBoundary );

    bool displayAirBoundary = i_dataBlock.inputValue(ia_displayAirBoundary, &stat).asBool();
    CHECK_MSTATUS(stat);
    m_beaker->setDisplayAirBoundary( displayAirBoundary );
}


MStatus WmBunsenNode::connectionMade( const MPlug& i_plug, const MPlug& i_otherPlug, bool i_asSrc ) 
{     
    MStatus stat;
    MStatus retVal( MS::kUnknownParameter );

    // It would be so great if we could do this here but Maya loads objects in some random order
    // so we can't guarantee that the mesh node will have a mesh attached when it is connected
    // here.
    /*
    if ( i_plug == ia_collisionMeshes)
    {
        cerr << "Connecting collision mesh to dynamics node\n";

        size_t idx = i_plug.logicalIndex();
       
                 
       MObject collisionMeshDataObj = otherPlug.node( &stat );
        CHECK_MSTATUS( stat );
        MFnDependencyNode collisionMeshNodeFn( collisionMeshDataObj );
        CollisionMeshNode *collisionMeshNode = (CollisionMeshNode*)collisionMeshNodeFn.userNode();
        
        CollisionMeshData* collisionMeshData = collisionMeshNode->collisionMeshData();
        
        m_beaker->addCollisionMesh( idx );
    }

    */
    
    return retVal;
}

MStatus WmBunsenNode::connectionBroken( const MPlug& i_plug, const MPlug& i_otherPlug, bool i_asSrc )
{
    MStatus retVal( MS::kUnknownParameter );

    if ( i_plug == ia_collisionMeshes )
    {
        int idx = i_plug.logicalIndex();
        m_beaker->removeCollisionMesh( idx );
    }

    return MStatus::kUnknownParameter;
}

bool WmBunsenNode::isBounded() const
{ 
    return false;
}

void* WmBunsenNode::creator()
{
    return new WmBunsenNode();
}

/*static */ MStatus WmBunsenNode::addNumericAttribute( MObject& i_attribute, MString i_longName, 
    MString i_shortName, MFnNumericData::Type i_type, double i_defaultValue, bool i_isInput,
    bool i_isArray )
{
    // Creates a numeric attribute with default attributes
    MStatus stat = MS::kSuccess;

    MFnNumericAttribute nAttr;
    i_attribute = nAttr.create( i_longName, i_shortName, i_type, i_defaultValue, &stat );
    if ( !stat )
    {
        cerr << "Failed to create attribute " << i_longName << endl;
        return stat;
    }
    if ( i_isInput )
        nAttr.setWritable( true );
    else
        nAttr.setWritable( false );
    
    if ( i_isArray )
        nAttr.setArray( true );

    stat = addAttribute( i_attribute );
    if ( !stat ) { stat.perror( "addAttribute " + i_longName ); return stat; }

    return stat;
}

MStatus WmBunsenNode::initialize()
{ 
    MStatus stat;
    
    {
        MFnEnumAttribute enumAttrFn;
        ia_solverType = enumAttrFn.create( "solverType", "sot", (short) RodTimeStepper::SYM_IMPL_EULER, & stat );
        CHECK_MSTATUS( stat );
        enumAttrFn.addField( "Implicit Euler",   (short) RodTimeStepper::IMPL_EULER );
        enumAttrFn.addField( "Symmetric Implicit Euler",  (short) RodTimeStepper::SYM_IMPL_EULER );
        enumAttrFn.addField( "Statics",   (short) RodTimeStepper::STATICS );
        enumAttrFn.setKeyable( false );
        enumAttrFn.setStorable( true );
        enumAttrFn.setWritable( true );
        enumAttrFn.setReadable( true );
        stat = addAttribute( ia_solverType );
        CHECK_MSTATUS( stat );
    }

    {
        MFnUnitAttribute uAttr;
        ia_time = uAttr.create( "time", "t", MTime( 0.0 ), &stat );
        if ( !stat) 
        {
            stat.perror("create ia_time attribute");
            return stat;
        }
        CHECK_MSTATUS( uAttr.setWritable(true) );
        CHECK_MSTATUS( uAttr.setConnectable(true) );
        CHECK_MSTATUS( uAttr.setStorable(false) );
        stat = addAttribute( ia_time );
        if ( !stat ) { stat.perror( "addAttribute ia_time" ); return stat; }
    }
    
    {
        MFnNumericAttribute nAttr;
    	ia_startTime = nAttr.create( "startTime", "stt", MFnNumericData::kDouble, 1.0, &stat );
        if ( !stat ) 
        {
            stat.perror( "create aStartTime attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( true );
        nAttr.setKeyable( true );  
        stat = addAttribute( ia_startTime );
        if ( !stat ) { stat.perror( "addAttribute ia_startTime" ); return stat; }
    }

   {
	MFnNumericAttribute nAttr;
    	ia_fps = nAttr.create( "framesPerSecond", "fps", MFnNumericData::kDouble, 24.0, &stat );
        if ( !stat ) 
        {
            stat.perror( "create fps attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );  
        stat = addAttribute( ia_fps );
        if ( !stat ) { stat.perror( "addAttribute ia_fps" ); return stat; }
    }

    {
        MFnNumericAttribute nAttr;
    	ia_substeps = nAttr.create( "subSteps", "sub", MFnNumericData::kInt, 10, &stat );
        if ( !stat ) 
        {
            stat.perror( "create substeps attribute");
            return stat;
        }
        nAttr.setWritable( true );
        //nAttr.setReadable( false );
        nAttr.setKeyable( true );  
        stat = addAttribute( ia_substeps );
        if ( !stat ) { stat.perror( "addAttribute ia_substeps" ); return stat; }
    }

    {
	MFnNumericAttribute nAttr;
    	ia_maxDt = nAttr.create( "maxDt", "mdt", MFnNumericData::kDouble, 0.01, &stat );
        if ( !stat ) 
        {
            stat.perror( "create maxDt attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );  
        stat = addAttribute( ia_maxDt );
        if ( !stat ) { stat.perror( "addAttribute ia_maxDt" ); return stat; }
    }

   {
	MFnNumericAttribute nAttr;
    	ia_maxIter = nAttr.create( "maxIter", "mitr", MFnNumericData::kInt, 100 , &stat );
        if ( !stat ) 
        {
            stat.perror( "create maxIter attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );  
        stat = addAttribute( ia_maxIter );
        if ( !stat ) { stat.perror( "addAttribute ia_maxIter" ); return stat; }
    }
    
    {
        MFnNumericAttribute nAttr;
        ia_rodsNodes = nAttr.create( "rodsNodes", "rod", MFnNumericData::kBoolean, true, &stat );
        CHECK_MSTATUS (stat );
        CHECK_MSTATUS( nAttr.setWritable( true ) );
        CHECK_MSTATUS( nAttr.setReadable( false ) );
        CHECK_MSTATUS( nAttr.setConnectable( true ) );
        CHECK_MSTATUS( nAttr.setArray( true ) );
        //CHECK_MSTATUS( tAttr.setStorable( false ) );
        stat = addAttribute( ia_rodsNodes );
        if ( !stat ) { stat.perror( "addAttribute ia_rodsNodes" ); return stat;}
    }
    
    {
        MFnNumericAttribute nAttr;
        oa_simStepTaken = nAttr.create( "simStepTaken", "sst", MFnNumericData::kBoolean, true, &stat );
        CHECK_MSTATUS (stat );
        CHECK_MSTATUS( nAttr.setWritable( false ) );
        CHECK_MSTATUS( nAttr.setReadable( true ) );
        CHECK_MSTATUS( nAttr.setConnectable( true ) );
        stat = addAttribute( oa_simStepTaken );
        if ( !stat ) { stat.perror( "addAttribute oa_simulatedRods" ); return stat; }
    }
    
    {
        MFnNumericAttribute nAttr;
        ia_numberOfThreads = nAttr.create( "numberOfThreads", "nut", MFnNumericData::kInt, 1, &stat );
        CHECK_MSTATUS (stat );
        CHECK_MSTATUS( nAttr.setWritable( true ) );
        CHECK_MSTATUS( nAttr.setReadable( false ) );
        CHECK_MSTATUS( nAttr.setConnectable( true ) );
        stat = addAttribute( ia_numberOfThreads );
        if ( !stat ) { stat.perror( "addAttribute ia_numberOfThreads" ); return stat; }
    }
    
    {
        MFnNumericAttribute nAttr;
        ia_gravity = nAttr.create( "gravity", "gr", MFnNumericData::k3Double, 0, &stat );
        if ( !stat ) {
            stat.perror( "create gravity attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( true );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_gravity );
        if (!stat) { stat.perror( "addAttribute ia_gravity" ); return stat; }
    }
 
   {
        MFnNumericAttribute nAttr;
        ia_solver = nAttr.create( "solver", "sol", MFnNumericData::kInt, 1, &stat );
        if ( !stat ) {
            stat.perror( "create solver attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_solver );
        if (!stat) { stat.perror( "addAttribute ia_solver" ); return stat; }
    }
 
    {
        MFnNumericAttribute	nAttr;
        ca_syncAttrs = nAttr.create( "syncAttrs", "sya", MFnNumericData::kDouble, 1.0, &stat );
        if ( !stat) 
        {
            stat.perror( "create ca_syncAttrs attribute" );
            return stat;
        }
        nAttr.setWritable( false );
        nAttr.setReadable( true );
        nAttr.setConnectable( true );
        nAttr.setKeyable( false );  
        stat = addAttribute( ca_syncAttrs );
        if (!stat) { stat.perror( "addAttribute ca_syncAttrs" ); return stat; }
	}
    
    {
        MFnNumericAttribute nAttr;
        ia_collisionMeshes = nAttr.create( "collisionMeshes", "com", MFnNumericData::kBoolean, false, &stat);
        if (!stat) {
            stat.perror("create ia_collisionMeshes attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        nAttr.setDisconnectBehavior( MFnAttribute::kDelete );
        nAttr.setArray( true );
        stat = addAttribute( ia_collisionMeshes );
        if (!stat) { stat.perror( "addAttribute ia_collisionMeshes" ); return stat; }
    }

    {
    	MFnMessageAttribute msgAttr;
    	ia_msgConstraintNodes = msgAttr.create( "ia_msgConstraintNodes", "cnm", &stat  );
        if (!stat) {
            stat.perror("create ia_msgConstraintNodes attribute");
            return stat;
        }
    	msgAttr.setArray( true );
        stat = addAttribute( ia_msgConstraintNodes );
        if (!stat) { stat.perror( "addAttribute ia_msgConstraintNodes" ); return stat; }
    }

    {
        MFnNumericAttribute nAttr;
        ia_plasticDeformations = nAttr.create( "plasticDeformations", "pde", MFnNumericData::kBoolean, false, &stat);
        if (!stat) {
            stat.perror("create ia_plasticDeformations attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_plasticDeformations );
        if (!stat) { stat.perror( "addAttribute ia_plasticDeformations" ); return stat; }
    }

    {
        MFnNumericAttribute nAttr;
        ia_collisionsEnabled = nAttr.create( "objectCollisionsEnabled", "oce", MFnNumericData::kBoolean, true, &stat);
        if (!stat) {
            stat.perror("create ia_collisionsEnabled attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_collisionsEnabled );
        if (!stat) { stat.perror( "addAttribute ia_collisionsEnabled" ); return stat; }
    }

    {
        MFnNumericAttribute nAttr;
        ia_selfCollisionPenaltyForces = nAttr.create( "selfCollisionPenaltyForces", "scp", MFnNumericData::kBoolean, false, &stat);
        if (!stat) {
            stat.perror("create ia_selfCollisionPenaltyForces attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_selfCollisionPenaltyForces );
        if (!stat) { stat.perror( "addAttribute ia_selfCollisionPenaltyForces" ); return stat; }
    }

    {
        MFnNumericAttribute nAttr;
        ia_isClumpingEnabled = nAttr.create( "clumping", "clu", MFnNumericData::kBoolean, false, &stat);
        if (!stat) {
            stat.perror("create ia_isClumpingEnabled attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_isClumpingEnabled );
        if (!stat) { stat.perror( "addAttribute ia_isClumpingEnabled" ); return stat; }
    }

    {
        MFnNumericAttribute nAttr;
        ia_clumpingCoefficient = nAttr.create( "clumpingCoefficient", "clc", MFnNumericData::kDouble, 0.3, &stat );
        if ( !stat ) {
            stat.perror("create ia_clumpingCoefficient attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_clumpingCoefficient );
        if (!stat) { stat.perror( "addAttribute ia_clumpingCoefficient" ); return stat; }
    }

    {
        MFnNumericAttribute nAttr;
        ia_fullSelfCollisions = nAttr.create( "fullSelfCollisions", "fsc", MFnNumericData::kBoolean, false, &stat);
        if (!stat) {
            stat.perror("create ia_fullSelfCollisions attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_fullSelfCollisions );
        if (!stat) { stat.perror( "addAttribute ia_fullSelfCollisions" ); return stat; }
    }

    {
        MFnNumericAttribute nAttr;
        ia_fullSelfCollisionIterations = nAttr.create( "fullSelfCollisionIterations", "fsci", MFnNumericData::kInt, 10, &stat);
        if (!stat) {
            stat.perror("create ia_fullSelfCollisionIterations attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_fullSelfCollisionIterations );
        if (!stat) { stat.perror( "addAttribute ia_fullSelfCollisionIterations" ); return stat; }
    }
    
    {
        MFnNumericAttribute nAttr;
        ia_fullSelfCollisionCOR = nAttr.create( "fullSelfCollisionCOR", "fscc", MFnNumericData::kDouble, 0.1, &stat);
        if (!stat) {
            stat.perror("create ia_fullSelfCollisions attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_fullSelfCollisionCOR );
        if (!stat) { stat.perror( "addAttribute ia_fullSelfCollisionCOR" ); return stat; }
    }
    
    {
        MFnTypedAttribute tAttr;
        MFnStringData fnStringData;
        MObject defaultString = fnStringData.create("");
        ia_timingsFile = tAttr.create( "timingsFile", "tif", MFnData::kString, defaultString, &stat );
        if ( !stat ) 
        {
            stat.perror( "create ia_timingsFile attribute" );
            return stat;
        }
        tAttr.setWritable( true );
        tAttr.setReadable( false );
        tAttr.setConnectable( true );
        stat = addAttribute( ia_timingsFile );
        if (!stat) { stat.perror( "addAttribute ia_timingsFile" ); return stat; }
    }
    {
        MFnNumericAttribute nAttr;
        ia_timingEnabled = nAttr.create( "timingEnabled", "tie", MFnNumericData::kBoolean, false, &stat );
        if (!stat) {
            stat.perror("create ia_timingEnabled attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_timingEnabled );
        if (!stat) { stat.perror( "addAttribute ia_timingEnabled" ); return stat; }
    }
    
    {
        MFnNumericAttribute nAttr;
        ia_enabled = nAttr.create( "enabled", "en", MFnNumericData::kBoolean, true, &stat );
        if (!stat) {
            stat.perror("create ia_enabled attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_enabled );
        if (!stat) { stat.perror( "addAttribute ia_enabled" ); return stat; }
    }
    stat = attributeAffects( ia_enabled, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_enabled->ca_syncAttrs" ); return stat; }

    {
        MFnNumericAttribute nAttr;
        ia_drawSubSteppedVertices = nAttr.create( "drawSubSteppedVertices", "dsv", MFnNumericData::kBoolean, false, &stat );
        if (!stat) {
            stat.perror("create ia_drawSubSteppedVertices attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_drawSubSteppedVertices );
        if (!stat) { stat.perror( "addAttribute ia_drawSubSteppedVertices" ); return stat; }
    }
    stat = attributeAffects( ia_drawSubSteppedVertices, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_enabled->ca_syncAttrs" ); return stat; }
    
    stat = attributeAffects( ia_time, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_time->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_startTime, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_startTimer->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_fps, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_fps->ca_syncAttrs" );return stat;}
    stat = attributeAffects( ia_substeps, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_substeps->ca_syncAttrs" );return stat;}
    stat = attributeAffects( ia_maxDt, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_maxDt->ca_syncAttrs" );return stat;}
    stat = attributeAffects( ia_maxIter, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_maxIter->ca_syncAttrs" );return stat;}
    stat = attributeAffects( ia_rodsNodes, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_rodsNodes->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_gravity, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_gravity->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_solver, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_rodsNodes->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_numberOfThreads, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_numberOfThreads->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_collisionMeshes, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_collisionMeshes->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_collisionsEnabled, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_collisionsEnabled->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_plasticDeformations, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_plasticDeformations->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_selfCollisionPenaltyForces, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_selfCollisionPenaltyForces->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_fullSelfCollisions, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_fullSelfCollisions->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_fullSelfCollisionIterations, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_fullCollisionIterations->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_fullSelfCollisionCOR, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_fullSelfCollisionCOR->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_timingsFile, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_timingsFile->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_timingEnabled, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_timingEnabled->ca_syncAttrs" ); return stat; }

    // Clumping
    stat = attributeAffects( ia_isClumpingEnabled, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_isClumpingEnabled->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_clumpingCoefficient, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_clumpingCoefficient->ca_syncAttrs" ); return stat; }
    
    stat = attributeAffects( ia_time, oa_simStepTaken );
    if (!stat) { stat.perror( "attributeAffects ia_time->oa_simulatedRods" ); return stat; }

    stat = attributeAffects( ia_solverType, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_solverType->ca_syncAttrs" ); return stat; }

    ////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // Volumetric Collisions
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////
    addNumericAttribute( ia_volumetricCollisions, "volumetricCollisions", "voc", MFnNumericData::kBoolean, false, true );
    stat = attributeAffects( ia_volumetricCollisions, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_volumetricCollisions->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_gridDX, "gridDX", "gdx", MFnNumericData::kDouble, 1.0, true );
    stat = attributeAffects( ia_gridDX, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_gridDX->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_targetEdgeDensity, "targetEdgeDensity", "ted", MFnNumericData::kDouble, 100.0, true );
    stat = attributeAffects( ia_targetEdgeDensity, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_targetEdgeDensity->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_volumetricRadius, "volumetricRadius", "vra", MFnNumericData::kDouble, 1.0, true );
    stat = attributeAffects( ia_volumetricRadius, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_volumetricRadius->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_flip, "flip", "flp", MFnNumericData::kDouble, 0.95, true );
    stat = attributeAffects( ia_flip, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_flip->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_slip, "slip", "slp", MFnNumericData::kDouble, .1, true );
    stat = attributeAffects( ia_slip, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_slip->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_separationCondition, "separationCondition", "vsc", MFnNumericData::k3Double, -1.0, true );
    stat = attributeAffects( ia_separationCondition, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_seperationCondition->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_displayGrid, "displayGrid", "dgr", MFnNumericData::kBoolean, false, true );
    stat = attributeAffects( ia_displayGrid, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_displayGrid->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_displayAirBoundary, "displayAirBoundary", "dab", MFnNumericData::kBoolean, false, true );
    stat = attributeAffects( ia_displayAirBoundary, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_displayAirBoundary->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_displayCollisionBoundary, "displayCollisionBoundary", "dcb", MFnNumericData::kBoolean, false, true );
    stat = attributeAffects( ia_displayCollisionBoundary, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_displayCollisionBoundary->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_displayGridVelocitiesMultiplier, "displayGridVelocitiesMultiplier", "dvm", MFnNumericData::kDouble, 0.0, true );
    stat = attributeAffects( ia_displayGridVelocitiesMultiplier, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_displayGridVelocitiesMultiplier->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_maxDisplayDensity, "maxDisplayDensity", "mdd", MFnNumericData::kDouble, 0.0, true );
    stat = attributeAffects( ia_maxDisplayDensity, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_maxDisplayDensity->ca_syncAttrs" ); return stat; }


    ////////////////////////////////////////////////////////////////////////////////////////////////
    //
    // XML File attributes
    //
    ////////////////////////////////////////////////////////////////////////////////////////////////

    {
        MFnTypedAttribute tAttr;
        MFnStringData fnStringData;
        MObject defaultString = fnStringData.create("");
        ia_XMLFilePath = tAttr.create( "XMLFilePath", "xfp", MFnData::kString, defaultString, &stat );
        if ( !stat )
        {
            stat.perror( "create ia_XMLFilePath attribute" );
            return stat;
        }
        tAttr.setWritable( true );
        tAttr.setReadable( false );
        tAttr.setConnectable( true );
        stat = addAttribute( ia_XMLFilePath );
        if (!stat) { stat.perror( "addAttribute ia_XMLFilePath" ); return stat; }
    }
    stat = attributeAffects( ia_XMLFilePath, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_XMLFilePath->oa_rodsChanged" ); return stat; }
    
    {
        MFnNumericAttribute nAttr;
        ia_writeToXMLFile = nAttr.create( "writeToXMLFile", "wtx", MFnNumericData::kBoolean, false, &stat );
        if (!stat) {
            stat.perror("create ia_drawSubSteppedVertices attribute");
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_writeToXMLFile );
        if (!stat) { stat.perror( "addAttribute ia_writeToXMLFile" ); return stat; }
    }
    stat = attributeAffects( ia_writeToXMLFile, ca_syncAttrs );
    if ( !stat ) { stat.perror( "attributeAffects ia_writeToXMLFile->oa_rodsChanged" ); return stat; }
    

    
    return MS::kSuccess;
}
