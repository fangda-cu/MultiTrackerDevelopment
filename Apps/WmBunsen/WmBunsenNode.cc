#include "WmBunsenNode.hh"

MTypeId WmBunsenNode::typeID( 0x001135, 0x18 ); 
MString WmBunsenNode::typeName( "wmBunsenNode" );
MObject WmBunsenNode::ca_syncAttrs;
MObject WmBunsenNode::ia_time;
MObject WmBunsenNode::ia_fps;
MObject WmBunsenNode::ia_maxDt;
MObject WmBunsenNode::ia_maxIter;
MObject WmBunsenNode::ia_startTime;
MObject WmBunsenNode::ia_rodsNodes;
MObject WmBunsenNode::ia_gravity;
MObject WmBunsenNode::ia_numberOfThreads;
MObject WmBunsenNode::ia_solver;
MObject WmBunsenNode::oa_simStepTaken;

WmBunsenNode::WmBunsenNode() : m_initialised( false ), m_beaker( NULL )
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

void WmBunsenNode::createRodDataFromRodNodes( MDataBlock& i_dataBlock, 
    ObjectControllerBase::SolverLibrary solverLibrary )
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
                WmBunsenRodNode* wmBunsenRodNode = ( WmBunsenRodNode* )rodNodeFn.userNode();
                
                // Since the rod node is purely there to fill in data that comes from its inputs
                // and attributes, we don't let it deal with memory allocation. This node is in 
                // charge of all that.
                m_beaker->createSpaceForRods( r, wmBunsenRodNode->numberOfRods() );
                
                wmBunsenRodNode->initialiseRodData( m_beaker->rodData( r ) );

                // Now the rod node has used initialised the undeformed positions for the rods
                // it owns. Since we are resetting the sim we need to actually now create the
                // rods and add them to the world.
                
                m_beaker->createRods( r, solverLibrary );
            }
            else
                CHECK_MSTATUS( stat );
        }
    }
}

MStatus WmBunsenNode::compute( const MPlug& i_plug, MDataBlock& i_dataBlock ) 
{
    MStatus stat;
    
   // cerr << "WmBunsenNode::compute plug = " << i_plug.name() << endl;
	
    if ( i_plug == ca_syncAttrs )
    {
        m_previousTime = m_currentTime;
        m_currentTime = i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
        CHECK_MSTATUS( stat );

        m_startTime = i_dataBlock.inputValue( ia_startTime, &stat ).asDouble();
        CHECK_MSTATUS( stat );

        m_framedt = 1.0 / i_dataBlock.inputValue(ia_fps, &stat).asDouble();
        CHECK_MSTATUS( stat );

        m_beaker->setDt(i_dataBlock.inputValue(ia_maxDt, &stat ).asDouble());
        CHECK_MSTATUS( stat );	

        const double3 &gravity = i_dataBlock.inputValue( ia_gravity, &stat ).asDouble3();
        CHECK_MSTATUS( stat );
        m_beaker->setGravity( Vec3d( gravity[0], gravity[1], gravity[2] ) );

    	m_beaker->setMaxIter(i_dataBlock.inputValue(ia_maxIter, &stat).asInt());
        CHECK_MSTATUS( stat );
    
        int numberOfThreads = i_dataBlock.inputValue( ia_numberOfThreads, &stat ).asInt();
        CHECK_MSTATUS( stat );
        
        int solver = i_dataBlock.inputValue( ia_solver, &stat ).asInt();
        
        if ( m_currentTime == m_startTime )
        {
            m_beaker->resetEverything();
            
            ObjectControllerBase::SolverLibrary solverLibrary;
            if ( solver == 0 )
                solverLibrary = ObjectControllerBase::MKL_SOLVER;
            else
                solverLibrary = ObjectControllerBase::PETSC_SOLVER;
            
            createRodDataFromRodNodes( i_dataBlock, solverLibrary );
        }

        pullOnAllRodNodes( i_dataBlock );
        
        if ( m_currentTime > m_previousTime ) 
        {
            // take a step of size 1.0/24.0
            m_beaker->takeTimeStep( numberOfThreads, m_framedt ); 
        }
    
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
    
    m_beaker->draw();

	// What did this line do? it was here from the devkit example. Is it to with point colouring
	//view.setDrawColor ( WmBunsenNode );

	glPopAttrib();
	i_view.endGL();
}

bool WmBunsenNode::isBounded() const
{ 
    return false;
}

void* WmBunsenNode::creator()
{
    return new WmBunsenNode();
}

MStatus WmBunsenNode::initialize()
{ 
    MStatus stat;
    
    {
        MFnUnitAttribute	uAttr;
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
        nAttr.setReadable( false );
        nAttr.setKeyable( true );  
        stat = addAttribute( ia_startTime );
        if ( !stat ) { stat.perror( "addAttribute ia_startTime" ); return stat; }
    }

    {
	MFnNumericAttribute nAttr;
    	ia_fps = nAttr.create( "framesPerSecond", "fps", MFnNumericData::kDouble, 24.0, &stat );
        if ( !stat ) 
        {
            stat.perror( "create framesPerSecond attribute");
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
        nAttr.setReadable( false );
        nAttr.setConnectable( true );
        stat = addAttribute( ia_gravity );
        if (!stat) { stat.perror( "addAttribute ia_gravity" ); return stat; }
    }
 
   {
        MFnNumericAttribute nAttr;
        ia_solver = nAttr.create( "solver", "sol", MFnNumericData::kInt, 0, &stat );
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

    stat = attributeAffects( ia_time, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_time->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_startTime, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_startTimer->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_fps, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_fps->ca_syncAttrs" );return stat;}
    stat = attributeAffects( ia_maxDt, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_maxDt->ca_syncAttrs" );return stat;}
    stat = attributeAffects( ia_maxIter, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_maxIter->ca_syncAttrs" );return stat;}
    stat = attributeAffects( ia_rodsNodes, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_rodsNodes->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_gravity, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_rodsNodes->ca_syncAttrs" ); return stat; }
	stat = attributeAffects( ia_solver, ca_syncAttrs );
    if (!stat) { stat.perror( "attributeAffects ia_rodsNodes->ca_syncAttrs" ); return stat; }
	stat = attributeAffects( ia_numberOfThreads, ca_syncAttrs );
	if (!stat) { stat.perror( "attributeAffects ia_numberOfThreads->ca_syncAttrs" ); return stat; }
    
    stat = attributeAffects( ia_time, oa_simStepTaken );
    if (!stat) { stat.perror( "attributeAffects ia_time->oa_simulatedRods" ); return stat; }
    
    return MS::kSuccess;
}
