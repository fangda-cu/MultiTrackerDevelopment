#include "WmBunsenRodNode.hh"

///////////////////////////////////////////////////////////////////////////////////
// WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING
// WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING
// The below typeID is NOT A VALID WETA ID
// WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING
// WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING
///////////////////////////////////////////////////////////////////////////////////

MTypeId WmBunsenRodNode::ia_typeID( 0x80007 ); 
MString WmBunsenRodNode::ia_typeName( "wmBunsenRodNode" );
MObject WmBunsenRodNode::ia_time;
MObject WmBunsenRodNode::ia_startTime;
MObject WmBunsenRodNode::ia_nurbsCurves;

MObject WmBunsenRodNode::ca_syncAttrs;
MObject WmBunsenRodNode::oa_rodsChanged;

WmBunsenRodNode::WmBunsenRodNode() : m_initialised( false ), mx_rodData( NULL ), mx_world( NULL )
{
}

WmBunsenRodNode::~WmBunsenRodNode()
{
}

void WmBunsenRodNode::initialiseRodData( vector<RodData*>* i_rodData )
{
    mx_rodData = i_rodData;
    
    updateRodDataFromInputs();
    
    // Since we're initialising the rod data that means the rod is about to be created. In which
    // case we need to set the current vertex positions since they will not get set by the
    // simulation until the user moves forward a frame.
    
    size_t numRods = mx_rodData->size();
    for ( size_t r=0; r<numRods; r++ )
    {
        size_t numVertices = (*mx_rodData)[ r ]->undeformedVertexPositions.size();
        (*mx_rodData)[ r ]->initialVertexPositions.resize( numVertices );
        for ( size_t v=0; v<numVertices; v++ )
        {
            (*mx_rodData)[ r ]->initialVertexPositions[ v ] = (*mx_rodData)[ r ]->undeformedVertexPositions[ v ];
        }
    }
}

void WmBunsenRodNode::updateRodDataFromInputs()
{    
    // just update a single imaginary rod to test...
 
    //mx_rodData->resize( 1 );
    
    // Create the data for the rods here, it is deleted by beaker when it is finished with it.
    // That's because this class has no control over when it is deleted and cleaning things up in
    // Beaker. The user may kill us in Maya at any time so Beaker/WmBunsenNode take responsibility
    // for cleaning up.
    
    (*mx_rodData)[0]->rodOptions.numVertices = 50;
    (*mx_rodData)[0]->rodOptions.YoungsModulus = 1000.0;
    (*mx_rodData)[0]->rodOptions.ShearModulus = 375.0;
    (*mx_rodData)[0]->rodOptions.radiusA = 0.5;
    (*mx_rodData)[0]->rodOptions.radiusB = 1.0;
    std::string frame = "time";
    if (frame == "time") (*mx_rodData)[0]->rodOptions.refFrame = ElasticRod::TimeParallel;
    else if (frame == "space") (*mx_rodData)[0]->rodOptions.refFrame = ElasticRod::SpaceParallel;

    Scalar radius = 20.0;

    //std::vector<Vec3d> vertices, undeformed;
    (*mx_rodData)[ 0 ]->undeformedVertexPositions.clear();
    for (int i = 0; i < (*mx_rodData)[0]->rodOptions.numVertices; ++i) 
    {
        (*mx_rodData)[ 0 ]->undeformedVertexPositions.push_back( Vec3d(radius * cos(i * M_PI / ((*mx_rodData)[0]->rodOptions.numVertices - 1)),
                                 radius * sin(i * M_PI / ((*mx_rodData)[0]->rodOptions.numVertices - 1)),
                                 0) );
    }
}

MStatus WmBunsenRodNode::compute( const MPlug& i_plug, MDataBlock& i_dataBlock ) 
{
    MStatus stat;
	
    if ( i_plug == ca_syncAttrs )
    {
        m_previousTime = m_currentTime;
        m_currentTime = i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
		CHECK_MSTATUS( stat );
        
	    m_startTime = i_dataBlock.inputValue( ia_startTime, &stat ).asDouble();
		CHECK_MSTATUS( stat );
		
        updateRodDataFromInputs();
        
		/*if ( m_currentTime == m_startTime )
        {
			// reinitialise
		}
		else if ( m_currentTime > m_previousTime ) 
        {
            if ( m_initialised )
   			{
                // take a step
            }
    	}*/

		MDataHandle outputData = i_dataBlock.outputValue ( ca_syncAttrs, &stat );
		if ( !stat )
        {
			stat.perror("WmBunsenRodNode::compute get ca_syncAttrs");
			return stat;
		}
		
        // We don't even need to put anything in the output handle as nothing uses it.
        // Just tell Maya it's clean so it doesn't repeatedly evaluate it.

		stat = i_dataBlock.setClean( i_plug );
		if ( !stat )
        {
			stat.perror("WmBunsenRodNode::compute setClean");
			return stat;
		}
	} 
    else if (  i_plug == oa_rodsChanged )
    {
        if ( mx_rodData != NULL )
        {
            // The Bunsen node is asking us to update the rod data in Beaker for our inputs
            // so do so here.....
        }
    }
    else
    {
		return MS::kUnknownParameter;
	}

	return MS::kSuccess;
}

void WmBunsenRodNode::draw( M3dView& i_view, const MDagPath& i_path,
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
		stat.perror( "WmBunsenRodNode::draw getting ca_syncAttrs" );
		return;
	}

	i_view.beginGL(); 
	glPushAttrib( GL_CURRENT_BIT | GL_POINT_BIT | GL_LINE_BIT );
    
	// draw dynamic Hair

	// What did this line do? it was here from the devkit example. Is it to with point colouring
	//view.setDrawColor ( WmBunsenRodNode );

	glPopAttrib();
	i_view.endGL();
}

bool WmBunsenRodNode::isBounded() const
{ 
	return false;
}

void* WmBunsenRodNode::creator()
{
	return new WmBunsenRodNode();
}

MStatus WmBunsenRodNode::initialize()
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
        MFnNumericAttribute	nAttr;
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
        MFnTypedAttribute   tAttr;  
        ia_nurbsCurves = tAttr.create( "nurbsCurves", "nc",
                                       MFnNurbsCurveData::kNurbsCurve, &stat );
        CHECK_MSTATUS( stat );
        CHECK_MSTATUS( tAttr.setArray(true) );
        CHECK_MSTATUS( tAttr.setReadable(false) );
        CHECK_MSTATUS( tAttr.setWritable(true) );
        CHECK_MSTATUS( tAttr.setUsesArrayDataBuilder(true) );
        stat = addAttribute( ia_nurbsCurves );
        CHECK_MSTATUS( stat );
    }
    
    {
        MFnNumericAttribute	 nAttr;
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
        oa_rodsChanged = nAttr.create( "rodsChanged", "rch", MFnNumericData::kBoolean, true , &stat );
        CHECK_MSTATUS ( stat );
        CHECK_MSTATUS( nAttr.setWritable( false ) );
        CHECK_MSTATUS( nAttr.setReadable( true ) );
        CHECK_MSTATUS( nAttr.setConnectable( true ) );
        stat = addAttribute( oa_rodsChanged );
        if ( !stat ) { stat.perror( "addAttribute oa_rodsChanged" ); return stat;}
    }
    
	stat = attributeAffects( ia_time, ca_syncAttrs );
	if ( !stat ) { stat.perror( "attributeAffects ia_time->ca_syncAttrs" ); return stat; }
	stat = attributeAffects( ia_startTime, ca_syncAttrs );
	if ( !stat ) { stat.perror( "attributeAffects ia_startTimer->ca_syncAttrs" ); return stat; }

    stat = attributeAffects( ia_nurbsCurves, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_nurbsCurves->oa_rodsChanged" ); return stat; }
    
	return MS::kSuccess;
}
