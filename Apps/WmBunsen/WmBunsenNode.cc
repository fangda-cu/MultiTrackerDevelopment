#include "WmBunsenNode.hh"

///////////////////////////////////////////////////////////////////////////////////
// WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING
// WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING
// The below typeID is NOT A VALID WETA ID
// WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING
// WARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNINGWARNING
///////////////////////////////////////////////////////////////////////////////////

MTypeId WmBunsenNode::ia_typeID( 0x80006 ); 
MString WmBunsenNode::ia_typeName( "wmBunsenNode" );
MObject WmBunsenNode::ca_syncAttrs;
MObject WmBunsenNode::ia_time;
MObject WmBunsenNode::ia_startTime;

WmBunsenNode::WmBunsenNode() : m_initialised( false ), m_beaker( NULL )
{
    m_beaker = new Beaker();
}

WmBunsenNode::~WmBunsenNode()
{
    if ( m_beaker != NULL )
        delete m_beaker;
}

MStatus WmBunsenNode::compute( const MPlug& i_plug, MDataBlock& i_dataBlock ) 
{
    MStatus stat;
	
    if ( i_plug == ca_syncAttrs )
    {
        m_previousTime = m_currentTime;
        m_currentTime = i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
		CHECK_MSTATUS( stat );
        
	    m_startTime = i_dataBlock.inputValue( ia_startTime, &stat ).asDouble();
		CHECK_MSTATUS( stat );
		
		if ( m_currentTime == m_startTime )
        {
			// reinitialise
		}
		else if ( m_currentTime > m_previousTime ) 
        {
            if ( m_initialised )
   			{
                // take a step
                m_beaker->takeTimeStep();
            }
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
	} else
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
	glPushAttrib( GL_CURRENT_BIT | GL_POINT_BIT | GL_LINE_BIT );
    
    if ( !m_initialised )
    {
        RodOptions opts;
        opts.numVertices = 50;
        opts.YoungsModulus = 1000.0;
        opts.ShearModulus = 375.0;
        opts.radiusA = 0.5;
        opts.radiusB = 1.0;
        std::string frame = "time";
        if (frame == "time") opts.refFrame = ElasticRod::TimeParallel;
        else if (frame == "space") opts.refFrame = ElasticRod::SpaceParallel;

        Scalar radius = 20.0;

        std::vector<Vec3d> vertices, undeformed;
        for (int i = 0; i < opts.numVertices; ++i) 
        {
            vertices.push_back(Vec3d(radius * cos(i * M_PI / (opts.numVertices - 1)),
                                     radius * sin(i * M_PI / (opts.numVertices - 1)),
                                     0));
        }
        
        m_beaker->addRod( vertices,
                          vertices,
                          opts );
        
        m_initialised = true;
    }
    else
        m_beaker->display();
    
	// draw dynamic Hair

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

	return MS::kSuccess;
}
