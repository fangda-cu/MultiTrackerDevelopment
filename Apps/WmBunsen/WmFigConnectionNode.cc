#include "WmFigConnectionNode.hh"
#include "WmFigRodNode.hh"

#include <maya/MPlugArray.h>

using namespace BASim;

// Required by Maya to identify node
/* static */ MTypeId WmFigConnectionNode::typeID ( 0x001135, 0xDA ); 
/* static */ MString WmFigConnectionNode::typeName( "wmFigConnectionNode" );

// 
/* static */ MObject WmFigConnectionNode::ia_time;
/* static */ MObject WmFigConnectionNode::ia_startTime;

/* static */ MObject WmFigConnectionNode::ia_rodNumber;
/* static */ MObject WmFigConnectionNode::ia_edgeNumber;
/* static */ MObject WmFigConnectionNode::ia_transformMatrix;

/* static */ MObject WmFigConnectionNode::ia_controllingEdge;
/* static */ MObject WmFigConnectionNode::oa_outTransformMatrix;

/* static */ MObject WmFigConnectionNode::ia_rodEdgeTransforms;
/* static */ MObject WmFigConnectionNode::oa_edgeTransform;

// We output the material frame to the rod node rather than having it connect directly because
// if this node is deleted we want the frame connection to be deleted too. So this node 
// routes it through.
///* static */ MObject WmFigConnectionNode::oa_materialFrame;

WmFigConnectionNode::WmFigConnectionNode() : m_startTime( 1 ), m_currentTime( 1 ), 
    m_previousTime( 1 ), m_controlledRodIndex( 0 ), m_controlledEdgeIndex( 0 )
{
}

WmFigConnectionNode::~WmFigConnectionNode()
{
}

void WmFigConnectionNode::getControlledRodInfo( unsigned int& o_rodIndex, unsigned int& o_edgeIndex, 
                                                EdgeTransform& o_edgeTransform )
{
    o_rodIndex = m_controlledRodIndex;
    o_edgeIndex = m_controlledEdgeIndex;
    o_edgeTransform = m_edgeTransform;
}


MStatus WmFigConnectionNode::compute( const MPlug& i_plug, MDataBlock& i_dataBlock ) 
{
    MStatus stat;
    
    //cerr << "WmFigConnectionNode::compute plug = " << i_plug.name() << endl;
	
    /*if ( i_plug == oa_materialFrame )
    {
        m_previousTime = m_currentTime;
        m_currentTime = i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
        CHECK_MSTATUS( stat );
        
        m_startTime = i_dataBlock.inputValue( ia_startTime, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        
        m_controlledEdgeIndex = i_dataBlock.inputValue( ia_edgeNumber, &stat ).asInt();
        CHECK_MSTATUS( stat );

        m_controlledRodIndex = i_dataBlock.inputValue( ia_rodNumber, &stat ).asInt();
        CHECK_MSTATUS( stat );

        m_controllingEdge = i_dataBlock.inputValue( ia_controllingEdge, &stat ).asBool();
        CHECK_MSTATUS( stat );
        
        MDataHandle materialFrameH = i_dataBlock.outputValue( oa_materialFrame, &stat );
        CHECK_MSTATUS( stat );
                
        materialFrameH.set( transformMatrix );
        
        materialFrameH.setClean();
        i_dataBlock.setClean( i_plug );
    }
    else*/ if ( i_plug == oa_outTransformMatrix )
    {
        // If the output transform matrix is being asked for then
        m_controlledEdgeIndex = i_dataBlock.inputValue( ia_edgeNumber, &stat ).asInt();
        CHECK_MSTATUS( stat );

        m_controlledRodIndex = i_dataBlock.inputValue( ia_rodNumber, &stat ).asInt();
        CHECK_MSTATUS( stat );

        // pull this so maya knows we care but it has no data, we need to go mining for it
        // on the node itself
        i_dataBlock.inputValue( ia_rodEdgeTransforms, &stat ).asBool();
        CHECK_MSTATUS( stat );

        MPlug rodEdgePlug( thisMObject(), ia_rodEdgeTransforms );
        
        MPlugArray rodsPlugs;
        MMatrix outMatrix;
        if ( rodEdgePlug.connectedTo( rodsPlugs, true, false, &stat ) )
        {
            MPlug rodEdgePlug = rodsPlugs[0];
            MObject rodNodeObj = rodEdgePlug.node( &stat );
            CHECK_MSTATUS( stat );
            
            MFnDependencyNode rodNodeDepFn( rodNodeObj );
            WmFigRodNode* rodNode = static_cast<WmFigRodNode*>( rodNodeDepFn.userNode() );
            
            if ( rodNode != NULL )
            {
                outMatrix = rodNode->getRodEdgeMatrix( m_controlledRodIndex, m_controlledEdgeIndex );
            }
            else
                outMatrix.setToIdentity();
        }
        else
            outMatrix.setToIdentity();
        
        MDataHandle outMatrixH = i_dataBlock.outputValue( oa_outTransformMatrix, &stat );
        CHECK_MSTATUS( stat );
        outMatrixH.set( outMatrix );
        
        outMatrixH.setClean();
        i_dataBlock.setClean( i_plug );
    }
    else if ( i_plug == oa_edgeTransform ) 
    {
        // We don't actually send the data to the rod node as we'd also have to send the rod number
        // and index so the rod node follows this connection and asks us for the data. We need to
        // get the inputs to make sure the data is correct when we're asked for it.
        
        m_controlledRodIndex = i_dataBlock.inputValue( ia_rodNumber ).asInt();
        m_controlledEdgeIndex = i_dataBlock.inputValue( ia_edgeNumber ).asInt();
        m_inputTransformMatrix = i_dataBlock.inputValue( ia_transformMatrix ).asMatrix();

        // Work out the material frame to send to the rod node when we are asked for it.
        m_edgeTransform.materialFrame.m1 = Vec3d( m_inputTransformMatrix( 0, 0 ),
                                m_inputTransformMatrix( 0, 1 ), m_inputTransformMatrix( 0, 2 ) );
        m_edgeTransform.materialFrame.m2 = Vec3d( m_inputTransformMatrix( 1, 0 ),
                                m_inputTransformMatrix( 1, 1 ), m_inputTransformMatrix( 1, 2 ) );
        m_edgeTransform.materialFrame.m3 = Vec3d( m_inputTransformMatrix( 2, 0 ),
                                m_inputTransformMatrix( 2, 1 ), m_inputTransformMatrix( 2, 2 ) );
        
        m_edgeTransform.position = Vec3d( m_inputTransformMatrix( 3, 0 ),
                                m_inputTransformMatrix( 3, 1 ), m_inputTransformMatrix( 3, 2 ) );
        
        i_dataBlock.outputValue( oa_edgeTransform ).setClean();
        i_dataBlock.setClean( i_plug );
    }
    else
    {
		return MS::kUnknownParameter;
	}

	return MS::kSuccess;
}

void WmFigConnectionNode::draw( M3dView& i_view, const MDagPath& i_path,
                            M3dView::DisplayStyle i_style,
                            M3dView::DisplayStatus i_status )
{ 
	MStatus stat;
	MObject thisNode = thisMObject();

	
	i_view.beginGL(); 
	glPushAttrib( GL_CURRENT_BIT | GL_POINT_BIT | GL_LINE_BIT );
       
    glPopAttrib();
	i_view.endGL();
}

MStatus WmFigConnectionNode::connectionMade( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc ) 
{    
    MStatus stat;
    MStatus retVal(MS::kUnknownParameter );
    
    return retVal;
}

MStatus WmFigConnectionNode::connectionBroken( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc )
{
    MStatus stat;
    MStatus retVal( MS::kUnknownParameter );
    
    return retVal;
}

bool WmFigConnectionNode::isBounded() const
{ 
	return false;
}

void* WmFigConnectionNode::creator()
{
	return new WmFigConnectionNode();
}

/*static */ MStatus WmFigConnectionNode::addNumericAttribute( MObject& i_attribute, MString i_longName, 
    MString i_shortName, MFnNumericData::Type i_type, double i_defaultValue, bool i_isInput, bool i_isOutput )
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
    if ( i_isOutput )
        nAttr.setReadable( true );
    
    stat = addAttribute( i_attribute );
    if ( !stat ) { stat.perror( "addAttribute " + i_longName ); return stat; }

    return stat;
}

/* static */ MStatus WmFigConnectionNode::initialize()
{ 
    MStatus stat;

 /*   {
        MFnMatrixAttribute mAttr;
        oa_materialFrame = mAttr.create( "materialFrame", "mf", MFnMatrixAttribute::kDouble, &stat );
        if ( !stat ) 
        {
            stat.perror("create oa_materialFrame attribute");
            return stat;
        }
        CHECK_MSTATUS( mAttr.setWritable( false ) );
        CHECK_MSTATUS( mAttr.setReadable( true ) );
        stat = addAttribute( oa_materialFrame );
        if ( !stat ) { stat.perror( "addAttribute oa_materialFrame" ); return stat; }
    }*/
    
    {
        MFnMatrixAttribute mAttr;
        ia_transformMatrix = mAttr.create( "transformMatrix", "tm", MFnMatrixAttribute::kDouble, &stat );
        if ( !stat ) 
        {
            stat.perror("create ia_transformMatrix attribute");
            return stat;
        }
        CHECK_MSTATUS( mAttr.setWritable( true ) );
        CHECK_MSTATUS( mAttr.setReadable( false ) );
        stat = addAttribute( ia_transformMatrix );
        if ( !stat ) { stat.perror( "addAttribute ia_transformMatrix" ); return stat; }
    }
    
   /*{
        MFnUnitAttribute	uAttr;
        ia_time = uAttr.create( "time", "t", MTime( 0.0 ), &stat );
        if ( !stat ) 
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
	//stat = attributeAffects( ia_time, oa_materialFrame );
	//if ( !stat ) { stat.perror( "attributeAffects ia_time->oa_materialFrame" ); return stat; }

    addNumericAttribute( ia_startTime, "startTime", "stt", MFnNumericData::kDouble, 1.0, true, false );
	stat = attributeAffects( ia_startTime, oa_materialFrame );
	if ( !stat ) { stat.perror( "attributeAffects ia_startTime->oa_materialFrame" ); return stat; }

    addNumericAttribute( ia_controllingEdge, "controllingEdge", "coe", MFnNumericData::kBoolean, false, true, false );
	//stat = attributeAffects( ia_controllingEdge, oa_materialFrame );
	//if ( !stat ) { stat.perror( "attributeAffects ia_controllingEdge->oa_materialFrame" ); return stat; }*/
    
    addNumericAttribute( ia_rodNumber, "rodNumber", "rn", MFnNumericData::kInt, 0, true, true );
	//stat = attributeAffects( ia_rodNumber, oa_materialFrame );
	//if ( !stat ) { stat.perror( "attributeAffects ia_rodNumber->oa_materialFrame" ); return stat; }
    
    addNumericAttribute( ia_edgeNumber, "edgeNumber", "en", MFnNumericData::kInt, 0, true, true );
	//stat = attributeAffects( ia_edgeNumber, oa_materialFrame );
	//if ( !stat ) { stat.perror( "attributeAffects ia_edgeNumber->oa_materialFrame" ); return stat; }
 
    {
        MFnMatrixAttribute mAttr;
        oa_outTransformMatrix = mAttr.create( "outTransformMatrix", "otm", MFnMatrixAttribute::kDouble, &stat );
        if ( !stat ) 
        {
            stat.perror("create oa_outTransformMatrix attribute");
            return stat;
        }
        CHECK_MSTATUS( mAttr.setWritable( false ) );
        CHECK_MSTATUS( mAttr.setReadable( true ) );
        stat = addAttribute( oa_outTransformMatrix );
        if ( !stat ) { stat.perror( "addAttribute oa_outTransformMatrix" ); return stat; }
    }
    
    addNumericAttribute( ia_rodEdgeTransforms, "rodEdgeTransforms", "ret", MFnNumericData::kBoolean, false, true, false );
	stat = attributeAffects( ia_rodEdgeTransforms, oa_outTransformMatrix );
	if ( !stat ) { stat.perror( "attributeAffects ia_rodEdgeTransforms->oa_outTransformMatrix" ); return stat; }
    stat = attributeAffects( ia_rodNumber, oa_outTransformMatrix );
	if ( !stat ) { stat.perror( "attributeAffects ia_rodEdgeTransforms->oa_outTransformMatrix" ); return stat; }
    stat = attributeAffects( ia_edgeNumber, oa_outTransformMatrix );
	if ( !stat ) { stat.perror( "attributeAffects ia_rodEdgeTransforms->oa_outTransformMatrix" ); return stat; }
 
    addNumericAttribute( oa_edgeTransform, "edgeTransform", "oet", MFnNumericData::kBoolean, false, false, true );
	stat = attributeAffects( ia_rodNumber, oa_edgeTransform );
    stat = attributeAffects( ia_edgeNumber, oa_edgeTransform );
    stat = attributeAffects( ia_transformMatrix, oa_edgeTransform );
	
	return MS::kSuccess;
}
