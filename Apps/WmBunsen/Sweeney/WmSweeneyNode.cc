#include "WmSweeneyNode.hh"
#include "WmFigConnectionNode.hh"

#include <maya/MFnMatrixAttribute.h>
#include <maya/MPlugArray.h>

#include <sys/stat.h>

using namespace BASim;


// Required by Maya to identify the node
/* static */ MTypeId WmSweeneyNode::typeID ( 0x001135, 0xF6  );
/* static */ MString WmSweeneyNode::typeName( "wmSweeneyNode" );

// Input attributes
/* static */ MObject WmSweeneyNode::ia_time;
/* static */ MObject WmSweeneyNode::ia_startTime;

// Barbershop specific inputs
/*static*/ MObject WmSweeneyNode::ia_strandVertices;
/*static*/ MObject WmSweeneyNode::ia_verticesPerStrand;

// Sync attributes
/* static */ MObject WmSweeneyNode::ca_sync;

WmSweeneyNode::WmSweeneyNode() : m_rodManager( NULL )
{    
}

WmSweeneyNode::~WmSweeneyNode()
{ 
}

MStatus WmSweeneyNode::compute( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    MStatus status;

    if ( i_plug == ca_sync )
    {
        m_currentTime = i_dataBlock.inputValue( ia_time ).asTime().value();
        m_startTime = i_dataBlock.inputValue( ia_startTime ).asDouble();
        
        MObject strandVerticesObj = i_dataBlock.inputValue( ia_strandVertices ).data();
        MFnVectorArrayData strandVerticesArrayData( strandVerticesObj, &status );
        CHECK_MSTATUS( status );
        
        m_strandVertices = strandVerticesArrayData.array( &status );
        CHECK_MSTATUS( status );
        
        m_numberOfVerticesPerStrand = i_dataBlock.inputValue( ia_verticesPerStrand ).asInt();        
        
        // Input has changed so rebuild rods.
        if ( m_currentTime == m_startTime )
        {
            initialiseRodFromBarberShopInput();            
        }
        else
        {
            if ( m_rodManager == NULL )
            {
                MGlobal::displayError( "Please rewind to start time to initialise simulation." );
            }
            
            m_rodManager->takeStep();
        }
        
        i_dataBlock.setClean( i_plug );
    }
    else
    {
        return MS::kUnknownParameter;
    }

    return MS::kSuccess;
}


void WmSweeneyNode::initialiseRodFromBarberShopInput()
{
    //TODO: Reset the manager and remove all rods before adding more
    delete m_rodManager;
    
    m_rodManager = new WmSweeneyRodManager();

    // Create one rod for each barbershop strand. Ignore the strand shape or length but do
    // take its initial direction as a follicle angle
    size_t currentVertexIndex = 0;
    
    // This should be a user parameter or the spacing should be?
    size_t verticesPerRod = 10;
    
    while ( currentVertexIndex < ( m_strandVertices.length() - 2 ) )
    {
        MVector direction = m_strandVertices[ currentVertexIndex + 1 ] 
                            - m_strandVertices[ currentVertexIndex ];
        direction.normalize();
        
        MVector currentVertex( m_strandVertices[ currentVertexIndex ] );
        
        vector< BASim::Vec3d > vertices;
        vertices.resize( verticesPerRod );
    
        for ( size_t v = 0; v < verticesPerRod; ++v )
        {
            vertices[ v ] = BASim::Vec3d( currentVertex.x, currentVertex.y, currentVertex.z );
            
            // this should be scaled so the rod ends up being the length the user asks for
            currentVertex += direction;
        }
        
        m_rodManager->addRod( vertices, m_currentTime );
        
        currentVertexIndex += m_numberOfVerticesPerStrand;
    }
    
    m_rodManager->initialiseSimulation( 1 / 24.0, 1.0 );
}


void WmSweeneyNode::draw( M3dView& i_view, const MDagPath& i_path,
                            M3dView::DisplayStyle i_style,
                            M3dView::DisplayStatus i_status )
{
    MStatus stat;
    MObject thisNode = thisMObject();

    MPlug syncPlug( thisNode, ca_sync );
    double d;
    stat = syncPlug.getValue( d );
    if ( !stat )
    {
        stat.perror( "WmSweeneyNode::draw getting ca_sync" );
        return;
    }


    i_view.beginGL();
    //glPushAttrib( GL_CURRENT_BIT | GL_POINT_BIT | GL_LINE_BIT );
	glPushAttrib( GL_ALL_ATTRIB_BITS );

	// draw dynamic Hair
    m_rodManager->drawAllRods();

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

/* static */ MStatus WmSweeneyNode::initialize()
{
    MStatus status;

    addNumericAttribute( ca_sync, "sync", "syn", MFnNumericData::kBoolean, false, false );

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
    status = attributeAffects( ia_time, ca_sync );
	if ( !status ) { status.perror( "attributeAffects ia_time->ca_sync" ); return status; }    
    
	addNumericAttribute( ia_startTime, "startTime", "stt", MFnNumericData::kDouble, 1.0, true );
	status = attributeAffects( ia_startTime, ca_sync );
	if ( !status ) { status.perror( "attributeAffects ia_startTime->ca_sync" ); return status; }
    
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
    
    addNumericAttribute( ia_verticesPerStrand, "verticesPerStrand", "vps", MFnNumericData::kInt, 12, true );
	status = attributeAffects( ia_verticesPerStrand, ca_sync );
	if ( !status ) { status.perror( "attributeAffects verticesPerStrand->ca_sync" ); return status; }

    return MS::kSuccess;

}



