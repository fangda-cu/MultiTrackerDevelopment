#include "WmSweeneyNode.hh"
#include "WmFigConnectionNode.hh"
#include "../WmBunsenCollisionMeshNode.hh"

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

// Hair Property Attributes
/* static */ MObject WmSweeneyNode::ia_length;
/* static */ MObject WmSweeneyNode::ia_edgeLength;
/* static */ MObject WmSweeneyNode::ia_verticesPerRod;

// Barbershop specific inputs
/*static*/ MObject WmSweeneyNode::ia_strandVertices;
/*static*/ MObject WmSweeneyNode::ia_verticesPerStrand;

// Sync attributes
/* static */ MObject WmSweeneyNode::ca_timeSync;
/* static */ MObject WmSweeneyNode::ca_rodPropertiesSync;

// Collision meshes
/* static */ MObject WmSweeneyNode::ia_collisionMeshes;
 
WmSweeneyNode::WmSweeneyNode() : m_rodManager( NULL )
{    
}

WmSweeneyNode::~WmSweeneyNode()
{ 
}

MStatus WmSweeneyNode::compute( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    MStatus status;

    cerr << "WmSweeneyNode::compute() with i_plug = " << i_plug.name() << endl;

    if ( i_plug == ca_rodPropertiesSync )
    {
        // Hair properties
        m_length = i_dataBlock.inputValue( ia_length ).asDouble();
        m_edgeLength = i_dataBlock.inputValue( ia_edgeLength ).asDouble();
        m_verticesPerRod = i_dataBlock.inputValue( ia_verticesPerRod ).asInt();
        
        cerr << "m_verticesPerRod = " << m_verticesPerRod << endl;
        
        MObject strandVerticesObj = i_dataBlock.inputValue( ia_strandVertices ).data();
        MFnVectorArrayData strandVerticesArrayData( strandVerticesObj, &status );
        CHECK_MSTATUS( status );
        
        m_strandVertices = strandVerticesArrayData.array( &status );
        CHECK_MSTATUS( status );
        
        m_numberOfVerticesPerStrand = i_dataBlock.inputValue( ia_verticesPerStrand ).asInt();        
        
        // If we've not yet created all the rods then create them
        if ( m_rodManager == NULL )
        {
            initialiseRodFromBarberShopInput( i_dataBlock );
        }
        else
        {
            // If the rods exist then just update their undeformed configuration but keep running
            // the simulation.
            
            // TODO: Add code to update undeformed configuration 
            // for just now, recreate the rod
            initialiseRodFromBarberShopInput( i_dataBlock );
        }
        
        i_dataBlock.setClean( i_plug );
    }
    else if ( i_plug == ca_timeSync )
    {
        m_currentTime = i_dataBlock.inputValue( ia_time ).asTime().value();
        m_startTime = i_dataBlock.inputValue( ia_startTime ).asDouble();

        if ( m_rodManager == NULL )
        {
        //    MGlobal::displayError( "Please rewind to start time to initialise simulation." );            
        }
        else if ( m_currentTime == m_startTime )
        {            
            initialiseRodFromBarberShopInput( i_dataBlock );
        }
        else
        {
            updateCollisionMeshes( i_dataBlock );
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
                m_rodManager->addCollisionMesh( triangleMesh, figMeshController );
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
    
    vector< BASim::Vec3d > vertices;
       
    for ( unsigned int inputStrandNumber = 0; inputStrandNumber < numberOfStrands; ++inputStrandNumber )
    {
        MVector direction = m_strandVertices[ currentVertexIndex + 1 ] 
                            - m_strandVertices[ currentVertexIndex ];
        direction.normalize();
                
        constructRodVertices( vertices, direction, m_strandVertices[ currentVertexIndex ] );
                    
        m_rodManager->addRod( vertices, m_startTime );
        
        cerr << "Creating rod at time " << m_startTime << endl;
        
        currentVertexIndex += m_numberOfVerticesPerStrand;        
    }
  
/*
    // For debugging add rods with specific vertices to see problems in Beaker vs Sweeney
    vertices.resize( 15 );
    vertices[ 0 ] = BASim::Vec3d( 1.1468563, 4.7516446, 5.1037321 );
    vertices[ 1 ] = BASim::Vec3d( 1.3113695, 5.4118051, 5.836619 );
    vertices[ 2 ] = BASim::Vec3d( 1.4758841, 6.0719659, 6.5695057);
    vertices[ 3 ] = BASim::Vec3d( 1.6404004, 6.7321264, 7.3023923);
    vertices[ 4 ] = BASim::Vec3d( 1.8049186, 7.3922864, 8.0352787);
    vertices[ 5 ] = BASim::Vec3d( 1.9694386, 8.0524464, 8.768165);
    vertices[ 6 ] = BASim::Vec3d( 2.1339601, 8.7126063, 9.5010505);
    vertices[ 7 ] = BASim::Vec3d( 2.2984837, 9.372766, 10.233936);
    vertices[ 8 ] = BASim::Vec3d( 2.4630091, 10.032925, 10.966822);
    vertices[ 9 ] = BASim::Vec3d( 2.6275362, 10.693084, 11.699707);
    vertices[ 10 ] = BASim::Vec3d( 2.7920654, 11.353243, 12.432592);
    vertices[ 11 ] = BASim::Vec3d( 2.9565963, 12.013402, 13.165476);
    vertices[ 12 ] = BASim::Vec3d( 3.1211287, 12.67356, 13.898361);
    vertices[ 13 ] = BASim::Vec3d( 3.2856628, 13.333719, 14.631245);
    vertices[ 14 ] = BASim::Vec3d( 3.4501989, 13.993877, 15.364129);
    
    m_rodManager->addRod( vertices, m_startTime );
    
    vertices[ 0 ] = BASim::Vec3d( 5.2480741, -0.48017403, 4.6308284 );
    vertices[ 1 ] = BASim::Vec3d( 5.9961056, -0.54642217, 5.2911769 );
    vertices[ 2 ] = BASim::Vec3d( 6.7441311, -0.61267031, 5.951532 );
    vertices[ 3 ] = BASim::Vec3d( 7.4921505, -0.67891851, 6.6118941 );
    vertices[ 4 ] = BASim::Vec3d( 8.240163, -0.74516679, 7.2722641 );
    vertices[ 5 ] = BASim::Vec3d( 8.9881685, -0.8114151, 7.9326418 );
    vertices[ 6 ] = BASim::Vec3d( 9.7361675, -0.87766339, 8.5930272 );
    vertices[ 7 ] = BASim::Vec3d( 10.484159, -0.9439117, 9.2534205 );
    vertices[ 8 ] = BASim::Vec3d( 11.232145, -1.0101601, 9.9138215 );
    vertices[ 9 ] = BASim::Vec3d( 11.980123, -1.0764085, 10.57423 );
    vertices[ 10 ] = BASim::Vec3d( 12.728095, -1.1426569, 11.234647 );
    vertices[ 11 ] = BASim::Vec3d( 13.476059, -1.2089054, 11.895071 );
    vertices[ 12 ] = BASim::Vec3d( 14.224018, -1.2751539, 12.555502 );
    vertices[ 13 ] = BASim::Vec3d( 14.971969, -1.3414024, 13.215941 );
    vertices[ 14 ] = BASim::Vec3d( 15.719913, -1.4076509, 13.876389 );
    
    m_rodManager->addRod( vertices, m_startTime );*/
  
    cerr << "initialiseRodFromBarberShopInput() - About to initialise simulation\n";
    m_rodManager->initialiseSimulation( 1 / 24.0, m_startTime );
    cerr << "initialiseRodFromBarberShopInput() - Simulation initialised at time " << m_startTime << endl;
}

void WmSweeneyNode::constructRodVertices( vector< BASim::Vec3d >& o_rodVertices, const MVector& i_direction,
                                  const MVector& i_rootPosition )
{
    cerr << "constructRodVertices(): About to construct rod vertices\n";
    cerr << "constructRodVertices(): m_verticesPerRod = " << m_verticesPerRod << "\n";
    
    // Construct a list of vertices for a rod with its root at i_rootPosition and going in direction
    // i_direction
    
    o_rodVertices.clear();
    
    MVector edge = i_direction;
    edge.normalize();
    
    edge *= m_length / ( m_verticesPerRod - 1 );
    
    cerr << "constructRodVertices(): edgeLength = " << edge.length() << "\n";    
    
    MVector currentVertex( i_rootPosition );
    
    // For now we only do straight hair...
    // 
    for ( int v = 0; v < m_verticesPerRod; ++v )
    {
        o_rodVertices.push_back( BASim::Vec3d( currentVertex.x, currentVertex.y, currentVertex.z ) );
        
        currentVertex += edge;
    }
    
    cerr << "constructRodVertices(): Finished constructing rod vertices\n";    
}

void WmSweeneyNode::draw( M3dView& i_view, const MDagPath& i_path,
                            M3dView::DisplayStyle i_style,
                            M3dView::DisplayStatus i_status )
{
    MStatus status;
    
    // Pull on the sync plugs to cause compute() to be called if any 
    // of the rod properties or time has changed.
    double d;
    
   /* if ( m_rodManager == NULL )
    {    
        // If this is the case we probably just loaded, so make sure the 
        // rod properties are loaded and stored before we call timeSync to
        // initialised everything
        MPlug propertiesSyncPlug( thisMObject(), ca_rodPropertiesSync );
        propertiesSyncPlug.getValue( d );
    }*/
    
    MPlug timeSyncPlug( thisMObject(), ca_timeSync );
    timeSyncPlug.getValue( d );
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

    addNumericAttribute( ca_timeSync, "timeSync", "syn", MFnNumericData::kBoolean, false, false );
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
    status = attributeAffects( ia_time, ca_timeSync );
	if ( !status ) { status.perror( "attributeAffects ia_time->ca_timeSync" ); return status; }    
    
	addNumericAttribute( ia_startTime, "startTime", "stt", MFnNumericData::kDouble, 1.0, true );
	status = attributeAffects( ia_startTime, ca_timeSync );
	if ( !status ) { status.perror( "attributeAffects ia_startTime->ca_timeSync" ); return status; }

	addNumericAttribute( ia_length, "length", "len", MFnNumericData::kDouble, 10.0, true );
	status = attributeAffects( ia_length, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_length->ca_rodPropertiesSync" ); return status; }

	addNumericAttribute( ia_edgeLength, "edgeLength", "ele", MFnNumericData::kDouble, 1.0, true );
	status = attributeAffects( ia_edgeLength, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_edgeLength->ca_rodPropertiesSync" ); return status; }

    addNumericAttribute( ia_verticesPerRod, "verticesPerRod", "cpr", MFnNumericData::kInt, 10, true );
	status = attributeAffects( ia_verticesPerRod, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects ia_verticesPerRod->ca_rodPropertiesSync" ); return status; }
    
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
    status = attributeAffects( ia_collisionMeshes, ca_timeSync );
	if ( !status ) { status.perror( "attributeAffects ia_collisionMeshes->ca_timeSync" ); return status; }
    
        
    addNumericAttribute( ia_verticesPerStrand, "verticesPerStrand", "vps", MFnNumericData::kInt, 12, true );
	status = attributeAffects( ia_verticesPerStrand, ca_rodPropertiesSync );
	if ( !status ) { status.perror( "attributeAffects verticesPerStrand->ca_rodPropertiesSync" ); return status; }

    return MS::kSuccess;
}



