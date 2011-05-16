#include <math.h>           
#include <maya/MIOStream.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MTime.h>
#include <maya/MFnVectorArrayData.h>

#include "WmSweeneyShape.hh"         
#include "WmSweeneyShapeUI.hh"     
#include "WmSweeneyShapeIterator.hh"

/*static*/ MObject WmSweeneyShape::ia_strandVertices;
/*static*/ MObject WmSweeneyShape::ia_verticesPerStrand;
/*static*/ MObject WmSweeneyShape::ia_time;
/*static*/ MObject WmSweeneyShape::ca_sync;

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// Shape implementation
//
////////////////////////////////////////////////////////////////////////////////

MTypeId WmSweeneyShape::id( 0x001135, 0xF6 );
MString WmSweeneyShape::typeName( "wmSweeneyShape" );

WmSweeneyShape::WmSweeneyShape() 
{    
    initialiseRodFromBarberShopInput();
}

WmSweeneyShape::~WmSweeneyShape() {}

void* WmSweeneyShape::creator()
//
// Description
//
//    Called internally to create a new instance of the users MPx node.
//
{
	return new WmSweeneyShape();
}

MStatus WmSweeneyShape::initialize()
//
// Description
//
//    Attribute (static) initialization.
//    See api_macros.h.
//
{
    MStatus status;
    
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
        MFnNumericAttribute numericAttrFn;
        ia_verticesPerStrand = numericAttrFn.create( "verticesPerStrand", "vps", 
                                                     MFnNumericData::kInt, 12, &status );
        CHECK_MSTATUS( status );
        numericAttrFn.setWritable( true );
        numericAttrFn.setReadable( false );
        numericAttrFn.setConnectable( true );
        status = addAttribute( ia_verticesPerStrand );
        CHECK_MSTATUS( status );
    }
    
    {
        MFnUnitAttribute unitAttrFn;
        ia_time = unitAttrFn.create( "time", "t", MTime( 0.0 ), &status );
        CHECK_MSTATUS( status );
        unitAttrFn.setWritable( true );
        unitAttrFn.setConnectable( true );
        unitAttrFn.setStorable( false );
        status = addAttribute( ia_time );
        CHECK_MSTATUS( status );
    }
    
    {
        MFnNumericAttribute numericAttrFn;
        ca_sync = numericAttrFn.create( "sync", "snc", MFnNumericData::kBoolean, false, &status );
        CHECK_MSTATUS( status );
        numericAttrFn.setWritable( false );
        numericAttrFn.setReadable( true );
        numericAttrFn.setConnectable( true );
        numericAttrFn.setKeyable( false );
        status = addAttribute( ca_sync );
        CHECK_MSTATUS( status );
    }
    
    attributeAffects( ia_strandVertices, ca_sync );
    attributeAffects( ia_verticesPerStrand, ca_sync );
    attributeAffects( ia_time, ca_sync );
    
    return MS::kSuccess;
}

MStatus WmSweeneyShape::compute( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    cerr << "WmSweeneyShape called with plug " << i_plug.name() << endl;

    MStatus status;

    if( i_plug == ca_sync )
    {
        i_dataBlock.inputValue( ca_sync, &status );
        CHECK_MSTATUS( status );

        m_currentTime = i_dataBlock.inputValue( ia_time ).asTime().value();
        
        MObject strandVerticesObj = i_dataBlock.inputValue( ia_strandVertices ).data();
        MFnVectorArrayData strandVerticesArrayData( strandVerticesObj, &status );
        CHECK_MSTATUS( status );
        
        m_strandVertices = strandVerticesArrayData.array( &status );
        CHECK_MSTATUS( status );
        
        m_numberOfVerticesPerStrand = i_dataBlock.inputValue( ia_verticesPerStrand ).asInt();        
        
        // Input has changed so rebuild rods.
       //initialiseRodFromBarberShopInput();
        
        i_dataBlock.setClean( i_plug );
    }
    else
    {
        return MStatus::kUnknownParameter;
    }

    return MStatus::kSuccess;
}

void WmSweeneyShape::initialiseRodFromBarberShopInput()
{
    // TODO: Grab the root strand distribution from a Barbershop node that gets plugged in
    // 
    
    //TODO: Reset the manager and remove all rods before adding more
    
    vector< BASim::Vec3d > vertices;
    vertices.resize( 10 );
    
    for ( unsigned int x = 0; x < 5; ++x )
    {
        for ( size_t v = 0; v < 10; ++v )
        {
            vertices[ v ] = BASim::Vec3d( x, 0.0, ( double )v );
        }
        m_rodManager.addRod( vertices );
    }
    
    m_rodManager.initialiseSimulation( 1 / 24.0, 1.0 );
}

void WmSweeneyShape::draw()
{
    //MObject thisNode( thisMObject() );    
    //bool dummy = false;

    //MPlug syncPlug( thisNode, ca_sync );
    //syncPlug.getValue( dummy );
    
    m_rodManager.drawAllRods();
}

MPxGeometryIterator* WmSweeneyShape::geometryIteratorSetup(MObjectArray& componentList,
													MObject& components,
													bool forReadOnly )
//
// Description
//
//    Creates a geometry iterator compatible with his shape.
//
// Arguments
//
//    componentList - list of components to be iterated
//    components    - component to be iterator
//    forReadOnly   -
//
// Returns
//
//    An iterator for the components
//
{
	WmSweeneyShapeIterator * result = NULL;
	if ( components.isNull() ) 
	{
		result = new WmSweeneyShapeIterator( getControlPoints(), componentList );
	}
	else 
	{
		result = new WmSweeneyShapeIterator( getControlPoints(), components );
	}
	return result;
}

bool WmSweeneyShape::acceptsGeometryIterator( bool writeable )
//
// Description
//
//    Specifies that this shape can provide an iterator for getting/setting
//    control point values.
//
// Arguments
//
//    writable - maya asks for an iterator that can set points if this is true
//
{
	return true;
}

bool WmSweeneyShape::acceptsGeometryIterator( MObject&, bool writeable, bool forReadOnly )
//
// Description
//
//    Specifies that this shape can provide an iterator for getting/setting
//    control point values.
//
// Arguments
//
//    writable   - maya asks for an iterator that can set points if this is true
//    forReadOnly - maya asking for an iterator for querying only
//
{
	return true;
}
