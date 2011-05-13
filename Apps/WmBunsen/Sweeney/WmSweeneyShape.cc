#include <math.h>           
#include <maya/MIOStream.h>

#include "WmSweeneyShape.hh"         
#include "WmSweeneyShapeUI.hh"     
#include "WmSweeneyShapeIterator.hh"

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
    // FIXME: This should happen in the compute
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
    
    
    return MS::kSuccess;
}

void WmSweeneyShape::initialiseRodFromBarberShopInput()
{
    // TODO: Grab the root strand distribution from a Barbershop node that gets plugged in
    // 
    
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
