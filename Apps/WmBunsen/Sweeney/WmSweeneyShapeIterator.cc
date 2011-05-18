///////////////////////////////////////////////////////////////////////////////
//
// apiMeshIterator.cpp
//
///////////////////////////////////////////////////////////////////////////////

#include <maya/MVectorArray.h>

#include "WmSweeneyShapeIterator.hh"

WmSweeneyShapeIterator::WmSweeneyShapeIterator( void * geom, MObjectArray & comps )
	: MPxGeometryIterator( geom, comps ),
	geometry( (MVectorArray*)geom )
{
	reset();
}

WmSweeneyShapeIterator::WmSweeneyShapeIterator( void * geom, MObject & comps )
	: MPxGeometryIterator( geom, comps ),
	geometry( (MVectorArray*)geom )
{
	reset();
}

/* override */
void WmSweeneyShapeIterator::reset()
//
// Description
//
//  	
//   Resets the iterator to the start of the components so that another
//   pass over them may be made.
//
{
	MPxGeometryIterator::reset();
	setCurrentPoint( 0 );
	if ( NULL != geometry ) 
	{
		int maxVertex = geometry->length();
		setMaxPoints( maxVertex );
	}
}

/* override */
MPoint WmSweeneyShapeIterator::point() const
//
// Description
//
//    Returns the point for the current element in the iteration.
//    This is used by the transform tools for positioning the
//    manipulator in component mode. It is also used by deformers.	 
//
{
	MPoint pnt;
	if ( NULL != geometry ) 
	{
		pnt = (*geometry)[index()];
	}
	return pnt;
}

/* override */
void WmSweeneyShapeIterator::setPoint( const MPoint & pnt ) const
//
// Description
//
//    Set the point for the current element in the iteration.
//    This is used by deformers.	 
//
{
	if ( NULL != geometry ) 
	{
		(*geometry)[index()] = pnt;
	}
}

/* override */
int WmSweeneyShapeIterator::iteratorCount() const
{
//
// Description
//
//    Return the number of vertices in the iteration.
//    This is used by deformers such as smooth skinning
//
	return geometry->length();
	
}

/* override */
bool WmSweeneyShapeIterator::hasPoints() const
//
// Description
//
//    Returns true since the shape data has points.
//
{
	return true;
}
