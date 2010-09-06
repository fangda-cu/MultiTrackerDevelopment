//-
// ==========================================================================
// Copyright 1995,2006,2008 Autodesk, Inc. All rights reserved.
//
// Use of this software is subject to the terms of the Autodesk
// license agreement provided at the time of installation or download,
// or which otherwise accompanies this software in either electronic
// or hard copy form.
// ==========================================================================
//+

///////////////////////////////////////////////////////////////////////////////
//
// apiMeshIterator.cpp
//
///////////////////////////////////////////////////////////////////////////////

#include <maya/MVectorArray.h>

#include "WmFigaroRodShapeIterator.hh"

WmFigaroRodShapeIterator::WmFigaroRodShapeIterator( void * geom, MObjectArray & comps )
	: MPxGeometryIterator( geom, comps ),
	geometry( (MVectorArray*)geom )
{
	reset();
}

WmFigaroRodShapeIterator::WmFigaroRodShapeIterator( void * geom, MObject & comps )
	: MPxGeometryIterator( geom, comps ),
	geometry( (MVectorArray*)geom )
{
	reset();
}

/* override */
void WmFigaroRodShapeIterator::reset()
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
MPoint WmFigaroRodShapeIterator::point() const
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
void WmFigaroRodShapeIterator::setPoint( const MPoint & pnt ) const
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

    cerr << "Moving point\n";
}

/* override */
int WmFigaroRodShapeIterator::iteratorCount() const
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
bool WmFigaroRodShapeIterator::hasPoints() const
//
// Description
//
//    Returns true since the shape data has points.
//
{
	return true;
}
