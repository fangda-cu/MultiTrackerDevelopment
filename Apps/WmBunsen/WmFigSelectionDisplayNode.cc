/*
 * WmFigSelectionDisplayNode.cc
 *
 *  Created on: Oct 18, 2010
 *      Author: dgould
 */

#include "WmFigSelectionDisplayNode.hh"
#include "WmFigUtils.hh"
#include "WmFigRodNode.hh"
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnPluginData.h>
#include <maya/MGlobal.h>
#include <maya/MPlugArray.h>
#include <maya/MFnMesh.h>
#include <maya/MSelectionList.h>
#include <maya/MItMeshVertex.h>
#include <iterator>
#include <float.h>
#include <GL/glut.h>
#include <limits.h>

const MTypeId WmFigSelectionDisplayNode::TypeId( 0x001156, 0x63 ); // Registered in /vol/weta/src/Wcommon/MTypeId.txt
const MString WmFigSelectionDisplayNode::TypeName( "wmFigSelectionDisplayNode" );

WmFigSelectionDisplayNode::WmFigSelectionDisplayNode()
{
}

void WmFigSelectionDisplayNode::draw(
		M3dView &view, const MDagPath & path,
		M3dView::DisplayStyle displayStyle,
		M3dView::DisplayStatus displayStatus )
{
	MObject figRodNode = getFigRodNodeFromSelection();
	if( figRodNode.isNull() )
		return;

	MFnDependencyNode nodeFn( figRodNode );
    WmFigRodNode *figRod = static_cast<WmFigRodNode*>( nodeFn.userNode() );
    if( !figRod )
        return;

	view.beginGL();
	glPushAttrib( GL_CURRENT_BIT );

    WmFigRodGroup* rodGroup = figRod->rodGroup();
    //MGlobal::displayInfo( MString("# real rods: ") + rodGroup->numberOfRealRods() );
    if ( rodGroup->numberOfRealRods() == 0 )
        return;

    const int nRods = rodGroup->numberOfRods();
    //MGlobal::displayInfo( MString("# rods: ") + rodGroup->numberOfRealRods() );

    MColor clr( 1.0, 1.0, 0.0 );
	glColor3f( clr.r, clr.g, clr.b );

	glPointSize( 6 );
	glBegin( GL_POINTS );

	for( int ri = 0; ri < nRods; ri++ )
    {
		const int nVerts = rodGroup->elasticRod( ri )->nv();
        for( int vi = 0; vi < nVerts; vi++ )
        {
            if( selection.containsRodVertex( ri, vi ) ) {
            	const Vec3d p = rodGroup->elasticRod( ri )->getVertex( vi );
                glVertex3f( p[0], p[1], p[2] );
            }
            //MGlobal::displayInfo( MString( "rod: ") + ri + " vertex: " + vi + " " + p[0] + ", " + p[1] + ", " + p[2] );
        }
    }
    glEnd();

	glPopAttrib();
	view.endGL();
}

bool WmFigSelectionDisplayNode::isBounded() const
{
	return false;
}

void *WmFigSelectionDisplayNode::creator()
{
	return new WmFigSelectionDisplayNode();
}

MStatus WmFigSelectionDisplayNode::initialize()
{
	return MS::kSuccess;
}
