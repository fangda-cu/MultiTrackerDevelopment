/*
 * WmFigConstraintNode.cc
 *
 *  Created on: Oct 12, 2010
 *      Author: dgould
 */
#include "WmFigConstraintNode.hh"
#include "WmFigRodComponentList.hh"
#include "WmFigRodNode.hh"
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnMatrixAttribute.h>
#include <maya/MFnMessageAttribute.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MPlugArray.h>

const MTypeId WmFigConstraintNode::TypeId( 0x001156, 0x60 ); // Registered in /vol/weta/src/Wcommon/MTypeId.txt
const MString WmFigConstraintNode::TypeName( "wmFigConstraintNode" );

// Basic attributes
//
MObject WmFigConstraintNode::enable;
MObject WmFigConstraintNode::constraintType;
MObject WmFigConstraintNode::stiffness;
MObject WmFigConstraintNode::targetWorldPosition;
MObject WmFigConstraintNode::rodVertices;
MObject WmFigConstraintNode::figRodNodeMsg;

WmFigConstraintNode::WmFigConstraintNode()
{
}

WmFigConstraintNode::~WmFigConstraintNode()
{
}

void WmFigConstraintNode::draw( M3dView &view, const MDagPath &path, M3dView::DisplayStyle style, M3dView::DisplayStatus status )
{
	MFnDependencyNode figConstraintFn( thisMObject() );

	MObject figRodNode;
	MPlugArray srcPlugs;
	MPlug figRodNodeMsgPlug = figConstraintFn.findPlug( "figRodNodeMsg" );
	figRodNodeMsgPlug.connectedTo( srcPlugs, true, false );
	if( srcPlugs.length() )
	{
		figRodNode = srcPlugs[0].node();
	}

	MFnDependencyNode figRodNodeFn( figRodNode );
    WmFigRodNode *figRod = static_cast<WmFigRodNode*>( figRodNodeFn.userNode() );
    if( !figRod )
        return;

    MPlug rodVerticesPlug = figConstraintFn.findPlug( "rodVertices" );
    MString rodVerticesTxt( rodVerticesPlug.asString() );
    WmFigRodComponentList rodComponentList;
    rodComponentList.unserialise( rodVerticesTxt );

    MPlug targetWorldPositionPlug = figConstraintFn.findPlug( "targetWorldPosition" );
    MObject targetWorldPositionData = targetWorldPositionPlug.asMObject();
    MFnNumericData numericDataFn( targetWorldPositionData );
    MPoint targetWorldPosition;
    numericDataFn.getData3Double( targetWorldPosition.x, targetWorldPosition.y, targetWorldPosition.z );

    WmFigRodGroup* rodGroup = figRod->rodGroup();
    //MGlobal::displayInfo( MString("# real rods: ") + rodGroup->numberOfRealRods() );
    if ( rodGroup->numberOfRealRods() == 0 )
        return;

    view.beginGL();
	glPushAttrib( GL_CURRENT_BIT );

    GLboolean lineStipple[ 1 ];
    glGetBooleanv( GL_LINE_STIPPLE, lineStipple );

    glLineStipple(2, 0x0C0F);
    glEnable( GL_LINE_STIPPLE );

	MColor clr( 0.85, 0.4, 0.0 );
	glColor3f( clr.r, clr.g, clr.b );

	glPointSize( 6 );
	glBegin( GL_POINTS );
	glVertex3f( targetWorldPosition.x, targetWorldPosition.y, targetWorldPosition.z );
	glEnd();

	MIntArray rodIds;
	MIntArray rodVertexIds;
	unsigned int iVertex;
	unsigned int iRod;
	unsigned int rodId, vertexId;

	rodComponentList.getRodIds( rodIds );
	for( iRod=0; iRod < rodIds.length(); iRod++ ) {
		rodId = rodIds[iRod];

		rodComponentList.getRodVertexIds( rodId, rodVertexIds );
		for( iVertex=0; iVertex < rodVertexIds.length(); iVertex++ ) {
			vertexId = rodVertexIds[ iVertex ];

			const BASim::Vec3d p = rodGroup->elasticRod( rodId )->getVertex( vertexId );

			glBegin( GL_POINTS );
			glVertex3f( p[0], p[1], p[2] );
			glEnd();

			// Draw line from rod vertex to world space position
			glBegin( GL_LINE_LOOP );
			glVertex3f( p[0], p[1], p[2] );
			glVertex3f( targetWorldPosition.x, targetWorldPosition.y, targetWorldPosition.z );
			glEnd();
		}
	}

    if( lineStipple[0] )
    	glEnable( GL_LINE_STIPPLE );
    else
    	glDisable( GL_LINE_STIPPLE );

	glPopAttrib();
	view.endGL();
}

bool WmFigConstraintNode::isBounded() const
{
	return false;
}


MStatus WmFigConstraintNode::compute( const MPlug &plug, MDataBlock &data )
{
	MStatus stat;
	return MS::kUnknownParameter; // Haven't computed the plug, pass it back up the class hierarchy
}

void *WmFigConstraintNode::creator()
{
	return new WmFigConstraintNode();
}

MStatus WmFigConstraintNode::initialize()
{
	MStatus stat;

	// Create and initialize custom attributes
	//
	MFnEnumAttribute eAttr;
    MFnMessageAttribute msgAttr;
    MFnNumericAttribute nAttr;
    MFnTypedAttribute tAttr;

	enable = nAttr.create( "enable", "e", MFnNumericData::kBoolean, true );

	constraintType = eAttr.create( "constraintType", "ct", 1 );
    eAttr.addField( "Fixed", 0 );
    eAttr.addField( "Rest", 1 );

	stiffness = nAttr.create( "stiffness", "stf", MFnNumericData::kDouble, 50.0 );
	nAttr.setKeyable( true );

	targetWorldPosition = nAttr.create( "targetWorldPosition", "twp", MFnNumericData::k3Double );
	nAttr.setDefault( 1.0, 0.0, 0.0 );
	nAttr.setKeyable( true );

	rodVertices = tAttr.create( "rodVertices", "rvt", MFnData::kString );

	figRodNodeMsg = msgAttr.create( "figRodNodeMsg", "frm" );
    msgAttr.setHidden( true );

    addAttribute( enable );
    addAttribute( constraintType );
    addAttribute( stiffness );
    addAttribute( targetWorldPosition );
    addAttribute( rodVertices );
    addAttribute( figRodNodeMsg );

	return MS::kSuccess;
}
