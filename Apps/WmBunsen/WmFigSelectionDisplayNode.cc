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
//#include <maya/glut.h>

const MTypeId WmFigSelectionDisplayNode::TypeId( 0x001156, 0x63 ); // Registered in /vol/weta/src/Wcommon/MTypeId.txt
const MString WmFigSelectionDisplayNode::TypeName( "wmFigSelectionDisplayNode" );

//MObject WmFigSelectionDisplayNode::displayParameters;
//MObject WmFigSelectionDisplayNode::displayParametersNonConstant;
//MObject WmFigSelectionDisplayNode::displayFlagVerticesForExactMatch;
//MObject WmFigSelectionDisplayNode::displayParametersDetails;
//MObject WmFigSelectionDisplayNode::displayTraining;
//MObject WmFigSelectionDisplayNode::displayWhichTrainingSample;
//MObject WmFigSelectionDisplayNode::displayVectors;
//MObject WmFigSelectionDisplayNode::displayFontSizeChoice;


WmFigSelectionDisplayNode::WmFigSelectionDisplayNode()
{
}

#if 0
bool WmFigSelectionDisplayNode::isRodVertexSelected( const unsigned int rodId, const unsigned int vertexId )
{
	unsigned int key = rodId;

	RodVertexIndices::iterator it;
	it = rodVertexIndices.find( key );
	if( it != rodVertexIndices.end() ) { // Found matching rod
		MIntArray &rodVertices = it->second;
		if( vertexId >= 0 && vertexId < rodVertices.length() ) {
			return bool( rodVertices[ vertexId ] );
		}
	}

	return false;
}

void WmFigSelectionDisplayNode::setRodVertexSelected( const unsigned int rodId, const unsigned int vertexId, const bool isSelected )
{
	unsigned int key = rodId;

	RodVertexIndices::iterator it;
	it = rodVertexIndices.find( key );
	if( it == rodVertexIndices.end() ) { // Didn't find matching rod
		std::pair<RodVertexIndices::iterator, bool> result;
		result = rodVertexIndices.insert( std::pair<unsigned int,VertexIndices>( key, VertexIndices() ) );
		it = result.first;
	}

	MIntArray &rodVertices = it->second;
	if( vertexId >= rodVertices.length() )  {
		const unsigned int prevLength = rodVertices.length();

		// Resize
		rodVertices.setLength( vertexId+1 );

		// Set all the newly added values to 0 (not selected)
		unsigned int i;
		for( i=prevLength; i < rodVertices.length(); i++ )
			rodVertices[i] = 0; // Not selected
	}

	rodVertices[ vertexId ] = int( isSelected );
}


void WmFigSelectionDisplayNode::clearSelection()
{
	unsigned int vi;
	RodVertexIndices::iterator it;
	for( it=rodVertexIndices.begin(); it != rodVertexIndices.end(); it++ ) {
		MIntArray &rodVertices = (*it).second;
		for( vi=0; vi < rodVertices.length(); vi++ )
			rodVertices[vi] = 0; // Not selected
	}
}

void WmFigSelectionDisplayNode::getSelectedRodIds( MIntArray &selectedRodIds )
{
	selectedRodIds.clear();

	unsigned int vi;
	RodVertexIndices::iterator it;
	for( it=rodVertexIndices.begin(); it != rodVertexIndices.end(); it++ ) {
		MIntArray &rodVertices = (*it).second;
		for( vi=0; vi < rodVertices.length(); vi++ ) {
			if( rodVertices[vi] ) {
				selectedRodIds.append( (*it).first );
				break;
			}
		}
	}
}

void WmFigSelectionDisplayNode::getSelectedRodVertexIds( const unsigned int rodId, MIntArray &selectedRodVertexIds )
{
	selectedRodVertexIds.clear();

	unsigned int key = rodId;
	RodVertexIndices::iterator it;
	it = rodVertexIndices.find( key );
	if( it != rodVertexIndices.end() ) { // Found matching rod
		MIntArray &rodVertices = (*it).second;
		unsigned int vi;
		for( vi=0; vi < rodVertices.length(); vi++ ) {
			if( rodVertices[vi] )
				selectedRodVertexIds.append( vi );
		}
	}
}
#endif

void WmFigSelectionDisplayNode::draw(
		M3dView &view, const MDagPath & path,
		M3dView::DisplayStyle displayStyle,
		M3dView::DisplayStatus displayStatus )
{
//	MFnDependencyNode thisNodeFn( thisMObject() );

	MObject figRodNode = getFigRodNodeFromSelection();
	if( figRodNode.isNull() )
		return;

	MFnDependencyNode nodeFn( figRodNode );
    WmFigRodNode *figRod = static_cast<WmFigRodNode*>( nodeFn.userNode() );
    if( !figRod )
        return;

	view.beginGL();
	glPushAttrib( GL_CURRENT_BIT );

//	glColor3f( 1.0, 0.0, 0.0 );
//	glBegin( GL_LINE_LOOP );
//	glVertex3f( 0.0, 0.0, 0.0 );
//	glVertex3f( 10.0, 0.0, 0.0 );
//	glEnd();
//	glBegin( GL_LINE_LOOP );
//	glVertex3f( 0.0, 0.0, 0.0 );
//	glVertex3f( 0.0, 10.0, 0.0 );
//	glEnd();
//	glBegin( GL_LINE_LOOP );
//	glVertex3f( 0.0, 0.0, 0.0 );
//	glVertex3f( 0.0, 0.0, 10.0 );
//	glEnd();

    WmFigRodGroup* rodGroup = figRod->rodGroup();
    //MGlobal::displayInfo( MString("# real rods: ") + rodGroup->numberOfRealRods() );
    if ( rodGroup->numberOfRealRods() == 0 )
        return;

//	MPoint nearClipPt, farClipPt, pt;
//	view.viewToWorld( x0, y0, nearClipPt, farClipPt );

    const size_t nRods = rodGroup->numberOfRods();
    //MGlobal::displayInfo( MString("# rods: ") + rodGroup->numberOfRealRods() );

#if 1
    MColor clr( 1.0, 1.0, 0.0 );
	glColor3f( clr.r, clr.g, clr.b );
#else

#define LEAD_COLOR                      18      // green
#define ACTIVE_COLOR                    15      // white

	// Draw vertices in the lead colour
	//
	view.setDrawColor( LEAD_COLOR, M3dView::kActiveColors );
#endif

	glPointSize( 6 );
	glBegin( GL_POINTS );

	for( size_t ri = 0; ri < nRods; ri++ )
    {
		const int nVerts = rodGroup->elasticRod( ri )->nv();
        for( size_t vi = 0; vi < nVerts; vi++ )
        {
//            if( isRodVertexSelected( ri, vi ) ) {
            if( selection.containsRodVertex( ri, vi ) ) {
            	const Vec3d p = rodGroup->elasticRod( ri )->getVertex( vi );
                glVertex3f( p[0], p[1], p[2] );
            }
            //MGlobal::displayInfo( MString( "rod: ") + ri + " vertex: " + vi + " " + p[0] + ", " + p[1] + ", " + p[2] );
            //glVertex3f( 1.0, 0.0, 0.0 );
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
#if 0
	MFnNumericAttribute nAttr;
    MFnEnumAttribute eAttr;

	displayParameters = nAttr.create( "displayParameters", "dbp", MFnNumericData::kBoolean, false );
	displayParametersNonConstant = nAttr.create( "displayParametersNonConstant", "dbn", MFnNumericData::kBoolean, true );
    displayFlagVerticesForExactMatch = nAttr.create( "displayFlagVerticesForExactMatch", "dbe", MFnNumericData::kBoolean, true );
    displayParametersDetails = nAttr.create( "displayParametersDetails", "dpr", MFnNumericData::kBoolean, false );

    displayTraining = nAttr.create( "displayTraining", "dtr", MFnNumericData::kBoolean, false );
    displayWhichTrainingSample = nAttr.create( "displayWhichTrainingSample", "dwts", MFnNumericData::kInt, 1 );
    nAttr.setMin(1);

    displayVectors = nAttr.create( "displayVectors", "dvc", MFnNumericData::kBoolean, false );

    displayFontSizeChoice = eAttr.create( "displayFontSizeChoice", "dfsc", 0 );
    eAttr.addField( "10 point", 0 );
    eAttr.addField( "12 point", 1 );
    eAttr.addField( "18 point", 2 );

	MFnTypedAttribute tAttr;
	inWPSDDataRef = tAttr.create( "inWPSDDataRef", "iwpsdRef", WmWeightedPSDDataRef::TypeId, MObject::kNullObj );
    tAttr.setReadable(true);
    tAttr.setWritable(true);
    tAttr.setConnectable(true);

    addAttribute( displayParameters );
    addAttribute( displayParametersNonConstant );
    addAttribute( displayFlagVerticesForExactMatch );
    addAttribute( displayParametersDetails );

    addAttribute( displayTraining );
	addAttribute( displayWhichTrainingSample );

    addAttribute( displayVectors );

	addAttribute( displayFontSizeChoice );

    addAttribute( inWPSDDataRef );

    //doDrawing = nAttr.create( "doDrawing", "ddr", MFnNumericData::kBoolean, true );

//    displayMinMeanRatios = nAttr.create( "displayMinMeanRatios", "dmmr", MFnNumericData::kBoolean, true );
//	minMeanRatioThreshold = nAttr.create( "minMeanRatioThreshold", "mmt", MFnNumericData::kDouble, 0.3 );
//	nAttr.setMin(0.0);
//	nAttr.setMax(1.0);
//	displayCurvatureConstraints = nAttr.create( "displayCurvatureConstraints", "dcc", MFnNumericData::kBoolean, true );
//
//	addAttribute( doDrawing );

//    addAttribute( displayMinMeanRatios );
//    addAttribute( minMeanRatioThreshold );
//    addAttribute( displayCurvatureConstraints );
#endif

	return MS::kSuccess;
}
