/*
 * WmFigUtils.cc
 *
 *  Created on: Oct 18, 2010
 *      Author: dgould
 */
#include "WmFigUtils.hh"
#include "WmFigRodNode.hh"
#include <maya/MGlobal.h>
#include <maya/MSelectionList.h>
#include <maya/MDagPath.h>
#include <maya/MStatus.h>
#include <maya/MFnDagNode.h>

MObject getFigRodNodeFromSelection()
{
    MSelectionList selectionList;
    MGlobal::getActiveSelectionList( selectionList );

    MDagPath rodNodeDagPath;
    MObject rodNodeObj = MObject::kNullObj;
    bool foundRodNode = true;
    MStatus stat;

    if ( selectionList.length() == 1 )
    {
        selectionList.getDagPath( 0, rodNodeDagPath, rodNodeObj );

        MDagPath childPath;

        // Look for a child as the user probably selected the transform
        MFnDagNode rodDagNodeFn( rodNodeDagPath.child( 0, &stat ), &stat );
        if ( stat )
        {
            stat = rodDagNodeFn.getPath( childPath );
            CHECK_MSTATUS( stat );
            childPath.extendToShape();

            MFnDependencyNode nodeFn( childPath.node( &stat ) );
            CHECK_MSTATUS( stat );

            if ( nodeFn.typeName() == WmFigRodNode::typeName )
            {
                rodNodeObj = childPath.node();
            }
            else
                foundRodNode = false;
        }
        else // Perhaps no child as the user selected the shape node directly
        {
            MFnDependencyNode nodeFn( rodNodeDagPath.node( &stat ) );
            CHECK_MSTATUS( stat );

            if ( nodeFn.typeName() == WmFigRodNode::typeName )
            {
                //    rodNodeObj = nodeFn.node();
                rodNodeObj = rodNodeDagPath.node( &stat );
            }
            else
                foundRodNode = false;
        }
    }
    else
        foundRodNode = false;

    if ( !foundRodNode )
    	return MObject::kNullObj;
    else
    	return rodNodeObj;
}
