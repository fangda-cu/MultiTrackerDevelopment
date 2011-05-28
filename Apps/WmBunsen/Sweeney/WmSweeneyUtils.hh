#ifndef WMSWEENEYUTILS_HH_
#define WMSWEENEYUTILS_HH_

#include "WmSweeneyNode.hh"

// This is a generic place for utility functions that are needed across files and do not really
// live in any one class.

namespace sweeney {
namespace utils {

// This is inlined a as it's in the .hh and to avoid multiple definition issues. It should really 
// go live in a .cc
inline MStatus findSelectedSweeneyNodeAndRodManager( WmSweeneyNode* o_wmSweeneyNode, 
    WmSweeneyRodManager* o_wmSweenyRodManager )
{
    MStatus status;
    
    MSelectionList selectionList;
    MGlobal::getActiveSelectionList( selectionList );

    MObject sweeneyNodeObj = MObject::kNullObj;
    
    if ( selectionList.length() != 1 )
    {
        MGlobal::displayError( "Please select a single wmSweeney node" );
        return MStatus::kFailure;
    }

    MDagPath nodeDagPath;
    MObject nodeObject;
    selectionList.getDagPath( 0, nodeDagPath, nodeObject );
    
    // Look for a child as the user probably selected the transform
    MFnDagNode dagNodeFn( nodeDagPath.child( 0, &status ), &status );
    if ( status )
    {
        MDagPath childPath;
        
        status = dagNodeFn.getPath( childPath );
        CHECK_MSTATUS( status );
        childPath.extendToShape();

        MFnDependencyNode dependencyNodeFn( childPath.node( &status ) );
        CHECK_MSTATUS( status );

        if ( dependencyNodeFn.typeName() == WmSweeneyNode::typeName )
        {
            sweeneyNodeObj = childPath.node( &status );
            CHECK_MSTATUS( status );
        }       
    }
    else // Perhaps no child as the user selected the shape node directly
    {
        MFnDependencyNode dependencyNodeFn( nodeDagPath.node( &status ) );
        CHECK_MSTATUS( status );

        if ( dependencyNodeFn.typeName() == WmSweeneyNode::typeName )
        {
            sweeneyNodeObj = nodeDagPath.node( &status );
            CHECK_MSTATUS( status );
        }        
    }
    
    if ( sweeneyNodeObj != MObject::kNullObj )
    {        
        MFnDependencyNode sweeneyNodeDepFn( sweeneyNodeObj, &status );
        CHECK_MSTATUS( status );
        
        o_wmSweeneyNode = dynamic_cast< WmSweeneyNode* >( sweeneyNodeDepFn.userNode() );    
        o_wmSweenyRodManager = o_wmSweeneyNode->rodManager();
    
        return MStatus::kSuccess;
    }
    
    return MStatus::kFailure;
}


}
}

#endif
