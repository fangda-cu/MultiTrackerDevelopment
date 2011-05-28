#include "WmSweeneyUtils.hh"

#include <maya/MDagPath.h>

namespace sweeney {
namespace utils {

MStatus findASelectedNodeByTypeName( MString& i_typeName, MObject* o_selectedNodeObject )
{
    // Only returns the first selected node of the type it finds.
    
    MStatus status;
    
    MSelectionList selectionList;
    MGlobal::getActiveSelectionList( selectionList );

    *o_selectedNodeObject = MObject::kNullObj;
    
    MDagPath nodeDagPath;
    MObject nodeObject;
    
    for ( unsigned int s = 0; s < selectionList.length(); ++s )
    {
        selectionList.getDagPath( s, nodeDagPath, nodeObject );
        
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

            if ( dependencyNodeFn.typeName() == i_typeName )
            {
                *o_selectedNodeObject = childPath.node( &status );
                CHECK_MSTATUS( status );
                
                return MStatus::kSuccess;
            }       
        }
        else // Perhaps no child as the user selected the shape node directly
        {
            MFnDependencyNode dependencyNodeFn( nodeDagPath.node( &status ) );
            CHECK_MSTATUS( status );

            if ( dependencyNodeFn.typeName() == i_typeName )
            {
                *o_selectedNodeObject = nodeDagPath.node( &status );
                CHECK_MSTATUS( status );
                
                return MStatus::kSuccess;
            }        
        }
    }
    
    return MStatus::kFailure;
}
    
MStatus findSelectedSweeneyNodeAndRodManager( WmSweeneyNode* o_wmSweeneyNode, 
    WmSweeneyRodManager* o_wmSweenyRodManager )
{
    MStatus status;
    
    MObject sweeneyNodeObj = MObject::kNullObj;

    findASelectedNodeByTypeName( WmSweeneyNode::typeName, &sweeneyNodeObj );
    
    if ( sweeneyNodeObj != MObject::kNullObj )
    {        
        MFnDependencyNode sweeneyNodeDepFn( sweeneyNodeObj, &status );
        CHECK_MSTATUS( status );
        
        o_wmSweeneyNode = dynamic_cast< WmSweeneyNode* >( sweeneyNodeDepFn.userNode() );    
        o_wmSweenyRodManager = o_wmSweeneyNode->rodManager();
    
        return MStatus::kSuccess;
    }
    
    MGlobal::displayError( "Please select a wmSweeneyNode." );
    
    return MStatus::kFailure;
}

}
}
