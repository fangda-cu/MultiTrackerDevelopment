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
    status = MGlobal::getActiveSelectionList( selectionList );
    CHECK_MSTATUS( status );
    
    MObject sweeneyNodeObj;
    status = selectionList.getDependNode( 0, sweeneyNodeObj );

    // If the user didn't select any nodes then fail
    if ( status != MStatus::kSuccess )
    {
        MGlobal::displayError( "Please select a single wmSweeney node" );
        return MStatus::kFailure;
    }
    
    MFnDependencyNode sweeneyNodeDepFn( sweeneyNodeObj );
    if ( sweeneyNodeDepFn.typeId() != WmSweeneyNode::typeID )
    {
        MGlobal::displayError( "Please select a single wmSweeney node" );
        return MStatus::kFailure;
    }
    
    // Get the actual node so we can talk to it directly
    o_wmSweeneyNode = dynamic_cast< WmSweeneyNode* >( sweeneyNodeDepFn.userNode() );
    
    o_wmSweenyRodManager = o_wmSweeneyNode->rodManager();
    
    return MStatus::kSuccess;
}

}
}

#endif
