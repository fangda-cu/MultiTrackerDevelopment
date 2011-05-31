#include "WmSwAddRodToolCommand.hh"
#include "../../WmSweeneyUtils.hh"

#include <maya/MFnIntArrayData.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnNurbsCurveData.h>

using namespace std;

namespace sweeney {

// NOTE: Keep this type name lowercase, since it becomes a MEL command name.
// Yes it's annoying.
MString WmSwAddRodToolCommand::typeName( "wmSwAddRodToolCommand" );

WmSwAddRodToolCommand::WmSwAddRodToolCommand() : m_sweeneyNode( NULL ), 
    m_rodManager( NULL )
{
    setCommandString( typeName );
    
    utils::findSelectedSweeneyNodeAndRodManager( m_sweeneyNode, m_rodManager );
    
    MString meshNodeTypename( "mesh" );
    MObject meshNodeObject;
    utils::findASelectedNodeByTypeName( meshNodeTypename, &meshNodeObject, &m_meshNodeDagPath );
    
    if ( meshNodeObject == MObject::kNullObj )
    {
        MGlobal::displayError( "Please select a wmSweeneyNode and a poly mesh node to create a rod." );
    }
}

WmSwAddRodToolCommand::~WmSwAddRodToolCommand() 
{
}

void* WmSwAddRodToolCommand::creator() 
{
    return new WmSwAddRodToolCommand;
}

bool WmSwAddRodToolCommand::isUndoable() const 
{
    return true;
}

MStatus WmSwAddRodToolCommand::doIt( const MArgList& i_args ) 
{
    m_args = i_args;

    if ( editRunState() == kEditRunHasJustStarted )
    {
        if ( updateContextOptions() == MS::kFailure )
        {
            return MS::kFailure;
        }

    }

    return redoIt();
}

MStatus WmSwAddRodToolCommand::redoIt()
{
    //////////////////////////////////////////////////////
    // 
    // if we see that the edit run state is 'complete', we
    // know they have entered this function as the result
    // of doing a re-do in the interface. 
    // This is not an undo, but it is a redo after
    // an undo.  It seems weird this would be set this way
    // at this point, but it makes sense.. the edit run was
    // completed.
    // 
    //////////////////////////////////////////////////////

    MStatus status;
    if ( editRunState() == kEditRunIsComplete )
    {
        // As we don't change anything we have nothing to undo
        return MS::kSuccess;
    }
    
    //////////////////////////////////////////////////////
    // 
    // parse the args we got from the first call to doIt().
    // also, set up other values.
    // 
    //////////////////////////////////////////////////////

    if ( m_args.length() != (uint) WmSwAddRodToolCommand::expectedArgCount ) 
    {
        MGlobal::displayError( "WmSwAddRodToolCommand::redoIt() recieved "
                               "wrong number of args." );
        return MS::kFailure;
    }

    int xMouse = m_args.asInt( 0 );
    int yMouse = m_args.asInt( 1 );

    //////////////////////////////////////////////////////
    // 
    // Handle the click and make the command do the work
    // 
    //////////////////////////////////////////////////////
    

    //////////////////////////////////////////////////////
    // 
    // we are now officially in the edit run.
    // 
    //////////////////////////////////////////////////////
    
    editRunState() = kEditRunIsGoing;    

    return MS::kSuccess;
}

MStatus WmSwAddRodToolCommand::undoIt() 
{
    // We change nothing so we have nothing to undo
    
    return MS::kSuccess;
}


MStatus WmSwAddRodToolCommand::updateContextOptions()
{
    MStatus stat;
    MString currentCtx;
    MGlobal::executeCommand("currentCtx", currentCtx);

    // We currently have no options to update

    return stat;
}

void WmSwAddRodToolCommand::createRod()
{
    MStatus status;
    MFnMesh growthMesh( m_meshNodeDagPath, &status );
    CHECK_MSTATUS( status );
    
    // Use acceleration structures to work out where on the mesh was hit
     
    // Add the rod to the rod manager and store a reference to the rod so we can 
    // tell the rod manager to delete the rod later if we are asked to undo this
    // addition.
}

MStatus WmSwAddRodToolCommand::finalize()
{
    MStatus status;

    // Finally, create the rod!
    createRod();

    //////////////////////////////////////////////////////
    // 
    // set this, so in redoIt(), we know to do the undo
    // 
    //////////////////////////////////////////////////////

    editRunState() = kEditRunIsComplete;
    

    //////////////////////////////////////////////////////
    //
    // this is a dumb entry for the journal, but there isn't
    // a simple way to boil down all the effects of what
    // the user has been doing interactively up until the
    // mouse was released..
    //
    //////////////////////////////////////////////////////
    MArgList command;
    command.addArg( commandString() );

    return MPxToolCommand::finalize();
}

}
