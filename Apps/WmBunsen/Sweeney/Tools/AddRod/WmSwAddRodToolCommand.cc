#include "WmSwAddRodToolCommand.hh"

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


MStatus WmSwAddRodToolCommand::finalize()
{
    MStatus status;

    //////////////////////////////////////////////////////
    // 
    // set this, so in redoIt(), we know to restore the
    // cvs of all the nurbs to nurbsVerticesNew.
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
