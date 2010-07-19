#include "WmFigCombToolCommand.hh"

using namespace std;

// NOTE: Keep this type name lowercase, since it becomes a MEL command name.
// Yes it's annoying.
MString WmFigCombToolCommand::typeName( "wmFigCombToolCommand" );

WmFigCombToolCommand::WmFigCombToolCommand()
{
    setCommandString( typeName );
}

WmFigCombToolCommand::~WmFigCombToolCommand() 
{
}

void* WmFigCombToolCommand::creator() 
{
    return new WmFigCombToolCommand;
}

bool WmFigCombToolCommand::isUndoable() const 
{
    return true;
}

MStatus WmFigCombToolCommand::doIt( const MArgList& i_args ) 
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

MStatus WmFigCombToolCommand::redoIt()
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

    MStatus stat;
    if ( editRunState() == kEditRunIsComplete )
    {
        // As we don't change anything we have nothing to undo
        return MS::kSuccess;
    }

    // We actually do nothing in here because we need mouseUp to happen
    // before we can figure out all the rods that were selected.
    // Unfortunately we only get the args sent to us on mouseDown.
    // So all the work is done in the context and it simply tells us
    // what it finds and we return it to Maya as a result.
    
    //////////////////////////////////////////////////////
    // 
    // parse the args we got from the first call to doIt().
    // also, set up other values.
    // 
    //////////////////////////////////////////////////////

 /*   if ( m_args.length() != (uint) WmFigCombToolCommand::expectedArgCount ) 
    {
        MGlobal::displayError( "WmFigCombToolCommand::redoIt() recieved "
                               "wrong number of args." );
        return MS::kFailure;
    }

    int xMouse = m_args.asInt( 0 );
    int yMouse = m_args.asInt( 1 );*/
    
    //////////////////////////////////////////////////////
    // 
    // Handle the click and make the command do the work
    // 
    //////////////////////////////////////////////////////
    
    //MPoint clickCentre;
    //MVector clickNormal;
    
    //////////////////////////////////////////////////////
    // 
    // we are now officially in the edit run.
    // 
    //////////////////////////////////////////////////////
    
    editRunState() = kEditRunIsGoing;    

    return MS::kSuccess;
}

MStatus WmFigCombToolCommand::undoIt() 
{
    // We change nothing so we have nothing to undo
    
    return MS::kSuccess;
}


MStatus WmFigCombToolCommand::updateContextOptions()
{
    MStatus stat;
    MString currentCtx;
    MGlobal::executeCommand("currentCtx", currentCtx);

    // We currently have no options to update

    return stat;
}


MStatus WmFigCombToolCommand::finalize()
{
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
    
    //////////////////////////////////////////////////////
    // This is slightly weird because we run the base class 
    // finalize after this which means that the user will 
    // see the result before the command. But that's just
    // the interesting way Maya makes us do things...
    //////////////////////////////////////////////////////
    
    MStringArray results;
    results.setLength( m_selectedRods.size() );
    
    for ( size_t r=0; r<m_selectedRods.size(); r++ )
    {
     //   appendToResult( results[ r ] );
        results[ r ] = m_selectedRods[ r ];
    }
    
    setResult( results );
    
    return MPxToolCommand::finalize();
}

