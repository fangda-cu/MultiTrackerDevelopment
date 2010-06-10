#include "WmFigControllerToolCommand.hh"

using namespace std;

// NOTE: Keep this type name lowercase, since it becomes a MEL command name.
// Yes it's annoying.
MString WmFigControllerToolCommand::typeName( "wmFigControllerToolCommand" );

WmFigControllerToolCommand::WmFigControllerToolCommand()
{
    setCommandString( typeName );
}

WmFigControllerToolCommand::~WmFigControllerToolCommand() 
{
}

void* WmFigControllerToolCommand::creator() 
{
    return new WmFigControllerToolCommand;
}

bool WmFigControllerToolCommand::isUndoable() const 
{
    return true;
}

MStatus WmFigControllerToolCommand::doIt( const MArgList& i_args ) 
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

MStatus WmFigControllerToolCommand::redoIt()
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

    //////////////////////////////////////////////////////
    // 
    // parse the args we got from the first call to doIt().
    // also, set up other values.
    // 
    //////////////////////////////////////////////////////

    
    //////////////////////////////////////////////////////
    // 
    // Handle the click and make the command do the work
    // 
    //////////////////////////////////////////////////////

    if ( m_selectedRods.size() > 0 )
    {
        /*MGlobal::executeCommand( MString( "$rodNode = `ls -sl`; \
                                    $locator = `spaceLocator -p 0 0 0`; \
                                    select -add $rodNode; " ) + 
                                    "wmFigaro -aeo -ro " + m_selectedRods[ 0 ] + " -ed 0;" );*/

        MSelectionList selectionList;

        MGlobal::getActiveSelectionList( selectionList );
        cerr << "selection list has " << selectionList.length() << " entries\n";
        
        
    }
    
    //////////////////////////////////////////////////////
    // 
    // we are now officially in the edit run.
    // 
    //////////////////////////////////////////////////////
    
    editRunState() = kEditRunIsGoing;    

    return MS::kSuccess;
}

MStatus WmFigControllerToolCommand::undoIt() 
{
    // We change nothing so we have nothing to undo
    
    return MS::kSuccess;
}


MStatus WmFigControllerToolCommand::updateContextOptions()
{
    MStatus stat;
    MString currentCtx;
    MGlobal::executeCommand("currentCtx", currentCtx);

    // We currently have no options to update

    return stat;
}


MStatus WmFigControllerToolCommand::finalize()
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

