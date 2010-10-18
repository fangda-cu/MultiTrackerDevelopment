#include "WmFigSelectionToolCommand.hh"
#include "WmFigSelectionDisplayNode.hh"
#include "WmFigUtils.hh"
#include <maya/MFnDependencyNode.h>

using namespace std;

// NOTE: Keep this type name lowercase, since it becomes a MEL command name.
// Yes it's annoying.
MString WmFigSelectionToolCommand::typeName( "wmFigSelectionToolCommand" );

WmFigSelectionToolCommand::WmFigSelectionToolCommand()
{
    setCommandString( typeName );
    selectionDisplayNode = MObject::kNullObj;
}

WmFigSelectionToolCommand::~WmFigSelectionToolCommand() 
{
}

void* WmFigSelectionToolCommand::creator() 
{
    return new WmFigSelectionToolCommand;
}

bool WmFigSelectionToolCommand::isUndoable() const 
{
    return true;
}

MStatus WmFigSelectionToolCommand::doIt( const MArgList& i_args ) 
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

MStatus WmFigSelectionToolCommand::redoIt()
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

 /*   if ( m_args.length() != (uint) WmFigSelectionToolCommand::expectedArgCount ) 
    {
        MGlobal::displayError( "WmFigSelectionToolCommand::redoIt() recieved "
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

MStatus WmFigSelectionToolCommand::undoIt() 
{
    // We change nothing so we have nothing to undo
    
    return MS::kSuccess;
}


MStatus WmFigSelectionToolCommand::updateContextOptions()
{
    MStatus stat;
    MString currentCtx;
    MGlobal::executeCommand("currentCtx", currentCtx);

    // We currently have no options to update

    return stat;
}


MStatus WmFigSelectionToolCommand::finalize()
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
    
    MString result;

    if( !selectionDisplayNode.isNull() ) {
		MFnDependencyNode nodeFn( selectionDisplayNode );
		WmFigSelectionDisplayNode &displayNode = *static_cast<WmFigSelectionDisplayNode*>( nodeFn.userNode() );

//		MObject figRodNode = getFigRodNodeFromSelection();
//		nodeFn.setObject( figRodNode );

		displayNode.selection.serialise( result );
    }

#if 0
    MStringArray results;

    if( !selectionDisplayNode.isNull() ) {
		MFnDependencyNode nodeFn( selectionDisplayNode );
		WmFigSelectionDisplayNode &displayNode = *static_cast<WmFigSelectionDisplayNode*>( nodeFn.userNode() );

		MObject figRodNode = getFigRodNodeFromSelection();
		nodeFn.setObject( figRodNode );

	    MIntArray selectedRodIds;
	    MIntArray selectedRodVertexIds;
	    displayNode.getSelectedRodIds( selectedRodIds );

	    //MGlobal::displayInfo( MString("# selected rods: ") + selectedRodIds.length() );
	    unsigned int iVertex;
	    unsigned int iRod;
	    unsigned int rodId, vertexId;
	    for( iRod=0; iRod < selectedRodIds.length(); iRod++ ) {
	    	rodId = selectedRodIds[iRod];

	    	MString selItemTxt( nodeFn.name() );
	    	selItemTxt += MString(".rod[") + rodId + "]";

	    	// @@@ combine contiguous indices like Maya does

	    	selItemTxt += MString(".vtx[");
	    	displayNode.getSelectedRodVertexIds( rodId, selectedRodVertexIds );
	    	for( iVertex=0; iVertex < selectedRodVertexIds.length(); iVertex++ ) {
	    		vertexId = selectedRodVertexIds[ iVertex ];
	    		selItemTxt += vertexId;
	    		if( iVertex < selectedRodVertexIds.length()-1 )
	    			selItemTxt += " ";
	    	}
	    	selItemTxt += MString("]");

	    	results.append( selItemTxt );
	    }
    }
#endif

#if 0
    MStringArray results;
    
    if( m_selected.size() )
    {
    	results.setLength( (unsigned int)m_selected.size() );

    	WmFigSelections::iterator it=m_selected.begin();
    	unsigned int i;
    	for( i=0; it!=m_selected.end(); ++it, i++ )
    	{
    		WmFigSelectedItem &selItem = *it;
    		MFnDependencyNode nodeFn( selItem.figRodNode );
    		MString selItemTxt( nodeFn.name() );
    		if( selItem.rodId != -1 )
    			selItemTxt += MString(".rod[") + selItem.rodId + "]";
    		if( selItem.rodVertexId != -1 )
    			selItemTxt += MString(".vtx[") + selItem.rodVertexId + "]";
    	    if( selItem.rodU != DBL_MAX )
    	    	selItemTxt += MString(".u[") + selItem.rodU + "]";

    	    results[i] = selItemTxt;
    	}
    }
    else
    {
        results.setLength( (unsigned int) m_selectedRods.size() );

        for ( unsigned int r=0; r < m_selectedRods.size(); ++r )
        {
         //   appendToResult( results[ r ] );
            //results[r] = MString( (unsigned int)m_selectedRods[r] );
            results[r] = MString() + (int)m_selectedRods[r];
        }
    }
#endif
    
    setResult( result );

    //MString result;
    //setResult( "yeh baby" );

    return MPxToolCommand::finalize();
}

