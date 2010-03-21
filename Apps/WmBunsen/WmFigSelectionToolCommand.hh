#ifndef _WMFIGSELECTIONOOLCOMMAND_HH_
#define _WMFIGSELECTIONOOLCOMMAND_HH_

#include <maya/MIOStream.h>
#include <math.h>
#include <stdlib.h>
#include <maya/MString.h>
#include <maya/MGlobal.h>
#include <maya/M3dView.h>
#include <maya/MDagPath.h>
#include <maya/MPxNode.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MItSelectionList.h>
#include <maya/MSelectionList.h>
#include <maya/MPxContextCommand.h>
#include <maya/MPxContext.h>
#include <maya/MEvent.h>
#include <maya/MPointArray.h>
#include <maya/MFnMesh.h>
#include <maya/MQuaternion.h>
#include <maya/MArgList.h>
#include <maya/MPxToolCommand.h>

class WmFigSelectionToolCommand : public MPxToolCommand
{
public:

    WmFigSelectionToolCommand();
    virtual ~WmFigSelectionToolCommand();

    //////////////////////////////////////////////////////
    // 
    // some constants.
    // 
    //////////////////////////////////////////////////////

    static const int expectedArgCount = 2;

    //////////////////////////////////////////////////////
    // 
    // mainly inherited maya stuff.
    // 
    //////////////////////////////////////////////////////

    virtual MStatus doIt( const MArgList& i_args );
    virtual MStatus redoIt();
    virtual MStatus undoIt();
    virtual bool isUndoable() const;
    virtual MStatus finalize();

    static void* creator();
	
    static MString typeName;
    
    //////////////////////////////////////////////////////
    // 
    // This helps tool commands to figure out what state they
    // are in while editting.
    // 
    //////////////////////////////////////////////////////

    enum EditRunState
    {
        kEditRunHasJustStarted,
        kEditRunIsGoing,
        kEditRunIsComplete
    };

    inline EditRunState& editRunState()
    {
        return m_editRunState;
    }

private:
    EditRunState m_editRunState;

    //////////////////////////////////////////////////////
    // 
    // these args are only used for passing the argument list
    // from doIt() to redoIt().
    // 
    //////////////////////////////////////////////////////

    MArgList m_args;

    MStatus updateContextOptions();
};

#endif
