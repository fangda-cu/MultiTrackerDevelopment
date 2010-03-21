#ifndef _WMFIGSELECTIONCONTEXT_HH_
#define _WMFIGSELECTIONCONTEXT_HH_

#include "WmFigSelectionToolCommand.hh"

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
#include <maya/MItCurveCV.h>
#include <maya/MPxContext.h>
#include <maya/MEventMessage.h>

class WmFigSelectionContext : public MPxContext 
{
public:
    enum MouseStatus { kPress, kDrag, kRelease };
    
    WmFigSelectionContext();
    
    //////////////////////////////////////////////////////
    // 
    // Inherited.
    // 
    //////////////////////////////////////////////////////

    // 
    // ..from MPxContext..
    // 

    virtual void        toolOnSetup( MEvent & event );
    virtual void        toolOffCleanup();
    virtual MStatus     doPress( MEvent & event );
    virtual MStatus     doDrag( MEvent & event );
    virtual MStatus     doRelease( MEvent & event );
    virtual MStatus     doEnterRegion( MEvent & event );
    virtual void        getClassName( MString& name ) const;
    
    static WmFigSelectionContext* activeContext;
    static MString typeName;
    
    //////////////////////////////////////////////////////
    // 
    // context functions.
    // 
    //////////////////////////////////////////////////////

private:
    MStatus             drawMarqueeSelectBox() const;
    void                doTool( MEvent &event );

    MString             m_helpString;
    
    short m_xStartMouse, m_yStartMouse;
    short m_xMouse, m_yMouse;
    
    //////////////////////////////////////////////////////
    // 
    // Values sent to the tool command, and the tool command
    // itself
    // 
    //////////////////////////////////////////////////////

    WmFigSelectionToolCommand* m_toolCommand;
};

#endif
