#ifndef WMSWADDRODCONTEXT_HH_
#define WMSWADDRODCONTEXT_HH_

#include "WmSwAddRodToolCommand.hh"
#include "../../WmSweeneyNode.hh"
#include "../../WmSweeneyUtils.hh"

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
#include <maya/MPoint.h>
#include <maya/MPxContext.h>

namespace sweeney {

class WmSwAddRodContext : public MPxContext 
{
public:
    enum MouseStatus { kPress, kDrag, kRelease };
    
    WmSwAddRodContext();
    
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
    
    static WmSwAddRodContext* activeContext;
    static MString typeName;
    
    //////////////////////////////////////////////////////
    // 
    // context functions.
    // 
    //////////////////////////////////////////////////////
    
    void drawManipulator();
    void drawManipulatorGL();
        
private:
    MStatus findSelectedShrubbinNodeAndTree();
    void doTool( MEvent &event );
    void updateBranchDirections();
    
    MString m_helpString;
    
    short m_xStartMouse, m_yStartMouse;
    short m_xMouse, m_yMouse;
    
    //////////////////////////////////////////////////////
    // 
    // Values sent to the tool command, and the tool command
    // itself
    // 
    //////////////////////////////////////////////////////

    WmSwAddRodToolCommand* m_toolCommand;
    
    WmSweeneyNode* m_sweeneyNode;
    WmSweeneyRodManager* m_rodManager;
};

}

#endif
