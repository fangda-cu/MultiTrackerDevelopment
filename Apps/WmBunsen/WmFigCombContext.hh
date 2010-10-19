#ifndef _WMFIGCOMBCONTEXT_HH_
#define _WMFIGCOMBCONTEXT_HH_

#include "WmBunsenNode.hh"

#include "WmFigCombToolCommand.hh"
#include "WmFigRodNode.hh"

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
#include <maya/MFloatMatrix.h>
#include <maya/MFnCamera.h>

class WmFigCombContext : public MPxContext 
{
public:
    enum MouseStatus { kPress, kDrag, kRelease };
    
    WmFigCombContext();
    
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
    
    static WmFigCombContext* activeContext;
    static MString typeName;
    
    //////////////////////////////////////////////////////
    // 
    // context functions.
    // 
    //////////////////////////////////////////////////////
    
     /**
     * @brief Implement this method to return the current 3D projection matrix of the camera being used to draw the furSet.
     * @return 16 floats (right handed coordinate system)
     */
    const float* projectionMatrix();
    
    /**
     * @brief Implement this method to return the current 3D modelview matrix of the camera being used to draw the furSet.
     * @return 16 floats (right handed coordinate system)
     */
    const float* modelViewMatrix();

private:
    void closestPointOnLine( const MPoint& i_pPt, const MPoint& i_aPt, const MPoint& i_bPt, MPoint& o_xPt, bool i_segment = false );
    void applyForceToRods();
    GLint findRodsUsingOpenGLSelection( const double i_centreX, const double i_centreY,
            const double i_width, const double i_height,  WmFigRodNode* rodNode,
            vector<GLuint>& o_selectedRodIndices );
    
    bool searchForRodsIn2DScreenRectangle( vector<int>& o_rodIndices );
    bool findRodNode();
    void makeRodNodeEvaluate();

    MStatus drawMarqueeSelectBox() const;
    void doTool( MEvent &event );

    MString m_helpString;
    
    short m_xStartMouse, m_yStartMouse;
    short m_xMouse, m_yMouse;
    short m_xMouseStart, m_yMouseStart;
    
    //////////////////////////////////////////////////////
    //
    // Projection and modelview matrix, stored in case they
    // are needed when Maya has changed them to 2D for
    // selection
    //
    //////////////////////////////////////////////////////
    
    mutable MFloatMatrix m_projectionMatrix;
    mutable MFloatMatrix m_modelViewMatrix;
    
    //////////////////////////////////////////////////////
    // 
    // Values sent to the tool command, and the tool command
    // itself
    // 
    //////////////////////////////////////////////////////

    WmFigCombToolCommand* m_toolCommand;
    WmFigRodNode* m_figRodNode;
    MObject m_figRodNodeMObject;
};

#endif
