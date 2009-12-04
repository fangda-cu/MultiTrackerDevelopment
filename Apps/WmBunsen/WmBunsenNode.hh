#ifndef WMBUNSENNODE_HH_
#define WMBUNSENNODE_HH_

#include <maya/MPxLocatorNode.h> 
#include <maya/MString.h> 
#include <maya/MTypeId.h> 
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MMatrix.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MDataBlock.h>
#include <maya/MDataHandle.h>
#include <maya/MColor.h>
#include <maya/M3dView.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MFnNurbsSurfaceData.h>
#include <maya/MFnPointArrayData.h>
#include <maya/MFnMatrixData.h>
#include <maya/MIOStream.h>
#include <maya/MTime.h>
#include <maya/MGlobal.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MVectorArray.h>

#include "Beaker.hh"
#include "WmBunsenRodNode.hh"

class WmBunsenNode : public MPxLocatorNode 
{
public:
    WmBunsenNode();
    virtual	~WmBunsenNode();
    virtual MStatus compute( const MPlug& i_plug, MDataBlock& i_dataBlock );
	virtual void draw( M3dView& i_view, const MDagPath& i_path, 
                       M3dView::DisplayStyle i_style,
                       M3dView::DisplayStatus i_status );
    virtual bool isBounded() const;
    virtual MStatus connectionMade( const  MPlug & i_plug, const  MPlug & i_otherPlug, bool i_asSrc );
    virtual MStatus connectionBroken( const  MPlug & i_plug, const  MPlug & i_otherPlug, bool i_asSrc );

    static void* creator();
    static MStatus initialize();
    
    static MTypeId typeID;
    static MString typeName;
    static MObject ia_time;
    static MObject ia_maxDt;
    static MObject ia_maxIter; //maximum number of newton iterations
    static MObject ia_startTime;
    static MObject ia_fps;
    static MObject ia_substeps;
    static MObject ia_rodsNodes;
    static MObject ia_gravity;
    static MObject ia_numberOfThreads;
    static MObject ia_solver;
    static MObject ia_collisionsEnabled;
    
    static MObject ia_collisionMeshes;
    
    static MObject ca_syncAttrs;
    static MObject oa_simStepTaken;
    
private:
    void pullOnAllRodNodes( MDataBlock& i_dataBlock );
    void createRodDataFromRodNodes( MDataBlock& i_dataBlock, 
                                    ObjectControllerBase::SolverLibrary solverLibrary );
    void updateAllCollisionMeshes( MDataBlock& i_dataBlock );
    
    double m_currentTime;
    double m_previousTime;
    double m_startTime;
    double m_framedt;
    
    bool m_initialised;
    
    Beaker* m_beaker;
};

#endif
