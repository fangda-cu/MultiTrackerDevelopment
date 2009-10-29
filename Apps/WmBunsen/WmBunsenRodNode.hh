#ifndef WMBUNSENRODNODE_HH_
#define WMBUNSENRODNODE_HH_

#include <maya/MPxLocatorNode.h> 
#include <maya/MString.h> 
#include <maya/MTypeId.h> 
#include <maya/MPlug.h>
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
#include <maya/MFnNurbsCurveData.h>

#include "Beaker.hh"

class WmBunsenRodNode : public MPxLocatorNode 
{
public:
    WmBunsenRodNode();
    virtual	~WmBunsenRodNode();
    virtual MStatus compute( const MPlug& i_plug, MDataBlock& i_dataBlock );
	virtual void draw( M3dView& i_view, const MDagPath& i_path, 
                       M3dView::DisplayStyle i_style,
                       M3dView::DisplayStatus i_status );
    virtual bool isBounded() const;
    static void* creator();
    static MStatus initialize();

    static MTypeId ia_typeID;
    static MString ia_typeName;
    static MObject ia_time;
    static MObject ia_startTime;
    static MObject ia_nurbsCurves;

    static MObject ca_syncAttrs;
    static MObject oa_rodsChanged;

private:
    double m_currentTime;
    double m_previousTime;
    double m_startTime;
    
    bool m_initialised;
    
    Beaker* m_beaker;
};

#endif
