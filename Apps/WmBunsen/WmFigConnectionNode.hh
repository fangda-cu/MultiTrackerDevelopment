#ifndef WMFIGCONNECTIONODE_HH_
#define WMFIGCONNECTIONODE_HH_

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
#include <maya/MFnIntArrayData.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MVectorArray.h>
#include <maya/MFnNurbsCurveData.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnStringData.h>
#include <maya/MFloatArray.h>
#include <maya/MFnMatrixAttribute.h>

#include "Beaker.hh"

//static const int FILE_FORMAT_VERSION = 1;

class WmFigConnectionNode : public MPxLocatorNode 
{
public:
    WmFigConnectionNode();
    virtual	~WmFigConnectionNode();
    virtual MStatus compute( const MPlug& i_plug, MDataBlock& i_dataBlock );
	virtual void draw( M3dView& i_view, const MDagPath& i_path, 
                       M3dView::DisplayStyle i_style,
                       M3dView::DisplayStatus i_status );
    virtual bool isBounded() const;
    static void* creator();
    static MStatus initialize();
    virtual MStatus connectionMade( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc );
    virtual MStatus connectionBroken( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc );

    static MTypeId typeID;
    static MString typeName;
    static MObject ia_time;
    static MObject ia_startTime;
    
    // Attributes defining which rod and edge this connecter is attached to
    static MObject ia_rodNumber;
    static MObject ia_edgeNumber;
    static MObject ia_transformMatrix;
    
   // static MObject oa_materialFrame;
    
    static MObject ia_controllingEdge;
    static MObject ia_rodEdgeTransforms;
    static MObject oa_outTransformMatrix;
    static MObject oa_edgeTransform;

    static MStatus addNumericAttribute( MObject& i_attribute, MString i_longName, 
                                        MString i_shortName,
                                        MFnNumericData::Type i_type, double i_defaultValue,
                                        bool i_isInput, bool i_isOutput  );
    
    void getControlledRodInfo( unsigned int& o_rodIndex, unsigned int& o_edgeIndex, 
                               EdgeTransform& o_edgeTransform );
private:
    void calculateMaterialFrame();
    
    double m_currentTime;
    double m_previousTime;
    double m_startTime;
    bool m_controllingEdge;    
    
    unsigned int m_controlledRodIndex;
    unsigned int m_controlledEdgeIndex;
    EdgeTransform m_edgeTransform;
    
    MMatrix m_inputTransformMatrix;
};

#endif
