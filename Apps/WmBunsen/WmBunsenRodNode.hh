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

    void initialiseRodData( vector<RodData*>* mx_rodDataMap );
    
private:
    double m_currentTime;
    double m_previousTime;
    double m_startTime;
    
    bool m_initialised;
    
    /// This is a pointer to the rod data stored in the beaker class that this node is supplying
    /// input data for. If it is not NULL then the beaker node has created the rods and whenever
    /// the BunsenNode asks for our rod data we actually dump it into this pointer. For efficiency
    /// there is  no point in passing it to Bunsen then having it stick it into Beaker. So we
    /// use the connection to Bunsen as a flag to indicate it's time to update the data in this
    /// pointer.
    //vector<RodData*>* mx_rodData;
    
    // We take a pointer to the map not the individual vector of rodData for this node because
    // the map is likely to be resized and the data will move around. The only pointer that will
    // not move is the pointer to the map.
    RodDataMap* mx_rodDataMap;
    
    Beaker* mx_beaker;
};

#endif
