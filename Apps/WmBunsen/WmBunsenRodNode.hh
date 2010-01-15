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
#include <maya/MFnIntArrayData.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MVectorArray.h>
#include <maya/MFnNurbsCurveData.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnStringData.h>
#include <maya/MFloatArray.h>
#include <maya/MRampAttribute.h>

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
    virtual MStatus connectionMade( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc );
    virtual MStatus connectionBroken( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc );

    static MTypeId typeID;
    static MString typeName;
    static MObject ia_time;
    static MObject ia_startTime;
    static MObject ia_nurbsCurves;
    static MObject ia_fozzieVertices;
    static MObject ia_percentageOfFozzieStrands;
    
    // Rod options
    static MObject ia_cvsPerRod;
    static MObject ia_youngsModulus;
    static MObject ia_shearModulus;
    static MObject ia_density;
    static MObject ia_minorRadius;
    static MObject ia_majorRadius;
    static MObject ia_vertexSpacing;
    static MObject ia_hairSpray;
    static MObject ia_hairSprayScaleFactor;
    static MObject ia_massDamping;

    // Caching
    static MObject ia_cachePath;
    static MObject ia_cacheFrame;
    static MObject ia_readFromCache;
    
    static MObject ca_syncAttrs;
    static MObject oa_rodsChanged;
    static MObject ia_simStepTaken;
    static MObject ca_simulationSync;
    
    // Output attributes heading to Fozzie
    static MObject oa_simulatedVertices;
    static MObject oa_nonSimulatedVertices;
    static MObject oa_verticesInEachRod;
    
    // Returns the number of rods this node has input data for
    size_t numberOfRods();

    void initialiseRodData( vector<RodData*>* i_rodDataMap );
    void updateRodDataFromInputs();
    
    static MStatus addNumericAttribute( MObject& i_attribute, MString i_longName, 
                                        MString i_shortName,
                                        MFnNumericData::Type i_type, double i_defaultValue,
                                        bool i_isInput );
    
private:
    void writeRodDataToCacheFile( MString i_cachePath );
    void readRodDataFromCacheFile( MString i_cachePath );
    void updateHairsprayScales( MDataBlock& i_dataBlock );
    
    double m_currentTime;
    double m_previousTime;
    double m_startTime;
    double m_massDamping;
    
    bool m_initialised;
    
    /// This is a pointer to the rod data stored in the beaker class that this node is supplying
    /// input data for. If it is not NULL then the beaker node has created the rods and whenever
    /// the BunsenNode asks for our rod data we actually dump it into this pointer. For efficiency
    /// there is  no point in passing it to Bunsen then having it stick it into Beaker. So we
    /// use the connection to Bunsen as a flag to indicate it's time to update the data in this
    /// pointer.
    vector<RodData*>* mx_rodData;
    RodOptions m_rodOptions;
    
    World* mx_world;
    
    size_t m_numberOfInputCurves;
    
    int m_percentageOfFozzieStrands;
    
    // If we're overriding the number of cvs per rod then this will be not -1
    int m_cvsPerRod;
};

#endif
