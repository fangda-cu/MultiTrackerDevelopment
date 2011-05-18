#ifndef WMSWEENEYNODE_HH_
#define WMSWEENEYNODE_HH_

/**
  * @file Apps/WmFigaro/Sweeney/WmSweeneyNode.hh
  * @author Alasdair Coull (acoull@wetafx.co.nz)
  * @date 17-05-2011
  */

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
#include <maya/MFnEnumAttribute.h>
#include <maya/MFileIO.h>
#include <maya/MFnCompoundAttribute.h>

#include "WmSweeneyRodManager.hh"

class WmSweeneyNode : public MPxLocatorNode 
{
public:
    WmSweeneyNode();
    virtual	~WmSweeneyNode();
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
    
    // Rod properties
    static MObject ia_length;
    static MObject ia_edgeLength;
    static MObject ia_verticesPerRod;
    static MObject ia_rodRadius;
    static MObject ia_rodPitch;
    
    // Collision meshes
    static MObject ia_collisionMeshes;
    
    // Sync attributes to force compute() when inputs change
    static MObject ca_rodPropertiesSync;
    
    static MObject ia_strandVertices;
    static MObject ia_verticesPerStrand;
    
    static MStatus addNumericAttribute( MObject& i_attribute, MString i_longName, 
                                        MString i_shortName,
                                        MFnNumericData::Type i_type, double i_defaultValue,
                                        bool i_isInput = true, bool i_isArray = false );
                                        
    void constructRodVertices( std::vector< BASim::Vec3d >& o_rodVertices, const MVector& i_direction,
                       const MVector& i_rootPosition );
        
private:
    void initialiseRodFromBarberShopInput( MDataBlock& i_dataBlock );
    void initialiseCollisionMeshes( MDataBlock &i_data );
    void updateCollisionMeshes( MDataBlock& i_dataBlock );
    
    double m_currentTime;
    double m_previousTime;
    double m_startTime;

    // Hair properties
    double m_edgeLength;
    double m_length;
    double m_rodRadius;
    double m_rodPitch;
    int m_verticesPerRod;
    
    WmSweeneyRodManager* m_rodManager;
    MVectorArray m_strandVertices;
    unsigned int m_numberOfVerticesPerStrand;
};

#endif // WmSweeneyNode_HH_
