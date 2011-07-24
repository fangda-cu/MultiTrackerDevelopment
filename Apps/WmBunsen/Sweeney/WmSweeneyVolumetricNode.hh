#ifndef _WMSWEENEYVOLUMETRICNODE_HH_
#define _WMSWEENEYVOLUMETRICNODE_HH_

#ifdef WETA
#include <weta/Wfigaro/Core/TriangleMesh.hh>
#include <weta/Wfigaro/Render/TriangleMeshRenderer.hh>
#else
#include "BASim/src/Core/TriangleMesh.hh"
#include "BASim/src/Render/TriangleMeshRenderer.hh"
#endif

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
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnStringData.h>
#include <maya/MFnMesh.h>

using namespace BASim;

class WmSweeneyVolumetricNode : public MPxLocatorNode
{
public:
    WmSweeneyVolumetricNode();
	virtual	~WmSweeneyVolumetricNode();
    virtual MStatus compute( const MPlug& i_plug, MDataBlock& i_data );
	virtual void draw( M3dView & view, const MDagPath & path, 
                       M3dView::DisplayStyle style,
                       M3dView::DisplayStatus status );

	virtual bool isBounded() const;
	
	virtual MStatus connectionMade( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc );
	virtual MStatus connectionBroken( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc );


    virtual void postConstructor();
	
	static  void* creator();
	static  MStatus initialize();

	Eigen::Quaternion<double> getQuaternion( )
    {
	    return m_quaternion;
    }

    Vec3d getCenter( )
    {
        return m_center;
    }

    Vec3d getScale( )
    {
        return m_scale;
    }

	static MTypeId typeID;
    static MString typeName;
	static MObject ia_inMesh;
    static MObject ia_startTime;
    static MObject ia_time;
    static MObject oa_meshData;

    static MObject ia_meshTransform;

    static MObject ia_charge;
    
    
   void initialise( const unsigned int i_volumetricMeshIndex,
                     TriangleMesh** o_currentMesh = NULL );
        			
private:
    MStatus updateVolumetricMeshFromMayaMesh( MFnMesh &meshFn, bool forceReset = false,
                                             std::string filename = "" );
    double m_currentTime;
    double m_previousTime;
    double m_startTime;
    
    double m_charge;

    Eigen::Quaternion<double> m_quaternion;
    Vec3d m_center;
    Vec3d m_scale;

    bool m_enableCaching;
    MString m_cacheDirectory;
    TriangleMesh* m_nextMesh;        // Mesh at next full frame
    TriangleMesh* m_previousMesh;    // mesh at previous full frame
    TriangleMesh* m_currentMesh;     // Mesh at current in between frame in simulation substeps
    TriangleMeshRenderer* m_triangleMeshRenderer;
    bool m_meshConnected;

    unsigned int m_volumetricMeshIndex;
    
    Eigen::Matrix4f m_meshTransformMatrix;
};


#endif
