#ifndef WMFIGRODNODE_HH_
#define WMFIGRODNODE_HH_

/**
  * @file Apps/WmFigaro/WmFigRodNode.hh
  * @author Alasdair Coull (acoull@wetafx.co.nz)
  * @date 20-04-2010
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

#include <set>

#include "Beaker.hh"
#include "WmFigRodNurbsInput.hh"
#include "WmFigRodBarbInput.hh"
#include "WmFigRodFileInput.hh"
#include "WmFigRodGroup.hh"


static const int FILE_FORMAT_VERSION = 1;

typedef std::tr1::unordered_map<unsigned int, EdgeTransform> EdgeTransformMap;
typedef std::tr1::unordered_map<unsigned int, EdgeTransformMap > EdgeTransformRodMap;

class WmFigRodNode : public MPxLocatorNode 
{
public:
    WmFigRodNode();
    virtual	~WmFigRodNode();
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
    static MObject oa_nurbsCurves;
    static MObject ia_barberShopVertices;
    static MObject ia_percentageOfBarberShopStrands;    
    
    // Rod options
    static MObject ia_solverType;
    static MObject ia_cvsPerRod;
    static MObject ia_gravity;
    static MObject ia_youngsModulus;
    static MObject ia_shearModulus;
    static MObject ia_viscosity;
    static MObject ia_density;
    static MObject ia_minorRadius;
    static MObject ia_majorRadius;
    static MObject ia_vertexSpacing;
    static MObject ia_minimumRodLength;
    static MObject ia_hairSpray;
    static MObject ia_hairSprayScaleFactor;
    static MObject ia_massDamping;
    static MObject ia_drawMaterialFrames;
    static MObject ia_lockFirstEdgeToInput;
    static MObject ia_simulationSet;
    
    // Drawing 
    static MObject ia_userDefinedColors;
    static MObject ia_draw3DRod;
    static MObject ia_drawScale;
    static MObject ca_drawDataChanged;
    
    // Caching
    static MObject ia_cachePath;
    static MObject ia_cacheFrame;
    static MObject ia_readFromCache;
        // If the user changes ia_readFromCache this will evaluated and caused the node to 
        // stop doing anything until time==startTime
    static MObject ca_cachingHasBeenToggled; 
    
    static MObject ca_syncAttrs;
    static MObject oa_rodsChanged;
    static MObject ia_simStepTaken;
    static MObject ca_simulationSync;
    
    // Output attributes heading to Fozzie
    static MObject oa_simulatedVertices;
    static MObject oa_nonSimulatedVertices;
    static MObject oa_verticesInEachRod;
    
    // Frames from and too Barbershop
    static MObject oa_materialFrames;
    static MObject oa_undeformedMaterialFrames;
    static MObject ia_strandRootFrames;
    
    static MObject ia_edgeTransforms;
    static MObject oa_edgeTransforms;
    
    static MObject oa_numberOfRods;

    // Returns the number of rods this node has input data for
    //size_t numberOfRods();

   // void initialiseRodData( vector<RodData*>* i_rodDataMap );
    
    static MStatus addNumericAttribute( MObject& i_attribute, MString i_longName, 
                                        MString i_shortName,
                                        MFnNumericData::Type i_type, double i_defaultValue,
                                        bool i_isInput = true, bool i_isArray = false );
    
    ///////////////////////////////////////////////////////////
    //
    // Accessors
    //
    ///////////////////////////////////////////////////////////
    
    MMatrix getRodEdgeMatrix( size_t i_rod, size_t i_edge );
    
    /*vector<RodData*>* getRodData()
    {
        return mx_rodData;
    }*/

    WmFigRodGroup* rodGroup()
    {
        return &m_rodGroup;
    }
    
private:
    // Compute helper functions
    void compute_oa_simulatedVertices( const MPlug& i_plug, MDataBlock& i_dataBlock );
    void compute_ca_simulationSync( const MPlug& i_plug, MDataBlock i_dataBlock );
    void compute_oa_rodsChanged( const MPlug& i_plug, MDataBlock& i_dataBlock );
    void compute_oa_nonSimulatedVertices( const MPlug& i_plug, MDataBlock& i_dataBlock );
    void compute_oa_verticesInEachRod( const MPlug& i_plug, MDataBlock& i_dataBlock );
    void compute_oa_materialFrames( const MPlug& i_plug, MDataBlock& i_dataBlock );
    void compute_oa_undeformedMaterialFrames( const MPlug& i_plug, MDataBlock& i_dataBlock );
    void compute_oa_EdgeTransforms( const MPlug& i_plug, MDataBlock& i_dataBlock );
    void compute_ca_drawDataChanged( const MPlug& i_plug, MDataBlock& i_dataBlock );
    void readCacheRelatedInputs( MDataBlock& i_dataBlock );
    void writeCacheIfNeeded( MDataBlock& i_dataBlock );

    void updateControlledEdgeArrayFromInputs( MDataBlock& i_dataBlock );
    void updateOrInitialiseRodDataFromInputs( MDataBlock& i_dataBlock );
    
    void getStrandRootFrames( MDataBlock& i_dataBlock, vector<MaterialFrame>& o_strandRootFrames );
    void updateHairsprayScales( MDataBlock& i_dataBlock );
    void initialiseRodData( MDataBlock& i_dataBlock );
    void updateSimulationSet( MString i_simulationSetString );

    MString getCacheFilename( MDataBlock& i_dataBlock );
    
    void updateKinematicEdgesFromInput();

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
  //  vector<RodData*>* mx_rodData;
    WmFigRodGroup m_rodGroup;

    RodOptions m_rodOptions;
    
    World* mx_world;
    
    size_t m_numberOfInputCurves;
    
    double m_percentageOfBarberShopStrands;
    
    // If we're overriding the number of cvs per rod then this will be not -1
    int m_verticesPerRod;
    
    MString m_cacheFilename;
  
    bool m_lockFirstEdgeToInput;
    
    //vector<MaterialFrame> m_materialFrames;
    vector<MaterialFrame> m_strandRootFrames;
    
    // this is a map of maps. The first map is for the rod, the second is for the edge.
    EdgeTransformRodMap m_controlledEdgeTransforms;
    
    std::tr1::unordered_map< unsigned int, Vec3d > m_rodColourMap;

    WmFigRodInputType* m_pRodInput;

    double m_vertexSpacing;
    double m_minimumRodLength;

    bool m_readFromCache;
    bool m_writeToCache;
    MString m_cachePath;
    RodTimeStepper::Method m_solverType;

    Vec3d m_gravity;
    std::set< size_t > m_simulationSet;

};

#endif // WMFIGRODNODE_HH_
