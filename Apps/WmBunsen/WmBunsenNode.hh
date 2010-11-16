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
#include <maya/MFnEnumAttribute.h>

#include "Beaker.hh"
#include "WmFigRodNode.hh"

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
    static MObject ia_solverType;
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
    static MObject ia_enabled;
    static MObject oa_dt;
    
    //Solver Tolerances
    static MObject ia_stol;
    static MObject ia_atol;
    static MObject ia_rtol;
    static MObject ia_inftol;

    static MObject ia_collisionMeshes;

    static MObject ia_msgConstraintNodes;

    static MObject ia_plasticDeformations;
    static MObject ia_isClumpingEnabled;
    static MObject ia_clumpingCoefficient;
    static MObject ia_selfCollisionPenaltyForces;
    static MObject ia_fullSelfCollisions;
    static MObject ia_fullSelfCollisionIterations;
    static MObject ia_fullSelfCollisionCOR;
    
    // Volumetric Collisions
    static MObject ia_volumetricCollisions;
    static MObject ia_gridDX;
    static MObject ia_targetEdgeDensity;
    static MObject ia_volumetricRadius;
    static MObject ia_flip;
    static MObject ia_slip;
    static MObject ia_separationCondition;
    static MObject ia_displayGrid;
    static MObject ia_displayAirBoundary;
    static MObject ia_displayCollisionBoundary;
    static MObject ia_displayGridVelocitiesMultiplier;
    static MObject ia_maxDisplayDensity;

    static MObject ca_syncAttrs;
    static MObject oa_simStepTaken;
    
    // File to store results of timings taken during simulation
    static MObject ia_timingsFile;
    static MObject ia_timingEnabled;

    // XML Output info
    static MObject ia_writeToXMLFile;
    static MObject ia_XMLFilePath;
    
    // Drawing
    static MObject ia_drawSubSteppedVertices;
    
private:
    void pullOnAllRodNodes( MDataBlock& i_dataBlock );
    void createRodDataFromRodNodes( MDataBlock& i_dataBlock );
    void updateAllCollisionMeshes( MDataBlock& i_dataBlock );
    void updateAllRodNodes( MDataBlock &i_dataBlock );
    void addRodsToWorld( MDataBlock& i_dataBlock );

    void addAllConstraints( MDataBlock &i_dataBlock );
    void updateAllConstraints( MDataBlock &i_dataBlock );

    static MStatus addNumericAttribute( MObject& i_attribute, MString i_longName, MString i_shortName, 
                                 MFnNumericData::Type i_type, double i_defaultValue, bool i_isInput,
                                 bool i_isArray = false );
    void getAllVolumetricCollisionAttributes( MDataBlock& i_dataBlock );


    double m_currentTime;
    double m_previousTime;
    double m_startTime;
    double m_framedt;
    
    bool m_initialised;
    
    bool m_enabled;

    Beaker* m_beaker;

    RodTimeStepper::Method m_solverType;

    bool m_writeXMLData;
    std::string m_xmlFilePath;
};

#endif
