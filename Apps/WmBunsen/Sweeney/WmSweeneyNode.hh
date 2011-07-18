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
#include <maya/MDagPath.h>
#include <maya/MColor.h>
#include <maya/M3dView.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnNurbsSurface.h>
#include <maya/MFnNurbsSurfaceData.h>
#include <maya/MFnPointArrayData.h>
#include <maya/MFnMatrixData.h>
#include <maya/MFnMesh.h>
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

#include <maya/MFnMatrixAttribute.h>
#include <maya/MPlugArray.h>
#include <maya/MQuaternion.h>

#include <sys/stat.h>

#include "WmSweeneyRodManager.hh"
#include "WmSweeneySubsetNode.hh"

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
    static MObject ia_rodsPerClump;
    static MObject ia_rodRadius;
    static MObject ia_rodAspectRatio;
    static MObject ia_rodRotation;
    static MObject ia_curlTightness;
    static MObject ia_curlRadius;
    static MObject ia_curlCount;
    static MObject ia_curlStart;
    static MObject ia_rodPitch;
    static MObject ia_fixCurlCount;
    static MObject ia_curlInXFrame;
    static MObject ia_preserveLengthVariation;
    static MObject ia_rodDamping;
    static MObject ia_rodCharge;
    static MObject ia_rodPower;
    static MObject ia_rodClumpSeparation;
    static MObject ia_rodClumpingRamp;

    // Sweeney Subset connection
    static MObject ia_fromSweeneySubsetNodes;

    // Collision meshes
    static MObject ia_collisionMeshes;
    
    //Solver Tolerances
    static MObject ia_stol;
    static MObject ia_atol;
    static MObject ia_rtol;
    static MObject ia_inftol;
    static MObject ia_numLineSearchIters;

    // Performance Tuning
      //GeneralParameters
    static MObject ia_enablePenaltyResponse;
    static MObject ia_implicitThickness;
    static MObject ia_implicitStiffness;
    static MObject ia_levelsetSubsampling;
    static MObject ia_inextensibilityThreshold;

      //Failuredetection
    static MObject ia_maxNumOfSolverIters;
    static MObject ia_maxNumOfCollisionIters;
    static MObject ia_enableExplosionDetection;
    static MObject ia_explosionDampening;
    static MObject ia_explosionThreshold;
    static MObject ia_stretchingThreshold;

     //FailureResponse
    static MObject ia_solverFailure;
    static MObject ia_collisionFailure;
    static MObject ia_explosionFailure;
    static MObject ia_stretchingFailure;
    static MObject ia_maxNumSolverSubsteps;
    static MObject ia_maxNumCollisionSubsteps;
    static MObject ia_maxNumExplosionSubsteps;
    static MObject ia_maxNumStretchingSubsteps;

    // Sync attributes to force compute() when inputs change
    static MObject ca_rodPropertiesSync;
    
    // Output to the guide curve deformer
    static MObject oa_simulatedNurbs;
    
    static MObject ia_strandVertices;
    static MObject ia_strandRootFrames;
    //static MObject ia_strandLengths;
    static MObject ia_verticesPerStrand;
    
    // Debug drawing
    static MObject ia_shouldDrawStrands;
    static MObject ia_shouldDrawRootFrames;
    static MObject ia_shouldDrawVelocity;
    static MObject ia_shouldDrawOnlySelected;

    static MStatus addNumericAttribute( MObject& i_attribute, MString i_longName, 
                                        MString i_shortName,
                                        MFnNumericData::Type i_type, double i_defaultValue,
                                        bool i_isInput = true, bool i_isArray = false );
                                        
    void constructRodVertices( std::vector< BASim::Vec3d >& o_rodVertices, const MVector& i_direction,
                       const MVector& i_rootPosition );

    MStatus subsetNodes( std::vector<WmSweeneySubsetNode*>& o_subsetNodes );
    WmSweeneyRodManager* rodManager();
        
private:
    void initialiseRodFromBarberShopInput( MDataBlock& i_dataBlock );
    void initialiseCollisionMeshes( MDataBlock &i_data );
    void updateCollisionMeshes( MDataBlock& i_dataBlock );
    void computeSubsetRodMapping( MDataBlock& i_dataBlock );
    void compute_oa_simulatedNurbs( const MPlug& i_plug, MDataBlock& i_dataBlock );

    void updateRods( );
    void updateRodParameters( const int& rodIdx, const bool& update_all_rods );

    // reformat all of these to take the updated value as a parameter
    void updateStrandLength( BASim::ElasticRod* current_rod, bool& update_rod,
            const BASim::Scalar updated_edge_length );
    void updateStrandCrossSection( BASim::ElasticRod* current_rod, bool& update_rod,
            const double i_rodRadius, const double i_rodAspectRatio );
    void updateStrandRotation( BASim::ElasticRod* current_rod, bool& update_rod,
            const double i_rodRotation );
    void updateStrandCurl( BASim::ElasticRod* current_rod,  bool& update_rod,
            BASim::Scalar curvature, BASim::Scalar torsion, const double i_curlStart,
            const double i_verticesPerRod, const bool i_curlInXFrame );
    // TODO(sainsley) : right now solver settings are on a global level
    // we will need to refactor this to a rod-level if we want finer
    // control of solver settings in subsets
    void updateSolverSettings( MDataBlock &i_dataBlock );

    // Helper functions for scalp mesh
    void getSurfaceTangent( MFnMesh& surface, BASim::Vec3d& surface_tan, int rodIdx, const BASim::Vec3d strand_tan );
    bool getScalpTangents( std::vector<BASim::Vec3d>& i_scalpTangents );
    void getMeshPath( MDagPath& meshPath, MStatus& status );
    void getSubsetRods( WmSweeneySubsetNode* subset, std::vector< BASim::ElasticRod* >& subset_rods );

    double m_currentTime;
    double m_previousTime;
    double m_startTime;

    // Hair properties
    double m_edgeLength;
    double m_length;
    double m_rodRadius;
    double m_rodAspectRatio;
    double m_rodRotation;
    double m_curlTightness;
    double m_curlRadius;
    double m_curlCount;
    double m_curlStart;
    double m_rodPitch;
    int m_verticesPerRod;
    int m_rodsPerClump;
    bool m_fixCurlCount;
    bool m_curlInXFrame;
    bool m_preserveLengthVariation;
    bool m_rodDamping;
    double m_rodCharge;
    double m_rodPower;
    double m_rodClumpSeparation;
    double m_rodClumpingRamp [];
    // Debug drawing
    bool m_shouldDrawStrands;
    bool m_shouldDrawRootFrames;
    bool m_shouldDrawVelocity;
    
    WmSweeneyRodManager* m_rodManager;
    MVectorArray m_strandVertices;
    MVectorArray m_strandRootFrames;

    // strand lengths from barbershop
    std::vector<double> m_strandLengths;

    std::vector<WmSweeneySubsetNode*> m_subsetNodes;

    unsigned int m_numberOfVerticesPerStrand;
};

#endif // WmSweeneyNode_HH_
