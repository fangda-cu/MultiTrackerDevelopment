#ifndef WMSWEENEYSUBSETNODE_HH_
#define WMSWEENEYSUBSETNODE_HH_

#include <maya/MPxLocatorNode.h>
#include <maya/M3dView.h>
#include <maya/MDagPath.h>


class WmSweeneySubsetNode : public MPxLocatorNode
{
  public:

    WmSweeneySubsetNode();
    
    virtual ~WmSweeneySubsetNode();

    
    static const MTypeId typeID;
    static const MString typeName;
        
    static void* creator();
    static MStatus initialize();
    

    virtual void postConstructor();
    
    virtual MStatus compute( 
        const MPlug& i_plug, 
        MDataBlock& i_dataBlock );

    virtual void draw( 
        M3dView& i_view, 
        const MDagPath& i_path, 
        M3dView::DisplayStyle i_style, 
        M3dView::DisplayStatus i_status );
    
    virtual bool isBounded() const;


    // Method to set the selected faces array
    void setScalpFaceIndices( const MIntArray i_indicies );
    // debug method for previous method
    void checkScalpFaceIndices();

    //////////////////////////////////////////////////////
    //
    // Accessor methods.
    //
    //////////////////////////////////////////////////////

    MIntArray getScalpFaceIndices( MDataBlock* i_dataBlock = NULL ) const;

    // Visibility
    bool isVisible( ) const;

    // Rod property accessors

    double getRodLength( MDataBlock* i_dataBlock = NULL ) const;
    double getRodRadius( MDataBlock* i_dataBlock = NULL ) const;
    double getRodAspectRatio( MDataBlock* i_dataBlock = NULL ) const;
    double getRodRotation( MDataBlock* i_dataBlock = NULL ) const;
    double getCurlTightness( MDataBlock* i_dataBlock = NULL ) const;
    double getCurlCount( MDataBlock* i_dataBlock = NULL ) const;
    double getCurlRadius( MDataBlock* i_dataBlock = NULL ) const;
    double getCurlStart( MDataBlock* i_dataBlock = NULL ) const;
    double getRodCharge( MDataBlock* i_dataBlock = NULL ) const;
    double getRodPower( MDataBlock* i_dataBlock = NULL ) const;
    double getRodClumpSeparation( MDataBlock* i_dataBlock = NULL ) const;

    int getVerticesPerRod( MDataBlock* i_dataBlock = NULL ) const;
    int getNumberOfClumps( MDataBlock* i_dataBlock = NULL ) const;

    bool getIsFixCurlCount( MDataBlock* i_dataBlock = NULL ) const;
    bool getIsCurlInXFrame( MDataBlock* i_dataBlock = NULL ) const;
    bool getIsPreserveLengthVariation( MDataBlock* i_dataBlock = NULL ) const;
    bool getIsRodDamping( MDataBlock* i_dataBlock = NULL ) const;

    // Solver settings accessor

    void getSolverSettings(
            double& i_stol, double& i_atol, double& i_rtol,
            double& i_inftol, int& numLineSearchIters,
            MDataBlock* i_dataBlock = NULL );


    //////////////////////////////////////////////////////
    // 
    // Attributes.
    // 
    //////////////////////////////////////////////////////
    
    // Scalp Face Indices
    static MObject ia_scalpFaceIndices;

    // Rod properties
    static MObject ia_length;
    static MObject ia_edgeLength;
    static MObject ia_verticesPerRod;
    static MObject ia_numberOfClumps;
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

    // Solver Tolerances
    static MObject ia_stol;
    static MObject ia_atol;
    static MObject ia_rtol;
    static MObject ia_inftol;
    static MObject ia_numLineSearchIters;

    // Connection to Sweeney Node
    static MObject oa_toSweeneyParentNode;


  private:

    static MStatus addNumericAttribute( MObject& i_attribute, MString i_longName,
                                                MString i_shortName,
                                                MFnNumericData::Type i_type, double i_defaultValue,
                                                bool i_isInput = true, bool i_isArray = false );

};

#endif // WMSWEENEYSUBSETNODE_H_