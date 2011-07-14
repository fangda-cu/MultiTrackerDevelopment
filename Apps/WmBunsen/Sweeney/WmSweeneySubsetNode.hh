#ifndef WMSWEENEYSUBSETNODE_HH_
#define WMSWEENEYSUBSETNODE_HH_


//////////////////////////////////////////////////////////
// 
// Headers.
// 
//////////////////////////////////////////////////////////

#include <maya/MPxLocatorNode.h>
#include <maya/M3dView.h>
#include <maya/MDagPath.h>

//////////////////////////////////////////////////////////
// 
// Class declaration.
// 
//////////////////////////////////////////////////////////

class WmSweeneySubsetNode : public MPxLocatorNode
{
  public:


    //////////////////////////////////////////////////////
    // 
    // Stuff we always need.
    // 
    //////////////////////////////////////////////////////

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
    void setScalpFaceIndices ( const MIntArray i_indicies );

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
