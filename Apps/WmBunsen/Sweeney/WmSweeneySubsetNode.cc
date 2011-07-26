
//////////////////////////////////////////////////////////
// 
// Headers.
// 
//////////////////////////////////////////////////////////

#include "WmSweeneySubsetNode.hh"

#include <maya/MFnNumericAttribute.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnEnumAttribute.h>
#include <maya/MFnDagNode.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MGlobal.h>
#include <maya/MPlug.h>
#include <maya/MBoundingBox.h>

const MTypeId WmSweeneySubsetNode::typeID( 0x001156, 0xA9 );
const MString WmSweeneySubsetNode::typeName( "wmSweeneySubsetNode" );


// Scalp Face Indices

/* static */ MObject WmSweeneySubsetNode::ia_scalpFaceIndices;

// Hair Property Attributes
/* static */ MObject WmSweeneySubsetNode::ia_length;
/* static */ MObject WmSweeneySubsetNode::ia_edgeLength;
/* static */ MObject WmSweeneySubsetNode::ia_verticesPerRod;
/* static */ MObject WmSweeneySubsetNode::ia_numberOfClumps;
/* static */ MObject WmSweeneySubsetNode::ia_rodRadius;
/* static */ MObject WmSweeneySubsetNode::ia_rodAspectRatio;
/* static */ MObject WmSweeneySubsetNode::ia_rodRotation;
/* static */ MObject WmSweeneySubsetNode::ia_curlTightness;
/* static */ MObject WmSweeneySubsetNode::ia_curlRadius;
/* static */ MObject WmSweeneySubsetNode::ia_curlCount;
/* static */ MObject WmSweeneySubsetNode::ia_curlStart;
/* static */ MObject WmSweeneySubsetNode::ia_rodPitch;
/* static */ MObject WmSweeneySubsetNode::ia_fixCurlCount;
/* static */ MObject WmSweeneySubsetNode::ia_curlInXFrame;
/* static */ MObject WmSweeneySubsetNode::ia_preserveLengthVariation;
/* static */ MObject WmSweeneySubsetNode::ia_rodDamping;
/* static */ MObject WmSweeneySubsetNode::ia_rodCharge;
/* static */ MObject WmSweeneySubsetNode::ia_rodPower;
/* static */ MObject WmSweeneySubsetNode::ia_rodClumpSeparation;

/* static */ MObject WmSweeneySubsetNode::ia_volumetricCharge;
/* static */ MObject WmSweeneySubsetNode::ia_volumetricScale;
/* static */ MObject WmSweeneySubsetNode::ia_drawGaussianVolume;

// Solver Tolerances
/* static */ MObject WmSweeneySubsetNode::ia_stol;
/* static */ MObject WmSweeneySubsetNode::ia_atol;
/* static */ MObject WmSweeneySubsetNode::ia_rtol;
/* static */ MObject WmSweeneySubsetNode::ia_inftol;
/* static */ MObject WmSweeneySubsetNode::ia_numLineSearchIters;

/* static */ MObject WmSweeneySubsetNode::oa_toSweeneyParentNode;

WmSweeneySubsetNode::WmSweeneySubsetNode()
{
}

/*virtual*/ void WmSweeneySubsetNode::postConstructor()
{
}

/*virtual*/ WmSweeneySubsetNode::~WmSweeneySubsetNode()
{
}

bool WmSweeneySubsetNode::isVisible( ) const
{
    MStatus status;

    MFnDagNode sweeneySubsetNode = MFnDagNode( thisMObject( ) );
    MFnDagNode sweeneySubsetTransformNode = MFnDagNode( sweeneySubsetNode.parent(0, &status ) );

    CHECK_MSTATUS( status );

    MPlug plug = sweeneySubsetTransformNode.findPlug( "visibility", true, &status );
    if ( !status )
    {
        status.perror( "cannot locate visibility plug for WmSweeneySubsetNode" );
        return false;
    }

    bool value = plug.asBool();

    plug = sweeneySubsetNode.findPlug( "visibility", true, &status );
    if ( !status )
    {
        status.perror( "cannot locate visibility plug for WmSweeneySubsetNode" );
        return false;
    }

    value = value && plug.asBool();

    return value;
}

void WmSweeneySubsetNode::setScalpFaceIndices( const MIntArray i_indices )
{
    cout << "WmSweeneySubsetNode::setScalpFaceIndices::setting indices for " << i_indices.length()
            << " faces " << endl;

    MPlug plug;
    MStatus status;

    MFnDagNode dagFn = MFnDagNode( thisMObject( ) );

    // grab face indices plug
    plug = dagFn.findPlug( ia_scalpFaceIndices, true, &status );
    if ( status != MStatus::kSuccess )
    {
        status.perror( "cannot locate scalp face indices from WmSweeneySubsetNode" );
        return;
    }

    MFnIntArrayData indicesDataFn;
    MObject indicesDataObj = indicesDataFn.create( & status );
    CHECK_MSTATUS( status );

    MIntArray indices = indicesDataFn.array( & status );
    CHECK_MSTATUS( status );

    indices = i_indices;

    status = plug.setValue( indicesDataObj );
    CHECK_MSTATUS( status );

    //checkScalpFaceIndices( );

}

void WmSweeneySubsetNode::checkScalpFaceIndices( )
{
    MPlug plug;
    MStatus status;

    MFnDagNode sweeneySubsetNode = MFnDagNode( thisMObject( ) );

    // grab face indices plug
    plug = sweeneySubsetNode.findPlug( ia_scalpFaceIndices, true, &status );
    if ( !status )
    {
        status.perror( "cannot locate scalp face indices from WmSweeneySubsetNode for indices check" );
        return;
    }

    MObject indicesDataObj = MObject::kNullObj;
    status = plug.getValue( indicesDataObj );
    CHECK_MSTATUS( status );

    MFnIntArrayData indicesDataFn( indicesDataObj, & status );
    CHECK_MSTATUS( status );

    MIntArray indices = indicesDataFn.array( & status );
    CHECK_MSTATUS( status );

    cout << "WmSweeneySubsetNode::checkScalpFaceIndices::begin face index list " << endl;
    for ( int i = 0; i < indices.length(); i++ )
    {
        cout << "face index: " << indices[ i ] << endl;
    }
    cout << "WmSweeneySubsetNode::checkScalpFaceIndices::end face index list " << endl;
}

MIntArray WmSweeneySubsetNode::getScalpFaceIndices( MDataBlock* i_dataBlock ) const
{

    MStatus status;

    MIntArray indices;

    MObject indicesDataObj;

    if ( i_dataBlock )
    {
        indicesDataObj = i_dataBlock->inputValue( ia_scalpFaceIndices, & status ).data();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_scalpFaceIndices );
        indicesDataObj = MObject::kNullObj;
        status = plug.getValue( indicesDataObj );
        CHECK_MSTATUS( status );
    }


    MFnIntArrayData indicesDataFn( indicesDataObj, & status );
    CHECK_MSTATUS( status );

    indices = indicesDataFn.array( & status );
    CHECK_MSTATUS( status );

    return indices;
}

// Rod property accessors

double WmSweeneySubsetNode::getRodLength( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    double value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_length, & status ).asDouble();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_length );
        value = plug.asDouble();
    }

    return value;
}

double WmSweeneySubsetNode::getRodRadius( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    double value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_rodRadius, & status ).asDouble();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_rodRadius );
        value = plug.asDouble();
    }

    return value;
}

double WmSweeneySubsetNode::getRodAspectRatio( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    double value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_rodAspectRatio, & status ).asDouble();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_rodAspectRatio );
        value = plug.asDouble();
    }

    return value;
}

double WmSweeneySubsetNode::getRodRotation( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    double value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_rodRotation, & status ).asDouble();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_rodRotation );
        value = plug.asDouble();
    }

    return value;
}

double WmSweeneySubsetNode::getCurlTightness( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    double value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_curlTightness, & status ).asDouble();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_curlTightness );
        value = plug.asDouble();
    }

    return value;
}

double WmSweeneySubsetNode::getCurlCount( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    double value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_curlCount, & status ).asDouble();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_curlCount );
        value = plug.asDouble();
    }

    return value;
}

double WmSweeneySubsetNode::getCurlRadius( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    double value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_curlRadius, & status ).asDouble();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_curlRadius );
        value = plug.asDouble();
    }

    return value;
}

double WmSweeneySubsetNode::getCurlStart( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    double value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_curlStart, & status ).asDouble();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_curlStart );
        value = plug.asDouble();
    }

    return value;
}

double WmSweeneySubsetNode::getRodCharge( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    double value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_rodCharge, & status ).asDouble();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_rodCharge );
        value = plug.asDouble();
    }

    return value;
}

double WmSweeneySubsetNode::getRodPower( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    double value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_rodPower, & status ).asDouble();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_rodPower );
        value = plug.asDouble();
    }

    return value;
}

double WmSweeneySubsetNode::getRodClumpSeparation( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    double value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_rodClumpSeparation, & status ).asDouble();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_rodClumpSeparation );
        value = plug.asDouble();
    }

    return value;
}

int WmSweeneySubsetNode::getVerticesPerRod( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    int value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_verticesPerRod, & status ).asInt();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_verticesPerRod );
        value = plug.asInt();
    }

    return value;
}

int WmSweeneySubsetNode::getNumberOfClumps( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    int value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_numberOfClumps, & status ).asInt();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_numberOfClumps );
        value = plug.asInt();
    }

    return value;
}

bool WmSweeneySubsetNode::getIsFixCurlCount( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    bool value;
    if ( i_dataBlock )
    {
        value = i_dataBlock->inputValue( ia_fixCurlCount, & status ).asBool();
        CHECK_MSTATUS( status );
    }
    else
    {
        MPlug plug( thisMObject(), ia_fixCurlCount );
        value = plug.asBool();
    }

    return value;
}

bool WmSweeneySubsetNode::getIsCurlInXFrame( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    bool value;
    if ( i_dataBlock )
    {
       value = i_dataBlock->inputValue( ia_curlInXFrame, & status ).asBool();
       CHECK_MSTATUS( status );
    }
    else
    {
       MPlug plug( thisMObject(), ia_curlInXFrame );
       value = plug.asBool();
    }

    return value;
}

bool WmSweeneySubsetNode::getIsPreserveLengthVariation( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    bool value;
    if ( i_dataBlock )
    {
       value = i_dataBlock->inputValue( ia_preserveLengthVariation, & status ).asBool();
       CHECK_MSTATUS( status );
    }
    else
    {
       MPlug plug( thisMObject(), ia_preserveLengthVariation );
       value = plug.asBool();
    }

    return value;
}

bool WmSweeneySubsetNode::getIsRodDamping( MDataBlock* i_dataBlock ) const
{
    MStatus status;

    bool value;
    if ( i_dataBlock )
    {
       value = i_dataBlock->inputValue( ia_rodDamping, & status ).asBool();
       CHECK_MSTATUS( status );
    }
    else
    {
       MPlug plug( thisMObject(), ia_rodDamping );
       value = plug.asBool();
    }

    return value;
}
    // Solver settings accessor

void WmSweeneySubsetNode::getSolverSettings(
        double& i_stol, double& i_atol, double& i_rtol,
        double& i_inftol, int& i_numLineSearchIters,
        MDataBlock* i_dataBlock )
{


    MStatus status;

    if ( i_dataBlock )
    {
        i_stol = i_dataBlock->inputValue( ia_stol, & status ).asDouble();
        CHECK_MSTATUS( status );
        i_atol = i_dataBlock->inputValue( ia_atol, & status ).asDouble();
        CHECK_MSTATUS( status );
        i_rtol = i_dataBlock->inputValue( ia_rtol, & status ).asDouble();
        CHECK_MSTATUS( status );
        i_inftol = i_dataBlock->inputValue( ia_inftol, & status ).asDouble();
        CHECK_MSTATUS( status );
        i_numLineSearchIters = i_dataBlock->inputValue( ia_numLineSearchIters, & status ).asInt();
        CHECK_MSTATUS( status );

    }
    else
    {
        MPlug s_plug( thisMObject(), ia_stol );
        i_stol = s_plug.asDouble();
        MPlug a_plug( thisMObject(), ia_atol );
        i_atol = a_plug.asDouble();
        MPlug r_plug( thisMObject(), ia_rtol );
        i_rtol = r_plug.asDouble();
        MPlug inf_plug( thisMObject(), ia_inftol );
        i_inftol = inf_plug.asDouble();
        MPlug nlsi_plug( thisMObject(), ia_numLineSearchIters );
        i_numLineSearchIters = nlsi_plug.asInt();
    }
}

/*static*/ void* WmSweeneySubsetNode::creator()
{
    return new WmSweeneySubsetNode();
}

/*virtual*/ bool WmSweeneySubsetNode::isBounded() const
{
    return false;
}

/*static*/ MStatus WmSweeneySubsetNode::initialize()
{   
    MStatus status;

    //
    // Sweeney node connection.
    //

    {
        MFnNumericAttribute numericAttr;
        oa_toSweeneyParentNode = numericAttr.create( "toSweeneyParentNode", "tsp",
            MFnNumericData::kBoolean, false, & status );
        CHECK_MSTATUS( status );
        numericAttr.setHidden( true );
        addAttribute( oa_toSweeneyParentNode );
    }

    {

        MFnTypedAttribute typedAttrFn;
        ia_scalpFaceIndices = typedAttrFn.create( "scalpFaceIndices", "sfi",
            MFnData::kIntArray, MObject::kNullObj, & status );
        status = addAttribute( ia_scalpFaceIndices );
        CHECK_MSTATUS( status );
    }

    // Rod properties
    {
        MFnNumericAttribute numericAttr;
        ia_length = numericAttr.create("length", "len", MFnNumericData::kDouble, 10.0, &status);
        CHECK_MSTATUS(status);
        CHECK_MSTATUS(numericAttr.setReadable(true));
        CHECK_MSTATUS(numericAttr.setWritable(true));
        CHECK_MSTATUS(numericAttr.setMin(1.0));
        CHECK_MSTATUS(numericAttr.setMax(100.0));
        status = addAttribute(ia_length);
        CHECK_MSTATUS(status);
    }

    {
        MFnNumericAttribute numericAttr;
        ia_rodRadius = numericAttr.create("rodRadius", "ror", MFnNumericData::kDouble, 0.005, &status);
        CHECK_MSTATUS(status);
        CHECK_MSTATUS(numericAttr.setReadable(true));
        CHECK_MSTATUS(numericAttr.setWritable(true));
        CHECK_MSTATUS(numericAttr.setMin(0.001));
        CHECK_MSTATUS(numericAttr.setMax(0.1));
        status = addAttribute(ia_rodRadius);
        CHECK_MSTATUS(status);
    }

    {
        MFnNumericAttribute numericAttr;
        ia_rodAspectRatio = numericAttr.create("rodAspectRatio", "roar", MFnNumericData::kDouble, 1.0, &status);
        CHECK_MSTATUS(status);
        CHECK_MSTATUS(numericAttr.setReadable(true));
        CHECK_MSTATUS(numericAttr.setWritable(true));
        CHECK_MSTATUS(numericAttr.setMin(0.1));
        CHECK_MSTATUS(numericAttr.setMax(10.0));
        status = addAttribute(ia_rodAspectRatio);
        CHECK_MSTATUS(status);
    }

    {
        MFnNumericAttribute numericAttr;
        ia_rodRotation = numericAttr.create("rodRotation", "rorot", MFnNumericData::kDouble, 0.0, &status);
        CHECK_MSTATUS(status);
        CHECK_MSTATUS(numericAttr.setReadable(true));
        CHECK_MSTATUS(numericAttr.setWritable(true));
        CHECK_MSTATUS(numericAttr.setMin(-1.0));
        CHECK_MSTATUS(numericAttr.setMax(1.0));
        status = addAttribute(ia_rodRotation);
        CHECK_MSTATUS(status);
    }

    addNumericAttribute(ia_rodCharge, "rodCharge", "rcg", MFnNumericData::kDouble, 0.0, true);

    addNumericAttribute(ia_rodPower, "rodPower", "rpw", MFnNumericData::kDouble, 1.0, true);

    addNumericAttribute(ia_rodClumpSeparation, "rodClumpSeparation", "rcdst", MFnNumericData::kDouble, 0.001, true);

    {
           MFnNumericAttribute numericAttr;
            ia_curlTightness = numericAttr.create( "globalCurlTightness", "crltight", MFnNumericData::kDouble, 0.0, &status );
            CHECK_MSTATUS( status );
            CHECK_MSTATUS( numericAttr.setReadable( true ) );
            CHECK_MSTATUS( numericAttr.setWritable( true ) );
            CHECK_MSTATUS( numericAttr.setMin( -2.0 ) );
            CHECK_MSTATUS( numericAttr.setMax( 2.0 ) );
            status = addAttribute( ia_curlTightness );
            CHECK_MSTATUS( status );

    }

    {
        MFnNumericAttribute numericAttr;
        ia_curlRadius = numericAttr.create( "curlRadius", "crlrad", MFnNumericData::kDouble, 0.0, &status );
        CHECK_MSTATUS( status );
        CHECK_MSTATUS( numericAttr.setReadable( true ) );
        CHECK_MSTATUS( numericAttr.setWritable( true ) );
        CHECK_MSTATUS( numericAttr.setMin( -5.0 ) );
        CHECK_MSTATUS( numericAttr.setMax( 5.0 ) );
        status = addAttribute( ia_curlRadius );
        CHECK_MSTATUS( status );

    }

    {
        MFnNumericAttribute numericAttr;
        ia_curlCount = numericAttr.create( "curlCount", "crlcnt", MFnNumericData::kDouble, 0.0, &status );
        CHECK_MSTATUS( status );
        CHECK_MSTATUS( numericAttr.setReadable( true ) );
        CHECK_MSTATUS( numericAttr.setWritable( true ) );
        CHECK_MSTATUS( numericAttr.setMin( 0.0 ) );
        CHECK_MSTATUS( numericAttr.setMax( 10.0 ) );
        status = addAttribute( ia_curlCount );
        CHECK_MSTATUS( status );

    }

    {
        MFnNumericAttribute numericAttr;
        ia_curlStart = numericAttr.create( "curlStart", "crlstrt", MFnNumericData::kDouble, 0.0, &status );
        CHECK_MSTATUS( status );
        CHECK_MSTATUS( numericAttr.setReadable( true ) );
        CHECK_MSTATUS( numericAttr.setWritable( true ) );
        CHECK_MSTATUS( numericAttr.setMin( 0.0 ) );
        CHECK_MSTATUS( numericAttr.setMax( 1.0 ) );
        status = addAttribute( ia_curlStart );
        CHECK_MSTATUS( status );
    }


    addNumericAttribute( ia_volumetricCharge, "volumetricForceCharge", "vc", MFnNumericData::kDouble, 0.0, true );

    addNumericAttribute( ia_volumetricScale, "volumetricForceScale", "vs", MFnNumericData::kDouble, 0.0, true );

    addNumericAttribute( ia_drawGaussianVolume, "drawGaussianVolume", "dgv", MFnNumericData::kBoolean, false, true );

    addNumericAttribute( ia_fixCurlCount, "fixCurlCount", "fixcurlc", MFnNumericData::kBoolean, false, true );

    addNumericAttribute( ia_curlInXFrame, "curlInXFrame", "curinx", MFnNumericData::kBoolean, true, true );

    addNumericAttribute( ia_preserveLengthVariation, "preserveLengthVariation", "plenvar", MFnNumericData::kBoolean, true, true );

    addNumericAttribute( ia_rodDamping, "rodDamping", "roddamp", MFnNumericData::kBoolean, true, true );

    addNumericAttribute( ia_rodPitch, "rodPitch", "rop", MFnNumericData::kDouble, 0.5, true );

    //Solver settings
    addNumericAttribute( ia_stol, "stol", "stl", MFnNumericData::kDouble, 99, true );

    addNumericAttribute( ia_atol, "atol", "atl", MFnNumericData::kDouble, 8, true );

    addNumericAttribute( ia_rtol, "rtol", "rtl", MFnNumericData::kDouble, 99, true );

    addNumericAttribute( ia_inftol, "inftol", "itl", MFnNumericData::kDouble, 8, true );

    addNumericAttribute( ia_numLineSearchIters, "numLineSearchIters", "nlsi", MFnNumericData::kInt, 2, true );

    // parameters set at the beginning of the simulation
    addNumericAttribute(ia_verticesPerRod, "verticesPerRod", "cpr", MFnNumericData::kInt, 10, true);

    addNumericAttribute(ia_numberOfClumps, "numberOfClumps", "rpc", MFnNumericData::kInt, 5, true);

    return MS::kSuccess;
}

/*virtual*/ MStatus WmSweeneySubsetNode::compute(
    const MPlug& i_plug, 
    MDataBlock& i_dataBlock )
{
    // compute logic handled in WmSweeneySubsetNode
    return MStatus::kSuccess;
}

/*virtual*/ void WmSweeneySubsetNode::draw(
    M3dView& i_view, 
    const MDagPath& i_path, 
    M3dView::DisplayStyle i_style, 
    M3dView::DisplayStatus i_status )
{
    // Don't draw anything.
}


/*static */MStatus WmSweeneySubsetNode::addNumericAttribute(MObject& i_attribute, MString i_longName, MString i_shortName,
        MFnNumericData::Type i_type, double i_defaultValue, bool i_isInput, bool i_isArray)
{
    // Creates a numeric attribute with default attributes
    MStatus stat = MS::kSuccess;

    MFnNumericAttribute nAttr;
    i_attribute = nAttr.create(i_longName, i_shortName, i_type, i_defaultValue, &stat);
    if (!stat)
    {
        cerr << "Failed to create attribute " << i_longName << endl;
        return stat;
    }
    if (i_isInput)
    {
        nAttr.setWritable(true);
    }
    else
    {
        nAttr.setWritable(false);
    }

    if (i_isArray)
    {
        nAttr.setArray(true);
    }

    stat = addAttribute(i_attribute);
    if (!stat)
    {
        stat.perror("addAttribute " + i_longName);
        return stat;
    }

    return stat;
}
