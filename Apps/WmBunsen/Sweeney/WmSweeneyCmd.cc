#include <inttypes.h>
#include <map>
#include <string>
#include <fstream>
#include <maya/MDagPath.h>
#include <maya/MDagPathArray.h>
#include <maya/MFnDagNode.h>
#include <maya/MItDag.h>
#include <maya/MFnDoubleArrayData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnVectorArrayData.h>
#include <maya/MFnIntArrayData.h>
#include <maya/MFnMeshData.h>
#include <maya/MGlobal.h>
#include <maya/MItMeshPolygon.h>
#include <maya/MItSelectionList.h>
#include <maya/MItDependencyNodes.h>
#include <maya/MPlug.h>
#include <maya/MPlugArray.h>
#include <maya/MItDependencyGraph.h>
#include <maya/MItMeshVertex.h>
#include <maya/MFnPointArrayData.h>
#include <maya/MFnStringArrayData.h>
#include <maya/MFnSingleIndexedComponent.h>
#include <maya/MQuaternion.h>
#include <maya/MEulerRotation.h>
#include <maya/MPxTransformationMatrix.h>
#include <maya/MPlug.h>

#include <values.h>
#include <maya/MProgressWindow.h>
#include <maya/MAnimControl.h>

#include "WmSweeneyCmd.hh"
#include "WmSweeneyNode.hh"
#include "WmSweeneyUtils.hh"
#include "../WmBunsenCollisionMeshNode.hh"

using namespace sweeney::utils;

// Statics
/* static */std::map<std::string, WmSweeneyHelp> WmSweeneyCmd::m_help;
/* static */MString WmSweeneyCmd::typeName( "wmSweeney" );
/* static */MStringArray WmSweeneyCmd::m_results;
//MObject WmSweeneyCmd::m_dynNode;

// CTOR

WmSweeneyCmd::WmSweeneyCmd() :
    m_undoable( false ), m_mArgDatabase( NULL ), m_mDagModifier( NULL ),
            m_selectedSweeneyNode( MObject::kNullObj )
/*,
 m_undo(NULL)*/
{
    m_results.clear();
}

// DTOR

WmSweeneyCmd::~WmSweeneyCmd()
{
    if ( m_mArgDatabase )
    {
        delete m_mArgDatabase;
    }

    if ( m_mDagModifier )
    {
        delete m_mDagModifier;
    }

    //   Wmaya::undo_t::deletelist( m_undo );
}

void* WmSweeneyCmd::creator()
{
    return new WmSweeneyCmd;
}

const char * const kCreateSweeneyNode( "-cs" );
const char * const kAddCollisionMeshes( "-acm" );
const char * const kCreateSweeneySubset( "-css" );
const char * const kSetSimulateAll( "-ssa" );
const char * const kCreateClumpsFromPelt( "-ccp" );

const char * const kHelp( "-h" );

MSyntax WmSweeneyCmd::syntaxCreator()
{
    m_help.clear();

    MSyntax mSyntax;

    p_AddFlag( mSyntax, kHelp, "-help", "Prints this information" );
    p_AddFlag( mSyntax, kCreateSweeneyNode, "-createSweeneyNode",
            "Creates rods from the selected Barbershop furset node." );
    p_AddFlag( mSyntax, kAddCollisionMeshes, "-addCollisionMesh",
            "Adds a a selected polygon mesh as a collision object for the selected Sweeney node." );
    p_AddFlag( mSyntax, kCreateSweeneySubset, "-createSweeneySubset",
            "Creates a new node for a subset of the scalp mesh with its own set of rod parameters" );
    p_AddFlag( mSyntax, kSetSimulateAll, "-simulateAllRods",
            "Applies slider settings only to all rods when simulation begins." );
    p_AddFlag( mSyntax, kCreateClumpsFromPelt, "-createClumpsFromPelt",
            "Creates clumps using centerlines from a wmPelt node." );

    mSyntax.setObjectType( MSyntax::kSelectionList );
    mSyntax.useSelectionAsDefault( true );
    mSyntax.enableEdit( true );

    return mSyntax;
}

MStatus WmSweeneyCmd::doIt( const MArgList &i_mArgList )
{
    MStatus mStatus( MS::kSuccess );
    MArgDatabase *mArgDatabase = new MArgDatabase( syntax(), i_mArgList, &mStatus );

    // Clean Up Old Dag Modifier
    if ( m_mDagModifier )
    {
        delete m_mDagModifier;
        m_mDagModifier = NULL;
    }

    // Clean Out SetAttr Undo Thing...
    //Wmaya::undo_t::deletelist( m_undo );

    if ( mArgDatabase )
    {
        if ( mStatus == MS::kSuccess )
        {
            // We first copy out all the default values
            m_mArgDatabase = mArgDatabase;

            // The command usually operates on whatever we have selected, so filter the selection
            // into the various categories opeations are going to care about
            MGlobal::getActiveSelectionList( m_undoSelectionList );

            MSelectionList opt_nodes;
            m_mArgDatabase->getObjects( opt_nodes );
            MGlobal::getActiveSelectionList( opt_nodes );

            getNodes( opt_nodes );

            // and perform any requested operations
            mStatus = redoIt();
        }
        else
            delete mArgDatabase;
    }

    return mStatus;
}

MStatus WmSweeneyCmd::redoIt()
{
    MStatus mStatus( MS::kSuccess );

    MString opt_fo;

    if ( m_mArgDatabase->isFlagSet( kHelp ) )
    {
        printHelp();
    }
    else
    {
        m_mDagModifier = new MDagModifier;
        m_undoable = true;

        if ( m_mArgDatabase->isFlagSet( kCreateSweeneyNode ) )
        {
            // Create
            createSweeneyNode();
        }
        if ( m_mArgDatabase->isFlagSet( kAddCollisionMeshes ) )
        {
            addCollisionMeshes();
        }
        if ( m_mArgDatabase->isFlagSet( kCreateSweeneySubset ) )
        {
            createSweeneySubsetNode();
        }
        if ( m_mArgDatabase->isFlagSet( kSetSimulateAll ) )
        {
            setSimulateAll();
        }
        if ( m_mArgDatabase->isFlagSet( kCreateClumpsFromPelt ) )
        {
            createClumpCenterLinesFromPelt();
        }

        setResult( m_results );
    }

    return mStatus;
}

MStatus WmSweeneyCmd::undoIt()
{
    return MS::kSuccess;
}

void WmSweeneyCmd::createSweeneyNode()
{
    MStatus stat;

    if ( m_barberShopNodeList.isEmpty() )
    {
        MGlobal::displayError( "Please select a Barbershop furset node to create rods from\n" );
        return;
    }

    // Create the Sweeney rods node
    MObject rodTObj; // Object for transform node
    MObject rodSObj; // Object for shape node
    MDagPath shapeDagPath;
    MObject pObj;
    MDagModifier dagModifier;
    MString rodShapeName = "";

    createDagNode( WmSweeneyNode::typeName.asChar(), WmSweeneyNode::typeName.asChar(), pObj,
            &rodTObj, &rodSObj, &dagModifier, rodShapeName );

    appendToResultString( rodShapeName );

    MDagPath rodDagPath;
    stat = MDagPath::getAPathTo( rodSObj, rodDagPath );
    CHECK_MSTATUS( stat );

    // Turn Off Inherits Transform On Transform Node as the rod node should never move from the origin
    MFnDependencyNode tFn( rodTObj );
    MPlug mPlug( tFn.findPlug( "inheritsTransform", true, &stat ) );
    CHECK_MSTATUS( stat );
    mPlug.setValue( false );

    // Lock the transform attributes because the object should not be moved
    // Transform first
    MPlug translatePlug( tFn.findPlug( "translate", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = translatePlug.setLocked( true );
    CHECK_MSTATUS( stat );
    MPlug scalePlug( tFn.findPlug( "scale", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = scalePlug.setLocked( true );
    CHECK_MSTATUS( stat );
    MPlug rotatePlug( tFn.findPlug( "rotate", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = rotatePlug.setLocked( true );
    CHECK_MSTATUS( stat );
    MPlug shearPlug( tFn.findPlug( "shear", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = shearPlug.setLocked( true );
    CHECK_MSTATUS( stat );

    // Shape now MFnDependencyNode sFn( rodSObj );
    MFnDependencyNode sFn( rodSObj );

    MPlug localTransPlug( sFn.findPlug( "localPosition", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = localTransPlug.setLocked( true );
    CHECK_MSTATUS( stat );
    MPlug localScalePlug( sFn.findPlug( "localScale", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = localScalePlug.setLocked( true );
    CHECK_MSTATUS( stat );

    MString timeStr( "connectAttr -f time1.outTime " + rodDagPath.fullPathName() + ".time" );
    dagModifier.commandToExecute( timeStr );
    dagModifier.doIt();

    cerr << "WmSweeneyCmd::createSweeneyNode - Connecting fur set\n";

    // Just take the first furset selected, selecting multiple sets is pointless anyway
    MDagPath dagPath;
    MObject component;
    m_barberShopNodeList.getDagPath( 0, dagPath, component );
    dagPath.extendToShape();

    MObject nodeObj = dagPath.node( &stat );
    CHECK_MSTATUS( stat );
    MFnDependencyNode barberShopNodeFn( nodeObj, &stat );

    MPlug barberShopStrandsPlug = barberShopNodeFn.findPlug( "strandVertices", true, &stat );
    CHECK_MSTATUS( stat );
    MPlug sweeneyVerticesPlug( sFn.findPlug( "strandVertices", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( barberShopStrandsPlug, sweeneyVerticesPlug );
    CHECK_MSTATUS( stat );

    // todo(sainsley) : grab root strand frame here "strandRootFrames" ?
    // Look at WMFigaroCmd
    MPlug barberShopRootFramesPlug = barberShopNodeFn.findPlug( "strandRootFrames", true, &stat );
    CHECK_MSTATUS( stat );
    MPlug sweeneyRootFramePlug( sFn.findPlug( "strandRootFrames", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( barberShopRootFramesPlug, sweeneyRootFramePlug );
    CHECK_MSTATUS( stat );

    MPlug numBarberShopCVsPlug = barberShopNodeFn.findPlug( "outCvsPerStrand", true, &stat );
    CHECK_MSTATUS( stat );
    MPlug numCVsPlug( sFn.findPlug( "verticesPerStrand", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( numBarberShopCVsPlug, numCVsPlug );
    CHECK_MSTATUS( stat );

    // set up clump vertex-power map
    stat = MGlobal::executeCommand( "select -r " + rodDagPath.fullPathName() + ";" );
    CHECK_MSTATUS( stat );

    // They want the curve to default like this. 0@0, 1@1, and also 1@0.1
    stat = MGlobal::executeCommand( MString( "setAttr -s 3 \".rodClumpingRamp[1:3]\"  " ) +
        "1.0 1.0 1.0 " +
        "0.1 1.0 1.0 " +
        "0.0 0.0 1.0 " +
        ";" );
    CHECK_MSTATUS( stat );

    stat = dagModifier.doIt();
    CHECK_MSTATUS( stat );
}

void WmSweeneyCmd::createSweeneySubsetNode()
{
    MStatus stat;

    // check that a sweeney node and a list of faces are selected
    if ( m_selectedSweeneyNode == MObject::kNullObj )
    {
        MGlobal::displayError( "Please select a wmSweeney node associated with these mesh faces." );
        return;
    }

    MIntArray faces;
    MString objectName;
    getSelectedComponents( kFaceComponent, objectName, faces );

    // check that at least one scalp face is selected
    if ( faces.length() < 1 )
    {
        MGlobal::displayError( "Please select a subset of faces from the BarberShop mesh." );
        return;
    }

    // Create the Sweeney subset node
    MObject rodSubsetTransformObj; // Object for transform node
    MObject rodSubsetShapeObj; // Object for shape node
    MDagPath shapeDagPath;
    MObject pObj;
    MDagModifier dagModifier;
    MString rodSubsetShapeName = "";

    createDagNode( WmSweeneySubsetNode::typeName.asChar(), WmSweeneySubsetNode::typeName.asChar(),
            pObj, &rodSubsetTransformObj, &rodSubsetShapeObj, &dagModifier, rodSubsetShapeName );

    appendToResultString( rodSubsetShapeName );

    MDagPath rodSubsetDagPath;
    stat = MDagPath::getAPathTo( rodSubsetShapeObj, rodSubsetDagPath );
    CHECK_MSTATUS( stat );

    // Turn Off Inherits Transform On Transform Node as the rod node should never move from the origin
    MFnDependencyNode tFn( rodSubsetTransformObj );
    MPlug mPlug( tFn.findPlug( "inheritsTransform", true, &stat ) );
    CHECK_MSTATUS( stat );
    mPlug.setValue( false );

    // Lock the transform attributes because the object should not be moved
    // Transform first
    MPlug translatePlug( tFn.findPlug( "translate", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = translatePlug.setLocked( true );
    CHECK_MSTATUS( stat );
    MPlug scalePlug( tFn.findPlug( "scale", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = scalePlug.setLocked( true );
    CHECK_MSTATUS( stat );
    MPlug rotatePlug( tFn.findPlug( "rotate", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = rotatePlug.setLocked( true );
    CHECK_MSTATUS( stat );
    MPlug shearPlug( tFn.findPlug( "shear", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = shearPlug.setLocked( true );
    CHECK_MSTATUS( stat );

    // Shape now MFnDependencyNode sFn( rodSObj );
    MFnDependencyNode sFn( rodSubsetShapeObj );

    MPlug localTransPlug( sFn.findPlug( "localPosition", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = localTransPlug.setLocked( true );
    CHECK_MSTATUS( stat );
    MPlug localScalePlug( sFn.findPlug( "localScale", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = localScalePlug.setLocked( true );
    CHECK_MSTATUS( stat );

    dagModifier.doIt();

    // Connect Sweeney node to subset node
    stat = dagModifier.connect( rodSubsetShapeObj, WmSweeneySubsetNode::oa_toSweeneyParentNode,
            m_selectedSweeneyNode, WmSweeneyNode::ia_fromSweeneySubsetNodes );
    CHECK_MSTATUS( stat );

    // Pass face indicies to subset node
    MFnDagNode dagNodeFn( rodSubsetShapeObj );

    WmSweeneySubsetNode* sweeneySubsetNode = ( WmSweeneySubsetNode* ) dagNodeFn.userNode();
    sweeneySubsetNode->setScalpFaceIndices( faces );

    stat = dagModifier.doIt();
    CHECK_MSTATUS( stat );
}

MStatus WmSweeneyCmd::createClumpCenterLinesFromPelt()
{
    // First let's check that the right nodes have been previously selected.
    if ( m_selectedSweeneyNode == MObject::kNullObj || m_selectedPeltNode == MObject::kNullObj )
    {
        MGlobal::displayError( "Please select a wmSweeney node and a wmPeltnode" );
        return MStatus::kFailure;
    }
    // TODO: allow the user to select the mesh directly if it already exists.

    // Get the mesh from the peltNode.
    MObject peltMesh;
    if ( MFnDependencyNode( m_selectedPeltNode ).findPlug( "meshOut" ).getValue( peltMesh )
            == MStatus::kFailure )
    {
        MGlobal::displayError(
                "The selected wmPeltnode doesn't have an output mesh. It should have been created with -mo" );
        return MStatus::kFailure;
    }
    // TODO: create the mesh if a mesh-less peltNode has been selected.

    // Extract the mesh centres.
    MItMeshPolygon polyIt( peltMesh );
    MPointArray centralArr;
    for ( MItMeshPolygon polyIt( peltMesh ); !polyIt.isDone(); polyIt.next() )
    {
        MPoint center = polyIt.center( MSpace::kWorld );
        std::cout << "Found centre point: " << center << '\n';
        centralArr.append( center );
    }

    MFnDependencyNode sweeneyFn( m_selectedSweeneyNode );
    WmSweeneyNode* sweeneyNode = ( WmSweeneyNode* ) sweeneyFn.userNode();
    sweeneyNode->createClumpCenterLinesFromPelt(centralArr);

    return MStatus::kSuccess;
}

void WmSweeneyCmd::appendToResultString( MString& i_resultString )
{
    MStatus stat;

    // FIXME:
    // I can't get appendToResult to work so I just store all the results and set them in one
    // go at the end of redoit

    //  MString s = currentStringResult( &stat );
    // cerr << "s = " << s << endl;

    /*  if ( !stat )
     {
     cerr << "setting!\n";
     setResult( i_resultString );
     }
     else*/
    {
        m_results.append( i_resultString );
    }

    //    MString s = currentStringResult( &stat );
    //  cerr << " after s = " << s << endl;
}

void WmSweeneyCmd::getNodes( MSelectionList i_opt_nodes )
{
    // Prepare all the selection lists that this command operates on

    MDagPath mDagPath;
    MObject mObj;
    MItDag dagit;
    MStatus stat;

    m_barberShopNodeList.clear();
    m_meshList.clear();
    m_sweeneyNodeList.clear();
    m_allOtherTransformNodesList.clear();

    for ( MItSelectionList sIt( i_opt_nodes ); !sIt.isDone(); sIt.next() )
    {
        sIt.getDagPath( mDagPath, mObj );

        cerr << "WmSweeneyCmd::getNodes() - selected Dag path = " << mDagPath.fullPathName()
                << endl;

        MFnDagNode dagFn( mDagPath.child( 0, &stat ), &stat );
        CHECK_MSTATUS( stat );
        MDagPath childPath;
        stat = dagFn.getPath( childPath );
        CHECK_MSTATUS( stat );
        childPath.extendToShape();

        //////////////////
        //
        // First check for Maya in built nodes

        if ( childPath.apiType() == MFn::kMesh )
        {
            mObj = childPath.node();
            m_meshList.add( childPath, mObj, false );
        }
        else
        {
            /////////////////////
            //
            // Now check for user created plugin nodes

            MFnDependencyNode nodeFn( childPath.node( &stat ) );
            CHECK_MSTATUS( stat );

            if ( nodeFn.typeName() == "wmBarbFurSetNode" )
            {
                mObj = childPath.node();
                m_barberShopNodeList.add( childPath, mObj, false );
            }
            else if ( nodeFn.typeName() == WmSweeneyNode::typeName )
            {
                mObj = childPath.node();
                m_selectedSweeneyNode = mObj;

                m_sweeneyNodeList.add( childPath, mObj, false );
            }
            else if ( nodeFn.typeName() == "wmPelt" )
            {
                mObj = childPath.node();
                m_selectedPeltNode = mObj;
                m_peltNodeList.add( childPath, mObj, false );
            }
            else
            {
                stat = m_allOtherTransformNodesList.add( mDagPath, mObj, false );
                CHECK_MSTATUS( stat );
            }
        }
    }
}

void WmSweeneyCmd::addCollisionMeshes()
{
    MStatus stat;

    if ( m_selectedSweeneyNode == MObject::kNullObj )
    {
        MGlobal::displayError( "Please select a wmSweeney node to connect the mesh to." );
        return;
    }

    MFnDependencyNode sweeneyNodeFn( m_selectedSweeneyNode, &stat );
    CHECK_MSTATUS( stat );

    MPlug sweeneyInputPlugArr = sweeneyNodeFn.findPlug( "collisionMeshes", true, &stat );
    CHECK_MSTATUS( stat );

    for ( unsigned int l = 0; l < m_meshList.length(); l++ )
    {
        MDagPath dagPath;
        MObject component;
        m_meshList.getDagPath( l, dagPath, component );
        dagPath.extendToShape();

        MObject nodeObj = dagPath.node( &stat );
        CHECK_MSTATUS( stat );
        MFnDependencyNode nodeFn( nodeObj, &stat );

        // Create a collisionMeshNode to pipe the mesh through on the way to the dynamics node

        MObject collisionMeshNodeTObj; // Object for transform node
        MObject collisionMeshNodeSObj; // Object for shape node
        MDagPath shapeDagPath;
        MObject pObj;
        MDagModifier dagModifier;
        MString collisionNodeShapeName = "";

        createDagNode( WmBunsenCollisionMeshNode::typeName.asChar(),
                WmBunsenCollisionMeshNode::typeName.asChar(), pObj, &collisionMeshNodeTObj,
                &collisionMeshNodeSObj, &dagModifier, collisionNodeShapeName );

        appendToResultString( collisionNodeShapeName );

        MDagPath collisionMeshNodeDagPath;
        stat = MDagPath::getAPathTo( collisionMeshNodeSObj, collisionMeshNodeDagPath );
        CHECK_MSTATUS( stat );

        MString
                timeStr(
                        "connectAttr -f time1.outTime " + collisionMeshNodeDagPath.fullPathName()
                                + ".time" );
        dagModifier.commandToExecute( timeStr );
        dagModifier.doIt();

        if ( nodeFn.typeName() == "mesh" )
        {
            MFnMesh meshFn( dagPath, &stat );
            CHECK_MSTATUS( stat );

            MPlug worldMeshPlug =
                    meshFn.findPlug( "worldMesh", true, &stat ).elementByLogicalIndex( 0, &stat );
            CHECK_MSTATUS( stat );

            MPlug collisionMeshNodeInPlug( collisionMeshNodeSObj,
                    WmBunsenCollisionMeshNode::ia_inMesh );

            stat = dagModifier.connect( worldMeshPlug, collisionMeshNodeInPlug );
            CHECK_MSTATUS( stat );
            stat = dagModifier.doIt();
            CHECK_MSTATUS( stat );

            // We need to track when the mesh transforms so we can do fast level set look ups
            // if the mesh is only rigidly moving.
            MPlug worldMatrixPlug =
                    meshFn.findPlug( "worldMatrix", true, &stat ).elementByLogicalIndex( 0, &stat );
            CHECK_MSTATUS( stat );

            MPlug collisionMeshNodeMatrixPlug( collisionMeshNodeSObj,
                    WmBunsenCollisionMeshNode::ia_meshTransform );

            stat = dagModifier.connect( worldMatrixPlug, collisionMeshNodeMatrixPlug );
            CHECK_MSTATUS( stat );
            stat = dagModifier.doIt();
            CHECK_MSTATUS( stat );

            MPlug collisionMeshNodeOutPlug( collisionMeshNodeSObj,
                    WmBunsenCollisionMeshNode::oa_meshData );

            unsigned int numElements = sweeneyInputPlugArr.numElements( &stat );
            CHECK_MSTATUS( stat );
            MPlug sweeneyMeshPlug = sweeneyInputPlugArr.elementByLogicalIndex( numElements );
            CHECK_MSTATUS( stat );
            if ( stat.error() )
            {
                MGlobal::displayError( stat.errorString() );
                return;
            }
            stat = dagModifier.connect( collisionMeshNodeOutPlug, sweeneyMeshPlug );
            CHECK_MSTATUS( stat );
            stat = dagModifier.doIt();
            CHECK_MSTATUS( stat );

            // Connect startTime so the rods reset when the Figaro node resets
            MPlug collisionNodeStartTimePlug( collisionMeshNodeSObj,
                    WmBunsenCollisionMeshNode::ia_startTime );
            CHECK_MSTATUS( stat );
            MPlug bunsenNodeStartTimePlug = sweeneyNodeFn.findPlug( "startTime", true, &stat );
            CHECK_MSTATUS( stat );
            stat = dagModifier.connect( bunsenNodeStartTimePlug, collisionNodeStartTimePlug );
            CHECK_MSTATUS( stat );
            stat = dagModifier.doIt();
            CHECK_MSTATUS( stat );
        }
    }
}

/*void WmSweeneyCmd::setSimulatedSubset( )
 {
 if ( m_selectedSweeneyNode == MObject::kNullObj )
 {
 MGlobal::displayError( "Please select a wmSweeney node associated with these mesh faces." );
 return;
 }

 MIntArray faces;
 MString objectName;
 getSelectedComponents( kFaceComponent, objectName, faces );

 MFnDagNode dagNodeFn( m_selectedSweeneyNode );

 WmSweeneyNode* sweeneyNode = (WmSweeneyNode*) dagNodeFn.userNode();

 sweeneyNode->setScalpSelection( faces );
 // now pass this list to the sweeney node

 //cout << " SELECTED FACES: " << endl;
 //for ( int i = 0; i < faces.length(); i++ )
 //{
 //  cout << faces[i] << endl;
 //}

 //std::map<int,int> scalpRodMap = getScalpRodMap();
 }*/

void WmSweeneyCmd::setSimulateAll()
{
    // tell rod manager to simulate all
}

void WmSweeneyCmd::p_AddFlag( MSyntax &i_mSyntax, const char * const i_shortName,
        const char * const i_longName, const char * const i_help,
        const MSyntax::MArgType i_argType1, const MSyntax::MArgType i_argType2,
        const MSyntax::MArgType i_argType3, const MSyntax::MArgType i_argType4,
        const MSyntax::MArgType i_argType5, const MSyntax::MArgType i_argType6 )
{
    CHECK_MSTATUS(
            i_mSyntax.addFlag( i_shortName, i_longName, i_argType1, i_argType2, i_argType3,
                    i_argType4, i_argType5, i_argType6 ) );
    WmSweeneyHelp &furHelp = m_help[i_shortName];
    furHelp.m_longName = i_longName;
    furHelp.m_help = i_help;

    if ( i_argType1 != MSyntax::kInvalidArgType && i_argType1 != MSyntax::kNoArg )
    {
        if ( i_argType1 == MSyntax::kLastArgType && furHelp.m_argTypes.size() )
        {
            furHelp.m_argTypes.push_back( furHelp.m_argTypes.back() );
        }
        else
        {
            furHelp.m_argTypes.push_back( i_argType1 );
        }
        if ( i_argType2 != MSyntax::kInvalidArgType && i_argType2 != MSyntax::kNoArg )
        {
            if ( i_argType2 == MSyntax::kLastArgType && furHelp.m_argTypes.size() )
            {
                furHelp.m_argTypes.push_back( furHelp.m_argTypes.back() );
            }
            else
            {
                furHelp.m_argTypes.push_back( i_argType2 );
            }
            if ( i_argType3 != MSyntax::kInvalidArgType && i_argType3 != MSyntax::kNoArg )
            {
                if ( i_argType3 == MSyntax::kLastArgType && furHelp.m_argTypes.size() )
                {
                    furHelp.m_argTypes.push_back( furHelp.m_argTypes.back() );
                }
                else
                {
                    furHelp.m_argTypes.push_back( i_argType3 );
                }
                if ( i_argType4 != MSyntax::kInvalidArgType && i_argType4 != MSyntax::kNoArg )
                {
                    if ( i_argType4 == MSyntax::kLastArgType && furHelp.m_argTypes.size() )
                    {
                        furHelp.m_argTypes.push_back( furHelp.m_argTypes.back() );
                    }
                    else
                    {
                        furHelp.m_argTypes.push_back( i_argType4 );
                    }
                    if ( i_argType5 != MSyntax::kInvalidArgType && i_argType5 != MSyntax::kNoArg )
                    {
                        if ( i_argType5 == MSyntax::kLastArgType && furHelp.m_argTypes.size() )
                        {
                            furHelp.m_argTypes.push_back( furHelp.m_argTypes.back() );
                        }
                        else
                        {
                            furHelp.m_argTypes.push_back( i_argType5 );
                        }
                        if ( i_argType6 != MSyntax::kInvalidArgType && i_argType6
                                != MSyntax::kNoArg )
                        {
                            if ( i_argType6 == MSyntax::kLastArgType && furHelp.m_argTypes.size() )
                            {
                                furHelp.m_argTypes.push_back( furHelp.m_argTypes.back() );
                            }
                            else
                            {
                                furHelp.m_argTypes.push_back( i_argType6 );
                            }
                        }
                    }
                }
            }
        }
    }
}

void WmSweeneyCmd::printHelp()
{
    /*  minfo <<
     "\n"
     "wmBunsen:"
     "\n"
     "DESCRIPTION:\n"
     "\n"
     "  Creates and edits WmBunsen nodes\n"
     "\n"
     "SWITCHES:\n" << "\n";

     // First Fing The Longest Short Name & Long Name

     const std::map<std::string, WmBunsenHelp>::const_iterator hEnd(m_help.end());
     std::map<std::string, WmBunsenHelp>::const_iterator hi;

     uint32_t sNameLen(0);

     for (hi = m_help.begin(); hi != hEnd; hi++) {
     if (hi->first.size() > sNameLen) {
     sNameLen = hi->first.size();
     }
     }

     sNameLen++;

     for (hi = m_help.begin(); hi != hEnd; hi++) {
     const uint32_t sLen(hi->first.size());
     for (uint32_t si = sLen; si < sNameLen; si++) {
     minfo << " ";
     }
     minfo << hi->first << " " << hi->second.m_longName;

     const std::vector<MSyntax::MArgType>::const_iterator aEnd(hi->second.m_argTypes.end());
     std::vector<MSyntax::MArgType>::const_iterator ai;

     for (ai = hi->second.m_argTypes.begin(); ai != aEnd; ai++) {
     minfo << " <";
     switch (*ai) {
     case MSyntax::kBoolean:
     minfo << "Boolean";
     break;
     case MSyntax::kLong:
     minfo << "Long";
     break;
     case MSyntax::kDouble:
     minfo << "Double";
     break;
     case MSyntax::kString:
     minfo << "String";
     break;
     case MSyntax::kUnsigned:
     minfo << "Unsigned";
     break;
     case MSyntax::kDistance:
     minfo << "Distance";
     break;
     case MSyntax::kAngle:
     minfo << "Angle";
     break;
     case MSyntax::kTime:
     minfo << "Time";
     break;
     case MSyntax::kSelectionItem:
     minfo << "SelectionItem";
     break;
     default:
     minfo << "Unknown";
     break;
     }
     minfo << ">";
     }

     minfo << "\n";
     for (uint32_t si = 0; si < sNameLen; si++) {
     minfo << " ";
     }
     minfo << " * " << hi->second.m_help << "\n";
     }

     minfo << endl;*/
}

MStatus WmSweeneyCmd::createDagNode( const char *transformName, const char *nodeType,
        MObject &parentObj, MObject *transformObjP, MObject *shapeObjP, MDagModifier *iDagModifier,
        MString& o_shapeName )
{
    MStatus rStatus = MS::kSuccess;

    MDagModifier tDagModifier;
    MDagModifier &mDagModifier = ( iDagModifier != NULL ) ? *iDagModifier : tDagModifier;

    MObject transformObj = mDagModifier.createNode( nodeType, MObject::kNullObj, &rStatus );

    if ( rStatus )
    {
        mDagModifier.doIt();

        if ( !transformObj.isNull() )
        {
            if ( transformObjP )
                *transformObjP = transformObj;

            if ( strcmp( transformName, nodeType ) )
            {
                rStatus = mDagModifier.renameNode( transformObj, transformName );
            }
            else
            {
                MDagPath dagPath;
                MFnDagNode( transformObj ).getPath( dagPath );
                dagPath.extendToShape();
                rStatus = mDagModifier.renameNode( transformObj,
                        MFnDependencyNode( dagPath.node() ).name() );
            }

            if ( rStatus )
            {
                mDagModifier.doIt();

                MDagPath shapeDagPath = MDagPath::getAPathTo( transformObj );
                const unsigned childCount = shapeDagPath.childCount();

                if ( childCount == 1 )
                {
                    shapeDagPath.extendToShape();
                    MObject shapeObj = shapeDagPath.node();

                    if ( shapeObjP )
                        *shapeObjP = shapeObj;

                    MFnTransform transformFn( transformObj );

                    const MString transformNameStr = transformFn.name();
                    const char *sp = transformNameStr.asChar();
                    const char *ep = sp + transformNameStr.length() - 1;
                    const char *cp = ep;

                    while ( cp >= sp && isdigit( *cp ) )
                        cp--;

                    MString shapeName;
                    if ( cp != sp )
                        shapeName = transformNameStr.substring( 0, int( cp - sp ) );
                    else
                        shapeName = transformNameStr;

                    shapeName += "Shape";

                    if ( cp != ep )
                        shapeName
                                += transformNameStr.substring( int( cp - sp + 1 ), int( ep - sp ) );

                    rStatus = mDagModifier.renameNode( shapeObj, shapeName );
                    if ( rStatus )
                    {
                        mDagModifier.doIt();
                        o_shapeName = shapeName;
                    }
                    else
                    {
                        cerr << "Couldn't Rename Shape To: \"" << shapeName << "\"" << endl;
                        o_shapeName = "";
                    }
                }
                else
                {
                    if ( shapeObjP )
                        *shapeObjP = MObject::kNullObj;

                    if ( childCount )
                        cerr << "Unknown Child Count == " << childCount << endl;
                }
            }
            else
            {
                cerr << "Couldn't Rename Transform Node To \"" << transformName << "\"" << endl;
            }

            if ( parentObj != MObject::kNullObj )
            {
                rStatus = mDagModifier.reparentNode( transformObj, parentObj );
                if ( rStatus )
                    mDagModifier.doIt();
                else
                    cerr << "Couldn't Reparent Node" << endl;
            }
        }
        else
        {
            cerr << "Couldn't Create Node Of \"" << nodeType << "\"" << endl;
            rStatus = MS::kFailure;
        }
    }
    else
    {
        cerr << "Couldn't Create Transform Node \"" << transformName << "\"" << endl;
    }

    return rStatus;
}
