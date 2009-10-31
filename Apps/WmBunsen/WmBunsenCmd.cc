#include "WmBunsenCmd.hh"

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

#include <values.h>
#include <maya/MProgressWindow.h>
#include <maya/MAnimControl.h>

// Statics
std::map<std::string, WmBunsenHelp> WmBunsenCmd::m_help;
MString WmBunsenCmd::typeName("wmBunsen");
//MObject WmBunsenCmd::m_dynNode;

// CTOR

WmBunsenCmd::WmBunsenCmd()
  : m_undoable(false),
    m_mArgDatabase(NULL),
    m_mDagModifier(NULL)/*,
    m_undo(NULL)*/
{
}

// DTOR

WmBunsenCmd::~WmBunsenCmd()
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

void* WmBunsenCmd::creator()
{
    return new WmBunsenCmd;
}

const char *const kCreateRods( "-cr" );
const char *const kCVsPerRod( "-cpr" );
const char *const kName( "-n" );
const char *const kHelp( "-h" );

MSyntax WmBunsenCmd::syntaxCreator()
{
    m_help.clear();

    MSyntax mSyntax;

    p_AddFlag( mSyntax, kHelp, "-help",
               "Prints this information" );
    p_AddFlag( mSyntax, kName, "-name",
               "Sets the name of the created WmBunsen node",
               MSyntax::kString );
    p_AddFlag( mSyntax, kCreateRods, "-createRods",
               "creates discrete elastic rods from the selected NURBS curves." );
    p_AddFlag( mSyntax, kCVsPerRod, "-cvsPerRod",
               "Sets the number of CVs to use for each created elastic rod.", MSyntax::kLong );
    
    mSyntax.setObjectType( MSyntax::kSelectionList );
    mSyntax.useSelectionAsDefault( true );
    mSyntax.enableEdit( true );

    return mSyntax;
}

MStatus WmBunsenCmd::doIt( const MArgList &i_mArgList )
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

            m_node_name = MString("wmBunsen");

            m_cvsPerRod = -1; // Use the number of cvs from the first curve unless the user specifies a number of cvs

            if ( m_mArgDatabase->isFlagSet( kName ) )
                m_mArgDatabase->getFlagArgument( kName, 0, m_node_name );
            if ( m_mArgDatabase->isFlagSet( kCreateRods ) )
                m_mArgDatabase->getFlagArgument( kCVsPerRod, 0, m_cvsPerRod );

            // The command usually operates on whatever we have selected, so filter the selection
            // into the various categories opeations are going to care about
            MGlobal::getActiveSelectionList( m_undoSelectionList );

            MSelectionList opt_nodes;
            m_mArgDatabase->getObjects( opt_nodes );
            MGlobal::getActiveSelectionList( opt_nodes );

            getNodes( opt_nodes );

            // Record the operation we're supposed to be performing, and perform setup
            /*if (m_mArgDatabase->isFlagSet( kHelp ) )
                m_op = Help;
            if (m_mArgDatabase->isFlagSet( kCreateRods ) )
                m_op = CreateRods;*/
            
            // and perform any requested operations

            mStatus = redoIt();
        }
        else
            delete mArgDatabase;
    }

    return mStatus;
}



void WmBunsenCmd::createWmBunsenRodNode( bool useNURBSInput )
{
    MStatus stat;
    size_t nCurves = 0;
    
    MObject wmBunsenNode;    

    createWmBunsenNode( wmBunsenNode );

    // Create the rods node
    MObject rodTObj;  // Object for transform node
    MObject rodSObj;  // Object for shape node
    MDagPath shapeDagPath;    
    MObject pObj;
    MDagModifier dagModifier;

    createDagNode( WmBunsenRodNode::ia_typeName.asChar(), WmBunsenRodNode::ia_typeName.asChar(), 
                   pObj, &rodTObj, &rodSObj, &dagModifier );
                          
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
    
    // Shape now    MFnDependencyNode sFn( rodSObj );
    MFnDependencyNode sFn( rodSObj );    
    MPlug localTransPlug( sFn.findPlug( "localPosition", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = localTransPlug.setLocked( true );
    CHECK_MSTATUS( stat );
    MPlug localScalePlug(sFn.findPlug( "localScale", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = localScalePlug.setLocked( true );
    CHECK_MSTATUS( stat );
        
    MString timeStr( "connectAttr -f time1.outTime " + rodDagPath.fullPathName() + ".time" );
    dagModifier.commandToExecute( timeStr );
    dagModifier.doIt();
    
    if ( useNURBSInput )
    {
        // Set the number of cvs for created rods, either -1 or a number the user has specified on the command line
      /*  cerr << "Setting cvs per rod to " << m_cvsPerRod << endl;
        MPlug cvsPlug( sFn.findPlug( "cvsPerRod", true, &stat ) );
        CHECK_MSTATUS( stat );
        cvsPlug.setValue( m_cvsPerRod );*/
        
        MPlug rodNodeInputPlugArr = sFn.findPlug( "nurbsCurves", &stat );
        CHECK_MSTATUS( stat );
        
        // if m_cvsPerRod = -1 then the user has not specified a number of cvs per rod so just take
        // the number of cvs from the first rod. If the user did specficy a number of cvs per rod
        // then all curves will be resampled to this number.
        size_t cvsPerCurve = 0;
        if ( m_cvsPerRod != -1 )
            cvsPerCurve = m_cvsPerRod;
    
        // Check if there are any NURBS curves connected as we can't create rods without curves    
        for ( uint l = 0; l < m_curve_list.length(); l++ ) 
        {
            MDagPath dagPath;
            MObject component;
            m_curve_list.getDagPath( l, dagPath, component );
            dagPath.extendToShape();
    
            MObject nodeObj = dagPath.node( &stat );        
            CHECK_MSTATUS( stat );
            MFnDependencyNode nodeFn( nodeObj,&stat );
            // Check if this is a NURBS curve
            if ( nodeFn.typeName() == "nurbsCurve" )
            {
                nCurves++;            
                MFnNurbsCurve curveFn( dagPath, &stat );
                CHECK_MSTATUS( stat );
                
                size_t nCVs = curveFn.numCVs( &stat );
                CHECK_MSTATUS( stat );
                
                // If we have not stored any curves then store the number of cvs from the 1st curve
                // and check that all curves have that many.
                if ( cvsPerCurve == 0 )
                    cvsPerCurve = nCVs;
                else 
                {
                    if ( nCVs != cvsPerCurve && m_cvsPerRod == -1 ) 
                    {
                        // This curve does not have the same number of cvs as the first curve and 
                        // the user has not specified reasmpling of curves so we must skip it. We
                        // can only handle rods with the same number of cvs.
    
                        MGlobal::displayWarning( "Selected curve does not have the same number of CVs as "
                                                 "first curve, skipping curve. Use dynoCmd -cr -cpr <num cvs> "
                                                 "to force curve resampling for curves with different numbers of cvs" );
                        continue;
                    }
                }                                    
                    
                MPlug curveShapePlug = curveFn.findPlug( "worldSpace", true, &stat ).elementByLogicalIndex(0, &stat );
                CHECK_MSTATUS( stat );
                
                MPlug rodPlug = rodNodeInputPlugArr.elementByLogicalIndex( rodNodeInputPlugArr.numElements(), &stat );
                CHECK_MSTATUS( stat );
                stat = dagModifier.connect( curveShapePlug, rodPlug );
                CHECK_MSTATUS( stat );
                stat = dagModifier.doIt();
                CHECK_MSTATUS( stat );
                
                // Experimentally create output NURBS curves and connect them up. This should be done
                // with the modifier.
                
                /*MDagPath curveTransformPath;
                MString cmd;
                m_curve_list.getDagPath( l, curveTransformPath, component );
                MString result;
                MStringArray results;
                cmd = "duplicate " +  curveTransformPath.fullPathName();
                cerr << cmd << endl;
                MGlobal::executeCommand( cmd, results );
                cerr << "result = " << results[0];
                cmd = "listRelatives -shapes " +  results[0];
                cerr << cmd << endl;
                MGlobal::executeCommand( cmd, results );
                cmd = "connectAttr " + rodDagPath.fullPathName() + ".NURBSCurves[" + l + "] "+ results[0] + ".create";
                cerr << cmd << endl;
                MGlobal::executeCommand( cmd );  */
            }
        }
    }

    // Connect up the rods node to the dynamics node
    MPlug rodNodeOutPlug( sFn.findPlug( "rodsChanged", true, &stat ) );
    CHECK_MSTATUS( stat );
    MFnDependencyNode bunsenDependNodeFn( wmBunsenNode );
    MPlug bunsenNodeRodArrayPlug( bunsenDependNodeFn.findPlug( "rodsNodes", true, &stat ) ); 
    CHECK_MSTATUS( stat );
    MPlug bunsenNodeRodPlug = bunsenNodeRodArrayPlug.elementByLogicalIndex( bunsenNodeRodArrayPlug.numElements(), &stat );
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( rodNodeOutPlug, bunsenNodeRodPlug );
    CHECK_MSTATUS( stat );
    
    stat = dagModifier.doIt();
    CHECK_MSTATUS( stat );
            
    //m_dynNode = rodSObj;

    // If we had no well formed curves to use for shaping then just return
    /*if ( nCurves == 0 ) {
        MGlobal::displayError( "Did not find any curves to create rods from. " );                               
        return;
    }*/
}

void WmBunsenCmd::createWmBunsenNode( MObject &o_wmBunsenNodeObj ) 
{    
 
    MStatus stat;
    
    MObject bunsenNodeTObj;  // Object for transform node
    MObject bunsenNodeSObj;  // Object for shape node
    MDagPath shapeDagPath;    
    MObject pObj;
    MDagModifier dagModifier;

    createDagNode( WmBunsenNode::ia_typeName.asChar(), WmBunsenNode::ia_typeName.asChar(), 
                   pObj, &bunsenNodeTObj, &bunsenNodeSObj, &dagModifier );
                          
    MDagPath bunsenNodeDagPath;
    stat = MDagPath::getAPathTo( bunsenNodeSObj, bunsenNodeDagPath );
    CHECK_MSTATUS( stat );
    
     // Turn Off Inherits Transform On Transform Node as the fur should never move from the origin
    MFnDependencyNode tFn( bunsenNodeTObj );
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
    
    // Shape now
    MFnDependencyNode sFn( bunsenNodeSObj );
    MPlug localTransPlug( sFn.findPlug( "localPosition", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = localTransPlug.setLocked( true );
    CHECK_MSTATUS( stat );
    MPlug localScalePlug( sFn.findPlug( "localScale", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = localScalePlug.setLocked( true );
    CHECK_MSTATUS( stat );
        
    MString timeStr( "connectAttr -f time1.outTime " + bunsenNodeDagPath.fullPathName() + ".time" );
    dagModifier.commandToExecute( timeStr );
    
    MString gravStr( "setAttr " + bunsenNodeDagPath.fullPathName() + ".gravity 0 -981 0" );
    dagModifier.commandToExecute( gravStr );
    
    dagModifier.doIt();

    o_wmBunsenNodeObj = bunsenNodeSObj;
}

void WmBunsenCmd::getNodes( MSelectionList i_opt_nodes )
{
    // Prepare all the selection lists that this command operates on

    MDagPath mDagPath;
    MObject mObj;
    MItDag dagit;
    MStatus stat;

    m_curve_list.clear();

    for (MItSelectionList sIt(i_opt_nodes); !sIt.isDone(); sIt.next())
    {
        sIt.getDagPath(mDagPath, mObj);
        
        MFnDagNode dagFn( mDagPath.child(0, &stat), &stat);
        CHECK_MSTATUS( stat );
        MDagPath childPath;
        stat = dagFn.getPath( childPath );
        CHECK_MSTATUS( stat );
        childPath.extendToShape();
        //MFnDependencyNode nodeFn(childPath.node(&stat), &stat);
        
         if(( childPath.apiType() == MFn::kNurbsCurve ) )
            {
                mObj = childPath.node();

                m_curve_list.add( childPath, mObj, true);

                //MFnDependencyNode curfn( childPath.node() );
            }
        
        // In this case we let the user select a root group, and then iterate through all
        // the curves below

//        dagit.reset(mDagPath, MItDag::kDepthFirst);

//        while(!dagit.isDone())
//        {
//            dagit.getPath(mDagPath);
//
//            if((mDagPath.apiType() == MFn::kNurbsCurve))
            //{
              //  mObj = mDagPath.node();

                //o_curveList.add(mDagPath, mObj, true);

                //MFnDependencyNode curfn(mDagPath.node());
            //}

//            dagit.next();
//        }
    }
}

MStatus WmBunsenCmd::redoIt()
{
    MStatus mStatus(MS::kSuccess);

    MString opt_fo;

    if (m_mArgDatabase->isFlagSet(kHelp))
    {
        printHelp();
    }
    else
    {
        m_mDagModifier = new MDagModifier;
        m_undoable = true;

        if ( m_mArgDatabase->isFlagSet( kCreateRods ) )
        {
            // Create
            createWmBunsenRodNode( true );
        }
        else
        {
            // Do other stuff such as return number of rods or whatever
        }
    }

    return mStatus;
}

MStatus WmBunsenCmd::undoIt()
{
    return MS::kSuccess;
}


void WmBunsenCmd::p_AddFlag(
                        MSyntax &i_mSyntax,
                        const char *const i_shortName,
                        const char *const i_longName,
                        const char *const i_help,
                        const MSyntax::MArgType i_argType1,
                        const MSyntax::MArgType i_argType2,
                        const MSyntax::MArgType i_argType3,
                        const MSyntax::MArgType i_argType4,
                        const MSyntax::MArgType i_argType5,
                        const MSyntax::MArgType i_argType6)
{
    i_mSyntax.addFlag(i_shortName, i_longName, i_argType1, i_argType2, i_argType3, i_argType4, i_argType5, i_argType6);
    WmBunsenHelp &furHelp = m_help[i_shortName];
    furHelp.m_longName = i_longName;
    furHelp.m_help = i_help;

    if (i_argType1 != MSyntax::kInvalidArgType && i_argType1 != MSyntax::kNoArg) {
        if (i_argType1 == MSyntax::kLastArgType && furHelp.m_argTypes.size()) {
            furHelp.m_argTypes.push_back(furHelp.m_argTypes.back());
        }
        else {
            furHelp.m_argTypes.push_back(i_argType1);
        }
        if (i_argType2 != MSyntax::kInvalidArgType && i_argType2 != MSyntax::kNoArg) {
            if (i_argType2 == MSyntax::kLastArgType && furHelp.m_argTypes.size()) {
                furHelp.m_argTypes.push_back(furHelp.m_argTypes.back());
            }
            else {
                furHelp.m_argTypes.push_back(i_argType2);
            }
            if (i_argType3 != MSyntax::kInvalidArgType && i_argType3 != MSyntax::kNoArg) {
                if (i_argType3 == MSyntax::kLastArgType && furHelp.m_argTypes.size()) {
                    furHelp.m_argTypes.push_back(furHelp.m_argTypes.back());
                }
                else {
                    furHelp.m_argTypes.push_back(i_argType3);
                }
                if (i_argType4 != MSyntax::kInvalidArgType && i_argType4 != MSyntax::kNoArg) {
                    if (i_argType4 == MSyntax::kLastArgType && furHelp.m_argTypes.size()) {
                        furHelp.m_argTypes.push_back(furHelp.m_argTypes.back());
                    }
                    else {
                        furHelp.m_argTypes.push_back(i_argType4);
                    }
                    if (i_argType5 != MSyntax::kInvalidArgType && i_argType5 != MSyntax::kNoArg) {
                        if (i_argType5 == MSyntax::kLastArgType && furHelp.m_argTypes.size()) {
                            furHelp.m_argTypes.push_back(furHelp.m_argTypes.back());
                        }
                        else {
                            furHelp.m_argTypes.push_back(i_argType5);
                        }
                        if (i_argType6 != MSyntax::kInvalidArgType && i_argType6 != MSyntax::kNoArg) {
                            if (i_argType6 == MSyntax::kLastArgType && furHelp.m_argTypes.size()) {
                                furHelp.m_argTypes.push_back(furHelp.m_argTypes.back());
                            }
                            else {
                                furHelp.m_argTypes.push_back(i_argType6);
                            }
                        }
                    }
                }
            }
        }
    }
}



void WmBunsenCmd::printHelp()
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

void WmBunsenCmd::getSelectedCurves( const MSelectionList &i_selectionList,
                                     MSelectionList &o_curveList)
{
    MDagPath mDagPath;
    MObject mObj;
    MItDag dagit;
    MStatus stat;

    o_curveList.clear();

    for (MItSelectionList sIt(i_selectionList); !sIt.isDone(); sIt.next())
    {
        sIt.getDagPath(mDagPath, mObj);
        
        MFnDagNode dagFn( mDagPath.child(0, &stat), &stat);
        CHECK_MSTATUS( stat );
        MDagPath childPath;
        stat = dagFn.getPath( childPath );
        CHECK_MSTATUS( stat );
        childPath.extendToShape();
        //MFnDependencyNode nodeFn(childPath.node(&stat), &stat);
        
         if(( childPath.apiType() == MFn::kNurbsCurve ) )
            {
                mObj = childPath.node();

                o_curveList.add( childPath, mObj, true);

                //MFnDependencyNode curfn( childPath.node() );
            }
        
        // In this case we let the user select a root group, and then iterate through all
        // the curves below

//        dagit.reset(mDagPath, MItDag::kDepthFirst);

//        while(!dagit.isDone())
//        {
//            dagit.getPath(mDagPath);
//
//            if((mDagPath.apiType() == MFn::kNurbsCurve))
            //{
              //  mObj = mDagPath.node();

                //o_curveList.add(mDagPath, mObj, true);

                //MFnDependencyNode curfn(mDagPath.node());
            //}

//            dagit.next();
//        }
    }
}

MStatus WmBunsenCmd::createDagNode( const char *transformName, const char *nodeType, 
    MObject &parentObj, MObject *transformObjP, MObject *shapeObjP, MDagModifier *iDagModifier )
{
    MStatus rStatus = MS::kSuccess;

    MDagModifier tDagModifier;
    MDagModifier &mDagModifier = (iDagModifier != NULL) ? *iDagModifier : tDagModifier;

    MObject transformObj = mDagModifier.createNode(nodeType, MObject::kNullObj, &rStatus);

      if (rStatus) {
        mDagModifier.doIt();
    
        if (!transformObj.isNull()) {
          if (transformObjP)
            *transformObjP = transformObj;
    
          if (strcmp(transformName, nodeType)) {
            rStatus = mDagModifier.renameNode(transformObj, transformName);
          } else {
            MDagPath dagPath;
            MFnDagNode(transformObj).getPath(dagPath);
            dagPath.extendToShape();
            rStatus = mDagModifier.renameNode(transformObj, MFnDependencyNode(dagPath.node()).name());
          }
    
          if (rStatus) {
            mDagModifier.doIt();
    
            MDagPath shapeDagPath = MDagPath::getAPathTo(transformObj);
            const unsigned childCount = shapeDagPath.childCount();
    
            if (childCount == 1) {
              shapeDagPath.extendToShape();
              MObject shapeObj = shapeDagPath.node();
    
              if (shapeObjP)
                *shapeObjP = shapeObj;
    
              MFnTransform transformFn(transformObj);
    
              const MString transformNameStr = transformFn.name();
              const char *sp = transformNameStr.asChar();
              const char *ep = sp + transformNameStr.length() - 1;
              const char *cp = ep;
    
              while (cp >= sp && isdigit(*cp))
                cp--;
    
              MString shapeName;
              if (cp != sp)
                shapeName = transformNameStr.substring(0, cp - sp);
              else
                shapeName = transformNameStr;
    
              shapeName += "Shape";
    
              if (cp != ep)
                shapeName += transformNameStr.substring(cp - sp + 1, ep - sp);
    
              rStatus = mDagModifier.renameNode(shapeObj, shapeName);
              if (rStatus)
                mDagModifier.doIt();
              else
                cerr << "Couldn't Rename Shape To: \"" << shapeName << "\"" << endl;
            }
            else {
              if (shapeObjP)
                *shapeObjP = MObject::kNullObj;
    
              if (childCount)
                cerr << "Unknown Child Count == " << childCount << endl;
            }
          }
          else {
            cerr << "Couldn't Rename Transform Node To \"" << transformName << "\"" << endl;
          }
    
          if (parentObj != MObject::kNullObj) {
            rStatus = mDagModifier.reparentNode(transformObj, parentObj);
            if (rStatus)
              mDagModifier.doIt();
            else
              cerr << "Couldn't Reparent Node" << endl;
          }
        }
        else {
          cerr << "Couldn't Create Node Of \"" << nodeType << "\"" << endl;
          rStatus = MS::kFailure;
        }
      }
      else {
        cerr << "Couldn't Create Transform Node \"" << transformName << "\"" << endl;
      }
    
    return rStatus;
}

