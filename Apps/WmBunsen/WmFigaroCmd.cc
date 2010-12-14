#include "WmFigaroCmd.hh"
#include "WmBunsenCollisionMeshNode.hh"
#include "WmFigConnectionNode.hh"
#include "WmFigaroRodShape/WmFigaroRodShape.hh"
#include "constraints/WmFigConstraintNode.hh"

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

// Statics
/* static */ std::map<std::string, WmBunsenHelp> WmFigaroCmd::m_help;
/* static */ MString WmFigaroCmd::typeName("wmFigaro");
/* static */ MStringArray WmFigaroCmd::m_results;
//MObject WmFigaroCmd::m_dynNode;

// CTOR

WmFigaroCmd::WmFigaroCmd()
  : m_undoable( false ),
    m_mArgDatabase( NULL ),
    m_mDagModifier( NULL ),
    m_selectedwmBunsenNode( MObject::kNullObj ),
    m_cacheFile( "" ),
    m_rodNumber( 0 ),
    m_edgeNumber( 0 ),
    m_color( 0, 0, 0 )/*,
    m_undo(NULL)*/
{
    m_results.clear();
}

// DTOR

WmFigaroCmd::~WmFigaroCmd()
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

void* WmFigaroCmd::creator()
{
    return new WmFigaroCmd;
}

const char *const kCreateRods( "-cr" );
const char *const kCreateRodsFromFozzie( "-crf" );
const char *const kCVsPerRod( "-cpr" );
const char *const kName( "-n" );
const char *const kAddCollisionMeshes( "-acm" );
const char *const kPreviewCache( "-prc" );
const char *const kCacheFile( "-caf" );
const char *const kAttachEdgeToObject( "-aeo" );
const char *const kEdge( "-ed" );
const char *const kRod( "-ro" );
const char *const kUserDefinedColour( "-udc" );
const char *const kDeleteDefinedColour( "-ddc" );
const char *const kCreateRodShapeNode( "-crs" );
const char *const kAttachParticles( "-aps" );

// Rod vertex constraints
//
const char *const kAddVertexConstraintFlag( "-avc");
const char *const kAddVertexConstraintLongFlag( "-addVertexConstraint");
const char *const kRodNodeNameFlag( "-rn");
const char *const kRodNodeNameLongFlag( "-rodNode");
const char *const kRodVerticesFlag( "-rv");
const char *const kRodVerticesLongFlag( "-rodVertices");
const char *const kTargetPositionFlag( "-tp");
const char *const kTargetPositionLongFlag( "-targetPosition");
const char *const kStiffnessFlag( "-st");
const char *const kStiffnessLongFlag( "-stiffness");
const char *const kVertexConstraintTypeFlag("-vct");
const char *const kVertexConstraintTypeLongFlag("-vertexConstraintType");

// Rod vertices
const char *const kRodVertexPositionFlag("-rvp");
const char *const kRodVertexPositionLongFlag("-rodVertexPosition");
const char *const kWorldSpaceFlag("-ws");
const char *const kWorldSpaceLongFlag("-worldSpace");

const char *const kHelp( "-h" );

MSyntax WmFigaroCmd::syntaxCreator()
{
    m_help.clear();

    MSyntax mSyntax;

    p_AddFlag( mSyntax, kHelp, "-help",
               "Prints this information" );
    p_AddFlag( mSyntax, kName, "-name",
               "Sets the name of the created WmFigaro node",
               MSyntax::kString );
    p_AddFlag( mSyntax, kCreateRods, "-createRods",
               "creates discrete elastic rods from the selected NURBS curves." );
    p_AddFlag( mSyntax, kCreateRodsFromFozzie, "-createRodsFromFozzie",
               "creates discrete elastic rods from the selected BarberShop node." );
    p_AddFlag( mSyntax, kCVsPerRod, "-cvsPerRod",
               "Sets the number of CVs to use for each created elastic rod.", MSyntax::kLong );
    p_AddFlag( mSyntax, kAddCollisionMeshes, "-addCollisionMesh",
               "Adds a a selected polygon mesh as a collision object for the selected Bunsen node." );
    
    p_AddFlag( mSyntax, kAttachEdgeToObject, "-attachEdgeToObject",
               "Creates a connection node to allow a rod edge to be controlled by the selected object." );
    p_AddFlag( mSyntax, kEdge, "-edgeNumber",
        "Edge to apply operation to.", MSyntax::kLong);
    p_AddFlag( mSyntax, kRod, "-rodNumber",
        "Rod to apply operation to.", MSyntax::kLong);
    
    p_AddFlag( mSyntax, kUserDefinedColour, "-userDefinedColour",
        "Colour to draw specified rods in.", MSyntax::kString);
    p_AddFlag( mSyntax, kDeleteDefinedColour, "-deleteUserDefinedColour",
        "Remove colour override from a rod.", MSyntax::kString);
    
    p_AddFlag( mSyntax, kPreviewCache, "-previewCache",
               "Creates wmFigaro and wmFigRodNodes and to load the specified cache file for preview." );
    p_AddFlag( mSyntax, kCacheFile, "-cacheFile",
        "Cache file to load and preview.", MSyntax::kString);
    p_AddFlag( mSyntax, kCreateRodShapeNode, "-createRodShape",
               "Creates a single rod as a shape node from a selected NURBS curve." );

    p_AddFlag( mSyntax, kAttachParticles, "-attachParticles",
               "Creates rods from a particle system." );


    // Rod vertex constraints
    //
    p_AddFlag( mSyntax,
    		   kAddVertexConstraintFlag,
    		   kAddVertexConstraintLongFlag,
    		   "add vertex constraint for rod ");

    p_AddFlag(  mSyntax,
    			kRodNodeNameFlag,
    			kRodNodeNameLongFlag,
    			"specify rod node name",
    			MSyntax::kString );

    p_AddFlag(  mSyntax,
    			kRodVerticesFlag,
    			kRodVerticesLongFlag,
    			"specify rod vertices",
    			MSyntax::kString );

    p_AddFlag(  mSyntax,
    			kTargetPositionFlag,
    			kTargetPositionLongFlag,
    			"specify target world position of the contraint",
    			MSyntax::kDouble,
    			MSyntax::kDouble,
    			MSyntax::kDouble );

    p_AddFlag(  mSyntax,
    			kStiffnessFlag,
    			kStiffnessLongFlag,
    			"specify the constraint's stiffness",
    			MSyntax::kDouble );

    p_AddFlag( mSyntax,
    		   kVertexConstraintTypeFlag,
    		   kVertexConstraintTypeLongFlag,
    		   "specify the vertex constraint type",
    		   MSyntax::kString );
    
    // Rod vertex position
    //
    p_AddFlag( mSyntax,
    			kRodVertexPositionFlag,
    			kRodVertexPositionLongFlag,
    		   "get the position of the rod vertex",
    		   MSyntax::kLong,
    		   MSyntax::kLong );
    p_AddFlag(  mSyntax,
    			kWorldSpaceFlag,
    			kWorldSpaceLongFlag,
    			"get position in world space");

    mSyntax.setObjectType( MSyntax::kSelectionList );
    mSyntax.useSelectionAsDefault( true );
    mSyntax.enableEdit( true );

    return mSyntax;
}

MStatus WmFigaroCmd::doIt( const MArgList &i_mArgList )
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

            m_node_name = MString("wmFigaro");

            m_cvsPerRod = -1; // Use the number of cvs from the first curve unless the user specifies a number of cvs

            if ( m_mArgDatabase->isFlagSet( kName ) )
                m_mArgDatabase->getFlagArgument( kName, 0, m_node_name );
            if ( m_mArgDatabase->isFlagSet( kCreateRods ) )
                m_mArgDatabase->getFlagArgument( kCVsPerRod, 0, m_cvsPerRod );
            if ( m_mArgDatabase->isFlagSet( kCacheFile ) )
            {
                m_mArgDatabase->getFlagArgument( kCacheFile, 0, m_cacheFile );
            }
            if ( m_mArgDatabase->isFlagSet( kEdge ) )
            {
                m_mArgDatabase->getFlagArgument( kEdge, 0, m_edgeNumber );
            }
            if ( m_mArgDatabase->isFlagSet( kRod ) )
            {
                m_mArgDatabase->getFlagArgument( kRod, 0, m_rodNumber );
            }
            if ( m_mArgDatabase->isFlagSet( kUserDefinedColour ) )
            {
                MString colourString;
                m_mArgDatabase->getFlagArgument( kUserDefinedColour, 0, colourString );
                
                // We pass the colour as a string because there is no int3 option
                // to commands. So we now need to split it up
                MStringArray subStrings;
                colourString.split( ' ', subStrings );
                
                if ( subStrings.length() != 3 )
                   displayError( "Please enter a colour with 3 components - r,g and b." );
               else
               {
                   m_color[ 0 ] = subStrings[ 0 ].asFloat();
                   m_color[ 1 ] = subStrings[ 1 ].asFloat();
                   m_color[ 2 ] = subStrings[ 2 ].asFloat();
               }
            }
            
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

MStatus WmFigaroCmd::redoIt()
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
            createWmFigRodNode( true );
        }
        if ( m_mArgDatabase->isFlagSet( kCreateRodsFromFozzie ) )
        {
            // Create
            createWmFigRodNode( false );
        }
        if ( m_mArgDatabase->isFlagSet( kAddCollisionMeshes ) )
        {
            addCollisionMeshes();
        }
        if ( m_mArgDatabase->isFlagSet( kPreviewCache ) )
        {
            createPreviewNodes();
        }
        if ( m_mArgDatabase->isFlagSet( kAttachEdgeToObject ) )
        {
            attatchEdgeToObject();
        }
        if ( m_mArgDatabase->isFlagSet( kUserDefinedColour ) )
        {
            setColorOfRod();
        }
        if ( m_mArgDatabase->isFlagSet( kDeleteDefinedColour ) )
        {
            setColorOfRod( true );
        }
        if ( m_mArgDatabase->isFlagSet( kCreateRodShapeNode ) )
        {
            createRodShapeNode();
        }
        if ( m_mArgDatabase->isFlagSet( kAttachParticles) )
        {
            attachParticles();
        }

        if( m_mArgDatabase->isFlagSet( kAddVertexConstraintFlag ))
        {
/*
	Example:

        	wmFigaro -addVertexConstraint -rodNode "wmFigRodNodeShape1" -rodVertices "1 rod[0].vtx[21]"; // The rodVertices string the string result from the "Select Rod Vertices Tool"

			Others:

			wmFigaro -addVertexConstraint -rodNode "wmFigRodNodeShape1" -rodIdVertices "1 rod[0].vtx[21]" -vertexConstraintType "fixed"

			- select wmFigRodNode then execute:
			wmFigaro -addVertexConstraint -rodVertices "1 rod[0].vtx[21]" -vertexConstraintType "fixed"

			World Position:

			// Set it explicitly
			wmFigaro -addVertexConstraint -rodNode "wmFigRodNodeShape1"  -rodVertices "1 rod[0].vtx[21]" -vertexConstraintType "fixed" -worldPosition 20 0 0

			// Use a locator
			- select the locator then execute:
			wmFigaro -addVertexConstraint -rodNode "wmFigRodNodeShape1"  -rodVertices "1 rod[0].vtx[21]" -vertexConstraintType "fixed"

*/
        	MStatus stat;
		    MSelectionList selList;
		    MGlobal::getActiveSelectionList( selList );

        	int constraintType = 1; // Default is rest constraint
			if( m_mArgDatabase->isFlagSet( kVertexConstraintTypeFlag ) )
			{
				MString constraintTypeTxt;
				m_mArgDatabase->getFlagArgument( kVertexConstraintTypeFlag, 0, constraintTypeTxt );

			    if( constraintTypeTxt == "fixed" )
			    	constraintType = 0;
			}

			MString rodVerticesTxt;
			if( m_mArgDatabase->isFlagSet( kRodVerticesFlag ) )
				m_mArgDatabase->getFlagArgument( kRodVerticesFlag, 0, rodVerticesTxt );

			MString figRodNodeName;
			if( m_mArgDatabase->isFlagSet( kRodNodeNameFlag ) )
				m_mArgDatabase->getFlagArgument( kRodNodeNameFlag, 0, figRodNodeName );

#if 0
			int rodId = -1;
			int vertexId = -1;

			if( rodVertexPath.length() )
			{
				// Parse the rod vertex path string, e.g. "|wmFigRodNode1.rod[0].vtx[21]"
				//
				MGlobal::executePythonCommand( MString("\"") + rodVertexPath.asChar() + "\".split('.')[0]", figRodNodeName );

				MString rodIdTxt;
				MGlobal::executePythonCommand( MString("\"") + rodVertexPath.asChar() + "\".split('.rod[')[-1].split(']')[0]", rodIdTxt );
				rodId = rodIdTxt.asInt();

				MString vertexIdTxt;
				MGlobal::executePythonCommand( MString("\"") + rodVertexPath.asChar() + "\".split('.vtx[')[-1].split(']')[0]", vertexIdTxt );
				vertexId = vertexIdTxt.asInt();

			}
			else
			{
				// Get the parts of the rod vertex path as separate arguments to the command
				//
				if( m_mArgDatabase->isFlagSet( kRodNodeNameFlag ) )
					m_mArgDatabase->getFlagArgument( kRodNodeNameFlag, 0, figRodNodeName );

				if( m_mArgDatabase->isFlagSet( kRodIdFlag ) )
					m_mArgDatabase->getFlagArgument( kRodIdFlag, 0, rodId );

				if(  m_mArgDatabase->isFlagSet( kVertIdFlag ) )
					m_mArgDatabase->getFlagArgument( kVertIdFlag, 0, vertexId );
			}
#endif

			// The rod's shape path wasn't given explicitly, get the first rod shape in the selection
			//
			if( !figRodNodeName.length() )
			{
			    unsigned int i;
				for( i=0; i < selList.length(); i++ )
			    {
			        MDagPath dagPath;
			        selList.getDagPath( i, dagPath );
			        dagPath.extendToShape();

			        MFnDependencyNode nodeFn( dagPath.node() );
			        if( nodeFn.typeId() == WmFigRodNode::typeID )
			        {
			        	figRodNodeName = dagPath.fullPathName();
			            break;
			        }
			    }
			}

			if( !figRodNodeName.length() )
			{
				MGlobal::displayInfo( "No rod shape path given or selected." );
				return MS::kFailure;
			}

			// Convert rodShapeTxt into MObject
			MSelectionList sel;
			MObject figRodNode;
			sel.add( figRodNodeName );
			sel.getDependNode( 0, figRodNode );

			if( figRodNode.isNull() )
			{
				MGlobal::displayInfo( "Rod shape: " + figRodNodeName + " doesn't exist." );
				return MS::kFailure;
			}

			MFnDependencyNode figRodNodeFn( figRodNode );

			// Get the wmFigaroNode node
			//
			MObject figNode;
			MPlugArray destPlugs;
			MPlug rodsChangedPlug = figRodNodeFn.findPlug( "rodsChanged" );
			rodsChangedPlug.connectedTo( destPlugs, false, true );
			if( destPlugs.length() )
			{
				figNode = destPlugs[0].node();
//				MFnDependencyNode nodeFn( figNode );
//				MGlobal::displayInfo( "FOUND: " + nodeFn.name() );
			}

			if( figNode.isNull() )
			{
				MGlobal::displayInfo( "No wmFigaroNode connected to the rod shape node: " + figRodNodeName );
				return MS::kFailure;
			}

			MFnDependencyNode figNodeFn( figNode );

			// Set the stiffness
			//
			double stiffness = 50.0;
			if( m_mArgDatabase->isFlagSet( kStiffnessFlag ) )
				m_mArgDatabase->getFlagArgument( kStiffnessFlag, 0, stiffness );


			// Create a wmFigConstraintNode
			//
			MString figConstraintNodeName;
			stat = MGlobal::executeCommand( MString("createNode \"") + WmFigConstraintNode::TypeName + "\"", figConstraintNodeName );
			if( !stat )
			{
				MGlobal::displayInfo( MString("Unable to create a node of type ") + WmFigConstraintNode::TypeName );
				return MS::kFailure;
			}

			MObject figConstraintNode;
			sel.clear();
			sel.add( figConstraintNodeName );
			sel.getDependNode( 0, figConstraintNode );

			MFnDependencyNode constraintNodeFn( figConstraintNode );
			MPlug constraintTypePlug = constraintNodeFn.findPlug( "constraintType" );
			constraintTypePlug.setValue( constraintType );

			// Given the worldPosition explicitly
			//
			MObject worldPositionLocatorNode;
			MPoint pos(0.0,0.0,0.0);
			if( m_mArgDatabase->isFlagSet( kTargetPositionFlag ) )
			{
				m_mArgDatabase->getFlagArgument( kTargetPositionFlag, 0, pos.x );
				m_mArgDatabase->getFlagArgument( kTargetPositionFlag, 1, pos.y );
				m_mArgDatabase->getFlagArgument( kTargetPositionFlag, 2, pos.z );
			}
			else // Is there a locator in the selection?
			{
			    unsigned int i;
				for( i=0; i < selList.length(); i++ )
			    {
			        MDagPath dagPath;
			        selList.getDagPath( i, dagPath );
			        dagPath.extendToShape();
			        if( dagPath.apiType() == MFn::kLocator )
			        {
			        	worldPositionLocatorNode = dagPath.node();
			        	break;
			        }
			    }
			}

			MPlug worldPositionPlug = constraintNodeFn.findPlug( "targetWorldPosition" );
			if( !worldPositionLocatorNode.isNull() )
			{

				MFnDependencyNode locatorFn( worldPositionLocatorNode );
				MPlug srcWorldPosPlug = locatorFn.findPlug( "worldPosition" );
				srcWorldPosPlug.selectAncestorLogicalIndex(0);

				MPlug srcWorldPosXPlug = srcWorldPosPlug.child(0);
				MPlug srcWorldPosYPlug = srcWorldPosPlug.child(1);
				MPlug srcWorldPosZPlug = srcWorldPosPlug.child(2);

				//srcWorldPosXPlug = srcWorldPosXPlug.elementByPhysicalIndex(0);
				//srcWorldPosXPlug.selectAncestorLogicalIndex( 0 );

				MPlug destWorldPosXPlug = constraintNodeFn.findPlug( "targetWorldPosition0" );
				MPlug destWorldPosYPlug = constraintNodeFn.findPlug( "targetWorldPosition1" );
				MPlug destWorldPosZPlug = constraintNodeFn.findPlug( "targetWorldPosition2" );

				MGlobal::executeCommand( "connectAttr -f " + srcWorldPosXPlug.info() + " " +  destWorldPosXPlug.info() );
				MGlobal::executeCommand( "connectAttr -f " + srcWorldPosYPlug.info() + " " +  destWorldPosYPlug.info() );
				MGlobal::executeCommand( "connectAttr -f " + srcWorldPosZPlug.info() + " " +  destWorldPosZPlug.info() );

				//MGlobal::displayInfo( MString("FOUND LOCATOR: ") + locatorFn.name() );
			}
			else
			{
				MFnNumericData numericData;
				MObject posData = numericData.create( MFnNumericData::k3Double );
				numericData.setData( pos.x, pos.y, pos.z );
				worldPositionPlug.setMObject( posData );
			}

//			MPlug rodIdPlug = constraintNodeFn.findPlug( "rodId" );
//		    rodIdPlug.setValue( rodId );
//
//		    MPlug vertexIdPlug = constraintNodeFn.findPlug( "vertexId" );
//		    vertexIdPlug.setValue( vertexId );

			MPlug rodVerticesPlug = constraintNodeFn.findPlug( "rodVertices" );
			rodVerticesPlug.setValue( rodVerticesTxt );
//			MFnStringData stringDataFn( rodVerticesPlug.asMObject() );
//			stringDataFn.create( rodVerticesTxt );
//			stringDataFn.set( rodVerticesTxt );

			MPlug stiffnessPlug = constraintNodeFn.findPlug( "stiffness" );
		    stiffnessPlug.setValue( stiffness );

		    MPlug figRodNodeMsgPlug = constraintNodeFn.findPlug( "figRodNodeMsg" );
			MGlobal::executeCommand( "connectAttr -f " + figRodNodeName + ".message " +  figRodNodeMsgPlug.info() );

			MGlobal::displayInfo( "Figaro node: " + figNodeFn.name() );

			// Register the constraint node with the wmFigNode node
			//
			MPlug msgConstraintNodesPlug = figNodeFn.findPlug( "ia_msgConstraintNodes" );
			int nConstraintNodes = msgConstraintNodesPlug.evaluateNumElements();

			MPlug msgConstrainNodePlug = msgConstraintNodesPlug.elementByLogicalIndex( nConstraintNodes ); // Add new element on the end
			MGlobal::executeCommand( "connectAttr -f " + figConstraintNodeName + ".message " + msgConstrainNodePlug.info() );

			MGlobal::displayInfo( MString("figConstraintNode: ") + figConstraintNodeName );
			MGlobal::displayInfo( figRodNodeName );
//			MGlobal::displayInfo( MString( "rodId: ") + rodId );
//			MGlobal::displayInfo( MString( "vertexId: ") + vertexId );
			MGlobal::displayInfo( MString( "rodVertices: ") + rodVerticesTxt );
			MGlobal::displayInfo( MString( "constraintType: ") + constraintType );
			MGlobal::displayInfo( MString( "worldPosition: ") + pos.x + ", " + pos.y + ", " + pos.z );
			MGlobal::displayInfo( MString( "stiffness: ") + stiffness );
        }
        
        // Get rod vertex position
        //
/*
		Example:

			wmFigaro -addVertexConstraint -rodNode "wmFigRodNodeShape1" -rodVertices "1 rod[0].vtx[21]"; // The rodVertices string the string result from the "Select Rod Vertices Tool"
*/
        if( m_mArgDatabase->isFlagSet( kRodVertexPositionFlag ) )
        {
        	int rodId;
        	int vertexId;
        	m_mArgDatabase->getFlagArgument( kRodVertexPositionFlag, 0, rodId );
        	m_mArgDatabase->getFlagArgument( kRodVertexPositionFlag, 1, vertexId );

        	bool inWorldSpace = false;
        	if( m_mArgDatabase->isFlagSet( kWorldSpaceFlag ) )
        		inWorldSpace = true;

			MString figRodNodeName;
			if( m_mArgDatabase->isFlagSet( kRodNodeNameFlag ) )
				m_mArgDatabase->getFlagArgument( kRodNodeNameFlag, 0, figRodNodeName );

			if( !figRodNodeName.length() )
			{
				MGlobal::displayInfo( "No rod shape name given." );
				return MS::kFailure;
			}

			// Convert figRodNodeName into MObject
			MSelectionList sel;
			MObject figRodNode;
			sel.add( figRodNodeName );
			sel.getDependNode( 0, figRodNode );

			if( figRodNode.isNull() )
			{
				MGlobal::displayInfo( "Rod shape: " + figRodNodeName + " doesn't exist." );
				return MS::kFailure;
			}

			MFnDependencyNode figRodNodeFn( figRodNode );
		    WmFigRodNode *figRod = static_cast<WmFigRodNode*>( figRodNodeFn.userNode() );
		    if( !figRod )
		    {
				MGlobal::displayInfo( "Unable to get WmFigRodNode." );
				return MS::kFailure;
		    }

		    WmFigRodGroup* rodGroup = figRod->rodGroup();
		    const int nRods = rodGroup->numberOfRods();

		    if( rodId < 0 || rodId >= nRods ) {
				MGlobal::displayInfo( MString("Rod index ") + rodId + " is invalid. The node contains " + nRods + " rods." );
				return MS::kFailure;
		    }

		    const int nVerts = rodGroup->elasticRod( rodId )->nv();

		    if( vertexId < 0 || vertexId >= nVerts ) {
				MGlobal::displayInfo( MString("Vertex index ") + vertexId + " is invalid. The rod contains " + nVerts + " vertices." );
				return MS::kFailure;
		    }

			const BASim::Vec3d p = rodGroup->elasticRod( rodId )->getVertex( vertexId );
			MPoint pt( p[0], p[1], p[2] );

			// Transform local space position to world space
			//
			if( inWorldSpace ) {
				MFnDagNode dagNodeFn( figRodNode );
				MObject transformNode = dagNodeFn.parent( 0 );
				dagNodeFn.setObject( transformNode );

				MDagPath transformDagPath;
				dagNodeFn.getPath( transformDagPath );

				MMatrix worldMtx = transformDagPath.inclusiveMatrix();
				pt *= worldMtx;
			}

			MDoubleArray coordinates( 3 );
			coordinates[0] = pt.x;
			coordinates[1] = pt.y;
			coordinates[2] = pt.z;
			setResult( coordinates );

			return MS::kSuccess;
        }


        setResult( m_results );
    }

    return mStatus;
}

MStatus WmFigaroCmd::undoIt()
{
    return MS::kSuccess;
}

void WmFigaroCmd::setColorOfRod( bool i_unset )
{
    MStatus stat;
    
    if ( m_figRodNodeList.length() != 1 )
    {
        displayError( "Please select the single wmFigRodNode that you wish to set the color of." );
        return;
    }
    
    MObject rodNodeObj;
    stat = m_figRodNodeList.getDependNode( 0, rodNodeObj);
    CHECK_MSTATUS( stat );
    
    MPlug colorPlugArray( rodNodeObj, WmFigRodNode::ia_userDefinedColors );
    
    MPlug colorPlug = colorPlugArray.elementByLogicalIndex( m_rodNumber, &stat );
    CHECK_MSTATUS( stat );
    
    MFnNumericData colorData;
    MObject colorObj = colorData.create( MFnNumericData::k3Double, &stat );
    CHECK_MSTATUS( stat );
 
    if ( i_unset )
        stat = colorData.setData3Double( -1, -1, -1 );
    else
        stat = colorData.setData3Double( m_color[ 0 ], m_color[ 1 ], m_color[ 2 ] );
    
    CHECK_MSTATUS( stat );
    
    stat = colorPlug.setValue( colorObj );
    CHECK_MSTATUS( stat );
}

void WmFigaroCmd::createWmFigRodNode( bool useNURBSInput, bool i_previewOnly,
                                                      MObject* o_rodNode )
{
    MStatus stat;
    size_t nCurves = 0;
    
    MObject wmBunsenNode = m_selectedwmBunsenNode;

    if ( useNURBSInput && m_nurbsCurveList.isEmpty() )
    {
        MGlobal::displayError( "Please select some NURBS curves to create rods from\n" );
        return;
    }
    else if ( !useNURBSInput && m_fozzieNodeList.isEmpty() && !i_previewOnly )
    {
        MGlobal::displayError( "Please select a Fozzie furset node to create rods from\n" );
        return;
    }
    
    if ( wmBunsenNode == MObject::kNullObj )
        createWmBunsenNode( wmBunsenNode );

    // Create the rods node
    MObject rodTObj;  // Object for transform node
    MObject rodSObj;  // Object for shape node
    MDagPath shapeDagPath;
    MObject pObj;
    MDagModifier dagModifier;
    MString rodShapeName = "";

    createDagNode( WmFigRodNode::typeName.asChar(), WmFigRodNode::typeName.asChar(),
                   pObj, &rodTObj, &rodSObj, &dagModifier, rodShapeName );

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
        for ( uint l = 0; l < m_nurbsCurveList.length(); l++ ) 
        {
            MDagPath dagPath;
            MObject component;
            m_nurbsCurveList.getDagPath( l, dagPath, component );
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
                m_nurbsCurveList.getDagPath( l, curveTransformPath, component );
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
    else if ( !i_previewOnly )
    {
        cerr << "Connecting fur set\n";
        // Just take the first furset selected, selecting multiple sets is pointless anyway
        MDagPath dagPath;
        MObject component;
        m_fozzieNodeList.getDagPath( 0, dagPath, component );
        dagPath.extendToShape();

        MObject nodeObj = dagPath.node( &stat );
        CHECK_MSTATUS( stat );
        MFnDependencyNode fozzieNodeFn( nodeObj,&stat );
        
        MPlug fozzieStrandsPlug = fozzieNodeFn.findPlug( "strandVertices", true, &stat );
        CHECK_MSTATUS( stat );
        
        MPlug fozzieVerticesPlug( sFn.findPlug( "fozzieVertices", true, &stat ) );
        CHECK_MSTATUS( stat );
        stat = dagModifier.connect( fozzieStrandsPlug, fozzieVerticesPlug );
        CHECK_MSTATUS( stat );
        
        MPlug numFozzieCVsPlug = fozzieNodeFn.findPlug( "outCvsPerStrand", true, &stat );
        CHECK_MSTATUS( stat );
        MPlug numCVsPlug( sFn.findPlug( "cvsPerRod", true, &stat ) );
        CHECK_MSTATUS( stat );
        stat = dagModifier.connect( numFozzieCVsPlug, numCVsPlug );
        CHECK_MSTATUS( stat );
        
        MPlug strandRootsBarbPlug = fozzieNodeFn.findPlug( "strandRootFrames", true, &stat );
        CHECK_MSTATUS( stat );
        MPlug strandRootPlug( sFn.findPlug( "strandRootFrames", true, &stat ) );
        CHECK_MSTATUS( stat );
        stat = dagModifier.connect( strandRootsBarbPlug, strandRootPlug );
        CHECK_MSTATUS( stat );
                                  
        stat = dagModifier.doIt();
        CHECK_MSTATUS( stat );
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
    
    // Connect up the rod output sim plug so we know when the sim moves forward in time
    MPlug rodNodeSimPlug( sFn.findPlug( "simStepTaken", true, &stat ) );
    CHECK_MSTATUS( stat );
    MPlug bunsenSimPlug = bunsenDependNodeFn.findPlug( "simStepTaken", true, &stat ); 
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( bunsenSimPlug, rodNodeSimPlug );
    CHECK_MSTATUS( stat );

    // Connect the solver type so that rod nodes know what type of solver to create
    // when building the Elastic rods
    MPlug rodNodeSolverPlug( sFn.findPlug( "solverType", true, &stat ) );
    CHECK_MSTATUS( stat );
    MPlug figaroSolverPlug = bunsenDependNodeFn.findPlug( "solverType", true, &stat ); 
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( figaroSolverPlug, rodNodeSolverPlug );
    CHECK_MSTATUS( stat );

    // Connect gravity so the rods inherit the global system gravity
    MPlug rodNodeGravityPlug( sFn.findPlug( "gravity", true, &stat ) );
    CHECK_MSTATUS( stat );
    MPlug figaroGravityPlug = bunsenDependNodeFn.findPlug( "gravity", true, &stat ); 
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( figaroGravityPlug, rodNodeGravityPlug );
    CHECK_MSTATUS( stat );
    
    // Connect startTime so the rods reset when the Figaro node resets
    MPlug rodStartTimePlug( sFn.findPlug( "startTime", true, &stat ) );
    CHECK_MSTATUS( stat );
    MPlug bunsenNodeStartTimePlug = bunsenDependNodeFn.findPlug( "startTime", true, &stat ); 
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( bunsenNodeStartTimePlug, rodStartTimePlug );
    CHECK_MSTATUS( stat );
    
    stat = dagModifier.doIt();
    CHECK_MSTATUS( stat );
    
    if ( o_rodNode != NULL )
        *o_rodNode = rodSObj;
    
    //m_dynNode = rodSObj;

    // If we had no well formed curves to use for shaping then just return
    /*if ( nCurves == 0 ) {
        MGlobal::displayError( "Did not find any curves to create rods from. " );                               
        return;
    }*/
}

void WmFigaroCmd::appendToResultString( MString& i_resultString )
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

void WmFigaroCmd::createWmBunsenNode( MObject &o_wmBunsenNodeObj ) 
{    
 
    MStatus stat;
    
    MObject bunsenNodeTObj;  // Object for transform node
    MObject bunsenNodeSObj;  // Object for shape node
    MDagPath shapeDagPath;
    MObject pObj;
    MDagModifier dagModifier;
    MString bunsenShapeName = "";

    createDagNode( WmBunsenNode::typeName.asChar(), WmBunsenNode::typeName.asChar(), 
                   pObj, &bunsenNodeTObj, &bunsenNodeSObj, &dagModifier, bunsenShapeName );
    
    appendToResultString( bunsenShapeName );
                          
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

void WmFigaroCmd::attachParticles()
{
    MStatus stat;

    MSelectionList selList;
    MGlobal::getActiveSelectionList( selList );
    // we're looking for a ParticleShape node
    bool foundParticleNode = false;
    MDagPath particleDagPath;
    for( int i=0; i < selList.length(); i++ )
    {
        selList.getDagPath( i, particleDagPath );
        particleDagPath.extendToShape();

        MFnDependencyNode nodeFn( particleDagPath.node() );
        if( nodeFn.typeId().id() == 1313882450 ) // typeid of nParticleShape; TODO: surely there's a better way?
        {
            foundParticleNode = true;
            break;
        }
    }

    if(!foundParticleNode)
    {
        MGlobal::displayError( "No particle node selected!\n" );
        return;
    }

    MObject wmBunsenNode = m_selectedwmBunsenNode;

    if ( wmBunsenNode == MObject::kNullObj )
        createWmBunsenNode( wmBunsenNode );

    // Create the rods node
    MObject rodTObj;  // Object for transform node
    MObject rodSObj;  // Object for shape node
    MDagPath shapeDagPath;
    MObject pObj;
    MDagModifier dagModifier;
    MString rodShapeName = "";

    createDagNode( WmFigRodNode::typeName.asChar(), WmFigRodNode::typeName.asChar(),
                   pObj, &rodTObj, &rodSObj, &dagModifier, rodShapeName );

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
    MPlug localScalePlug(sFn.findPlug( "localScale", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = localScalePlug.setLocked( true );
    CHECK_MSTATUS( stat );

    MString timeStr( "connectAttr -f time1.outTime " + rodDagPath.fullPathName() + ".time" );
    dagModifier.commandToExecute( timeStr );
    dagModifier.doIt();

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



    MPlug bunsenEnabledPlug = bunsenDependNodeFn.findPlug("enabled", &stat);
    CHECK_MSTATUS( stat );
    bunsenEnabledPlug.setBool(false);

    // Connect up the rod output sim plug so we know when the sim moves forward in time
    MPlug rodNodeSimPlug( sFn.findPlug( "simStepTaken", true, &stat ) );
    CHECK_MSTATUS( stat );
    MPlug bunsenSimPlug = bunsenDependNodeFn.findPlug( "simStepTaken", true, &stat );
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( bunsenSimPlug, rodNodeSimPlug );
    CHECK_MSTATUS( stat );

    // Connect the solver type so that rod nodes know what type of solver to create
    // when building the Elastic rods
    MPlug rodNodeSolverPlug( sFn.findPlug( "solverType", true, &stat ) );
    CHECK_MSTATUS( stat );
    MPlug figaroSolverPlug = bunsenDependNodeFn.findPlug( "solverType", true, &stat );
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( figaroSolverPlug, rodNodeSolverPlug );
    CHECK_MSTATUS( stat );

    // Connect gravity so the rods inherit the global system gravity
    MPlug rodNodeGravityPlug( sFn.findPlug( "gravity", true, &stat ) );
    CHECK_MSTATUS( stat );
    MPlug figaroGravityPlug = bunsenDependNodeFn.findPlug( "gravity", true, &stat );
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( figaroGravityPlug, rodNodeGravityPlug );
    CHECK_MSTATUS( stat );

    // Connect startTime so the rods reset when the Figaro node resets
    MPlug rodStartTimePlug( sFn.findPlug( "startTime", true, &stat ) );
    CHECK_MSTATUS( stat );
    MPlug bunsenNodeStartTimePlug = bunsenDependNodeFn.findPlug( "startTime", true, &stat );
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( bunsenNodeStartTimePlug, rodStartTimePlug );
    CHECK_MSTATUS( stat );

    stat = dagModifier.doIt();
    CHECK_MSTATUS( stat );

    // connect the partice node and the rod node
    MFnDependencyNode particleNodeFn( particleDagPath.node() );
    MString cmdStr("connectAttr -f " + particleNodeFn.name() + ".position " + rodDagPath.fullPathName() + ".particlePositions;");
    dagModifier.commandToExecute( cmdStr );
    dagModifier.doIt();

    cmdStr = MString("connectAttr -f " + particleNodeFn.name() + ".curveCounts " + rodDagPath.fullPathName() + ".perRodParticleCounts;");
    dagModifier.commandToExecute( cmdStr );
    dagModifier.doIt();
}


void WmFigaroCmd::createRodShapeNode()
{
    if ( m_nurbsCurveList.isEmpty() )
    {
        MGlobal::displayError( "Please select a NURBS curves to create a rod from.\n" );
        return;
    }

    MStatus stat;
    // Create the rods node
    MObject rodTObj;  // Object for transform node
    MObject rodSObj;  // Object for shape node
    MDagPath shapeDagPath;
    MObject pObj;
    MDagModifier dagModifier;
    MString rodShapeName = "";

    // Set the control points of the rod shape to be the vertices of the input NURBS curve
    MDagPath dagPath;
    MObject component;
    for(int curve_idx = 0; curve_idx < m_nurbsCurveList.length(); curve_idx++)
    {

        createDagNode( WmFigaroRodShape::typeName.asChar(), WmFigaroRodShape::typeName.asChar(),
                       pObj, &rodTObj, &rodSObj, &dagModifier, rodShapeName );

        appendToResultString( rodShapeName );

        MDagPath rodDagPath;
        stat = MDagPath::getAPathTo( rodSObj, rodDagPath );
        CHECK_MSTATUS( stat );


        //m_nurbsCurveList.getDagPath( 0, dagPath, component );
        m_nurbsCurveList.getDagPath( curve_idx, dagPath, component );

        dagPath.extendToShape();

        MFnNurbsCurve curveFn( dagPath, &stat );
        CHECK_MSTATUS( stat );

        MPointArray nurbsPoints;
        curveFn.getCVs( nurbsPoints, MSpace::kWorld );

        MFnDependencyNode rodShapeFn( rodSObj );

        // Now get the input attribute of the rod and set
        MPlug controlPointsPlugArr = rodShapeFn.findPlug( "controlPoints", true, &stat );
        CHECK_MSTATUS( stat );

        for ( int p=0; p<nurbsPoints.length(); ++p )
        {
            MPlug cvPlug = controlPointsPlugArr.elementByLogicalIndex( p, &stat );
            CHECK_MSTATUS( stat );
            MPlug xPlug = cvPlug.child( 0, &stat );
            CHECK_MSTATUS( stat );
            MPlug yPlug = cvPlug.child( 1, &stat );
            CHECK_MSTATUS( stat );
            MPlug zPlug = cvPlug.child( 2, &stat );
            CHECK_MSTATUS( stat );
            xPlug.setValue( nurbsPoints[ p ].x );
            yPlug.setValue( nurbsPoints[ p ].y );
            zPlug.setValue( nurbsPoints[ p ].z );
        }

        //now connect the transformation of the rod and the curve
    //    MFnDependencyNode rodTransformFn( rodTObj, &stat);
    //    CHECK_MSTATUS( stat );
    //    MFnDagNode curveTransformFn( dagPath.transform(), &stat );
    //    CHECK_MSTATUS( stat );
    //    dagModifier.connect( curveTransformFn.findPlug("translate"),
    //            rodTransformFn.findPlug( "translate"));
    //    dagModifier.connect( curveTransformFn.findPlug("rotate"),
    //            rodTransformFn.findPlug( "rotate"));
    //    dagModifier.connect( curveTransformFn.findPlug("scale"),
    //            rodTransformFn.findPlug( "scale"));
    //
    //    stat = dagModifier.doIt();
    //    CHECK_MSTATUS( stat );



        // Now connect the output from the rod to the NURBS curve
        MPlug cvsPlug = rodShapeFn.findPlug( "controlVertex", & stat );
        CHECK_MSTATUS( stat );

        for( unsigned int cv = 0; cv < nurbsPoints.length(); cv++ )
        {
            MPlug cvPlug = cvsPlug.elementByLogicalIndex( cv, &stat );
            CHECK_MSTATUS( stat );

            stat = dagModifier.connect( cvPlug, curveFn.findPlug( "controlPoints" ).elementByLogicalIndex( cv, & stat ) );
            CHECK_MSTATUS( stat );
        }

    }
    dagModifier.doIt();

    // Tell the node to initialise itself. Seems like pulling on a plug would be nice
    /*WmFigaroRodShape* wmFigaroRodShape = dynamic_cast< WmFigaroRodShape* >( rodShapeFn.userNode( &stat ) );
    CHECK_MSTATUS( stat );

    wmFigaroRodShape->resetSimulation();*/

}

void WmFigaroCmd::getNodes( MSelectionList i_opt_nodes )
{
    // Prepare all the selection lists that this command operates on

    MDagPath mDagPath;
    MObject mObj;
    MItDag dagit;
    MStatus stat;

    m_nurbsCurveList.clear();
    m_fozzieNodeList.clear();
    m_meshList.clear();
    m_figRodNodeList.clear();
    m_allOtherTransformNodesList.clear();

    for (MItSelectionList sIt(i_opt_nodes); !sIt.isDone(); sIt.next())
    {
        sIt.getDagPath(mDagPath, mObj);
        
        MFnDagNode dagFn( mDagPath.child(0, &stat), &stat);
        CHECK_MSTATUS( stat );
        MDagPath childPath;
        stat = dagFn.getPath( childPath );
        CHECK_MSTATUS( stat );
        childPath.extendToShape();

        //////////////////
        //
        // First check for Maya in built nodes
        
        if ( childPath.apiType() == MFn::kNurbsCurve )
        {
            mObj = childPath.node();
            
            m_nurbsCurveList.add( childPath, mObj, false);
        }
        else if ( childPath.apiType() == MFn::kMesh )
        {
            // First we add it to the list of all other objects since the user may just be wanting
            // to use this mesh to attach to ta rod.
            stat = m_allOtherTransformNodesList.add( mDagPath, mObj, false );
            CHECK_MSTATUS( stat );
            
            mObj = childPath.node();
            
            m_meshList.add( childPath, mObj, false);
        } 
        else
        {
            /////////////////////
            //
            // Now check for user created plugin nodes
            
            MFnDependencyNode nodeFn( childPath.node( &stat ) );
            CHECK_MSTATUS( stat );
            
            if ( nodeFn.typeName() == WmBunsenNode::typeName )
            {
                m_selectedwmBunsenNode = childPath.node();
            }
            else if ( nodeFn.typeName() == "wmBarbFurSetNode" )
            {
                mObj = childPath.node();
                m_fozzieNodeList.add( childPath, mObj, false );
            } 
            else if ( nodeFn.typeName() == WmFigRodNode::typeName )
            {
                mObj = childPath.node();
                m_figRodNodeList.add( childPath, mObj, false );
            } 
            else
            {
                stat = m_allOtherTransformNodesList.add( mDagPath, mObj, false );
                CHECK_MSTATUS( stat );
            }
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

void WmFigaroCmd::quaternionFromMatrix( MMatrix& a, MQuaternion& q ) 
{
  double trace = a( 0 ,0 ) + a( 1, 1 ) + a( 2, 2 );
  if( trace > 0 ) 
  {
    double s = 0.5f / sqrtf( trace + 1.0f );
    q.w = 0.25f / s;
    q.x = ( a(2,1) - a(1,2) ) * s;
    q.y = ( a(0,2) - a(2,0) ) * s;
    q.z = ( a(1,0) - a(0,1) ) * s;
  } else {
    if ( a(0,0) > a(1,1) && a(0,0) > a(2,2) ) {
      double s = 2.0f * sqrtf( 1.0f + a(0,0) - a(1,1) - a(2,2));
      q.w = (a(2,1) - a(1,2) ) / s;
      q.x = 0.25f * s;
      q.y = (a(0,1) + a(1,0) ) / s;
      q.z = (a(0,2) + a(2,0) ) / s;
    } else if (a(1,1) > a(2,2)) {
      double s = 2.0f * sqrtf( 1.0f + a(1,1) - a(0,0) - a(2,2));
      q.w = (a(0,2) - a(2,0) ) / s;
      q.x = (a(0,1) + a(1,0) ) / s;
      q.y = 0.25f * s;
      q.z = (a(1,2) + a(2,1) ) / s;
    } else {
      double s = 2.0f * sqrtf( 1.0f + a(2,2) - a(0,0) - a(1,1) );
      q.w = (a(1,0) - a(0,1) ) / s;
      q.x = (a(0,2) + a(2,0) ) / s;
      q.y = (a(1,2) + a(2,1) ) / s;
      q.z = 0.25f * s;
    }
  }
}


void WmFigaroCmd::attatchEdgeToObject()
{
    // Find the first selected object
    MDagPath mDagPath;
    MObject mObj;
    MItDag dagit;
    MStatus stat;

    if ( m_allOtherTransformNodesList.isEmpty() )
    {
        displayError( "Please select an object to attach to the rod edge." );
        return;
    }
    if ( m_figRodNodeList.isEmpty() )
    {
        displayError( "Please select a Figaro rod node to attach to the object." );
        return;
    }
    
    // Currently we just force it to be controlling edge 0 of the rod
    // Which edge an object controls is defined by the index of the array
    // it connects into on the rod node.

    MObject connectionNodeTObj;  // Object for transform node
    MObject connectionNodeSObj;  // Object for shape node
    MObject pObj;
    MString connectionShapeName;
    MDagModifier dagModifier;
            
    createDagNode( WmFigConnectionNode::typeName.asChar(), 
                   WmFigConnectionNode::typeName.asChar(), 
                   pObj, &connectionNodeTObj, &connectionNodeSObj, &dagModifier,
                   connectionShapeName );
    
    appendToResultString( connectionShapeName );
    
    MDagPath connnectionNodePath;
    stat = MDagPath::getAPathTo( connectionNodeSObj, connnectionNodePath );
    CHECK_MSTATUS( stat );

    /*
    MPlug connectionEdgePlug( connectionNodeSObj, WmFigConnectionNode::ia_edgeNumber );
    CHECK_MSTATUS( stat );
    stat = connectionEdgePlug.setValue( m_edgeNumber );
    CHECK_MSTATUS( stat );
    
    MPlug connectionRodPlug( connectionNodeSObj, WmFigConnectionNode::ia_rodNumber );
    CHECK_MSTATUS( stat );
    stat = connectionEdgePlug.setValue( m_rodNumber );
    CHECK_MSTATUS( stat );
    */
    
    // The above code appears to work but in Maya the attribute is never updated, the below works
    // but I am unsure why one does and the other doesn't.
    MFnDependencyNode connectionNodeFn( connectionNodeSObj, &stat );
    CHECK_MSTATUS( stat );
    MPlug edgeNumPlug = connectionNodeFn.findPlug( "edgeNumber", true, &stat );
    CHECK_MSTATUS( stat );
    edgeNumPlug.setValue( m_edgeNumber );
    MPlug rodNumPlug = connectionNodeFn.findPlug( "rodNumber", true, &stat );
    CHECK_MSTATUS( stat );
    rodNumPlug.setValue( m_rodNumber );
    
    MPlug connectionTransformPlug( connectionNodeSObj, WmFigConnectionNode::ia_transformMatrix );
    CHECK_MSTATUS( stat );
    
    MObject transformNodeObj;
    stat = m_allOtherTransformNodesList.getDependNode( 0, transformNodeObj );
    CHECK_MSTATUS( stat );
    
    MFnDependencyNode transformDepFn( transformNodeObj, &stat );
    CHECK_MSTATUS( stat );
    
    MPlug inputTransformPlug = transformDepFn.findPlug( "xformMatrix", true, &stat );
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( inputTransformPlug, connectionTransformPlug );
    CHECK_MSTATUS( stat );
    stat = dagModifier.doIt();
    CHECK_MSTATUS( stat );
    
    // Connect edge output connection to input on connection node so the edges can drive
    // meshes.
    
    MObject figRodNodeObj;
    stat = m_figRodNodeList.getDependNode( 0, figRodNodeObj );
    CHECK_MSTATUS( stat );
    
    MFnDependencyNode figRodNodeDepFn( figRodNodeObj, &stat );
    CHECK_MSTATUS( stat );
    
    MPlug outEdgeTransformPlug = figRodNodeDepFn.findPlug( "outEdgeTransforms", true, &stat );
    CHECK_MSTATUS( stat );
    
    MPlug connectionEdgeInputPlug( connectionNodeSObj, WmFigConnectionNode::ia_rodEdgeTransforms );
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( outEdgeTransformPlug, connectionEdgeInputPlug );
    CHECK_MSTATUS( stat );
    stat = dagModifier.doIt();
    CHECK_MSTATUS( stat );
    
    MPlug outputTransformPlug( connectionNodeSObj, WmFigConnectionNode::oa_outTransformMatrix );
    CHECK_MSTATUS( stat );
    
    MObject matrixObj;
    stat = outputTransformPlug.getValue( matrixObj );
    CHECK_MSTATUS( stat );
    
    MFnMatrixData matrixData( matrixObj, &stat );
    CHECK_MSTATUS( stat );
    
    // Set the input objects transform to match this data so that it lies on the edge before
    // it starts controlling the edge. Unfortunatley we can't set a matrix directly for some
    // reason we need to set a quaternion.
    
    MPxTransformationMatrix transformMatrix( matrixData.matrix() );
    MVector translation = transformMatrix.translation();
    
    // Now set the translation part of the transform
    MPlug xTranslatePlug = transformDepFn.findPlug( "translateX", true, &stat );
    CHECK_MSTATUS( stat );
    xTranslatePlug.setValue( translation[ 0 ] );
    MPlug yTranslatePlug = transformDepFn.findPlug( "translateY", true, &stat );
    CHECK_MSTATUS( stat );
    yTranslatePlug.setValue( translation[ 1 ] );
    MPlug zTranslatePlug = transformDepFn.findPlug( "translateZ", true, &stat );
    CHECK_MSTATUS( stat );
    zTranslatePlug.setValue( translation[ 2 ] );
    
    // And the rotation which is nastily in euler angles...
    
    MEulerRotation eulerRotation = transformMatrix.eulerRotation();
    
    MPlug xRotatePlug = transformDepFn.findPlug( "rotateX", true, &stat );
    CHECK_MSTATUS( stat );
    xRotatePlug.setValue( eulerRotation.x );
    MPlug yRotatePlug = transformDepFn.findPlug( "rotateY", true, &stat );
    CHECK_MSTATUS( stat );
    yRotatePlug.setValue( eulerRotation.y );
    MPlug zRotatePlug = transformDepFn.findPlug( "rotateZ", true, &stat );
    CHECK_MSTATUS( stat );
    zRotatePlug.setValue( eulerRotation.z );
    
    // Now the input object is sitting on the correct segment of the rod and oriented correctly.
    // It is ready to control that rod segment so connect the connection node to the rod node
    // so it knows the segment is being controlled.
    
    MPlug edgeTransformPlugArr = figRodNodeDepFn.findPlug( "edgeTransforms", true, &stat );
    CHECK_MSTATUS( stat );
    
    MPlug edgeTransformPlug = edgeTransformPlugArr.elementByLogicalIndex( edgeTransformPlugArr.numElements(), &stat );
    CHECK_MSTATUS( stat );
    
    MPlug connectionEdgeTransformPlug( connectionNodeSObj, WmFigConnectionNode::oa_edgeTransform );
    CHECK_MSTATUS( stat );
    stat = dagModifier.connect( connectionEdgeTransformPlug, edgeTransformPlug );
    CHECK_MSTATUS( stat );
    stat = dagModifier.doIt();
    CHECK_MSTATUS( stat );
}

void WmFigaroCmd::addCollisionMeshes()
{
    MStatus stat;
    
    if ( m_selectedwmBunsenNode == MObject::kNullObj )
    {
       MGlobal::displayError( "Please select a wmFigaro node to connect the mesh to." );
       return;
    }
    
    MFnDependencyNode bunsenNodeFn( m_selectedwmBunsenNode, &stat );
    CHECK_MSTATUS( stat );
    
    MPlug bunsenInputPlugArr = bunsenNodeFn.findPlug( "collisionMeshes", true, &stat );
    CHECK_MSTATUS( stat );
    
    for ( uint l=0; l<m_meshList.length(); l++ )
    {
        MDagPath dagPath;
        MObject component;
        m_meshList.getDagPath( l, dagPath, component );
        dagPath.extendToShape();

        MObject nodeObj = dagPath.node( &stat );
        CHECK_MSTATUS( stat );
        MFnDependencyNode nodeFn( nodeObj, &stat );

        // Create a collisionMeshNode to pipe the mesh through on the way to the dynamics node
    
        MObject collisionMeshNodeTObj;  // Object for transform node
        MObject collisionMeshNodeSObj;  // Object for shape node
        MDagPath shapeDagPath;    
        MObject pObj;
        MDagModifier dagModifier;
        MString collisionNodeShapeName = "";

        createDagNode( WmBunsenCollisionMeshNode::typeName.asChar(), 
                       WmBunsenCollisionMeshNode::typeName.asChar(), 
                       pObj, &collisionMeshNodeTObj, &collisionMeshNodeSObj, &dagModifier,
                       collisionNodeShapeName );
        
        appendToResultString( collisionNodeShapeName );
                          
        MDagPath collisionMeshNodeDagPath;
        stat = MDagPath::getAPathTo( collisionMeshNodeSObj, collisionMeshNodeDagPath );
        CHECK_MSTATUS( stat );

        MString timeStr( "connectAttr -f time1.outTime " + collisionMeshNodeDagPath.fullPathName() + ".time" );
        dagModifier.commandToExecute( timeStr );
        dagModifier.doIt();   
            
        if ( nodeFn.typeName() == "mesh" )
        {
            MFnMesh meshFn( dagPath, &stat );
            CHECK_MSTATUS( stat );

            MPlug worldMeshPlug = meshFn.findPlug( "worldMesh", true, &stat ).elementByLogicalIndex( 0, &stat );
            CHECK_MSTATUS( stat );

            MPlug collisionMeshNodeInPlug( collisionMeshNodeSObj, WmBunsenCollisionMeshNode::ia_inMesh );
            
            stat = dagModifier.connect( worldMeshPlug, collisionMeshNodeInPlug );
            CHECK_MSTATUS( stat );
            stat = dagModifier.doIt();
            CHECK_MSTATUS( stat );    

            MPlug collisionMeshNodeOutPlug( collisionMeshNodeSObj, WmBunsenCollisionMeshNode::oa_meshData );
            
            unsigned int numElements = bunsenInputPlugArr.numElements( &stat );
            CHECK_MSTATUS( stat );
            MPlug bunsenMeshPlug = bunsenInputPlugArr.elementByLogicalIndex( numElements );
            CHECK_MSTATUS( stat );
            if ( stat.error() )
            {
                MGlobal::displayError( stat.errorString() );
                return;
            }
            stat = dagModifier.connect( collisionMeshNodeOutPlug, bunsenMeshPlug );
            CHECK_MSTATUS( stat );
            stat = dagModifier.doIt();
            CHECK_MSTATUS( stat );
            
            // Connect startTime so the rods reset when the Figaro node resets
            MPlug collisionNodeStartTimePlug( collisionMeshNodeSObj, WmBunsenCollisionMeshNode::ia_startTime );
            CHECK_MSTATUS( stat );
            MPlug bunsenNodeStartTimePlug = bunsenNodeFn.findPlug( "startTime", true, &stat ); 
            CHECK_MSTATUS( stat );
            stat = dagModifier.connect( bunsenNodeStartTimePlug, collisionNodeStartTimePlug );
            CHECK_MSTATUS( stat );
            stat = dagModifier.doIt();
            CHECK_MSTATUS( stat );
        }
    }
}

void WmFigaroCmd::createPreviewNodes()
{
    MStatus stat;
    
    if ( m_cacheFile == "" )
    {
        MGlobal::displayError( "Please speicify a cache file using the -cacheFile flags" );
        return;
    }
    // Set this to null so the rod creation function builds us a bunsen node
    m_selectedwmBunsenNode = MObject::kNullObj;
    MObject rodNode;
    createWmFigRodNode( false, true, &rodNode );
    
    MFnDependencyNode rodNodeFn( rodNode );
    MPlug cachePathPlug( rodNodeFn.findPlug( "cachePath", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = cachePathPlug.setValue( m_cacheFile );
    CHECK_MSTATUS( stat );
    MPlug readFromCachePlug( rodNodeFn.findPlug( "readFromCache", true, &stat ) );
    CHECK_MSTATUS( stat );
    stat = readFromCachePlug.setValue( true );
    CHECK_MSTATUS( stat );
}

void WmFigaroCmd::p_AddFlag(
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
    CHECK_MSTATUS( i_mSyntax.addFlag(i_shortName, i_longName, i_argType1, 
        i_argType2, i_argType3, i_argType4, i_argType5, i_argType6) );
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



void WmFigaroCmd::printHelp()
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

void WmFigaroCmd::getSelectedCurves( const MSelectionList &i_selectionList,
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

MStatus WmFigaroCmd::createDagNode( const char *transformName, const char *nodeType, 
    MObject &parentObj, MObject *transformObjP, MObject *shapeObjP, MDagModifier *iDagModifier,
    MString& o_shapeName )
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
                shapeName = transformNameStr.substring(0, int(cp - sp));
              else
                shapeName = transformNameStr;
    
              shapeName += "Shape";
    
              if (cp != ep)
                shapeName += transformNameStr.substring(int(cp - sp + 1), int(ep - sp));
    
              rStatus = mDagModifier.renameNode(shapeObj, shapeName);
              if (rStatus)
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

void WmFigaroCmd::addVertexConstraint( int i_rodInd, double i_dis, int i_vertId, MPoint i_pos )
{
    //get selected rod node
    MObject rodNode = MObject::kNullObj;

    MStatus status;
    MSelectionList list;
    MGlobal::getActiveSelectionList( list );

    for( unsigned int l = 0; l < list.length(); l++ )
    {
        MDagPath path;
        list.getDagPath( l, path );
        path.extendToShape();

        MFnDependencyNode fn( path.node(), &status );
        CHECK_MSTATUS( status );

        if( fn.typeId() == WmFigRodNode::typeID )
        {
            rodNode = path.node();
            break;
        }
    }

    if( rodNode != MObject::kNullObj )
    {
        MFnDependencyNode rodFn( rodNode );
        WmFigRodNode *rodNode = static_cast< WmFigRodNode*>( rodFn.userNode());

        WmFigRodGroup*  rodGroup = rodNode->rodGroup();

        if( i_rodInd < rodGroup->numberOfRods() )
        {
            if( i_vertId < 0 )
                i_vertId = rodGroup->elasticRod( i_rodInd )->nv() -1 ;
            
            if( i_vertId < rodGroup->elasticRod( i_rodInd )->nv() )
            {
                BASim::Vec3d target_position = BASim::Vec3d( i_pos.x, i_pos.y, i_pos.z );
                rodGroup->collisionStepper( i_rodInd )->setVertexPositionPenalty(i_vertId, target_position, 50. , 1);
            }
            else
            {
                displayError( "Invalid input vertex id!" );
            }
        }
        else
        {
            displayError( "Invalid input rod id! " );
        }

    }
    else
    {
        displayError( "Please select a figaro rod node! " );
    }

}

