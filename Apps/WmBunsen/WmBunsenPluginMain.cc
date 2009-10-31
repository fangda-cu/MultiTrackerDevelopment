#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>

#include "WmBunsenRodNode.hh"
#include "WmBunsenNode.hh"
#include "WmBunsenCmd.hh"

MStatus initializePlugin( MObject obj )
{ 
	MStatus stat;
	MString errStr;
	char szVersion[1024];
  	sprintf( szVersion,"build %s %s",__DATE__,__TIME__ );

  	MFnPlugin plugin( obj,"Alasdair Coull", szVersion, "Any" );

    MGlobal::startErrorLogging("/tmp/MayaErrorLog.txt");
    cerr << "WmBunsen loading....Maya is logging errors to " << MGlobal::errorLogPathName() << endl;

    stat = plugin.registerNode( WmBunsenNode::ia_typeName, WmBunsenNode::ia_typeID,
                                WmBunsenNode::creator,
                                WmBunsenNode::initialize,
                                WmBunsenNode::kLocatorNode );
    if ( !stat )
    {
        stat.perror( "RegisterNode WmBunsenNode failed" );
        return stat;
    }
    
    stat = plugin.registerNode( WmBunsenRodNode::ia_typeName, WmBunsenRodNode::ia_typeID,
                                WmBunsenRodNode::creator,
                                WmBunsenRodNode::initialize,
                                WmBunsenRodNode::kLocatorNode );
    if ( !stat )
    {
        stat.perror( "RegisterNode WmBunsenRodNode failed" );
        return stat;
    }
    
    stat = plugin.registerCommand( WmBunsenCmd::typeName, WmBunsenCmd::creator, 
                                   WmBunsenCmd::syntaxCreator );
    if ( !stat ) {
        stat.perror( "registerCommand wmBunsen failed" );
        return stat;     
    }

    MGlobal::executeCommand( "source WmBunsen.mel", false );
    CHECK_MSTATUS( plugin.registerUI( "wmBunsenAddMainMenu", "wmBunsenRemoveMainMenu" ) );
    return stat;

    
    return MS::kSuccess;
}

MStatus uninitializePlugin( MObject obj)
{
	MStatus stat;
	MString errStr;
	MFnPlugin plugin( obj );
	
    //
    // Deregister nodes
	stat = plugin.deregisterNode( WmBunsenNode::ia_typeID );
    if( !stat ) 
    {
        stat.perror( "DeregisterNode WmBunsenNode failed" );
    }
    
    stat = plugin.deregisterNode( WmBunsenRodNode::ia_typeID );
    if( !stat ) 
    {
        stat.perror( "DeregisterNode WmBunsenRodNode failed" );
    }
    
    // Deregister custom commands
    stat = plugin.deregisterCommand( WmBunsenCmd::typeName );
    if (!stat) {
        stat.perror( "deregister command wmBunsen failed" );
    }

    MGlobal::stopErrorLogging();

	return stat;
}
