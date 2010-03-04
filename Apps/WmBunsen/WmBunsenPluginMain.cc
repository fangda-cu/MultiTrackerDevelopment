#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>

#include "WmBunsenRodNode.hh"
#include "WmBunsenNode.hh"
#include "WmBunsenCollisionMeshNode.hh"
#include "WmFigaroCmd.hh"
#include "WmFigConnectionNode.hh"

MStatus initializePlugin( MObject obj )
{ 
	MStatus stat;
	MString errStr;
	char szVersion[1024];
  	sprintf( szVersion,"build %s %s",__DATE__,__TIME__ );

  	MFnPlugin plugin( obj,"Alasdair Coull", szVersion, "Any" );

    MGlobal::startErrorLogging("/tmp/MayaErrorLog.txt");
    cerr << "WmFigaro loading....Maya is logging errors to " << MGlobal::errorLogPathName() << endl;

    stat = plugin.registerNode( WmBunsenNode::typeName, WmBunsenNode::typeID,
                                WmBunsenNode::creator,
                                WmBunsenNode::initialize,
                                WmBunsenNode::kLocatorNode );
    if ( !stat )
    {
        stat.perror( "RegisterNode WmFigaroNode failed" );
        return stat;
    }
    
    stat = plugin.registerNode( WmBunsenRodNode::typeName, WmBunsenRodNode::typeID,
                                WmBunsenRodNode::creator,
                                WmBunsenRodNode::initialize,
                                WmBunsenRodNode::kLocatorNode );
    if ( !stat )
    {
        stat.perror( "RegisterNode WmFigRodNode failed" );
        return stat;
    }

    stat = plugin.registerNode( WmBunsenCollisionMeshNode::typeName, WmBunsenCollisionMeshNode::typeId,
                                WmBunsenCollisionMeshNode::creator,
                                WmBunsenCollisionMeshNode::initialize,
                                WmBunsenCollisionMeshNode::kLocatorNode );
    if ( !stat )
    {
        stat.perror( "RegisterNode WmFigCollisionNode failed" );
        return stat;
    }

    stat = plugin.registerNode( WmFigConnectionNode::typeName, WmFigConnectionNode::typeID,
                                WmFigConnectionNode::creator,
                                WmFigConnectionNode::initialize,
                                WmFigConnectionNode::kLocatorNode );
    if ( !stat )
    {
        stat.perror( "RegisterNode WmFigConnectionNode failed" );
        return stat;
    }
    
    stat = plugin.registerCommand( WmFigaroCmd::typeName, WmFigaroCmd::creator, 
                                   WmFigaroCmd::syntaxCreator );
    if ( !stat ) {
        stat.perror( "registerCommand wmFigaro failed" );
        return stat;     
    }

    MGlobal::executeCommand( "source WmFigaro.mel", false );
    CHECK_MSTATUS( plugin.registerUI( "wmFigaroAddMainMenu", "wmFigaroRemoveMainMenu" ) );
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
	stat = plugin.deregisterNode( WmBunsenNode::typeID );
    if( !stat ) 
    {
        stat.perror( "DeregisterNode WmFigaroNode failed" );
    }
    
    stat = plugin.deregisterNode( WmBunsenRodNode::typeID );
    if( !stat ) 
    {
        stat.perror( "DeregisterNode WmFigRodNode failed" );
    }
    
    stat = plugin.deregisterNode( WmBunsenCollisionMeshNode::typeId );
    if( !stat ) 
    {
        stat.perror( "DeregisterNode WmFigCollisionNode failed" );
    }
    
    stat = plugin.deregisterNode( WmFigConnectionNode::typeID );
    if( !stat ) 
    {
        stat.perror( "DeregisterNode WmFigConnectionNode failed" );
    }
    
    // Deregister custom commands
    stat = plugin.deregisterCommand( WmFigaroCmd::typeName );
    if (!stat) {
        stat.perror( "deregister command wmFigaro failed" );
    }

    MGlobal::stopErrorLogging();

	return stat;
}
