#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>

#include "WmFigRodNode.hh"
#include "WmBunsenNode.hh"
#include "WmBunsenCollisionMeshNode.hh"
#include "WmFigaroCmd.hh"
#include "WmFigConnectionNode.hh"
#include "WmFigSelectionContext.hh"
#include "WmFigSelectionContextCommand.hh"
#include "WmFigSelectionToolCommand.hh"
#include "WmFigControllerContext.hh"
#include "WmFigControllerContextCommand.hh"
#include "WmFigControllerToolCommand.hh"
//#include "WmFigCombContext.hh"
//#include "WmFigCombContextCommand.hh"
//#include "WmFigCombToolCommand.hh"
#include "WmFigaroRodShape/WmFigaroRodShape.hh"
#include "WmFigaroRodShape/WmFigaroRodShapeUI.hh"
#include "WmFigaroRodShape/WmFigaroRodShapeIterator.hh"
//#include "constraints/WmFigConstraintNode.hh"
#include "WmFigSelectionDisplayNode.hh"

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
    
    stat = plugin.registerNode( WmFigRodNode::typeName, WmFigRodNode::typeID,
                                WmFigRodNode::creator,
                                WmFigRodNode::initialize,
                                WmFigRodNode::kLocatorNode );
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
/*
    stat = plugin.registerNode( WmFigConstraintNode::TypeName,
								WmFigConstraintNode::TypeId,
								WmFigConstraintNode::creator,
								WmFigConstraintNode::initialize,
								MPxNode::kLocatorNode);
    if ( !stat )
    {
        stat.perror( "RegisterNode WmFigConstraintNode failed" );
        return stat;
    }*/

    stat = plugin.registerNode( WmFigSelectionDisplayNode::TypeName,
    							WmFigSelectionDisplayNode::TypeId,
    							WmFigSelectionDisplayNode::creator,
    							WmFigSelectionDisplayNode::initialize,
    							WmFigConnectionNode::kLocatorNode);
    if ( !stat )
    {
        stat.perror( "RegisterNode WmFigSelectionDisplayNode failed" );
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
    
    if ( plugin.registerContextCommand( WmFigSelectionContext::typeName,
            WmFigSelectionContextCommand::creator,
            WmFigSelectionToolCommand::typeName,
            WmFigSelectionToolCommand::creator ) != MS::kSuccess )
    if ( !stat ) {
        stat.perror( "registerContextCommand WmFigSelectionContext failed" );
        return stat;     
    }

    if ( plugin.registerContextCommand( WmFigControllerContext::typeName,
            WmFigControllerContextCommand::creator,
            WmFigControllerToolCommand::typeName,
            WmFigControllerToolCommand::creator ) != MS::kSuccess )
    if ( !stat ) {
        stat.perror( "registerContextCommand WmFigControllerContext failed" );
        return stat;     
    }

    /*stat = plugin.registerContextCommand( WmFigCombContext::typeName,
            WmFigCombContextCommand::creator,
            WmFigCombToolCommand::typeName,
            WmFigCombToolCommand::creator );
    if ( !stat ) {
        stat.perror( "registerContextCommand WmFigCombContext failed" );
        return stat;     
    }*/

    stat = plugin.registerShape( WmFigaroRodShape::typeName, WmFigaroRodShape::id,
                                   &WmFigaroRodShape::creator,
                                   &WmFigaroRodShape::initialize,
                                   &WmFigaroRodShapeUI::creator );
    if ( !stat ) {
        stat.perror( "registerShape WmFigaroRodShape failed" );
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
    
    stat = plugin.deregisterNode( WmFigRodNode::typeID );
    if( !stat ) 
    {
        stat.perror( "DeregisterNode WmFigRodNode failed" );
    }
    
    stat = plugin.deregisterNode( WmBunsenCollisionMeshNode::typeId );
    if( !stat ) 
    {
        stat.perror( "DeregisterNode WmFigCollisionNode failed" );
    }
    /*
    stat = plugin.deregisterNode( WmFigConstraintNode::TypeId );
    if( !stat )
    {
        stat.perror( "DeregisterNode WmFigConstraintNode failed" );
    }*/

    stat = plugin.deregisterNode( WmFigSelectionDisplayNode::TypeId );
    if( !stat )
    {
        stat.perror( "DeregisterNode WmFigSelectionDisplayNode failed" );
    }

    stat = plugin.deregisterNode( WmFigConnectionNode::typeID );
    if( !stat ) 
    {
        stat.perror( "DeregisterNode WmFigConnectionNode failed" );
    }
    
    // Deregister custom commands
    stat = plugin.deregisterCommand( WmFigaroCmd::typeName );
    if (!stat) 
    {
        stat.perror( "deregister command wmFigaro failed" );
    }

    if ( plugin.deregisterContextCommand( WmFigSelectionContext::typeName, WmFigSelectionToolCommand::typeName ) != MS::kSuccess )
    {
        stat.perror( "deregister context command wmFigSelection failed" );
    }

    if ( plugin.deregisterContextCommand( WmFigControllerContext::typeName, WmFigControllerToolCommand::typeName ) != MS::kSuccess )
    {
        stat.perror( "deregister context command wmFigController failed" );
    }

    /*if ( plugin.deregisterContextCommand( WmFigCombContext::typeName, WmFigCombToolCommand::typeName ) != MS::kSuccess )
    {
        stat.perror( "deregister context command wmFigComb failed" );
    }*/

    if ( plugin.deregisterNode( WmFigaroRodShape::id ) != MS::kSuccess )
    {
        stat.perror( "deregisterNode WmFigaroRodShape failed" );
    }

    MGlobal::stopErrorLogging();

	return stat;
}
