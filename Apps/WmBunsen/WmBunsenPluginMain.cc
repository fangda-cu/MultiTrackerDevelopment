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

#include "Sweeney/WmSweeneyNode.hh"
#include "Sweeney/WmSweeneySubsetNode.hh"
#include "Sweeney/WmSweeneyVolumetricNode.hh"
#include "Sweeney/Tools/AddRod/WmSwAddRodContext.hh"
#include "Sweeney/Tools/AddRod/WmSwAddRodContextCommand.hh"
#include "Sweeney/Tools/AddRod/WmSwAddRodToolCommand.hh"
#include "Sweeney/WmSweeneyCmd.hh"

using namespace sweeney;

MStatus initializePlugin( MObject obj )
{
	MStatus stat;
	MString errStr;
	char szVersion[1024];
  	sprintf( szVersion,"build %s %s",__DATE__,__TIME__ );

  	MFnPlugin plugin( obj,"Alasdair Coull", szVersion, "Any" );

    MGlobal::startErrorLogging("/tmp/MayaProblemLog.txt");
    cerr << "WmFigaro loading....Maya is logging problems to " << MGlobal::errorLogPathName() << endl;

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
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // 
    // Sweeeny initialisation
    // 
    ///////////////////////////////////////////////////////////////////////////////////////////////

    MStatus status;
    
    status = plugin.registerNode( WmSweeneyNode::typeName, WmSweeneyNode::typeID,
                                WmSweeneyNode::creator,
                                WmSweeneyNode::initialize,
                                WmSweeneyNode::kLocatorNode );
    if ( !status )
    {
        status.perror( "RegisterNode WmSweeneyNode failed" );
        return status;
    }
    
    status = plugin.registerNode( WmSweeneySubsetNode::typeName, WmSweeneySubsetNode::typeID,
                                    WmSweeneySubsetNode::creator,
                                    WmSweeneySubsetNode::initialize,
                                    WmSweeneySubsetNode::kLocatorNode );
    if ( !status )
    {
        status.perror( "RegisterNode WmSweeneySubsetNode failed" );
    }

    status = plugin.registerNode( WmSweeneyVolumetricNode::typeName, WmSweeneyVolumetricNode::typeID,
                                        WmSweeneyVolumetricNode::creator,
                                        WmSweeneyVolumetricNode::initialize,
                                        WmSweeneyVolumetricNode::kLocatorNode );
    if ( !status )
    {
        status.perror( "RegisterNode WmSweeneyVolumetricNode failed" );
    }

    status = plugin.registerCommand( WmSweeneyCmd::typeName, WmSweeneyCmd::creator, 
                                   WmSweeneyCmd::syntaxCreator );
    if ( !status ) {
        status.perror( "registerCommand wmSweeney failed" );
        return status;     
    }
    
    if ( plugin.registerContextCommand( WmSwAddRodContext::typeName,
            WmSwAddRodContextCommand::creator,
            WmSwAddRodToolCommand::typeName,
            WmSwAddRodToolCommand::creator ) != MS::kSuccess )
    if ( !status ) 
    {
        status.perror( "registerContextCommand WmSwAddRodContext failed" );
        return status;
    }
    
    MGlobal::executeCommand( "source wmSweeney.mel", false );
    
    // This will trash the Figaro UI registration, so we just call the add and remove menu functions
    // directly.
    //CHECK_MSTATUS( plugin.registerUI( "wmSweeneyAddMainMenu", "wmSweeneyRemoveMainMenu" ) );
    MGlobal::executeCommand( "wmSweeneyAddMainMenu", false );

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

    stat = plugin.deregisterCommand( WmSweeneyCmd::typeName );
    if (!stat) 
    {
        stat.perror( "deregister command WmSweeneyCmd failed" );
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
    
    ///////////////////////////////////////////////////////////////////////////////////////////////
    // 
    // Sweeeny uninitialisation
    // 
    ///////////////////////////////////////////////////////////////////////////////////////////////

    MStatus status;

    if ( plugin.deregisterNode( WmSweeneyNode::typeID ) != MS::kSuccess )
    {
        status.perror( "DeregisterNode WmSweeneyNode failed" );
    }
    
    if ( plugin.deregisterNode( WmSweeneySubsetNode::typeID ) != MS::kSuccess )
    {
        status.perror( "DeregisterNode WmSweeneySubsetNode failed" );
    }

    if ( plugin.deregisterNode( WmSweeneyVolumetricNode::typeID ) != MS::kSuccess )
    {
        status.perror( "DeregisterNode WmSweeneyVolumetricNode failed" );
    }


    if ( plugin.deregisterContextCommand( WmSwAddRodContext::typeName, WmSwAddRodToolCommand::typeName ) != MS::kSuccess )
    {
        status.perror( "Deregister context command WmSwAddRodContext failed" );
    }
    
    // Again, explicitly remove the menu as we can't use the add/removeUI functionality due to
    // it being in the figaro code and I don't want to mix the two.
    MGlobal::executeCommand( "wmSweeneyRemoveMainMenu", false );
    
    MGlobal::stopErrorLogging();

	return stat;
}
