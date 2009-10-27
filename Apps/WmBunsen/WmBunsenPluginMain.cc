#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>

#include "WmBunsenNode.hh"

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
        stat.perror("DeregisterNode WmBunsenNode failed");        
    }

    MGlobal::stopErrorLogging();

	return stat;	
}
