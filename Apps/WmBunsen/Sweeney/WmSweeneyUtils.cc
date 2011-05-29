#include "WmSweeneyUtils.hh"

#include <maya/MDagPath.h>

namespace sweeney {
namespace utils {

MStatus findASelectedNodeByTypeName( MString& i_typeName, MObject* o_selectedNodeObject, 
    MDagPath* o_selectedNodDagpath )
{
    // Only returns the first selected node of the type it finds.
    
    MStatus status;
    
    MSelectionList selectionList;
    MGlobal::getActiveSelectionList( selectionList );

    *o_selectedNodeObject = MObject::kNullObj;
    
    MDagPath nodeDagPath;
    MObject nodeObject;
    
    for ( unsigned int s = 0; s < selectionList.length(); ++s )
    {
        selectionList.getDagPath( s, nodeDagPath, nodeObject );
        
        // Look for a child as the user probably selected the transform
        MFnDagNode dagNodeFn( nodeDagPath.child( 0, &status ), &status );
        if ( status )
        {
            MDagPath childPath;
            
            status = dagNodeFn.getPath( childPath );
            CHECK_MSTATUS( status );
            childPath.extendToShape();

            MFnDependencyNode dependencyNodeFn( childPath.node( &status ) );
            CHECK_MSTATUS( status );

            if ( dependencyNodeFn.typeName() == i_typeName )
            {
                *o_selectedNodeObject = childPath.node( &status );
                CHECK_MSTATUS( status );
                
                if ( o_selectedNodDagpath != NULL )
                {
                    *o_selectedNodDagpath = childPath;
                }
                
                return MStatus::kSuccess;
            }       
        }
        else // Perhaps no child as the user selected the shape node directly
        {
            MFnDependencyNode dependencyNodeFn( nodeDagPath.node( &status ) );
            CHECK_MSTATUS( status );

            if ( dependencyNodeFn.typeName() == i_typeName )
            {
                *o_selectedNodeObject = nodeDagPath.node( &status );
                CHECK_MSTATUS( status );
                
                if ( o_selectedNodDagpath != NULL )
                {
                    *o_selectedNodDagpath = nodeDagPath;
                }
                
                return MStatus::kSuccess;
            }        
        }
    }
    
    return MStatus::kFailure;
}
    
MStatus findSelectedSweeneyNodeAndRodManager( WmSweeneyNode* o_wmSweeneyNode, 
    WmSweeneyRodManager* o_wmSweenyRodManager, MDagPath* o_selectedNodDagpath )
{
    MStatus status;
    
    MObject sweeneyNodeObj = MObject::kNullObj;

    findASelectedNodeByTypeName( WmSweeneyNode::typeName, &sweeneyNodeObj );
    
    if ( sweeneyNodeObj != MObject::kNullObj )
    {        
        MFnDependencyNode sweeneyNodeDepFn( sweeneyNodeObj, &status );
        CHECK_MSTATUS( status );
        
        o_wmSweeneyNode = dynamic_cast< WmSweeneyNode* >( sweeneyNodeDepFn.userNode() );    
        o_wmSweenyRodManager = o_wmSweeneyNode->rodManager();
    
        return MStatus::kSuccess;
    }
    
    MGlobal::displayError( "Please select a wmSweeneyNode." );
    
    return MStatus::kFailure;
}

MStatus findPointOnMeshFrom2dScreenCoords( MFnMesh& i_meshFn, short i_x, short i_y, 
    MFloatPoint& o_position, MFloatVector& o_normal )
{
    MStatus status; 
    
    M3dView view = M3dView::active3dView();

    MPoint startPoint;
    MVector rayDirection;
    status = view.viewToWorld( i_x, i_y, startPoint, rayDirection );
    CHECK_MSTATUS( status );
    if ( status.error() )
    {
        return status;
    }
    
    MFloatPoint startFloatPoint( startPoint.x, startPoint.y, startPoint.z );
    MFloatPoint rayFloatDirection( rayDirection.x, rayDirection.y, rayDirection.z );
        
    status = findClosestRayMeshIntersection( i_meshFn, startFloatPoint, rayFloatDirection, 
                                             o_position, o_normal );
        
    return status;
}

MStatus findClosestRayMeshIntersection( MFnMesh& i_meshFn, const MFloatPoint& i_rayStart, const MFloatVector& i_rayVec,
    MFloatPoint& o_hit, MFloatVector& o_normal ) 
{
    MS status = MStatus::kSuccess;

    MMeshIsectAccelParams mmAccelParams = i_meshFn.autoUniformGridParams();

    float hitRayParam, hitBary1, hitBary2;
    MFloatPoint hitPoint;
    int hitFace, hitTriangle;

    bool intersected = i_meshFn.closestIntersection( i_rayStart, i_rayVec, NULL, NULL, false,
        MSpace::kWorld, 9999.9f, false, &mmAccelParams, hitPoint, &hitRayParam, &hitFace, 
        &hitTriangle, &hitBary1, &hitBary2, (float)1e-6,  & status );
    
    CHECK_MSTATUS( status );

    if ( intersected ) 
    {
        int vertList[ 3 ];
        status = i_meshFn.getPolygonTriangleVertices( hitFace, hitTriangle, vertList );
        CHECK_MSTATUS( status );

		MVector	normals[3];
		status = i_meshFn.getVertexNormal( vertList[ 0 ], normals[ 0 ], MSpace::kWorld );
		status = i_meshFn.getVertexNormal( vertList[ 1 ], normals[ 1 ], MSpace::kWorld );
		status = i_meshFn.getVertexNormal( vertList[ 2 ], normals[ 2 ], MSpace::kWorld );

		MVector n(((hitBary1)*normals[ 0 ]) + ((hitBary2)*normals[ 1 ]) + ((1-hitBary1-hitBary2)*normals[ 2 ]));		

        n.normalize();

        o_hit = hitPoint;
        o_normal =  n;

        return MStatus::kSuccess;
    }
    else 
    {
        o_hit = MFloatPoint( 999999.0f, 999999.0f, 999999.0f );
        o_normal = MVector( 0.0f, 0.0f, 0.0f );

        return MStatus::kFailure;
    }
}




}
}
