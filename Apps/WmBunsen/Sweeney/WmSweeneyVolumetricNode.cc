#include "WmSweeneyVolumetricNode.hh"

#include <sstream>

#include <maya/MFnMatrixAttribute.h>

using namespace std;

MTypeId WmSweeneyVolumetricNode::typeID( 0x001156, 0xAB );
MString WmSweeneyVolumetricNode::typeName( "wmSweeneyVolumetricNode" );
MObject WmSweeneyVolumetricNode::ia_time;
MObject WmSweeneyVolumetricNode::ia_startTime;
MObject WmSweeneyVolumetricNode::ia_inMesh;
MObject WmSweeneyVolumetricNode::ia_meshTransform;
MObject WmSweeneyVolumetricNode::oa_meshData;
    
WmSweeneyVolumetricNode::WmSweeneyVolumetricNode()
    : m_triangleMeshRenderer( NULL ), m_currentMesh( NULL ),
      m_previousMesh( NULL ), m_nextMesh( NULL )
{
}

WmSweeneyVolumetricNode::~WmSweeneyVolumetricNode()
{
    delete m_currentMesh;
    delete m_previousMesh;
    delete m_nextMesh;
    delete m_triangleMeshRenderer;
}

void WmSweeneyVolumetricNode::initialise( const unsigned int i_volumetricMeshIndex,
                                            TriangleMesh** o_currentMesh )
{
    MStatus status = MStatus::kSuccess;

    m_volumetricMeshIndex = i_volumetricMeshIndex;
        
    if ( m_currentMesh != NULL )
    {
        delete m_currentMesh;
    }
    if ( m_previousMesh != NULL )
    {
        delete m_previousMesh;
    }
    if ( m_nextMesh != NULL )
    {
        delete m_nextMesh;
    }
    if ( m_triangleMeshRenderer != NULL )
    {
        delete m_triangleMeshRenderer;
    }

    m_currentMesh = new TriangleMesh();
    m_previousMesh = new TriangleMesh();
    m_nextMesh = new TriangleMesh();

    // Now we need to ensure that the mesh variables have the mesh data from the input mesh.
    // ::compute() has already been called and we could initialise the meshes there but then we
    // have the problem of not knowing when the simulation has been reset and we need to delete
    // and rebuild everything. The only time we know that is when this function is called. 
    // So the safest thing to do is delete and rebuild everything here then call 
    // updateVolumetricMeshFromMayaMesh() to fill in the data *before* we pass it to Beaker.
    // The reason is that after this function is called, BunsenNode will call Beaker to add the
    // rods to the world at which point the BAGroomingStepper will keep track of the ndof of the
    // meshes. If we don't fill the mesh with the correct number of vertices before that happens
    // then it will screw up everything when it tries to do collision detection later.
    MDataBlock dataBlock = forceCache();
    MDataHandle meshH = dataBlock.inputValue( ia_inMesh, &status );
    if ( status.error() || meshH.type() != MFnData::kMesh )
    {
        cerr << "Problem trying to get mesh from input\n";    
        return;
    }

    MObject inMeshObj = meshH.asMeshTransformed();
    MFnMesh meshFn( inMeshObj );
    
    updateVolumetricMeshFromMayaMesh( meshFn );

    *o_currentMesh = m_currentMesh;

    m_triangleMeshRenderer = new TriangleMeshRenderer( *m_currentMesh );
}

MStatus WmSweeneyVolumetricNode::compute( const MPlug& i_plug, MDataBlock& i_data )
{
    MStatus stat;

   // cerr << "WmSweeneyVolumetricNode::compute() called with plug = " << i_plug.name() << endl;

    if ( i_plug == oa_meshData ) 
    {

        // Get transformation data
        MMatrix transMatrix = i_data.inputValue( ia_meshTransform, &stat ).asMatrix();
        CHECK_MSTATUS( stat );

        MTransformationMatrix transform = MTransformationMatrix( transMatrix );
        CHECK_MSTATUS( stat );

        stat = transform.getRotationQuaternion( m_quaternion.x(), m_quaternion.y(), m_quaternion.z(), m_quaternion.w() );

        stat = transform.getScale( &m_scale[0], MSpace::kTransform );
        CHECK_MSTATUS( stat );

        MVector Mcenter = transform.getTranslation( MSpace::kTransform, &stat );
        CHECK_MSTATUS( stat );
        m_center = Vec3d( Mcenter[0], Mcenter[1], Mcenter[2] );

        // Set up time
        m_previousTime = m_currentTime;
        m_currentTime = i_data.inputValue( ia_time, &stat).asTime().value();
            
        m_startTime = i_data.inputValue( ia_startTime, &stat ).asDouble();
        CHECK_MSTATUS( stat );

        if ( m_currentTime != m_previousTime || m_currentTime == m_startTime ) 
        {
            // We need to force a mesh update even if time hasn't changed but
            // currentTime == startTime. This is because if the scene loaded on start time then
            // the EZBake node won't have been initialised and the user is likely to click 
            // <- meaning go to start frame. This will initialise EZBake and the rods
            // but not the collision data unless we have this exception. 
            // This shouldn't cause extra work because this compute will only get
            // called if one of our inputs has changed and we need to do work anyway.
            // FIXME: Could we just get rid of the entire 'if' then?
            
            MDataHandle meshH = i_data.inputValue( ia_inMesh, &stat );
            if ( !stat.error() && meshH.type() == MFnData::kMesh )
            {
                MObject inMeshObj = meshH.asMeshTransformed();
                MFnMesh meshFn( inMeshObj );
   
                updateVolumetricMeshFromMayaMesh( meshFn );
            }
        }
        
        CHECK_MSTATUS( stat );
        
        MDataHandle outputData = i_data.outputValue ( oa_meshData, &stat );
        if (!stat) 
        {
            stat.perror( "WmSweeneyVolumetricNode::compute get oa_meshData" );
            return stat;
        }
        
        outputData.set( true );
        
        stat = i_data.setClean( i_plug );
        if (!stat) 
        {
            stat.perror( "WmSweeneyVolumetricNode::compute setClean" );
            return stat;
        }
    } 
    else 
    {
        return MS::kUnknownParameter;
    }

    return MS::kSuccess;
}

void WmSweeneyVolumetricNode::draw( M3dView & view, const MDagPath & path,
        M3dView::DisplayStyle style, M3dView::DisplayStatus status )
{
    //
    // Nothing to draw
    //
}

MStatus WmSweeneyVolumetricNode::updateVolumetricMeshFromMayaMesh( MFnMesh &i_meshFn,
    bool forceReset, std::string i_filename )
{    
    MStatus stat = MStatus::kSuccess;
    cerr << "Updating volumetric mesh from maya mesh\n";
    
    MPointArray points;
    i_meshFn.getPoints( points, MSpace::kWorld );
    
    if ( points.length() == 0 )
    {
        cerr << "WmSweeneyVolumetricNode::initialise ERROR, no points in mesh\n";
        return MStatus::kFailure;
    }
    
    if ( m_previousMesh == NULL )
    {
        cerr << "WmSweeneyVolumetricNode::updateCollisionMeshFromMayaMesh() - No mesh to update yet\n";
        return MStatus::kFailure;
    }

    if ( m_previousMesh->nv() == 0 )
    {
        cerr << "Mesh not initialised, initialising\n";
        
        MIntArray triangleCounts;
        MIntArray triangleVertexIndices;
    
        i_meshFn.getTriangles( triangleCounts, triangleVertexIndices );
    
        // TODO: Check if we can preallocate all the vertices/faces as this must be really slow
        // adding things one at a time.  
        for ( unsigned int v = 0; v < points.length(); ++v )
        {
            TriangleMesh::vertex_handle vertexHandle1 = m_currentMesh->addVertex();
            m_currentMesh->getVertex( vertexHandle1 ) = Vec3d( points[ v ].x, points[ v ].y, points[ v ].z );

            TriangleMesh::vertex_handle vertexHandle2 = m_previousMesh->addVertex();
            m_previousMesh->getVertex( vertexHandle2 ) = Vec3d( points[ v ].x, points[ v ].y, points[ v ].z );

            TriangleMesh::vertex_handle vertexHandle3 = m_nextMesh->addVertex();
            m_nextMesh->getVertex( vertexHandle3 ) = Vec3d( points[ v ].x, points[ v ].y, points[ v ].z );
        }
    
        for ( int i = 0; i < triangleVertexIndices.length(); i += 3 )
        {
            int index1 = triangleVertexIndices[ i ];
            int index2 = triangleVertexIndices[ i + 1 ];
            int index3 = triangleVertexIndices[ i + 2 ];
    
            TriangleMesh::vertex_handle vhx( index1 );
            TriangleMesh::vertex_handle vhy( index2 );
            TriangleMesh::vertex_handle vhz( index3 );
        
            // Generate three edges for the face
            m_currentMesh->addEdge( vhx, vhy ); m_currentMesh->addEdge( vhy, vhz ); m_currentMesh->addEdge( vhz, vhx );
            m_previousMesh->addEdge( vhx, vhy ); m_previousMesh->addEdge( vhy, vhz ); m_previousMesh->addEdge( vhz, vhx );
            m_nextMesh->addEdge( vhx, vhy ); m_nextMesh->addEdge( vhy, vhz ); m_nextMesh->addEdge( vhz, vhx );
        
            // Generate a triangular face
            m_currentMesh->addFace( vhx, vhy, vhz );
            m_previousMesh->addFace( vhx, vhy, vhz );
            m_nextMesh->addFace( vhx, vhy, vhz );
        }
        
    }
    else
    {
        cerr << "Updating pre-initialised mesh\n";    

        MPointArray points;
        i_meshFn.getPoints( points, MSpace::kWorld );
        vector<BASim::Vec3d> vPoints;
                
        if ( points.length() != m_previousMesh->nv() )
        {
            cerr << "WmSweeneyVolumetricNode::updateVolumetricMeshFromMayaMesh ERROR - number of points in maya mesh != points in BASim mesh!\n";
            return MStatus::kFailure;
        }

        // First make the previous mesh equal the next mesh as we're about to update next mesh
        for( TriangleMesh::vertex_iter previousMeshItr = m_previousMesh->vertices_begin(),
             nextMeshItr = m_nextMesh->vertices_begin();
             previousMeshItr != m_previousMesh->vertices_end(); ++previousMeshItr, ++nextMeshItr )
        {
            m_previousMesh->getVertex( *previousMeshItr ) = m_nextMesh->getVertex( *nextMeshItr );
        }
        
        // Now update the next mesh with the points of the Maya mesh
        unsigned int v = 0;
        for( TriangleMesh::vertex_iter itr = m_nextMesh->vertices_begin(); 
            itr != m_nextMesh->vertices_end(); ++itr )
        {
            m_nextMesh->getVertex( *itr ) = Vec3d( points[ v ].x, points[ v ].y, points[ v ].z );
            ++v;
        }
    }

    return stat;
}

MStatus WmSweeneyVolumetricNode::connectionMade( const  MPlug & i_plug, const  MPlug& i_otherPlug,
                                                   bool i_asSrc ) 
{    
    MStatus stat;

    return MS::kUnknownParameter;
}

MStatus WmSweeneyVolumetricNode::connectionBroken( const  MPlug& i_plug,
            const  MPlug& i_otherPlug, bool i_asSrc) 
{

    if ( i_plug == ia_inMesh)
    {
        // Get rid of all the mesh data so that Beaker knows we have nothing to provide.
        // Otherwise will get ghost collisions happening with whatever mesh positions were left
        // in the data.
        //m_collisionMeshData->clearAll();
        
        m_meshConnected = true;
    }    
                
    return MStatus::kUnknownParameter;
}

bool WmSweeneyVolumetricNode::isBounded() const
{ 
    return false;
}

void* WmSweeneyVolumetricNode::creator()
{
    return new WmSweeneyVolumetricNode();
}

void WmSweeneyVolumetricNode::postConstructor()
{
    setExistWithoutInConnections( false );
}

MStatus WmSweeneyVolumetricNode::initialize()
{ 
    MStatus    stat;
    
    {    
        MFnUnitAttribute    uAttr;
        ia_time = uAttr.create( "time", "t", MTime( 0.0 ), &stat );
        if (!stat) {
            stat.perror( "create ia_time attribute" );
            return stat;
        }
        uAttr.setWritable( true );
        uAttr.setConnectable( true );
        uAttr.setStorable( false );
        stat = addAttribute( ia_time );
        if (!stat) { stat.perror( "addAttribute time" ); return stat;}
    }
    
    {
        MFnNumericAttribute nAttr;
        ia_startTime = nAttr.create( "startTime", "stt", MFnNumericData::kDouble, 1.0, &stat );
        if (!stat) {
            stat.perror( "create ia_startTime attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );
        stat = addAttribute(ia_startTime);
        if (!stat) { stat.perror("addAttribute startTime"); return stat;}
    }

    {
      MStatus stat = MStatus::kSuccess;
      MFnTypedAttribute inMeshFn;
      ia_inMesh = inMeshFn.create( "inMesh", "in", MFnData::kMesh, MObject::kNullObj, &stat );
      if (!stat)
      {
         stat.perror( "create ia_inMesh attribute") ;
         return stat;
      }
      inMeshFn.setWritable( true );
      inMeshFn.setReadable( true );
      inMeshFn.setConnectable( true );
      inMeshFn.setDisconnectBehavior( MFnAttribute::kReset );
      inMeshFn.setStorable( false );
      inMeshFn.setArray( false );
      addAttribute( ia_inMesh );
      if (!stat)
      {
         stat.perror("addAttribute ia_inMesh");
         return stat;
      }
    }

    {
        MFnMatrixAttribute nAttr;
        ia_meshTransform = nAttr.create( "meshTransform", "mtr", MFnMatrixAttribute::kDouble, &stat );
        if (!stat) 
        {
            stat.perror( "create ia_meshTransform attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );
        stat = addAttribute( ia_meshTransform );
        if(!stat){ stat.perror("addAttribute ia_meshTransform"); return stat;}
    }

    {
        MFnNumericAttribute nAttr;
        oa_meshData = nAttr.create("meshData", "md", MFnNumericData::kBoolean, false, &stat);
        if (!stat) {
            stat.perror("create oa_meshData attribute");
            return stat;
        }
        nAttr.setWritable(false);
        nAttr.setReadable(true);
        nAttr.setConnectable(true);
        // Tell Maya that this node is dependant on world space transformations of the node
        nAttr.setWorldSpace( true );
        stat = addAttribute( oa_meshData );
        if (!stat) { stat.perror("addAttribute oa_meshData"); return stat;}
    }
    
    stat = attributeAffects( ia_time, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_time->oa_meshData" ); return stat;}
    stat = attributeAffects(ia_startTime, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_startTime->oa_meshData" ); return stat;}
    stat = attributeAffects( ia_inMesh, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_time->oa_meshData" ); return stat;}
    stat = attributeAffects( ia_meshTransform, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_meshTransform->oa_meshData" ); return stat;}
    
    return MS::kSuccess;
}


