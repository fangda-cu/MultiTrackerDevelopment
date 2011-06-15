#include "WmBunsenCollisionMeshNode.hh"

#include <sstream>

using namespace std;

MTypeId WmBunsenCollisionMeshNode::typeId( 0x001135, 0x1C );
MString WmBunsenCollisionMeshNode::typeName( "wmFigCollisionNode" );
MObject WmBunsenCollisionMeshNode::ia_time;
MObject WmBunsenCollisionMeshNode::ia_startTime;
MObject WmBunsenCollisionMeshNode::ia_inMesh;
MObject WmBunsenCollisionMeshNode::oa_meshData;
MObject WmBunsenCollisionMeshNode::ia_friction;
MObject WmBunsenCollisionMeshNode::ia_thickness;
MObject WmBunsenCollisionMeshNode::ia_separationStrength;
MObject WmBunsenCollisionMeshNode::ia_damping;
MObject WmBunsenCollisionMeshNode::ia_coefficientOfRestitution;
MObject WmBunsenCollisionMeshNode::ia_fullCollisions;
MObject WmBunsenCollisionMeshNode::ia_drawCollisionData;
MObject WmBunsenCollisionMeshNode::ia_createLevelSet;
MObject WmBunsenCollisionMeshNode::ia_levelSetCellSize;
MObject WmBunsenCollisionMeshNode::ia_drawLevelSet;
    
WmBunsenCollisionMeshNode::WmBunsenCollisionMeshNode()
    : m_friction( 0.0 ), m_thickness( 1.0 ), m_fullCollisions( false ),
    m_beaker( NULL ), m_meshController( NULL ), m_triangleMeshRenderer( NULL ),
    m_currentMesh( NULL ), m_previousMesh( NULL ), m_nextMesh( NULL )
{
}

WmBunsenCollisionMeshNode::~WmBunsenCollisionMeshNode() 
{
    delete m_currentMesh;
    delete m_previousMesh;
    delete m_nextMesh;
    delete m_meshController;
    delete m_triangleMeshRenderer;
}

void WmBunsenCollisionMeshNode::initialise( Beaker* i_beaker, const unsigned int i_collisionMeshIndex,
                                            TriangleMesh** o_currentMesh,  
                                            WmFigMeshController** o_meshController)
{
    MStatus status = MStatus::kSuccess;

    m_beaker = i_beaker;
    m_collisionMeshIndex = i_collisionMeshIndex;
        
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
    if ( m_meshController != NULL )
    {
        delete m_meshController;
    }

    m_currentMesh = new TriangleMesh();
    m_previousMesh = new TriangleMesh();
    m_nextMesh = new TriangleMesh();

    // NOTE: We're passing in bogus dt values but they are never used, the BARodStepper overrides
    // them.... Yes it's not clear but I'm working the Maya code into Breannan's BASimulator structure.
    m_meshController = new WmFigMeshController( m_currentMesh, m_previousMesh, m_nextMesh, 
                                                m_startTime, 1.0/24.0);

    // Now we need to ensure that the mesh variables have the mesh data from the input mesh.
    // ::compute() has already been called and we could initialise the meshes there but then we
    // have the problem of not knowing when the simulation has been reset and we need to delete
    // and rebuild everything. The only time we know that is when this function is called. 
    // So the safest thing to do is delete and rebuild everything here then call 
    // updateCollisionMeshFromMayaMesh() to fill in the data *before* we pass it to Beaker. 
    // The reason is that after this function is called, BunsenNode will call Beaker to add the
    // rods to the world at which point the BARodStepper will keep track of the ndof of the 
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

    // Grab the level set details as we need them for creating the level set as this is
    // called by the SweeneyNode and we're circumventing the compute() evaluation to speed
    // things up and to ensure Sweeney can set up the sim in the right order
    // TODO: Re-evaluate if we really need to call this function from Sweeney or can just
    // cleverly pull on attributes in the right order so that it all happens in compute()
    // thus making things a lot simpler.
    float levelsetCellSize = dataBlock.inputValue( ia_levelSetCellSize, &status ).asDouble();
    CHECK_MSTATUS( status );
    bool drawLevelSet = dataBlock.inputValue( ia_drawLevelSet, &status ).asBool();
    CHECK_MSTATUS( status );
   
    m_meshController->setLevelSetCellSize( levelsetCellSize );
    m_meshController->drawLevelSet( drawLevelSet );
    
    updateCollisionMeshFromMayaMesh( meshFn );

    if ( m_beaker != NULL )
    {
        // If we are being controlled by beaker then initialise it
        m_beaker->initialiseCollisionMesh( m_currentMesh, m_meshController->currentLevelSet(), 
                                           m_meshController, m_collisionMeshIndex );         
    }
    else
    {
        // If this is a Sweeney sim then just pass back the objects it needs to build the
        // BARodStepper object.
        *o_currentMesh = m_currentMesh;
        *o_meshController = m_meshController;       
    }

    m_triangleMeshRenderer = new TriangleMeshRenderer( *m_currentMesh );
}


MStatus WmBunsenCollisionMeshNode::compute( const MPlug& i_plug, MDataBlock& i_data ) 
{
    MStatus stat;

   // cerr << "WmBunsenCollisionMeshNode::compute() called with plug = " << i_plug.name() << endl;

    if ( i_plug == oa_meshData ) 
    {
        m_previousTime = m_currentTime;
        m_currentTime = i_data.inputValue( ia_time, &stat).asTime().value();
            
        m_startTime = i_data.inputValue( ia_startTime, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        
        // Level Set details
        float levelsetCellSize = i_data.inputValue( ia_levelSetCellSize, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        bool createLevelSet = i_data.inputValue( ia_createLevelSet, &stat ).asBool();
        CHECK_MSTATUS( stat );
        bool drawLevelSet = i_data.inputValue( ia_drawLevelSet, &stat ).asBool();
        CHECK_MSTATUS( stat );
        
        if ( m_meshController != NULL )
        {
         //   m_meshController->createLevelSet( createLevelSet );
            m_meshController->setLevelSetCellSize( levelsetCellSize );
            m_meshController->drawLevelSet( drawLevelSet );
        }

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
   
                updateCollisionMeshFromMayaMesh( meshFn );
            }
        }

        m_drawCollisionData = i_data.inputValue( ia_drawCollisionData, &stat ).asBool();
        CHECK_MSTATUS( stat );
        
//        m_collisionMeshData->setLevelsetDx( m_levelsetDx );

        /*m_friction = i_data.inputValue( ia_friction, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_collisionMeshData->setFriction( m_friction );

        m_thickness = i_data.inputValue( ia_thickness, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_collisionMeshData->setThickness( m_thickness );
        
        m_fullCollisions = i_data.inputValue( ia_fullCollisions, &stat ).asBool();
        CHECK_MSTATUS( stat );
        m_collisionMeshData->setFullCollisions( m_fullCollisions );
        
        Real separationStrength = i_data.inputValue( ia_separationStrength, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_collisionMeshData->setSeparationStrength( separationStrength );
        
        Real damping = i_data.inputValue( ia_damping, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_collisionMeshData->setDamping( damping );

        Real coefficientOfRestitution = i_data.inputValue( ia_coefficientOfRestitution, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_collisionMeshData->setCoefficientOfRestitution( coefficientOfRestitution );*/
        
        MDataHandle outputData = i_data.outputValue ( oa_meshData, &stat );
        if (!stat) 
        {
            stat.perror( "wmBunsenCollisionMeshNode::compute get oa_meshData" );
            return stat;
        }
        
        outputData.set( true );
        
        stat = i_data.setClean( i_plug );
        if (!stat) 
        {
            stat.perror( "wmBunsenCollisionMeshNode::compute setClean" );
            return stat;
        }
    } 
    else 
    {
        return MS::kUnknownParameter;
    }

    return MS::kSuccess;
}

void WmBunsenCollisionMeshNode::draw( M3dView & view, const MDagPath & path,
        M3dView::DisplayStyle style, M3dView::DisplayStatus status )
{     
    view.beginGL(); 
    
    if ( m_meshController )
    {
        m_meshController->draw();
    }
    
    if ( m_drawCollisionData )
    {        
        // TODO: I think we could use one of the BASim render classes to do this        
        glPointSize( 5.0 );

        glBegin( GL_POINTS );
        for( TriangleMesh::vertex_iter itr = m_currentMesh->vertices_begin(); 
             itr != m_currentMesh->vertices_end(); ++itr )
        {
            Vec3d p = m_currentMesh->getVertex( *itr );
            glVertex3d( p[ 0 ], p[ 1 ], p[ 2 ] );
        }
        glEnd();
    
        glPointSize( 1.0 );
    }
    
    view.endGL();
}

MStatus WmBunsenCollisionMeshNode::updateCollisionMeshFromMayaMesh( MFnMesh &i_meshFn, 
    bool forceReset, std::string i_filename )
{    
    MStatus stat = MStatus::kSuccess;
    
    cerr << "Updating collision mesh from maya mesh\n";
    
    MPointArray points;
    i_meshFn.getPoints( points, MSpace::kWorld );
    
    if ( points.length() == 0 )
    {
        cerr << "WmBunsenCollisionMeshNode::initialise ERROR, no points in mesh\n";
        return MStatus::kFailure;
    }
    
    if ( m_previousMesh == NULL )
    {
        cerr << "WmBunsenCollisionMeshNode::updateCollisionMeshFromMayaMesh() - No mesh to update yet\n";
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
        
        // Send the controller a list of indices for the level set code
        // This is a bit messy, all the level set code needs streamlined 
        // - Alasdair
        std::vector< unsigned int > indices;
        indices.resize( triangleVertexIndices.length() );
        
        for ( size_t t = 0; t < indices.size(); ++ t )
        {
            indices[ t ] = triangleVertexIndices[ (unsigned int)t ];
        }
        
        cerr << "Setup meshes with " << triangleVertexIndices.length() << " indices." << endl;
        m_meshController->setTriangleIndices( indices );
    }
    else
    {
        cerr << "Updating pre-initialised mesh\n";    

        MPointArray points;
        i_meshFn.getPoints( points, MSpace::kWorld );
        vector<BASim::Vec3d> vPoints;
                
        if ( points.length() != m_previousMesh->nv() )
        {
            cerr << "WmBunsenCollisionMeshNode::updateCollisionMeshFromMayaMesh ERROR - number of points in maya mesh != points in BASim mesh!\n";
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
    
    if ( m_meshController != NULL )
    {
        double simTime = m_currentTime;
        m_meshController->updateNextMayaTime( simTime );
    }

    return stat;
}

MStatus WmBunsenCollisionMeshNode::connectionMade( const  MPlug & i_plug, const  MPlug& i_otherPlug, 
                                                   bool i_asSrc ) 
{    
    MStatus stat;

/*    if( i_plug == ia_inMesh )
    {
        MObject meshObj;
        stat = i_plug.getValue( meshObj );
        CHECK_MSTATUS( stat );
        MFnMesh meshFn( meshObj, &stat );
        CHECK_MSTATUS( stat );

       // updateCollisionMeshFromMayaMesh( meshFn, true );
    }*/

    return MS::kUnknownParameter;
}

MStatus WmBunsenCollisionMeshNode::connectionBroken( const  MPlug& i_plug, 
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

bool WmBunsenCollisionMeshNode::isBounded() const
{ 
    return false;
}

void* WmBunsenCollisionMeshNode::creator()
{
    return new WmBunsenCollisionMeshNode();
}

void WmBunsenCollisionMeshNode::postConstructor()
{
    setExistWithoutInConnections( false );
}

MStatus WmBunsenCollisionMeshNode::initialize()
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
        MFnNumericAttribute nAttr;
        ia_levelSetCellSize = nAttr.create( "levelSetCellSize", "lsc", MFnNumericData::kDouble, 1.0, &stat );
        if (!stat) 
        {
            stat.perror( "create ia_levelSetCellSize attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );
        stat = addAttribute( ia_levelSetCellSize );
        if(!stat){ stat.perror("addAttribute ia_levelSetCellSize"); return stat;}
    }

    {
        MFnNumericAttribute nAttr;
        ia_createLevelSet = nAttr.create( "createLevelSet", "cls", MFnNumericData::kBoolean, false, &stat);
        if (!stat) {
            stat.perror( "create ia_createLevelSet attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );
        stat = addAttribute( ia_createLevelSet );
        if(!stat){ stat.perror("addAttribute ia_createLevelSet"); return stat;}
    }

    {
        MFnNumericAttribute nAttr;
        ia_drawLevelSet = nAttr.create( "drawLevelSet", "dls", MFnNumericData::kBoolean, false, &stat );
        if (!stat) {
            stat.perror( "create ia_drawLevelSet attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );
        stat = addAttribute( ia_drawLevelSet );
        if(!stat){ stat.perror("addAttribute ia_drawLevelSet"); return stat;}
    }

    {
        MFnNumericAttribute nAttr;
        ia_friction = nAttr.create("friction", "fri", MFnNumericData::kDouble, 0.4, &stat);
        if (!stat) {
            stat.perror( "create ia_friction attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );
        stat = addAttribute(ia_friction);
        if (!stat) { stat.perror("addAttribute friction"); return stat;}
    }

    {
        MFnNumericAttribute nAttr;
        ia_separationStrength = nAttr.create("separationStrength", "ss", MFnNumericData::kDouble, 40.0, &stat);
        if (!stat) {
            stat.perror( "create ia_separationStrength attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );
        stat = addAttribute(ia_separationStrength);
        if (!stat) { stat.perror("addAttribute separationStrength"); return stat;}
    }

    {
        MFnNumericAttribute nAttr;
        ia_damping = nAttr.create("damping", "dmp", MFnNumericData::kDouble, 0.7, &stat);
        if (!stat) {
            stat.perror( "create ia_damping attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );
        stat = addAttribute(ia_damping);
        if (!stat) { stat.perror("addAttribute damping"); return stat;}
    }

    {
        MFnNumericAttribute nAttr;
        ia_coefficientOfRestitution = nAttr.create("coefficientOfRestitution", "cor", MFnNumericData::kDouble, 0.1, &stat);
        if (!stat) {
            stat.perror( "create ia_coefficientOfRestitution attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );
        stat = addAttribute(ia_coefficientOfRestitution);
        if (!stat) { stat.perror("addAttribute coefficientOfRestitution"); return stat;}
    }

    {
        MFnNumericAttribute nAttr;
        ia_thickness = nAttr.create("thickness", "thk", MFnNumericData::kDouble, 0.1, &stat);
        if (!stat) {
            stat.perror( "create ia_thickness attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );
        stat = addAttribute(ia_thickness);
        if (!stat) { stat.perror("addAttribute thickness"); return stat;}
    }

    {
        MFnNumericAttribute nAttr;
        ia_fullCollisions = nAttr.create("edgeCollisions", "ec", MFnNumericData::kBoolean, false, &stat);
        if (!stat) {
            stat.perror( "create ia_fullCollisions attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );
        stat = addAttribute(ia_fullCollisions);
        if (!stat) { stat.perror("addAttribute fullCollisions"); return stat;}
    }

    {
        MFnNumericAttribute nAttr;
        ia_drawCollisionData = nAttr.create("drawCollisionData", "dcd", MFnNumericData::kBoolean, false, &stat);
        if (!stat) {
            stat.perror( "create ia_drawCollisionData attribute" );
            return stat;
        }
        nAttr.setWritable( true );
        nAttr.setReadable( false );
        nAttr.setKeyable( true );
        stat = addAttribute(ia_drawCollisionData);
        if (!stat) { stat.perror("addAttribute drawCollisionData"); return stat;}
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
        stat = addAttribute( oa_meshData );
        if (!stat) { stat.perror("addAttribute oa_meshData"); return stat;}
    }

    stat = attributeAffects( ia_time, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_time->oa_meshData" ); return stat;}
    stat = attributeAffects(ia_startTime, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_startTime->oa_meshData" ); return stat;}
    stat = attributeAffects( ia_inMesh, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_time->oa_meshData" ); return stat;}
    stat = attributeAffects( ia_friction, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_friction->oa_meshData" ); return stat;}
    stat = attributeAffects( ia_separationStrength, oa_meshData );
    if (!stat) { stat.perror("attributeAffects ia_separationStrength->oa_meshData"); return stat;}
    stat = attributeAffects( ia_damping, oa_meshData );
    if (!stat) { stat.perror("attributeAffects ia_damping->oa_meshData"); return stat;}
    stat = attributeAffects( ia_coefficientOfRestitution, oa_meshData );
    if (!stat) { stat.perror("attributeAffects ia_coefficientOfRestitution->oa_meshData"); return stat;}
    stat = attributeAffects( ia_fullCollisions, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_fullCollisions->oa_meshData" ); return stat;}
    stat = attributeAffects( ia_drawCollisionData, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_drawCollisionData->oa_meshData" ); return stat;}
    stat = attributeAffects( ia_levelSetCellSize, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_levelSetCellSize->oa_meshData" ); return stat;}
    stat = attributeAffects( ia_createLevelSet, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_createLevelSet->oa_meshData" ); return stat;}
    stat = attributeAffects( ia_drawLevelSet, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_drawLevelSet->oa_meshData" ); return stat;}
    
    return MS::kSuccess;
}


