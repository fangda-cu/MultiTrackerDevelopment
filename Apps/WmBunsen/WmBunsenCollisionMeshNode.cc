#include "WmBunsenCollisionMeshNode.hh"

#include <sstream>

using namespace std;

MTypeId WmBunsenCollisionMeshNode::typeId( 0x001135, 0x1C );
MString WmBunsenCollisionMeshNode::typeName( "wmFigCollisionNode" );
MObject WmBunsenCollisionMeshNode::ia_time;
MObject WmBunsenCollisionMeshNode::ia_startTime;
MObject WmBunsenCollisionMeshNode::ia_inMesh;
MObject WmBunsenCollisionMeshNode::oa_meshData;
MObject WmBunsenCollisionMeshNode::ia_levelsetDx;
MObject WmBunsenCollisionMeshNode::ia_friction;
MObject WmBunsenCollisionMeshNode::ia_thickness;
MObject WmBunsenCollisionMeshNode::ia_separationStrength;
MObject WmBunsenCollisionMeshNode::ia_damping;
MObject WmBunsenCollisionMeshNode::ia_coefficientOfRestitution;
MObject WmBunsenCollisionMeshNode::ia_fullCollisions;
MObject WmBunsenCollisionMeshNode::ia_drawCollisionData;

WmBunsenCollisionMeshNode::WmBunsenCollisionMeshNode()
    : m_levelsetDx( 0.0 ), m_friction( 0.0 ), m_thickness( 1.0 ), m_fullCollisions( false ),
    m_beaker( NULL ), m_meshController( NULL ), m_triangleMeshRenderer( NULL )
{
}

WmBunsenCollisionMeshNode::~WmBunsenCollisionMeshNode() 
{
    delete m_meshController;
    delete m_triangleMeshRenderer;
}

void WmBunsenCollisionMeshNode::initialise( Beaker* i_beaker, const unsigned int i_collisionMeshIndex )
{
    MStatus status = MStatus::kSuccess;

    m_beaker = i_beaker;
    m_collisionMeshIndex = i_collisionMeshIndex;
        
    m_meshController = new WmFigMeshController( &m_triangleMesh, m_startTime, ( 1.0 / 24.0 ) / 10.0 );
    m_beaker->initialiseCollisionMesh( &m_triangleMesh, m_meshController, m_collisionMeshIndex );    

    m_triangleMeshRenderer = new TriangleMeshRenderer( m_triangleMesh );
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

        m_levelsetDx = i_data.inputValue( ia_levelsetDx, &stat ).asDouble();
        CHECK_MSTATUS( stat );

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
    if ( m_drawCollisionData )
    {
        view.beginGL(); 
        
        // TODO: I think we could use one of the BASim render classes to do this
        
        glPointSize( 5.0 );

        glBegin( GL_POINTS );
        for( TriangleMesh::vertex_iter itr = m_triangleMesh.vertices_begin(); 
             itr != m_triangleMesh.vertices_end(); ++itr )
        {
            Vec3d p = m_triangleMesh.getVertex(*itr);
            glVertex3d( p[ 0 ], p[ 1 ], p[ 2 ] );
        }
        glEnd();

        view.endGL();

        glPointSize( 1.0 );
    }
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
    
    if ( m_triangleMesh.nv() == 0 )
    {
        cerr << "Mesh not initialised, initialising\n";
        
        MIntArray triangleCounts;
        MIntArray triangleVertexIndices;
    
        i_meshFn.getTriangles( triangleCounts, triangleVertexIndices );
    
        // TODO: Check if we can preallocate all the vertices/faces as this must be really slow
        // adding things one at a time.  
        for ( unsigned int v = 0; v < points.length(); ++v )
        {
            TriangleMesh::vertex_handle vertexHandle = m_triangleMesh.addVertex();
            m_triangleMesh.getVertex( vertexHandle ) = Vec3d( points[ v ].x, points[ v ].y, points[ v ].z );
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
            m_triangleMesh.addEdge( vhx, vhy );
            m_triangleMesh.addEdge( vhy, vhz );
            m_triangleMesh.addEdge( vhz, vhx );
        
            // Generate a triangular face
            m_triangleMesh.addFace( vhx, vhy, vhz );
        }
    }
    else
    {
        cerr << "Updating pre-initialised mesh\n";    

        MPointArray points;
        i_meshFn.getPoints( points, MSpace::kWorld );
        vector<BASim::Vec3d> vPoints;
                
        if ( points.length() != m_triangleMesh.nv() )
        {
            cerr << "WmBunsenCollisionMeshNode::updateCollisionMeshFromMayaMesh ERROR - number of points in maya mesh != points in BASim mesh!\n";
            return MStatus::kFailure;
        }
        
        unsigned int v = 0;
        for( TriangleMesh::vertex_iter itr = m_triangleMesh.vertices_begin(); 
            itr != m_triangleMesh.vertices_end(); ++itr )
        {
            m_triangleMesh.getVertex( *itr ) = Vec3d( points[ v ].x, points[ v ].y, points[ v ].z );
            ++v;
        }
    }
    
    if ( m_meshController != NULL )
    {
        double simTime = m_currentTime - m_startTime;
        m_meshController->updateInterpolatedMeshes( simTime );
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
        ia_startTime = nAttr.create("startTime", "stt", MFnNumericData::kDouble, 1.0, &stat);
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
    ia_levelsetDx = nAttr.create("levelsetDx", "ldx", MFnNumericData::kDouble, 0.0, &stat);
    if (!stat) {
        stat.perror( "create ia_levelsetDx attribute" );
        return stat;
    }
    nAttr.setWritable( true );
    nAttr.setReadable( false );
    nAttr.setKeyable( true );
    stat = addAttribute(ia_levelsetDx);
    if(!stat){ stat.perror("addAttribute levelsetDx"); return stat;}
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
    stat = attributeAffects( ia_levelsetDx, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_levelsetDx->oa_meshData"); return stat;}
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
    
    return MS::kSuccess;
}


