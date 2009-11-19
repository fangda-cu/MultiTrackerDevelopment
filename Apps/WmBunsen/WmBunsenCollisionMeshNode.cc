#include "WmBunsenCollisionMeshNode.hh"

#include <sstream>

MTypeId WmBunsenCollisionMeshNode::typeId( 0x001135, 0x1C );
MString WmBunsenCollisionMeshNode::typeName( "wmBunsenCollisionMeshNode" );
MObject WmBunsenCollisionMeshNode::ia_time;
MObject WmBunsenCollisionMeshNode::ia_startTime;
MObject WmBunsenCollisionMeshNode::ia_inMesh;
MObject WmBunsenCollisionMeshNode::oa_meshData;
MObject WmBunsenCollisionMeshNode::ia_levelsetDx;
MObject WmBunsenCollisionMeshNode::ia_friction;
MObject WmBunsenCollisionMeshNode::ia_thickness;
MObject WmBunsenCollisionMeshNode::ia_separationStrength;
MObject WmBunsenCollisionMeshNode::ia_coefficientOfRestitution;
MObject WmBunsenCollisionMeshNode::ia_fullCollisions;
MObject WmBunsenCollisionMeshNode::ia_drawCollisionData;

WmBunsenCollisionMeshNode::WmBunsenCollisionMeshNode()
    : m_levelsetDx(0.0), m_friction(0.0), m_thickness(1.0), m_fullCollisions(false), 
    m_drawCollisionData(false), m_collisionMeshData( NULL )
{
    m_collisionMeshData = new CollisionMeshData;
}

WmBunsenCollisionMeshNode::~WmBunsenCollisionMeshNode() 
{
    if ( m_collisionMeshData != NULL )
        delete m_collisionMeshData;
}

MStatus WmBunsenCollisionMeshNode::compute( const MPlug& i_plug, MDataBlock& i_data ) 
{
	MStatus stat;

    //cerr << "WmBunsenCollisionMeshNode::compute() called with plug = " << plug.name() << endl;

	if ( i_plug == oa_meshData ) 
    {
        m_previousTime = m_currentTime;
	    m_currentTime = i_data.inputValue( ia_time, &stat).asTime().value();
			
        m_startTime = i_data.inputValue( ia_startTime, &stat ).asDouble();
        CHECK_MSTATUS( stat );

        if ( m_currentTime != m_previousTime ) 
        {
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
//        m_collisionMeshData->setLevelsetDx( m_levelsetDx );

        m_friction = i_data.inputValue( ia_friction, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_collisionMeshData->setFriction( m_friction );

        m_thickness = i_data.inputValue( ia_thickness, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_collisionMeshData->setThickness( m_thickness );
        
        m_fullCollisions = i_data.inputValue( ia_fullCollisions, &stat ).asBool();
        CHECK_MSTATUS( stat );
        m_collisionMeshData->setFullCollisions( m_fullCollisions );
        
        m_drawCollisionData = i_data.inputValue( ia_drawCollisionData, &stat ).asBool();
        CHECK_MSTATUS( stat );
        
        Real separationStrength = i_data.inputValue( ia_separationStrength, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_collisionMeshData->setSeparationStrength( separationStrength );
        
        Real coefficientOfRestitution = i_data.inputValue( ia_coefficientOfRestitution, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_collisionMeshData->setCoefficientOfRestitution( coefficientOfRestitution );
		
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
	MStatus stat;
	MObject thisNode = thisMObject();

	view.beginGL(); 
	glPushAttrib(GL_CURRENT_BIT | GL_POINT_BIT | GL_LINE_BIT);

    if ( m_drawCollisionData && m_collisionMeshData )
        m_collisionMeshData->draw();

	glPopAttrib();
	view.endGL();
}

MStatus WmBunsenCollisionMeshNode::updateCollisionMeshFromMayaMesh( MFnMesh &i_meshFn, 
    bool forceReset, std::string i_filename )
{    
    MStatus stat;

    MPointArray points;
    i_meshFn.getPoints( points, MSpace::kWorld );
    
    vector<BASim::Vec3d> vPoints;
    
    // We don't want anything in BASim knowing about Maya so convert all the points
    // into a BASim format before passing them in.
    unsigned int numPoints = points.length();
    vPoints.resize( numPoints );
    for ( unsigned int p=0; p<numPoints; p++ )
    {
        vPoints[p][0] = points[p][0];
        vPoints[p][1] = points[p][1];
        vPoints[p][2] = points[p][2];
    }
        
    if ( m_currentTime == m_startTime || forceReset )
        m_collisionMeshData->reset( vPoints );
    else
        m_collisionMeshData->update( vPoints, "", (int)m_currentTime );

    return stat;
}

MStatus WmBunsenCollisionMeshNode::connectionMade( const  MPlug & i_plug, const  MPlug& i_otherPlug, 
                                                   bool i_asSrc ) 
{    
    MStatus stat;

    if( i_plug == ia_inMesh )
    {
        MObject meshObj;
        stat = i_plug.getValue( meshObj );
        CHECK_MSTATUS( stat );
        MFnMesh meshFn( meshObj, &stat );
        CHECK_MSTATUS( stat );
       
        m_collisionMeshData->currPositions.resize( meshFn.numVertices() );
        m_collisionMeshData->prevPositions.resize( meshFn.numVertices() );
        m_collisionMeshData->oldPositions.resize( meshFn.numVertices() );
        m_collisionMeshData->newPositions.resize( meshFn.numVertices() );
        m_collisionMeshData->velocities.resize( meshFn.numVertices() );

        MIntArray triangleCounts;
        MIntArray triangleInds;
        meshFn.getTriangles( triangleCounts, triangleInds );
    
        m_collisionMeshData->triangleIndices.resize( triangleInds.length() );
        
        for ( size_t i=0; i<triangleInds.length(); ++i )
            m_collisionMeshData->triangleIndices[i] = triangleInds[i];

        updateCollisionMeshFromMayaMesh( meshFn, true );
    }

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
        m_collisionMeshData->clearAll();
        
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
	MStatus	stat;
	
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
        ia_friction = nAttr.create("friction", "fri", MFnNumericData::kDouble, 0.0, &stat);
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
        ia_separationStrength = nAttr.create("separationStrength", "ss", MFnNumericData::kDouble, 10000.0, &stat);
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
        ia_thickness = nAttr.create("thickness", "thk", MFnNumericData::kDouble, 1.0, &stat);
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
        ia_fullCollisions = nAttr.create("fullCollisions", "fc", MFnNumericData::kBoolean, false, &stat);
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
    stat = attributeAffects( ia_coefficientOfRestitution, oa_meshData );
    if (!stat) { stat.perror("attributeAffects ia_coefficientOfRestitution->oa_meshData"); return stat;}
    stat = attributeAffects( ia_fullCollisions, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_fullCollisions->oa_meshData" ); return stat;}
    stat = attributeAffects( ia_drawCollisionData, oa_meshData );
    if (!stat) { stat.perror( "attributeAffects ia_drawCollisionData->oa_meshData" ); return stat;}
    
	return MS::kSuccess;
}


