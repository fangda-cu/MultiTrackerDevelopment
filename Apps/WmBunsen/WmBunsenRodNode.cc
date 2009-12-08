#include "WmBunsenRodNode.hh"

using namespace BASim;

// Required by Maya to identify node
/* static */ MTypeId WmBunsenRodNode::typeID ( 0x001135, 0x19 ); 
/* static */ MString WmBunsenRodNode::typeName( "wmBunsenRodNode" );

// 
/* static */ MObject WmBunsenRodNode::ia_time;
/* static */ MObject WmBunsenRodNode::ia_startTime;
/* static */ MObject WmBunsenRodNode::ia_nurbsCurves;
/* static */ MObject WmBunsenRodNode::ia_fozzieVertices;

// User adjustable rod Options
/* static */ MObject WmBunsenRodNode::ia_cvsPerRod;
/* static */ MObject WmBunsenRodNode::ia_youngsModulus;
/* static */ MObject WmBunsenRodNode::ia_shearModulus;
/* static */ MObject WmBunsenRodNode::ia_density;
/* static */ MObject WmBunsenRodNode::ia_minorRadius;
/* static */ MObject WmBunsenRodNode::ia_majorRadius;
/* static */ MObject WmBunsenRodNode::ia_vertexSpacing;

// Disk cacheing
/* static */ MObject WmBunsenRodNode::ia_cachePath;
/* static */ MObject WmBunsenRodNode::ia_cacheFrame;
/* static */ MObject WmBunsenRodNode::ia_readFromCache;

// Output and cached attributes
/* static */ MObject WmBunsenRodNode::ca_syncAttrs;
/* static */ MObject WmBunsenRodNode::oa_rodsChanged;
/* static */ MObject WmBunsenRodNode::ia_simStepTaken;
/* static */ MObject WmBunsenRodNode::ca_simulationSync;

/* static */ MObject WmBunsenRodNode::oa_simulatedVertices;
/* static */ MObject WmBunsenRodNode::oa_nonSimulatedVertices;
/* static */ MObject WmBunsenRodNode::oa_verticesInEachRod;

WmBunsenRodNode::WmBunsenRodNode() : m_initialised( false ), mx_rodData( NULL ), mx_world( NULL ),
                                     m_numberOfInputCurves( 0 ), m_cvsPerRod( -1 )
{
    m_rodOptions.YoungsModulus = 100000.0;
    m_rodOptions.ShearModulus = 375.0;
    m_rodOptions.density = 0.01;
    m_rodOptions.radiusA = 0.1;
    m_rodOptions.radiusB = 0.1;
}

WmBunsenRodNode::~WmBunsenRodNode()
{
}

void WmBunsenRodNode::initialiseRodData( vector<RodData*>* i_rodData )
{
    mx_rodData = i_rodData;
    MStatus stat;
 
    // We may get called from not inside a compute by the WmBunsenNode initialising all its
    // data. So build our own datablock to use.
    MDataBlock dataBlock = MPxNode::forceCache();

    ///////////////////////////////////////////////////////////////////////////////
    // If there is a Fozzie node connected then it takes control over any connected
    // NURBS curves.
    ///////////////////////////////////////////////////////////////////////////////    
    bool fozzieNodeConnected;
    MPlug fozziePlug( thisMObject(), ia_fozzieVertices );
    if ( fozziePlug.isConnected() )
    {
        MDataHandle verticesH = dataBlock.inputValue( ia_fozzieVertices, &stat );
        CHECK_MSTATUS( stat );
        MFnVectorArrayData verticesData( verticesH.data(), &stat );
        CHECK_MSTATUS( stat );
        MVectorArray vertices = verticesData.array( &stat );
    
        unsigned int numStrands = vertices.length() / m_cvsPerRod;
            
        size_t inputVertexIndex = 0;
        for ( unsigned int i = 0; i < numStrands; i++ ) 
        {
            (*mx_rodData)[i]->rodOptions = m_rodOptions;
            
            // Override the number of cvs at the moment to match the input curve
            (*mx_rodData)[i]->rodOptions.numVertices = m_cvsPerRod;
            
            // Make sure we have enough space to store the date for each CV. This should only
            // ever cause a resize when we are called by initialiseRodData().
            (*mx_rodData)[ i ]->undeformedVertexPositions.resize( m_cvsPerRod );
            (*mx_rodData)[ i ]->initialVertexPositions.resize( m_cvsPerRod );
            (*mx_rodData)[ i ]->prevVertexPositions.resize( m_cvsPerRod );
            (*mx_rodData)[ i ]->nextVertexPositions.resize( m_cvsPerRod );
            
            std::string frame = "time";
            if ( frame == "time" )
                (*mx_rodData)[ i ]->rodOptions.refFrame = BASim::ElasticRod::TimeParallel;
            else if (frame == "space")
                (*mx_rodData)[ i ]->rodOptions.refFrame = BASim::ElasticRod::SpaceParallel;

            for ( int c=0; c<m_cvsPerRod; c++ ) 
            {
                MVector cv = vertices[ inputVertexIndex ];
                (*mx_rodData)[ i ]->undeformedVertexPositions[ c ]  = Vec3d( cv.x, cv.y, cv.z );
                inputVertexIndex++;
            }
        }
    }
    else
    {
        MDataHandle inputCurveH;
        MObject inputCurveObj;    
        MArrayDataHandle inArrayH = dataBlock.inputArrayValue( ia_nurbsCurves, &stat );
        CHECK_MSTATUS(stat);
       
        size_t numCurvesConnected = inArrayH.elementCount();
    
        for ( unsigned int i = 0; i < numCurvesConnected; i++ )
        {
            inArrayH.jumpToElement( i );
            inputCurveH = inArrayH.inputValue( &stat );
            CHECK_MSTATUS( stat );
    
            // Use asNurbsCurveTransformed to get the curve data as we
            // want it in world space.
            inputCurveObj = inputCurveH.asNurbsCurveTransformed();
            MFnNurbsCurve inCurveFn( inputCurveObj );
            
            MPoint cv;
            int nCVs = 0;
            //if ( m_vertexSpacing == 0.0 )
                nCVs = inCurveFn.numCVs();
            /*else
            {
                nCVs = inCurveFn.length() / m_vertexSpacing + 1;
            }*/
           
            (*mx_rodData)[i]->rodOptions = m_rodOptions;
            
            // Override the number of cvs at the moment to match the input curve
            (*mx_rodData)[i]->rodOptions.numVertices = nCVs;
            
            // Make sure we have enough space to store the date for each CV. This should only
            // ever cause a resize when we are called by initialiseRodData().
            (*mx_rodData)[ i ]->undeformedVertexPositions.resize( nCVs );
            (*mx_rodData)[ i ]->initialVertexPositions.resize( nCVs );
            (*mx_rodData)[ i ]->prevVertexPositions.resize( nCVs );
            (*mx_rodData)[ i ]->nextVertexPositions.resize( nCVs );
            
            std::string frame = "time";
            if ( frame == "time" )
                (*mx_rodData)[ i ]->rodOptions.refFrame = BASim::ElasticRod::TimeParallel;
            else if (frame == "space")
                (*mx_rodData)[ i ]->rodOptions.refFrame = BASim::ElasticRod::SpaceParallel;
        
            for ( int c = 0; c < (*mx_rodData)[ i ]->rodOptions.numVertices; ++c )
            {
                MPoint cv;
               // stat = inCurveFn.getCV( c,cv,MSpace::kWorld );
                stat = inCurveFn.getCV( c,cv,MSpace::kObject );
                CHECK_MSTATUS( stat );
                
                Vec3d inputCurveVertex( cv.x, cv.y, cv.z );
                (*mx_rodData)[ i ]->undeformedVertexPositions[ c ] = inputCurveVertex; 
            }
        }
    }
    
    // Since we're initialising the rod data that means the rod is about to be created. In which
    // case we need to set the current vertex positions since they will not get set by the
    // simulation until the user moves forward a frame.
    
    size_t numRods = mx_rodData->size();
    for ( size_t r=0; r<numRods; r++ )
    {
        size_t numVertices = (*mx_rodData)[ r ]->undeformedVertexPositions.size();
        (*mx_rodData)[ r ]->initialVertexPositions.resize( numVertices );
        for ( size_t v=0; v<numVertices; v++ )
        {
            (*mx_rodData)[ r ]->initialVertexPositions[ v ] = (*mx_rodData)[ r ]->undeformedVertexPositions[ v ];
            (*mx_rodData)[ r ]->prevVertexPositions[ v ] = (*mx_rodData)[ r ]->undeformedVertexPositions[ v ];
            (*mx_rodData)[ r ]->nextVertexPositions[ v ] = (*mx_rodData)[ r ]->undeformedVertexPositions[ v ];
        }
    }
}

void WmBunsenRodNode::updateRodDataFromInputs()
{    
    if ( mx_rodData == NULL )
    {
        MGlobal::displayError( "Please rewind simulation to initialise\n" );
        return;
    }
    
    MStatus stat;
 
    // We may get called from not inside a compute by the WmBunsenNode initialising all its
    // data. So build our own datablock to use.
    MDataBlock dataBlock = MPxNode::forceCache();
    
    bool fozzieNodeConnected;
    MPlug fozziePlug( thisMObject(), ia_fozzieVertices );
    if ( fozziePlug.isConnected() )
    {
        MDataHandle verticesH = dataBlock.inputValue( ia_fozzieVertices, &stat );
        CHECK_MSTATUS( stat );
        MFnVectorArrayData verticesData( verticesH.data(), &stat );
        CHECK_MSTATUS( stat );
        MVectorArray vertices = verticesData.array( &stat );
    
        size_t numStrands = vertices.length() / m_cvsPerRod;
        
        if ( mx_rodData->size() != numStrands )
        {
            MGlobal::displayError( "Number of rods does not equal number of Fozzie strands, rewind simulation to reset" );
            return;
        }
            
        size_t inputVertexIndex = 0;
        for ( size_t i = 0; i < numStrands; i++ ) 
        {
            BASim::ElasticRod* rod = (*mx_rodData)[ i ]->rod;
            if ( rod != NULL )
            {
                rod->setRadius( m_rodOptions.radiusA, m_rodOptions.radiusB );
                rod->setYoungsModulus( m_rodOptions.YoungsModulus );
                rod->setShearModulus( m_rodOptions.ShearModulus );
                rod->setDensity(m_rodOptions.density);
                
                for ( int c = 0; c < m_cvsPerRod ; ++c ) 
                {
                    MVector cv = vertices[ inputVertexIndex ];
            
                    Vec3d inputCurveVertex( cv.x, cv.y, cv.z );
                    (*mx_rodData)[ i ]->prevVertexPositions[ c ] = (*mx_rodData)[ i ]->nextVertexPositions[ c ];
                    (*mx_rodData)[ i ]->nextVertexPositions[ c ] = inputCurveVertex;
                    
                    inputVertexIndex++;
                }
            }
        }
    }
    else
    {
        if ( mx_rodData->size() != m_numberOfInputCurves )
        {
            MGlobal::displayError( "Number of rods does not equal number of input curves, rewind simulation to reset" );
            return;
        }
        
        MDataHandle inputCurveH;
        MObject inputCurveObj;
        MPoint pt;
        MFnNurbsCurveData dataCreator;
        MMatrix mat;
        MPoint rootP;
        MArrayDataHandle inArrayH = dataBlock.inputArrayValue( ia_nurbsCurves, &stat );
        CHECK_MSTATUS(stat);
       
        size_t numCurvesConnected = inArrayH.elementCount();
    
        for ( unsigned int i = 0; i < numCurvesConnected; i++ )
        {
            inArrayH.jumpToElement( i );
            inputCurveH = inArrayH.inputValue( &stat );
            CHECK_MSTATUS( stat );
    
            // Use asNurbsCurveTransformed to get the curve data as we
            // want it in world space.
            inputCurveObj = inputCurveH.asNurbsCurveTransformed();
            MFnNurbsCurve inCurveFn( inputCurveObj );
            
            MPoint cv;
            int nCVs = inCurveFn.numCVs();
           
            BASim::ElasticRod* rod = (*mx_rodData)[ i ]->rod;
            if ( rod != NULL )
            {
                rod->setRadius( m_rodOptions.radiusA, m_rodOptions.radiusB );
                rod->setYoungsModulus( m_rodOptions.YoungsModulus );
                rod->setShearModulus( m_rodOptions.ShearModulus );
                rod->setDensity(m_rodOptions.density);
                
                int numVertices = rod->nv();
                for ( int c = 0; c < numVertices ; ++c ) 
                {
                    MPoint cv;
                   // stat = inCurveFn.getCV( c,cv,MSpace::kWorld );
                    stat = inCurveFn.getCV( c,cv,MSpace::kObject );
                    CHECK_MSTATUS( stat );
    
                    Vec3d inputCurveVertex( cv.x, cv.y, cv.z );
                    (*mx_rodData)[ i ]->prevVertexPositions[ c ] = (*mx_rodData)[ i ]->nextVertexPositions[ c ];
                    (*mx_rodData)[ i ]->nextVertexPositions[ c ] = inputCurveVertex;
                }
            }
        }
    }
}

MStatus WmBunsenRodNode::compute( const MPlug& i_plug, MDataBlock& i_dataBlock ) 
{
    MStatus stat;
    
    //cerr << "WmBunsenRodNode::compute plug = " << i_plug.name() << endl;
	
    if ( i_plug == oa_rodsChanged )
    {
        m_previousTime = m_currentTime;
        m_currentTime = i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
		CHECK_MSTATUS( stat );
        
	    m_startTime = i_dataBlock.inputValue( ia_startTime, &stat ).asDouble();
		CHECK_MSTATUS( stat );
        
        m_rodOptions.YoungsModulus = i_dataBlock.inputValue( ia_youngsModulus, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_rodOptions.ShearModulus = i_dataBlock.inputValue( ia_shearModulus, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_rodOptions.density = i_dataBlock.inputValue( ia_density, &stat).asDouble();
        CHECK_MSTATUS( stat );
        m_rodOptions.radiusA = i_dataBlock.inputValue( ia_minorRadius, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_rodOptions.radiusB = i_dataBlock.inputValue( ia_majorRadius, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        MString cachePath = i_dataBlock.inputValue( ia_cachePath, &stat ).asString();
        CHECK_MSTATUS( stat );
        //m_vertexSpacing = i_dataBlock.inputValue( ia_vertexSpacing, &stat ).asDouble();
        //CHECK_MSTATUS( stat );
        m_cvsPerRod = i_dataBlock.inputValue( ia_cvsPerRod, &stat ).asInt();
        CHECK_MSTATUS( stat );
		
        bool readFromCache = i_dataBlock.inputValue( ia_readFromCache, &stat ).asDouble();
        if ( readFromCache )
        {
            // Make sure that every rod is disabled before we read from the cache file
            if ( mx_rodData != NULL )
            {
                size_t numRods = mx_rodData->size();
                for ( size_t r=0; r<numRods; r++ )
                {
                    (*mx_rodData)[r]->stepper->setEnabled( false );
                }
                
                readRodDataFromCacheFile( cachePath );
            }
        }
        else if ( mx_rodData != NULL )
        {
            // The Bunsen node is asking us to update the rod data in Beaker for our inputs
            // so do so here.....
            
            // Make sure that every rod is disabled before we read from the cache file
            size_t numRods = mx_rodData->size();
            for ( size_t r=0; r<numRods; r++ )
                (*mx_rodData)[r]->stepper->setEnabled( true );
            
            updateRodDataFromInputs();
        }
        
        stat = i_dataBlock.setClean( i_plug );
		if ( !stat )
        {
			stat.perror("WmBunsenRodNode::compute setClean");
			return stat;
		}
    }
    else if ( i_plug == ca_simulationSync )
    {
        // If we're here then the Bunsen node has moved forward in time and we should pull the
        // simulated data from mx_rodData.
        
        // Later when we add nurbs curve outputs we should do everything we do in here, so maybe
        // make this if i_plug == ca_simulationSync || i_plug == oa_nurbsCurves
        
        i_dataBlock.inputValue( ia_simStepTaken, &stat ).asBool();
        CHECK_MSTATUS( stat );
        
        bool cacheFrame = i_dataBlock.inputValue( ia_cacheFrame, &stat ).asBool();
        CHECK_MSTATUS( stat );
		
        MString cachePath = i_dataBlock.inputValue( ia_cachePath, &stat ).asString();
        CHECK_MSTATUS( stat );
        
        bool readFromCache = i_dataBlock.inputValue( ia_readFromCache, &stat ).asDouble();
        
        if ( cacheFrame )
        {
            if ( !readFromCache )
                writeRodDataToCacheFile( cachePath );
            else
                MGlobal::displayWarning( "Ignoring request to cache frame as reading from cache!\n" );
        }
        
        stat = i_dataBlock.setClean( i_plug );
		if ( !stat )
        {
			stat.perror("WmBunsenRodNode::compute setClean");
			return stat;
		}
    }
    else if ( i_plug == oa_simulatedVertices )
    {        
        // We're here because a Fozzie fur set's wire deformer is asking for information
        // on how to deform the hair.
        
        // First pull all the inputs to make sure we're up to date.
        i_dataBlock.inputValue( ca_simulationSync, &stat ).asBool();
        CHECK_MSTATUS( stat );
        i_dataBlock.inputValue( oa_rodsChanged, &stat ).asBool();
        CHECK_MSTATUS( stat );

        // The above may have been clean so just make sure we actually read time
        i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
		CHECK_MSTATUS( stat );

        MDataHandle simulatedVerticesH = i_dataBlock.outputValue( oa_simulatedVertices, &stat );
        CHECK_MSTATUS( stat );
        MFnVectorArrayData simulatedVerticesArrayData;
        MStatus stat2=simulatedVerticesH.set( simulatedVerticesArrayData.create( &stat ) );
        CHECK_MSTATUS( stat );
        CHECK_MSTATUS( stat2 );
        
        MVectorArray simulatedVerticesArray = simulatedVerticesArrayData.array( &stat );
        CHECK_MSTATUS( stat );
       
        if ( mx_rodData != NULL )
        {            
            size_t numRods = mx_rodData->size();
            unsigned int idx = 0;
            cerr << "TO FOZZIE: numRods = " << numRods << endl;
            for ( size_t r = 0; r < numRods; r++ )
            {
                unsigned int verticesInRod = (*mx_rodData)[ r ]->rod->nv();
                simulatedVerticesArray.setLength( (unsigned int) ( simulatedVerticesArray.length() + verticesInRod ) );
                
                for ( unsigned int v = 0; v < verticesInRod; v++ )
                {
                    Vec3d pos = (*mx_rodData)[ r ]->rod->getVertex( v );
                    simulatedVerticesArray[ idx ] = MVector( pos[0], pos[1], pos[2] );
                    idx++;
                }
            }
        }
        
        cerr << "TO FOZZIE: simulated vertices array is " << simulatedVerticesArray.length() << " elements long " << endl;
        
        if ( simulatedVerticesArray.length() > 5 )
            for ( unsigned int p=0; p<5; p++ )
                cerr << "simulated verts array, point " << p << simulatedVerticesArray[p] << endl;
        
        simulatedVerticesH.setClean();
        i_dataBlock.setClean( i_plug );
        
        /*MPlug sPlug(thisMObject(), oa_simulatedVertices);
        MObject nodeAttr;
        stat = sPlug.getValue( nodeAttr );
        CHECK_MSTATUS( stat );

        MFnPointArrayData ptFn( nodeAttr, &stat);
        CHECK_MSTATUS( stat );

        MPointArray simulatedPts = ptFn.array( &stat );
        CHECK_MSTATUS( stat );

        cout<<"TO FOZZIE: Get Connected simpulated array : "<<simulatedPts.length()<<endl;*/
    }
    else if ( i_plug == oa_nonSimulatedVertices )
    {
        // First pull all the inputs to make sure we're up to date.
        i_dataBlock.inputValue( ca_simulationSync, &stat ).asBool();
        CHECK_MSTATUS( stat );
        i_dataBlock.inputValue( oa_rodsChanged, &stat ).asBool();
        CHECK_MSTATUS( stat );

        // The above may have been clean so just make sure we actually read time
        i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
		CHECK_MSTATUS( stat );

        MDataHandle nonSimulatedVerticesH = i_dataBlock.outputValue( oa_nonSimulatedVertices, &stat );
        CHECK_MSTATUS( stat );        
        MFnVectorArrayData nonSimulatedVerticesArrayData;
        nonSimulatedVerticesH.set( nonSimulatedVerticesArrayData.create( &stat ) );
        CHECK_MSTATUS( stat );

        MVectorArray nonSimulatedVerticesArray = nonSimulatedVerticesArrayData.array( &stat );
        CHECK_MSTATUS( stat );

        if ( mx_rodData != NULL )
        {
            size_t numRods = mx_rodData->size();            
            
            unsigned int idx = 0;
            
            cerr << "TO FOZZIE: numRods = " << numRods << endl;
            for ( size_t r = 0; r < numRods; r++ )
            {
                unsigned int verticesInRod = (*mx_rodData)[ r ]->rod->nv();
                nonSimulatedVerticesArray.setLength( (unsigned int) ( nonSimulatedVerticesArray.length() + verticesInRod ) );
                
                for ( unsigned int v = 0; v < verticesInRod; v++ )
                {
                    Vec3d pos = (*mx_rodData)[ r ]->nextVertexPositions[ v ];
                    nonSimulatedVerticesArray[ idx ] = MVector( pos[0], pos[1], pos[2] );
                    idx++;
                }
            }
            cerr << "TO FOZZIE: NON simulated vertices array is " << nonSimulatedVerticesArray.length() << " elements long " << endl;
        }
        
        nonSimulatedVerticesH.setClean();
        i_dataBlock.setClean( i_plug );
    }
    else if ( i_plug ==  oa_verticesInEachRod )
    {
        // First pull all the inputs to make sure we're up to date.
        i_dataBlock.inputValue( ca_simulationSync, &stat ).asBool();
        CHECK_MSTATUS( stat );
        i_dataBlock.inputValue( oa_rodsChanged, &stat ).asBool();
        CHECK_MSTATUS( stat );

        // The above may have been clean so just make sure we actually read time
        i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
		CHECK_MSTATUS( stat );

        MDataHandle verticesPerRodH = i_dataBlock.outputValue( oa_verticesInEachRod, &stat );
        CHECK_MSTATUS( stat );        
        MFnIntArrayData verticesPerRodArrayData;
        verticesPerRodH.set( verticesPerRodArrayData.create( &stat ) );
        CHECK_MSTATUS( stat );
        
        if ( mx_rodData != NULL )
        {
            MIntArray verticesPerRodArray = verticesPerRodArrayData.array( &stat );
            size_t numRods = mx_rodData->size();
            verticesPerRodArray.setLength( numRods );
            unsigned int idx = 0;
            
            cerr << "TO FOZZIE: numRods = " << numRods << endl;
            for ( size_t r = 0; r < numRods; r++ )
            {
                unsigned int verticesInRod = (*mx_rodData)[ r ]->rod->nv();
                verticesPerRodArray[ r ] = verticesInRod;
            }
        }
        
        verticesPerRodH.setClean();        
        i_dataBlock.setClean( i_plug );
    }
    else
    {
		return MS::kUnknownParameter;
	}

	return MS::kSuccess;
}

size_t WmBunsenRodNode::numberOfRods()
{
    // Fozzie nodes take precidence over nurbs curve inputs. We currently do not handle
    // both on the one node ( although this would probably be easy to do ).
    
    MStatus stat;
    MDataBlock dataBlock = MPxNode::forceCache();
    
    bool fozzieNodeConnected;
    MPlug fozziePlug( thisMObject(), ia_fozzieVertices );
    if ( fozziePlug.isConnected() )
    {
        MDataHandle verticesH = dataBlock.inputValue( ia_fozzieVertices, &stat );
        CHECK_MSTATUS( stat );
        MFnVectorArrayData verticesData( verticesH.data(), &stat );
        CHECK_MSTATUS( stat );
        MVectorArray vertices = verticesData.array( &stat );
    
        return ( vertices.length() / m_cvsPerRod );
    }
    else
    {
        return m_numberOfInputCurves;
    }
}

void WmBunsenRodNode::writeRodDataToCacheFile( MString i_cachePath )
{
    long magicNumber = 2051978;
    
    if ( mx_rodData == NULL )
        return;     // No rod data then nothing to put in the file
    
    if ( i_cachePath == "" )
    {
        MGlobal::displayError( "Empty cache path, not writing anything to disk." );
        return;
    }
    
    MString fileName = i_cachePath + "." + m_currentTime + ".bun";
    
    cerr << "Writing sim data to disk in file: '" << fileName << "'\n";
    
    FILE *fp;
    fp = fopen( fileName.asChar(), "w" );
    if ( fp == NULL )
    {
        MGlobal::displayError( MString( "Problem opening file " + fileName + " for writing." ) );
        return;
    }
    
    // write data from mx_rodData into the file....
    
    // FIXME: create a proper cache file format. This one is bad, it at least needs a proper header.
    size_t numRods = mx_rodData->size();
    fwrite( &numRods, sizeof( size_t ), 1, fp );
  
    for ( size_t r=0; r<numRods; r++ )
    {
        // Write number of vertices in this rod
        BASim::ElasticRod* rod = (*mx_rodData)[ r ]->rod;
     
        if ( rod == NULL )
        {
            cerr << "WTF, rod is NULL \n";
            continue;
        }
        
        size_t numVertices = rod->nv();
        fwrite( &numVertices, sizeof( size_t ), 1, fp );
      
        for ( size_t v=0; v<numVertices; v++ )
        {
            double pos[3];
            
            // Wonder if its safe to write Vec3ds. Need to check what's in them.
            // Really should package all this and write it as one.
            Vec3d vertex = rod->getVertex( v );
            pos[ 0 ] = vertex.x();
            pos[ 1 ] = vertex.y();
            pos[ 2 ] = vertex.z();
            
            fwrite( &pos[0], sizeof( double ), 3, fp );
                
            // We should write velocities and anything else we need into the file too so that
            // we can restart a simulation from a cached file.
        }
    }
    
    fclose ( fp );
}

void WmBunsenRodNode::readRodDataFromCacheFile( MString i_cachePath )
{
    long magicNumber = 2051978;
    
    if ( mx_rodData == NULL )
        return;     // No rod data, nowhere to store rod data
    
    if ( i_cachePath == "" )
    {
        MGlobal::displayError( "Empty cache path, can't read in rod data!" );
        return;
    }
    
    MString fileName = i_cachePath + "." + m_currentTime + ".bun";
    
    cerr << "Reading sim data from disk in file: '" << fileName << "'\n";
    
    FILE *fp;
    fp = fopen( fileName.asChar(), "r" );
    if ( fp == NULL )
    {
        MGlobal::displayError( MString( "Problem opening file " + fileName + " for reading." ) );
        return;
    }
    
    // FIXME: create a proper cache class which I can call just to read or write.
    size_t numRods;
    fread( &numRods, sizeof( size_t ), 1, fp );
    
    if ( numRods != mx_rodData->size() )
    {
        cerr << "Wrong number of rods in file, have not implemented changing it yet\n";
        return;
    }
    
    for ( size_t r=0; r<numRods; r++ )
    {
        BASim::ElasticRod* rod = (*mx_rodData)[ r ]->rod;
        if ( rod == NULL )
        {
            cerr << "WTF, rod is NULL \n";
            continue;
        }
        
        size_t numVertices;
        fread( &numVertices, sizeof( size_t ), 1, fp );
        
        if ( (int)numVertices != rod->nv() )
        {
            // FIXME: We really should be able to do this. The users should be able to just stick any cache file
            // on the attribute and we should allocate rodData for the correct number.
            // Maybe we need a cache node pluggin in to this node?
            MGlobal::displayError( "Can't read in cached rod data as number of vertices do not match rods in scene\n" );
            continue;
        }
        
        for ( size_t v=0; v<numVertices; v++ )
        {
            double pos[3];
            
            // Wonder if its safe to write Vec3ds. Need to check what's in them.
            // Really should package all this and write it as one.
            fread( &pos[0], sizeof( double ), 3, fp );
            
            Vec3d vertex( pos[0], pos[1], pos[2] );
            rod->setVertex( v, vertex );
        }
    }
    
    fclose ( fp );
}

void WmBunsenRodNode::draw( M3dView& i_view, const MDagPath& i_path,
                            M3dView::DisplayStyle i_style,
                            M3dView::DisplayStatus i_status )
{ 
	MStatus stat;
	MObject thisNode = thisMObject();

	MPlug syncPlug( thisNode, ca_simulationSync );
	double d; 
	stat = syncPlug.getValue( d );
	if ( !stat ) 
    {
		stat.perror( "WmBunsenRodNode::draw getting ca_simulationSync" );
		return;
	}
   
	i_view.beginGL(); 
	glPushAttrib( GL_CURRENT_BIT | GL_POINT_BIT | GL_LINE_BIT );
    
	// draw dynamic Hair

	// What did this line do? it was here from the devkit example. Is it to with point colouring
	//view.setDrawColor ( WmBunsenRodNode );

	glPopAttrib();
	i_view.endGL();
}

MStatus WmBunsenRodNode::connectionMade( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc ) 
{    
    MStatus stat;
    MStatus retVal(MS::kUnknownParameter );

    if( plug == ia_nurbsCurves )
    {
        m_numberOfInputCurves++;
    }

    return retVal;
}

MStatus WmBunsenRodNode::connectionBroken( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc )
{
    MStatus stat;
    MStatus retVal( MS::kUnknownParameter );

    if( plug == ia_nurbsCurves )
    {
        m_numberOfInputCurves--;
    }

    return retVal;
}


bool WmBunsenRodNode::isBounded() const
{ 
	return false;
}

void* WmBunsenRodNode::creator()
{
	return new WmBunsenRodNode();
}

/*static */ MStatus WmBunsenRodNode::addNumericAttribute( MObject& i_attribute, MString i_longName, 
    MString i_shortName, MFnNumericData::Type i_type, double i_defaultValue, bool i_isInput )
{
    // Creates a numeric attribute with default attributes
    MStatus stat = MS::kSuccess;
    
    MFnNumericAttribute nAttr;
    i_attribute = nAttr.create( i_longName, i_shortName, i_type, i_defaultValue, &stat );
    if ( !stat )
    {
        cerr << "Failed to create attribute " << i_longName << endl;
        return stat;
    }
    if ( i_isInput )
        nAttr.setWritable( true );
    else
        nAttr.setWritable( false );
    
    stat = addAttribute( i_attribute );
    if ( !stat ) { stat.perror( "addAttribute " + i_longName ); return stat; }

    return stat;
}

/* static */ MStatus WmBunsenRodNode::initialize()
{ 
    MStatus stat;

   // addNumericAttribute( ca_syncAttrs, "syncAttrs", "sya", MFnNumericData::kDouble, 1.0, false );
    addNumericAttribute( oa_rodsChanged, "rodsChanged", "rch", MFnNumericData::kBoolean, true, false );
    addNumericAttribute( ca_simulationSync, "simulationSync", "sis", MFnNumericData::kBoolean, false, false );
    
    addNumericAttribute( ia_simStepTaken, "simStepTaken", "sst", MFnNumericData::kBoolean, false, true );
	stat = attributeAffects( ia_simStepTaken, ca_simulationSync );
	if ( !stat ) { stat.perror( "attributeAffects ia_simStepTaken->ca_syncAttrs" ); return stat; }
    
    {
        MFnUnitAttribute	uAttr;
        ia_time = uAttr.create( "time", "t", MTime( 0.0 ), &stat );
        if ( !stat) 
        {
            stat.perror("create ia_time attribute");
            return stat;
        }
        CHECK_MSTATUS( uAttr.setWritable(true) );
        CHECK_MSTATUS( uAttr.setConnectable(true) );
        CHECK_MSTATUS( uAttr.setStorable(false) );
        stat = addAttribute( ia_time );
        if ( !stat ) { stat.perror( "addAttribute ia_time" ); return stat; }
    }
	stat = attributeAffects( ia_time, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_time->ca_syncAttrs" ); return stat; }
    
    addNumericAttribute( ia_startTime, "startTime", "stt", MFnNumericData::kDouble, 1.0, true );
	stat = attributeAffects( ia_startTime, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_startTime->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_youngsModulus, "youngsModulus", "ymo", MFnNumericData::kDouble, 100000.0, true );
    stat = attributeAffects( ia_youngsModulus, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_youngsModulus->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_shearModulus, "shearModulus", "shm", MFnNumericData::kDouble, 375.0, true );
    stat = attributeAffects( ia_shearModulus, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_shearModulus->ca_syncAttrs" ); return stat; }

	addNumericAttribute( ia_density, "density", "dns", MFnNumericData::kDouble, 0.01, true);
	stat = attributeAffects(ia_density, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_density->ca_syncAttrs" ); return stat; }
    
    addNumericAttribute( ia_minorRadius, "minorRadius", "mir", MFnNumericData::kDouble, 0.1, true );
    stat = attributeAffects( ia_minorRadius, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_minorRadius->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_majorRadius, "majorRadius", "mar", MFnNumericData::kDouble, 0.1, true );
    stat = attributeAffects( ia_majorRadius, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_majorRadius->ca_syncAttrs" ); return stat; }
    
    addNumericAttribute( ia_vertexSpacing, "vertexSpacing", "rvs", MFnNumericData::kDouble, 0.0, true );
    stat = attributeAffects( ia_vertexSpacing, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_majorRadius->ca_syncAttrs" ); return stat; }
 
    addNumericAttribute( ia_cvsPerRod, "cvsPerRod", "cvr", MFnNumericData::kInt, -1, true );
    stat = attributeAffects( ia_cvsPerRod, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_cvsPerRod->ca_syncAttrs" ); return stat; }
    
    addNumericAttribute( ia_cacheFrame, "cacheFrame", "caf", MFnNumericData::kBoolean, false, true );
    stat = attributeAffects( ia_cacheFrame, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_cacheFrame->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_cacheFrame, ca_simulationSync );
	if ( !stat ) { stat.perror( "attributeAffects ia_cacheFrame->ca_syncAttrs" ); return stat; }
        
    addNumericAttribute( ia_readFromCache, "readFromCache", "rfc", MFnNumericData::kBoolean, false, true );
    stat = attributeAffects( ia_readFromCache, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_readFromCache->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_cacheFrame, ca_simulationSync );
	if ( !stat ) { stat.perror( "attributeAffects ia_cacheFrame->ca_syncAttrs" ); return stat; }
    
    {
        MFnTypedAttribute tAttr;
        MFnStringData fnStringData;
        MObject defaultString = fnStringData.create("");
        ia_cachePath = tAttr.create( "cachePath", "cap", MFnData::kString, defaultString, &stat );
        if ( !stat ) 
        {
            stat.perror( "create cachePath attribute" );
            return stat;
        }
        tAttr.setWritable( true );
        tAttr.setReadable( false );
        tAttr.setConnectable( true );
        stat = addAttribute( ia_cachePath );
        if (!stat) { stat.perror( "addAttribute cachePath" ); return stat; }
    }
    stat = attributeAffects( ia_cachePath, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_cachePath->ca_simulationSync" ); return stat; }
    stat = attributeAffects( ia_cachePath, ca_simulationSync );
	if ( !stat ) { stat.perror( "attributeAffects ia_cachePath->ca_simulationSync" ); return stat; }
    
    {
        MFnTypedAttribute   tAttr;  
        ia_nurbsCurves = tAttr.create( "nurbsCurves", "nc",
                                       MFnNurbsCurveData::kNurbsCurve, &stat );
        CHECK_MSTATUS( stat );
        CHECK_MSTATUS( tAttr.setArray(true) );
        CHECK_MSTATUS( tAttr.setReadable(false) );
        CHECK_MSTATUS( tAttr.setWritable(true) );
        CHECK_MSTATUS( tAttr.setDisconnectBehavior( MFnAttribute::kDelete ) );
        CHECK_MSTATUS( tAttr.setUsesArrayDataBuilder(true) );
        stat = addAttribute( ia_nurbsCurves );
        CHECK_MSTATUS( stat );
    }
    stat = attributeAffects( ia_nurbsCurves, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_nurbsCurves->oa_rodsChanged" ); return stat; }
    
    {
        MFnTypedAttribute tAttr;  
        ia_fozzieVertices = tAttr.create( "fozzieVertices", "fov",
                                          MFnData::kVectorArray, &stat );
        CHECK_MSTATUS( stat );
        CHECK_MSTATUS( tAttr.setReadable( false ) );
        CHECK_MSTATUS( tAttr.setWritable( true ) );
        stat = addAttribute( ia_fozzieVertices );
        CHECK_MSTATUS( stat );
    }
    stat = attributeAffects( ia_fozzieVertices, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_fozzieVertices->oa_rodsChanged" ); return stat; }
    
    // Outputs to plug back into Fozzie, should really be a compound attribute
    //
    {
        MFnTypedAttribute tAttr;
        oa_simulatedVertices = tAttr.create( "simulatedVertices", "sve",
                                             MFnData::kVectorArray, &stat );
        CHECK_MSTATUS( stat );
        CHECK_MSTATUS( tAttr.setReadable( true ) );
        CHECK_MSTATUS( tAttr.setWritable( false ) );
        stat = addAttribute( oa_simulatedVertices );
        CHECK_MSTATUS( stat );
    }
    
    {
        MFnTypedAttribute tAttr;  
        oa_nonSimulatedVertices = tAttr.create( "nonSimulatedVertices", "nsv",
                                             MFnData::kVectorArray, &stat );
        CHECK_MSTATUS( stat );
        CHECK_MSTATUS( tAttr.setReadable( true ) );
        CHECK_MSTATUS( tAttr.setWritable( false ) );
        stat = addAttribute( oa_nonSimulatedVertices );
        CHECK_MSTATUS( stat );
    }
    
    {
        MFnTypedAttribute tAttr;  
        oa_verticesInEachRod = tAttr.create( "verticesInEachRod", "ver",
                                             MFnData::kIntArray, &stat );
        CHECK_MSTATUS( stat );
        CHECK_MSTATUS( tAttr.setReadable( true ) );
        CHECK_MSTATUS( tAttr.setWritable( false ) );
        stat = addAttribute( oa_verticesInEachRod );
        CHECK_MSTATUS( stat );
    }
    
    stat = attributeAffects( ia_time, oa_verticesInEachRod );
	stat = attributeAffects( ia_time, oa_nonSimulatedVertices );
	stat = attributeAffects( ia_time, oa_simulatedVertices );
    
	return MS::kSuccess;
}
