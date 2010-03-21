#include "WmBunsenRodNode.hh"
#include "WmFigConnectionNode.hh"

#include <maya/MFnMatrixAttribute.h>
#include <maya/MPlugArray.h>

using namespace BASim;

// Required by Maya to identify node
/* static */ MTypeId WmBunsenRodNode::typeID ( 0x001135, 0x19 );
/* static */ MString WmBunsenRodNode::typeName( "wmFigRodNode" );

//
/* static */ MObject WmBunsenRodNode::ia_time;
/* static */ MObject WmBunsenRodNode::ia_startTime;
/* static */ MObject WmBunsenRodNode::ia_nurbsCurves;
/* static */ MObject WmBunsenRodNode::oa_nurbsCurves;
/* static */ MObject WmBunsenRodNode::ia_fozzieVertices;
/* static */ MObject WmBunsenRodNode::ia_percentageOfFozzieStrands;

// User adjustable rod Options
/* static */ MObject WmBunsenRodNode::ia_cvsPerRod;
/* static */ MObject WmBunsenRodNode::ia_youngsModulus;
/* static */ MObject WmBunsenRodNode::ia_shearModulus;
/* static */ MObject WmBunsenRodNode::ia_viscosity;
/* static */ MObject WmBunsenRodNode::ia_density;
/* static */ MObject WmBunsenRodNode::ia_minorRadius;
/* static */ MObject WmBunsenRodNode::ia_majorRadius;
/* static */ MObject WmBunsenRodNode::ia_vertexSpacing;
/* static */ MObject WmBunsenRodNode::ia_hairSpray;
/* static */ MObject WmBunsenRodNode::ia_hairSprayScaleFactor;
/* static */ MObject WmBunsenRodNode::ia_massDamping;
/* static */ MObject WmBunsenRodNode::ia_drawMaterialFrames;
/* static */ MObject WmBunsenRodNode::ia_lockFirstEdgeToInput;
/* static */ MObject WmBunsenRodNode::ia_userDefinedColors;

// Disk caching
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

/* static */ MObject WmBunsenRodNode::oa_materialFrames;
/* static */ MObject WmBunsenRodNode::oa_undeformedMaterialFrames;
/* static */ MObject WmBunsenRodNode::ia_strandRootFrames;

/* static */ MObject WmBunsenRodNode::oa_numberOfRods;


// For being controlled by external objects and for controlling external objects
/* static */ MObject WmBunsenRodNode::ia_edgeTransforms;
/* static */ MObject WmBunsenRodNode::oa_edgeTransforms;

WmBunsenRodNode::WmBunsenRodNode() : m_massDamping( 10 ), m_initialised( false ),
                                     mx_rodData( NULL ), mx_world( NULL ),
                                     m_numberOfInputCurves( 0 ), m_percentageOfFozzieStrands( 100 ),
                                     m_cvsPerRod( -1 ),
                                     m_cachePath( "" ), m_cacheFilename( "" )
{
    m_rodOptions.YoungsModulus = 1000.0; /* megapascal */
    m_rodOptions.ShearModulus = 340.0;   /* megapascal */
    m_rodOptions.viscosity = 10.0;       /* poise */
    m_rodOptions.density = 1.3;          /* grams per cubic centimeter */
    m_rodOptions.radiusA = 0.05;         /* millimeter */
    m_rodOptions.radiusB = 0.05;         /* millimeter */
    m_strandRootFrames.clear();
    m_controlledEdgeTransforms.clear();
}

WmBunsenRodNode::~WmBunsenRodNode()
{
}

MMatrix WmBunsenRodNode::getRodEdgeMatrix( size_t i_rod, size_t i_edge )
{
    MMatrix identMatrix;
    identMatrix.setToIdentity();
    
    if ( mx_rodData == NULL )
        return identMatrix;
    
    if ( i_rod >= mx_rodData->size() )
        return identMatrix;
    
    if ( (*mx_rodData)[ i_rod ]->rod == NULL )
        return identMatrix;
    
    if ( i_edge >= (*mx_rodData)[ i_rod ]->rod->ne() )
        return identMatrix;
    
    Vec3d edgePos = ( (*mx_rodData)[ i_rod ]->rod->getVertex( i_edge ) + 
                      (*mx_rodData)[ i_rod ]->rod->getVertex( i_edge + 1 ) ) / 2.0;
    
    MMatrix edgeMatrix;
    edgeMatrix( 3, 0 ) = edgePos[ 0 ]; edgeMatrix( 3, 1 ) = edgePos[ 1 ]; edgeMatrix( 3, 2 ) = edgePos[ 2 ];
    
    Vec3d edge = (*mx_rodData)[ i_rod ]->rod->getEdge( i_edge );
    edgeMatrix( 0, 0 ) = edge[ 0 ]; edgeMatrix( 0, 1 ) = edge[ 1 ]; edgeMatrix( 0, 2 ) = edge[ 2 ];
    
    Vec3d material1 = (*mx_rodData)[ i_rod ]->rod->getMaterial1( i_edge );
    edgeMatrix( 1, 0 ) = material1[ 0 ]; edgeMatrix( 1, 1 ) = material1[ 1 ]; edgeMatrix( 1, 2 ) = material1[ 2 ];
    
    Vec3d material2 = (*mx_rodData)[ i_rod ]->rod->getMaterial2( i_edge );
    edgeMatrix( 2, 0 ) = material2[ 0 ]; edgeMatrix( 2, 1 ) = material2[ 1 ]; edgeMatrix( 2, 2 ) = material2[ 2 ];
    
    return edgeMatrix;
}

FILE* WmBunsenRodNode::readNumberOfRodsFromFile( const MString i_cacheFilename, size_t& o_numRodsInFile,
    bool closeFileAfterReading )
{
    FILE *fp = NULL;

    if ( i_cacheFilename == "" )
    {
        MGlobal::displayError( "Empty cache path, can't read in rod data!" );
        return NULL;
    }

    fp = fopen( i_cacheFilename.asChar(), "r" );
    if ( fp == NULL )
    {
        MGlobal::displayError( MString( "Problem opening file " + i_cacheFilename + " for reading." ) );
        return NULL;
    }

/*    int fileFormat;
    fread( &fileFormat, sizeof( int ), 1, fp );
    if ( fileFormat != FILE_FORMAT_VERSION )
    {
        Mglobal::displayError( "Unsupported cache file format version\n" );
        fclose( fp );
        fp = NULL
        o_numRodsInFile = 0;

    }*/

    size_t ret = fread( &o_numRodsInFile, sizeof( size_t ), 1, fp );

    if ( closeFileAfterReading )
    {
        fclose( fp );
        fp = NULL;
    }

    return fp;
}

bool WmBunsenRodNode::readDataFromRodCacheFile( const MString i_cacheFilename, size_t& o_numRodsInFile,
    vector<vector<Vec3d> >& o_rodVertices, vector<vector<Vec3d> >& o_unsimulatedRodVertices )
{

    FILE* fp = readNumberOfRodsFromFile( i_cacheFilename, o_numRodsInFile, false );
    if ( fp == NULL )
        return false;

    o_rodVertices.resize( o_numRodsInFile );
    o_unsimulatedRodVertices.resize( o_numRodsInFile );

    for ( size_t r=0; r<o_numRodsInFile; r++ )
    {
        size_t numVertices;
        size_t ret = fread( &numVertices, sizeof( size_t ), 1, fp );

        o_rodVertices[ r ].resize( numVertices );
        o_unsimulatedRodVertices[ r ].resize( numVertices );

        for ( size_t v=0; v<numVertices; v++ )
        {
            double pos[3];

            ret = fread( &pos[0], sizeof( double ), 3, fp );
            o_rodVertices[ r ][ v ] = Vec3d( pos[0], pos[1], pos[2] );

            ret = fread( &pos[0], sizeof( double ), 3, fp );
            o_unsimulatedRodVertices[ r ][ v ] = Vec3d( pos[0], pos[1], pos[2] );
        }
    }

    fclose( fp );
    return true;
}

MString WmBunsenRodNode::getCacheFilename( MDataBlock& i_dataBlock )
{
    MStatus stat;

    m_cachePath = i_dataBlock.inputValue( ia_cachePath, &stat ).asString();
    CHECK_MSTATUS( stat );

    if ( m_cachePath == "" )
        return "";

    m_currentTime = i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
    CHECK_MSTATUS( stat );

    m_cacheFilename = ( m_cachePath + "." + m_currentTime + ".fig" );

    return m_cacheFilename;
}

/**
    This function initialises the RodData vector to the correct size and the initial
    rod data. It is called when we are about to recreate all the rods so time has just been
    set back to 'startTime'.
*/
void WmBunsenRodNode::initialiseRodData( vector<RodData*>* i_rodData )
{
    mx_rodData = i_rodData;
    MStatus stat;

    // We may get called from _outside_ a compute by the WmBunsenNode initialising all its
    // data and it may be inside a compute and ours may not have been called
    // to get all this data yet. So build our own datablock to use... Some may say this
    // should be done with maya DG connections...

    MDataBlock dataBlock = MPxNode::forceCache();

    bool fozzieNodeConnected;
    MPlug fozziePlug( thisMObject(), ia_fozzieVertices );
    bool readFromCache = dataBlock.inputValue( ia_readFromCache, &stat ).asBool();

    // Reading from the cache takes precidence over anything else connected in
    if ( readFromCache )
    {
        MString cacheFilename = getCacheFilename( dataBlock );

        if ( mx_rodData == NULL )
            return;     // No rod data, nowhere to store rod data

        size_t numRodsInFile;
        vector<vector<Vec3d> > rodVertices;
        vector<vector<Vec3d> > unsimulatedRodVertices;

        if ( !readDataFromRodCacheFile( cacheFilename, numRodsInFile, rodVertices, unsimulatedRodVertices ) )
            return;

        if ( numRodsInFile != mx_rodData->size() )
        {
            MGlobal::displayError( "Wrong number of rods in file, have not implemented changing it yet" );
            MGlobal::displayError( MString( "Expected " ) +  (double)mx_rodData->size() + ", got " + (double)numRodsInFile );
            return;
        }

        for ( size_t i=0; i<numRodsInFile; i++ )
        {
            (*mx_rodData)[i]->rodOptions = m_rodOptions;

            std::string frame = "time";
            if ( frame == "time" )
                (*mx_rodData)[ i ]->rodOptions.refFrame = BASim::ElasticRod::TimeParallel;
            else if (frame == "space")
                (*mx_rodData)[ i ]->rodOptions.refFrame = BASim::ElasticRod::SpaceParallel;

            size_t numVertices = rodVertices[ i ].size();
            (*mx_rodData)[i]->rodOptions.numVertices = (int)numVertices;

            // Make sure we have enough space to store the date for each CV. This should only
            // ever cause a resize when we are called by initialiseRodData().
            (*mx_rodData)[ i ]->undeformedVertexPositions.resize( numVertices );
            (*mx_rodData)[ i ]->initialVertexPositions.resize( numVertices );
            (*mx_rodData)[ i ]->prevVertexPositions.resize( numVertices );
            (*mx_rodData)[ i ]->currVertexPositions.resize( numVertices );
            (*mx_rodData)[ i ]->nextVertexPositions.resize( numVertices );

            (*mx_rodData)[ i ]->undeformedMaterialFrame.resize( numVertices - 1 );

            for ( size_t v=0; v<numVertices; v++ )
            {
                (*mx_rodData)[ i ]->undeformedVertexPositions[ v ]  = rodVertices[ i ][ v ];

                // FIXME: Why do i set the next position to be the unsimulated position,
                // this seems wrong!
                (*mx_rodData)[ i ]->nextVertexPositions[ v ] = unsimulatedRodVertices[ i ][ v ];
            }
        }
    }
    else if ( fozziePlug.isConnected() )
    {
        ///////////////////////////////////////////////////////////////////////////////
        // If there is a Fozzie node connected then it takes control over any connected
        // NURBS curves.
        ///////////////////////////////////////////////////////////////////////////////
        MDataHandle verticesH = dataBlock.inputValue( ia_fozzieVertices, &stat );
        CHECK_MSTATUS( stat );
        MFnVectorArrayData verticesData( verticesH.data(), &stat );
        CHECK_MSTATUS( stat );
        MVectorArray vertices = verticesData.array( &stat );

        unsigned int numStrands = vertices.length() / m_cvsPerRod;
        // The user may have requested that we use less than the total number
        // of Fozzie strands, so scale by that %
        numStrands *= (m_percentageOfFozzieStrands/100.0);

        // store the material frames coming from barbershop
        vector<MaterialFrame> strandRootFrames;
        getStrandRootFrames( dataBlock, strandRootFrames );

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
            (*mx_rodData)[ i ]->currVertexPositions.resize( m_cvsPerRod );
            (*mx_rodData)[ i ]->nextVertexPositions.resize( m_cvsPerRod );

            (*mx_rodData)[ i ]->undeformedMaterialFrame.resize( m_cvsPerRod - 1 );

            std::string frame = "time";
            if ( frame == "time" )
                (*mx_rodData)[ i ]->rodOptions.refFrame = BASim::ElasticRod::TimeParallel;
            else if (frame == "space")
                (*mx_rodData)[ i ]->rodOptions.refFrame = BASim::ElasticRod::SpaceParallel;

            for ( int c=0; c<m_cvsPerRod; c++ )
            {
                MVector cv = vertices[ (int)inputVertexIndex ];
                (*mx_rodData)[ i ]->undeformedVertexPositions[ c ]  = Vec3d( cv.x, cv.y, cv.z );
                inputVertexIndex++;
            }

            // We need to add edge data so that each from the first edge will be locked to the input curve
            if ( strandRootFrames.size() > i && m_lockFirstEdgeToInput )
            {
                // The strand root frames may not have been connected for some reason so this don't
                // rely on having data

                // rod is probably null at this point but rodData will deal with it nicely
                (*mx_rodData)[ i ]->addKinematicEdge( 0, (*mx_rodData)[ i ]->rod, &strandRootFrames[ i ] );
            }
            else
            {
                (*mx_rodData)[ i ]->addKinematicEdge( 0 );
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
            (*mx_rodData)[ i ]->currVertexPositions.resize( nCVs );
            (*mx_rodData)[ i ]->nextVertexPositions.resize( nCVs );

            (*mx_rodData)[ i ]->undeformedMaterialFrame.resize( nCVs - 1 );

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

            if ( m_lockFirstEdgeToInput )
            {
                // It doesnt matter that the edge data is garbage because we don't use the frames
                // from the nurb curve (since there are none!) We just add it so that Beaker
                // knows it isn't to be simulated. Since we didn't pass in a material frame
                // it will know not to use frame locking.
                (*mx_rodData)[ i ]->addKinematicEdge( 0 );
            }
            
            EdgeTransformRodMap::iterator rodIt = m_controlledEdgeTransforms.find( i );
            if ( rodIt != m_controlledEdgeTransforms.end() )
            {
                for ( EdgeTransformMap::iterator edgeIt = m_controlledEdgeTransforms[ i ].begin();
                      edgeIt != m_controlledEdgeTransforms[ i ].end();
                      edgeIt++ )
                {
                    // Now add in ::compute()
                    //(*mx_rodData)[ i ]->addKinematicEdge( edgeIt->first, (*mx_rodData)[ i ]->rod, &(edgeIt->second.materialFrame) );
                    
                    // Set the next positions of these vertices to be wherever the input controller
                    // says they should be.
                    Vec3d position =  edgeIt->second.position;
                    if ( (*mx_rodData)[ i ]->rod != NULL )
                    {
                        double length = (*mx_rodData)[ i ]->rod->getEdge( edgeIt->first ).norm();
                        Vec3d edge = edgeIt->second.materialFrame.m1;
                        edge.normalize();
                        Vec3d start = position - ( edge * length / 2.0 );
                        Vec3d end = position + ( edge * length / 2.0 );
                        (*mx_rodData)[ i ]->undeformedVertexPositions[ edgeIt->first ] = start;
                        (*mx_rodData)[ i ]->undeformedVertexPositions[ edgeIt->first + 1 ] = end;                        
                    }
                }
            }
        }
    }

    // Set things that are the same no matter what input is being used.
    size_t numRods = mx_rodData->size();
    for ( size_t r=0; r<numRods; r++ )
    {
        // Set mass damping for this rod
        (*mx_rodData)[ r ]->massDamping = m_massDamping;

        // Since we're initialising the rod data that means the rod is about to be created. In which
        // case we need to set the current vertex positions since they will not get set by the
        // simulation until the user moves forward a frame.
        size_t numVertices = (*mx_rodData)[ r ]->undeformedVertexPositions.size();
        (*mx_rodData)[ r ]->initialVertexPositions.resize( numVertices );
        for ( size_t v=0; v<numVertices; v++ )
        {
            (*mx_rodData)[ r ]->initialVertexPositions[ v ] = (*mx_rodData)[ r ]->undeformedVertexPositions[ v ];
            (*mx_rodData)[ r ]->prevVertexPositions[ v ] = (*mx_rodData)[ r ]->undeformedVertexPositions[ v ];
            (*mx_rodData)[ r ]->currVertexPositions[ v ] = (*mx_rodData)[ r ]->undeformedVertexPositions[ v ];
            (*mx_rodData)[ r ]->nextVertexPositions[ v ] = (*mx_rodData)[ r ]->undeformedVertexPositions[ v ];
        }
    }

    // We need to make sure we have the spline attr data for the rods since compute may not have been called yet
    updateHairsprayScales( dataBlock );
}

void WmBunsenRodNode::readRodDataFromCacheFile()
{
    if ( mx_rodData == NULL )
        return;     // No rod data, nowhere to store rod data

    size_t numRodsInFile;
    vector<vector<Vec3d> > rodVertices;
    vector<vector<Vec3d> > unsimulatedRodVertices;

    if ( !readDataFromRodCacheFile( m_cacheFilename, numRodsInFile, rodVertices,
                                    unsimulatedRodVertices ) )
        return;

    if ( numRodsInFile != mx_rodData->size() )
    {
        MGlobal::displayError( "Rewind simulation to reset before cache file can be read" );
        return;
    }

    for ( size_t r=0; r<numRodsInFile; r++ )
    {
        BASim::ElasticRod* rod = (*mx_rodData)[ r ]->rod;
        if ( rod == NULL )
        {
            MGlobal::displayError( "Rewind simulation to reset before cache file can be read" );
            continue;
        }

        size_t numVertices = rodVertices[ r ].size();

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
                rod->setVertex( (int)v, rodVertices[ r ][ v ] );

                // FIXME: Why do i set the next position to be the unsimulated position,
                // this seems wrong!
                (*mx_rodData)[ r ]->nextVertexPositions[ v ] = unsimulatedRodVertices[ r ][ v ];
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

        // The user may have requested that we use less than the total number
        // of Fozzie strands, so scale by that %
        numStrands *= (m_percentageOfFozzieStrands/100.0);

        if ( mx_rodData->size() != numStrands )
        {
            MGlobal::displayError( "Number of rods does not equal number of Fozzie strands, rewind simulation to reset" );
            return;
        }

        // store the material frames coming from barbershop
        vector<MaterialFrame> strandRootFrames;
        getStrandRootFrames( dataBlock, strandRootFrames );

        size_t inputVertexIndex = 0;
        for ( size_t i = 0; i < numStrands; i++ )
        {
            BASim::ElasticRod* rod = (*mx_rodData)[ i ]->rod;
            if ( rod != NULL )
            {
                rod->setRadius( m_rodOptions.radiusA, m_rodOptions.radiusB );
                rod->setYoungsModulus( m_rodOptions.YoungsModulus );
                rod->setShearModulus( m_rodOptions.ShearModulus );
                rod->setViscosity( m_rodOptions.viscosity );
                rod->setDensity(m_rodOptions.density);

                for ( int c = 0; c < m_cvsPerRod ; ++c )
                {
                    MVector cv = vertices[ (int)inputVertexIndex ];

                    Vec3d inputCurveVertex( cv.x, cv.y, cv.z );
                    (*mx_rodData)[ i ]->prevVertexPositions[ c ] = (*mx_rodData)[ i ]->nextVertexPositions[ c ];
                    // Set the current position to be the prev as it will be moved forward in substeps by
                    // Beaker as it takes simulation steps.
                    (*mx_rodData)[ i ]->currVertexPositions[ c ] = (*mx_rodData)[ i ]->nextVertexPositions[ c ];
                    (*mx_rodData)[ i ]->nextVertexPositions[ c ] = inputCurveVertex;

                    inputVertexIndex++;
                }

                // The strand root frames may not have been connected for some reason so this don't
                // rely on having data
                if ( strandRootFrames.size() > i && m_lockFirstEdgeToInput )
                {
                    if ( m_currentTime == m_startTime )
                        (*mx_rodData)[ i ]->resetKinematicEdge( 0, (*mx_rodData)[ i ]->rod, strandRootFrames[ i ] );
                    else
                        (*mx_rodData)[ i ]->updateKinematicEdge( 0, strandRootFrames[ i ] );
                }
                else // remove the entry in the map
                {
                    // FIXME:
                    // We haven't released a version of barbershop that gives the frames so
                    // dont remove the edge just yet, leave it and we'll get the first edge
                    // locked we just won't track rotation properly without the frames.

                    //(*mx_rodData)[ i ]->removeKinematicEdge( 0 );
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
                rod->setViscosity( m_rodOptions.viscosity );
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
                    // Set the current position to be the prev as it will be moved forward in substeps by
                    // Beaker as it takes simulation steps.
                    (*mx_rodData)[ i ]->currVertexPositions[ c ] = (*mx_rodData)[ i ]->nextVertexPositions[ c ];
                    (*mx_rodData)[ i ]->nextVertexPositions[ c ] = inputCurveVertex;
                }
             
                EdgeTransformRodMap::iterator rodIt = m_controlledEdgeTransforms.find( i );
                if ( rodIt != m_controlledEdgeTransforms.end() )
                {
                    for ( EdgeTransformMap::iterator edgeIt = m_controlledEdgeTransforms[ i ].begin();
                          edgeIt != m_controlledEdgeTransforms[ i ].end();
                          edgeIt++ )
                    {
                        if ( m_currentTime == m_startTime )
                            (*mx_rodData)[ i ]->resetKinematicEdge( edgeIt->first, (*mx_rodData)[ i ]->rod, (edgeIt->second.materialFrame) );
                        else
                            (*mx_rodData)[ i ]->updateKinematicEdge( edgeIt->first, (edgeIt->second.materialFrame) );
                        
                        // Set the next positions of these vertices to be wherever the input controller
                        // says they should be.
                        Vec3d position =  edgeIt->second.position;
                        if ( (*mx_rodData)[ i ]->rod != NULL )
                        {
                            double length = (*mx_rodData)[ i ]->rod->getEdge( edgeIt->first ).norm();
                            Vec3d edge = edgeIt->second.materialFrame.m1;
                            edge.normalize();
                            Vec3d start = position - ( edge * length / 2.0 );
                            Vec3d end = position + ( edge * length / 2.0 );
                            (*mx_rodData)[ i ]->nextVertexPositions[ edgeIt->first ] = start;
                            (*mx_rodData)[ i ]->nextVertexPositions[ edgeIt->first + 1 ] = end;
                        }
                    }
                }
            }
        }
    }
}

void WmBunsenRodNode::updateHairsprayScales( MDataBlock& i_dataBlock )
{
    if ( mx_rodData == NULL )
        return;

    // FIXME:
    // This is crazy, each rod does not need a copy of this data. It should
    // be stored once.

    MStatus stat;
    MIntArray indexes;
    MFloatArray positions;
    MFloatArray values;
    MIntArray interps;

    Scalar hairSprayScaleFactor = i_dataBlock.inputValue( ia_hairSprayScaleFactor, &stat ).asDouble();
    CHECK_MSTATUS(stat);

    MPlug hsPlug( thisMObject(), ia_hairSpray);
    MRampAttribute weightRamp(hsPlug, &stat);
    CHECK_MSTATUS( stat );
    weightRamp.getEntries(indexes, positions, values, interps, &stat);
    CHECK_MSTATUS( stat );

    // FIXME:
    // doing this so many times is insanely inefficient...
    // Although it will let the user paint per strand attraction values.

    for ( size_t s=0; s<(*mx_rodData).size(); s++ )
    {
        (*mx_rodData)[s]->hairSprayScaleFactor = hairSprayScaleFactor;
        (*mx_rodData)[s]->forceWeightMap.clear();

        if ( indexes.length() > 0 )
        {
            for(unsigned int i=0; i<indexes.length(); i++)
            {
                (*mx_rodData)[s]->forceWeightMap[positions[i]] = std::pair<float, int16_t>(values[i], interps[i]);
            }
        }
        else
        {
          // There is a weird bug where it only lets me add one point sometimes so default to 4
          // points then i can move them.

            positions.append( 0.0 );
            values.append( 0 );
            interps.append( 2 );

            (*mx_rodData)[s]->forceWeightMap[positions[0]] = std::pair< float, int16_t>( 0.0, 0 );
            weightRamp.addEntries( positions, values, interps, &stat );

            positions.append( 0.25 );
            values.append( 0. );
            interps.append( 2 );

            (*mx_rodData)[s]->forceWeightMap[positions[1]] = std::pair< float, int16_t>( 0.25, 0 );
            weightRamp.addEntries( positions, values, interps, &stat );

            positions.append( 0.5 );
            values.append( 0 );
            interps.append( 2 );

            (*mx_rodData)[s]->forceWeightMap[positions[2]] = std::pair< float, int16_t>( 0.5, 0 );
            weightRamp.addEntries( positions, values, interps, &stat );

            positions.append( 0.75 );
            values.append( 0 );
            interps.append( 2 );

            (*mx_rodData)[s]->forceWeightMap[positions[3]] = std::pair< float, int16_t>( 0.75, 0 );
            weightRamp.addEntries( positions, values, interps, &stat );
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

        m_rodOptions.YoungsModulus = 1e7 * i_dataBlock.inputValue( ia_youngsModulus, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_rodOptions.ShearModulus = 1e7 * i_dataBlock.inputValue( ia_shearModulus, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_rodOptions.viscosity = i_dataBlock.inputValue( ia_viscosity, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_rodOptions.density = i_dataBlock.inputValue( ia_density, &stat).asDouble();
        CHECK_MSTATUS( stat );
        m_rodOptions.radiusA = 1e-1 * i_dataBlock.inputValue( ia_minorRadius, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_rodOptions.radiusB = 1e-1 * i_dataBlock.inputValue( ia_majorRadius, &stat ).asDouble();
        CHECK_MSTATUS( stat );
        m_massDamping = i_dataBlock.inputValue( ia_massDamping, &stat ).asDouble();
        CHECK_MSTATUS( stat );

        m_lockFirstEdgeToInput = i_dataBlock.inputValue( ia_lockFirstEdgeToInput, &stat ).asBool();
        CHECK_MSTATUS( stat );

        // This will get and store the cache filename for the current frame
        getCacheFilename( i_dataBlock );

        //m_vertexSpacing = i_dataBlock.inputValue( ia_vertexSpacing, &stat ).asDouble();
        //CHECK_MSTATUS( stat );
        m_cvsPerRod = i_dataBlock.inputValue( ia_cvsPerRod, &stat ).asInt();
        CHECK_MSTATUS( stat );

        m_percentageOfFozzieStrands = i_dataBlock.inputValue( ia_percentageOfFozzieStrands, &stat ).asInt();
        CHECK_MSTATUS( stat );
		
        //////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////
        
        MDataHandle inputCurveH;
        MObject inputCurveObj;    
        MArrayDataHandle inEdgeArrayH = i_dataBlock.inputArrayValue( ia_edgeTransforms, &stat );
        CHECK_MSTATUS(stat);
       
        size_t numNodesConnected = inEdgeArrayH.elementCount();
    
        for ( unsigned int e = 0; e < numNodesConnected; e++ )
        {
            // Ask for the input on them all to make sure the connection nodes evaluate
            inEdgeArrayH.jumpToElement( e );
            inEdgeArrayH.inputValue( &stat );
            CHECK_MSTATUS( stat );
        }
        
        // I don't think we need the above since the below will pull on all the attrs. Is the stuff
        // below even safe in a compute() since we're not using the data block? I can't see how
        // it can cause any evaluation that the above wouldn't.
        
        // First get rid of all the constraints as some of them may have been removed.
        // We could do this in connectionBroken but there is always weirdness when loading
        // and tracking stuff in connection made/broken.
        
        // FIXME: This is pretty darn innefficient, deleting then recreating
        // Although it does let us create and delete them at any frame in the simulation
        
        if ( mx_rodData != NULL )
        {            
            for ( EdgeTransformRodMap::iterator rodIt = m_controlledEdgeTransforms.begin();
                rodIt != m_controlledEdgeTransforms.end(); rodIt++ )
            {
                for ( EdgeTransformMap::iterator edgeIt = rodIt->second.begin();
                      edgeIt != rodIt->second.end();
                      edgeIt++ )
                {
                    //if ( (*mx_rodData)[ rodIt->first ]->rod != NULL )
                        (*mx_rodData)[ rodIt->first ]->removeKinematicEdge( edgeIt->first );
                }
                m_controlledEdgeTransforms[ rodIt->first ].clear();
            }
            m_controlledEdgeTransforms.clear();
        }
        
        MPlug myEdgePlugs( thisMObject(), ia_edgeTransforms );
        if ( myEdgePlugs.isArray() )
        {
            for ( unsigned int p=0; p<myEdgePlugs.numElements(); p++ )
            {
                MPlug edgePlug = myEdgePlugs.elementByPhysicalIndex( p, &stat );
                CHECK_MSTATUS( stat );
                
                MPlugArray connectionNodePlugs;
                if ( edgePlug.connectedTo( connectionNodePlugs, true, false, &stat ) )
                {
                    MPlug connectionNodePlug = connectionNodePlugs[0];
                    MObject connectionNodeObj = connectionNodePlug.node( &stat );
                    CHECK_MSTATUS( stat );
                    
                    MFnDependencyNode connectionNodeDepFn( connectionNodeObj );
                    WmFigConnectionNode* connectionNode = static_cast<WmFigConnectionNode*>( connectionNodeDepFn.userNode() );
                    
                    if ( connectionNode != NULL )
                    {
                        unsigned int rod, edge;
                        EdgeTransform edgeTransform;
                        connectionNode->getControlledRodInfo( rod, edge, edgeTransform );
                        
                        m_controlledEdgeTransforms[ rod ][ edge ] = edgeTransform;
                        
                        (*mx_rodData)[ rod ]->addKinematicEdge( edge, (*mx_rodData)[ rod ]->rod, &(edgeTransform.materialFrame) );
                    }
                }
            }
        }
        inEdgeArrayH.setClean();
        
        //////////////////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////////
        
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

                readRodDataFromCacheFile();
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
            updateHairsprayScales( i_dataBlock );
        }
        
        
        stat = i_dataBlock.setClean( i_plug );
        if ( !stat )
        {
            stat.perror("WmBunsenRodNode::compute setClean");
            return stat;
        }
    }
    else if ( i_plug == ca_simulationSync || i_plug == oa_numberOfRods )
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
                writeRodDataToCacheFile();
            else
                MGlobal::displayWarning( "Ignoring request to cache frame as reading from cache!\n" );
        }
        
        MDataHandle numRodsH = i_dataBlock.outputValue( oa_numberOfRods, &stat);
        CHECK_MSTATUS( stat );
        if ( mx_rodData != NULL )
            numRodsH.set( (int)mx_rodData->size() );
        else
            numRodsH.set( 0 );
        
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
            verticesPerRodArray.setLength( (unsigned int)numRods );
            unsigned int idx = 0;

            cerr << "TO FOZZIE: numRods = " << numRods << endl;
            for ( size_t r = 0; r < numRods; r++ )
            {
                unsigned int verticesInRod = (*mx_rodData)[ r ]->rod->nv();
                verticesPerRodArray[ (int)r ] = verticesInRod;
            }
        }

        verticesPerRodH.setClean();
        i_dataBlock.setClean( i_plug );
    }
    else if ( i_plug == oa_materialFrames )
    {
        //FIXME: All the barbershop output attributes can be done at the same time and then set
        // all clean. We don't need to duplicate all this code.
        // As soon as we have the pipeline working end to end then go back and refactor all this.

        // First pull all the inputs to make sure we're up to date.
        i_dataBlock.inputValue( ca_simulationSync, &stat ).asBool();
        CHECK_MSTATUS( stat );
        i_dataBlock.inputValue( oa_rodsChanged, &stat ).asBool();
        CHECK_MSTATUS( stat );

        // The above may have been clean so just make sure we actually read time
        i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
        CHECK_MSTATUS( stat );

        MDataHandle materialFramesH = i_dataBlock.outputValue( oa_materialFrames, &stat );
        CHECK_MSTATUS( stat );
        MFnVectorArrayData materialFramesArrayData;
        MStatus stat2=materialFramesH.set( materialFramesArrayData.create( &stat ) );
        CHECK_MSTATUS( stat );
        CHECK_MSTATUS( stat2 );

        MVectorArray materialFramesArray = materialFramesArrayData.array( &stat );
        CHECK_MSTATUS( stat );

        if ( mx_rodData != NULL )
        {
            size_t numRods = mx_rodData->size();
            unsigned int idx = 0;
            for ( size_t r = 0; r < numRods; r++ )
            {
                ElasticRod* rod = (*mx_rodData)[ r ]->rod;
                if ( rod != NULL )
                {
                    unsigned int edgesInRod = rod->ne();
                    materialFramesArray.setLength( (unsigned int) ( materialFramesArray.length() + edgesInRod*3 ) );

                    for ( unsigned int e = 0; e < edgesInRod; e++ )
                    {
                        Vec3d m1 =  rod->getMaterial1( e );
                        Vec3d m2 =  rod->getMaterial2( e );
                        Vec3d m3 =  rod->getEdge( e );
                        m3.normalize();

                        materialFramesArray[ idx ] = MVector( m1[0], m1[1], m1[2] );
                        idx++;
                        materialFramesArray[ idx ] = MVector( m2[0], m2[1], m2[2] );
                        idx++;
                        materialFramesArray[ idx ] = MVector( m3[0], m3[1], m3[2] );
                        idx++;
                    }
                }
            }
        }

        materialFramesH.setClean();
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
    else if ( i_plug == oa_undeformedMaterialFrames )
    {
        // First pull all the inputs to make sure we're up to date.
        i_dataBlock.inputValue( ca_simulationSync, &stat ).asBool();
        CHECK_MSTATUS( stat );
        i_dataBlock.inputValue( oa_rodsChanged, &stat ).asBool();
        CHECK_MSTATUS( stat );

        // The above may have been clean so just make sure we actually read time
        double time = i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
        CHECK_MSTATUS( stat );
        double startTime = i_dataBlock.inputValue( ia_startTime, &stat ).asDouble();
        CHECK_MSTATUS( stat );

        vector<MaterialFrame> strandRootFrames;

        if ( stat == MS::kSuccess )
            getStrandRootFrames( i_dataBlock, strandRootFrames );
        else
        {
            // Make sure there is no data in this so we don't try and use it below.
            strandRootFrames.clear();
        }

        MDataHandle undeformedMaterialFramesH = i_dataBlock.outputValue( oa_undeformedMaterialFrames, &stat );
        CHECK_MSTATUS( stat );
        MFnVectorArrayData undeformedMaterialFramesData;
        MStatus stat2=undeformedMaterialFramesH.set( undeformedMaterialFramesData.create( &stat ) );
        CHECK_MSTATUS( stat );
        CHECK_MSTATUS( stat2 );

        MVectorArray undeformedMaterialFramesArray = undeformedMaterialFramesData.array( &stat );
        CHECK_MSTATUS( stat );

        // If time==startTime then we need to store the strand root frames as these
        // are the strand root frames in the groom pose
        if ( time == startTime )
        {
            m_strandRootFrames.resize( strandRootFrames.size() );
            m_strandRootFrames = strandRootFrames;

            // since time==startTime and the rods should have been initialised before we
            // got here then we can safely store the material frames as the frames
            // in the groom pose

            if ( mx_rodData != NULL )
            {
                size_t numRods = mx_rodData->size();
                unsigned int idx = 0;
                for ( size_t r = 0; r < numRods; r++ )
                {
                    ElasticRod* rod = (*mx_rodData)[ r ]->rod;
                    if ( rod != NULL )
                    {
                        undeformedMaterialFramesArray.setLength( (unsigned int) ( undeformedMaterialFramesArray.length() + rod->ne()*3 ) );

                        for ( size_t e=0; e<(size_t)rod->ne(); e++ )
                        {
                            // currently we only store the undeformed frame for the first vertex
                            Vec3d m1 =  rod->getMaterial1( e );
                            Vec3d m2 =  rod->getMaterial2( e );
                            Vec3d m3 =  rod->getEdge( e );
                            m3.normalize();

                            (*mx_rodData)[ r ]->undeformedMaterialFrame[ e ].m1 = m1;
                            (*mx_rodData)[ r ]->undeformedMaterialFrame[ e ].m2 = m2;
                            (*mx_rodData)[ r ]->undeformedMaterialFrame[ e ].m3 = m3;

                            undeformedMaterialFramesArray[ idx ] = MVector( m1[0], m1[1], m1[2] );
                            idx++;
                            undeformedMaterialFramesArray[ idx ] = MVector( m2[0], m2[1], m2[2] );
                            idx++;
                            undeformedMaterialFramesArray[ idx ] = MVector( m3[0], m3[1], m3[2] );
                            idx++;
                        }
                    }
                }
            }
        }
        else
        {
            // We're not at startTime so we must be simulating and need to work out what the
            // unsimulated material frames would be at this point based on the input curve.

            if ( mx_rodData != NULL )
            {
                size_t numRods = mx_rodData->size();

                // if we have no strand root frames for this frame then we can't do anything.
                // We use < rather than == because we can choose to only use a percentage
                // of the input strands from barbershop in which case rods is < strandRootFrames
                if  ( strandRootFrames.size() <= numRods )
                {
                    MGlobal::displayError( "No strand root frames to use in deformation" );
                }
                else
                {
                    unsigned int idx = 0;
                    for ( size_t r = 0; r < numRods; r++ )
                    {
                        ElasticRod* rod = (*mx_rodData)[ r ]->rod;
                        if ( rod != NULL )
                        {
                        /*  m_strandRootFrames is what we create the first matrix from then we create
                            the second matrix from the strand root frames we were given this frame.
                            The rod frames should not be in either matrix!!!!*/

                            Vec3d im1 = m_strandRootFrames[r].m1;
                            Vec3d im2 = m_strandRootFrames[r].m2;
                            Vec3d im3 = m_strandRootFrames[r].m3;

                            Vec3d cm1 = strandRootFrames[r].m1;
                            Vec3d cm2 = strandRootFrames[r].m2;
                            Vec3d cm3 = strandRootFrames[r].m3;

                            double dim[4][4] = {{ im1[0], im1[1], im1[2], 0.0 },
                                               { im2[0], im2[1], im2[2], 0.0 },
                                               { im3[0], im3[1], im3[2], 0.0 },
                                               {    0.0,    0.0,    0.0, 1.0 }};
                            MMatrix im( dim );

                            double dcm[4][4] = {{ cm1[0], cm1[1], cm1[2], 0.0 },
                                                { cm2[0], cm2[1], cm2[2], 0.0 },
                                                { cm3[0], cm3[1], cm3[2], 0.0 },
                                                {    0.0,    0.0,    0.0, 1.0 }};
                            MMatrix cm( dcm );

                            undeformedMaterialFramesArray.setLength( (unsigned int) ( undeformedMaterialFramesArray.length() + rod->ne()*3 ) );

                            for ( size_t e=0; e<(size_t)rod->ne(); e++ )
                            {
                                // currently we only store the undeformed frame for the first vertex
                                Vec3d m1 =  (*mx_rodData)[ r ]->undeformedMaterialFrame[ r ].m1;
                                MVector mayaM1( m1[0], m1[1], m1[2] );
                                Vec3d m2 =   (*mx_rodData)[ r ]->undeformedMaterialFrame[ r ].m2;
                                MVector mayaM2( m2[0], m2[1], m2[2] );
                                Vec3d m3 = (*mx_rodData)[ r ]->undeformedMaterialFrame[ r ].m3;
                                MVector mayaM3( m3[0], m3[1], m3[2] );

                                // remove initial transform and apply current...
                                mayaM1 = mayaM1 * im.inverse() * cm;
                                mayaM2 = mayaM2 * im.inverse() * cm;
                                mayaM3 = mayaM3 * im.inverse() * cm;

                                undeformedMaterialFramesArray[ idx ] = mayaM1;
                                idx++;
                                undeformedMaterialFramesArray[ idx ] = mayaM2;
                                idx++;
                                undeformedMaterialFramesArray[ idx ] = mayaM3;
                                idx++;
                            }
                        }
                    }
                }
            }

        }
        undeformedMaterialFramesH.setClean();
        i_dataBlock.setClean( i_plug );
    }
    else if ( i_plug == oa_edgeTransforms )
    {
        i_dataBlock.inputValue( ia_edgeTransforms, &stat );
        CHECK_MSTATUS( stat );
        i_dataBlock.inputValue( ia_simStepTaken, &stat );
        CHECK_MSTATUS( stat );
        
        // Due to the complexities of the number of vertices changing on a rod we
        // dont actually pass any info here, we merely signal that the output has changed.
        // Then it is up to the connection node to grab the correct data from mx_rodData.
        MDataHandle edgeTransforms = i_dataBlock.outputValue( oa_edgeTransforms, &stat );
        CHECK_MSTATUS( stat );
        
        edgeTransforms.setClean();
        i_dataBlock.setClean( i_plug );
    }
    else
    {
		return MS::kUnknownParameter;
	}

	return MS::kSuccess;
}

void WmBunsenRodNode::getStrandRootFrames( MDataBlock& i_dataBlock, vector<MaterialFrame>& o_strandRootFrames )
{
    MStatus stat;

    MDataHandle strandRootFramesH = i_dataBlock.inputValue( ia_strandRootFrames, &stat );
    CHECK_MSTATUS( stat );
    MFnVectorArrayData strandRootFrameVecData( strandRootFramesH.data(), &stat );
    
    if ( !stat )    // Strand root frames are not connected so just leave
        return;
        
    CHECK_MSTATUS( stat );

    MVectorArray strandRootFrameVec;
    strandRootFrameVec = strandRootFrameVecData.array( &stat );

    o_strandRootFrames.resize( strandRootFrameVec.length()/3 );
    MVector v;
    size_t idx = 0;
    for ( size_t rIdx=0; rIdx<o_strandRootFrames.size(); rIdx++ )
    {
        v = strandRootFrameVec[idx++];
        o_strandRootFrames[rIdx].m1 = Vec3d( v[0], v[1], v[2] );
        v = strandRootFrameVec[idx++];
        o_strandRootFrames[rIdx].m2 = Vec3d( v[0], v[1], v[2] );
        v = strandRootFrameVec[idx++];
        o_strandRootFrames[rIdx].m3 = Vec3d( v[0], v[1], v[2] );
    }
}

size_t WmBunsenRodNode::numberOfRods()
{
    // Fozzie nodes take precidence over nurbs curve inputs. We currently do not handle
    // both on the one node ( although this would probably be easy to do ).

    MStatus stat;
    MDataBlock dataBlock = MPxNode::forceCache();

    bool fozzieNodeConnected;
    MPlug fozziePlug( thisMObject(), ia_fozzieVertices );
    bool readFromCache = dataBlock.inputValue( ia_readFromCache, &stat ).asDouble();

    if ( readFromCache )
    {
        MString fileName = getCacheFilename( dataBlock );

        size_t numRods = 0;
        readNumberOfRodsFromFile( fileName, numRods );

        return numRods;

    }
    else if ( fozziePlug.isConnected() )
    {
        MDataHandle verticesH = dataBlock.inputValue( ia_fozzieVertices, &stat );
        CHECK_MSTATUS( stat );
        MFnVectorArrayData verticesData( verticesH.data(), &stat );
        CHECK_MSTATUS( stat );
        MVectorArray vertices = verticesData.array( &stat );

        return ( vertices.length() / m_cvsPerRod ) * (m_percentageOfFozzieStrands/100.0);
    }
    else
    {
        return m_numberOfInputCurves;
    }
}

void WmBunsenRodNode::writeRodDataToCacheFile()
{
    if ( mx_rodData == NULL )
        return;     // No rod data then nothing to put in the file

    if ( m_cacheFilename == "" )
    {
        MGlobal::displayError( "Empty cache path, not writing anything to disk." );
        return;
    }

    cerr << "Writing sim data to disk in file: '" << m_cacheFilename << "'\n";

    FILE *fp;
    fp = fopen( m_cacheFilename.asChar(), "w" );
    if ( fp == NULL )
    {
        MGlobal::displayError( MString( "Problem opening file " + m_cacheFilename + " for writing." ) );
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
            Vec3d vertex = rod->getVertex( (unsigned int)v );
            pos[ 0 ] = vertex.x();
            pos[ 1 ] = vertex.y();
            pos[ 2 ] = vertex.z();

            fwrite( &pos[0], sizeof( double ), 3, fp );

            // Now store the unsimulated positions so that we can send them to Barbershop
            // when we want to simulate with it.
            vertex = (*mx_rodData)[ r ]->nextVertexPositions[ v ];
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

    if ( mx_rodData != NULL )
    {
        for ( size_t r=0; r<mx_rodData->size(); r++ )
            (*mx_rodData)[ r ]->rodRenderer->render();
    }

    MPlug drawMaterialFramesPlug( thisNode, ia_drawMaterialFrames );
	bool draw;
	drawMaterialFramesPlug.getValue( draw );
	
    /*if ( draw && mx_rodData != NULL )
    {
        size_t numRods = mx_rodData->size();
        unsigned int idx = 0;
        for ( size_t r = 0; r < numRods; r++ )
        {
            ElasticRod* rod = (*mx_rodData)[ r ]->rod;
            if ( rod != NULL )
            {
                unsigned int edgesInRod = rod->ne();

                glBegin( GL_LINES );
                for ( unsigned int e = 0; e < edgesInRod; e++ )
                {
                    Vec3d m1 =  rod->getMaterial1( e );
                    Vec3d m2 =  rod->getMaterial2( e );
                    Vec3d m3 =  rod->getEdge( e );
                    m3.normalize();

                    Vec3d p = rod->getVertex( e );
                    Vec3d p1 = rod->getVertex( e + 1 );
                    p = ( p + p1 ) / 2.0;

                    glColor3d(1,0,0);
                    glVertex3d( p[0], p[1], p[2] );
                    glVertex3d( p[0] + m1[0], p[1] + m1[1], p[2] + m1[2] );

                    glColor3d(0,1,0);
                    glVertex3d( p[0], p[1], p[2] );
                    glVertex3d( p[0] + m2[0], p[1] + m2[1], p[2] + m2[2] );

                    glColor3d(0,0,1);
                    glVertex3d( p[0], p[1], p[2] );
                    glVertex3d( p[0] + m3[0], p[1] + m3[1], p[2] + m3[2] );
                }
                glEnd();

                if ( m_strandRootFrames.size() > r )
                {
                    // Temporary drawing of frames from barbershop
                    glLineWidth( 5.0 );
                    glBegin( GL_LINES );

                    Vec3d m1 =  m_strandRootFrames[ r ].m1;
                    Vec3d m2 =  m_strandRootFrames[ r ].m2;
                    Vec3d m3 =  m_strandRootFrames[ r ].m3;

                    Vec3d p = rod->getVertex( 0 );

                    glColor3d(1,0,0);
                    glVertex3d( p[0], p[1], p[2] );
                    glVertex3d( p[0] + m1[0], p[1] + m1[1], p[2] + m1[2] );

                    glColor3d(0,1,0);
                    glVertex3d( p[0], p[1], p[2] );
                    glVertex3d( p[0] + m2[0], p[1] + m2[1], p[2] + m2[2] );

                    glColor3d(0,0,1);
                    glVertex3d( p[0], p[1], p[2] );
                    glVertex3d( p[0] + m3[0], p[1] + m3[1], p[2] + m3[2] );
                    glEnd();
                    glLineWidth( 1.0 );
                }

                glBegin( GL_LINES );
                for ( unsigned int e = 0; e < edgesInRod; e++ )
                {
                    Vec3d m1 = (*mx_rodData)[ r ]->undeformedMaterialFrame[ e ].m1;
                    Vec3d m2 = (*mx_rodData)[ r ]->undeformedMaterialFrame[ e ].m2;
                    Vec3d m3 = (*mx_rodData)[ r ]->undeformedMaterialFrame[ e ].m3;

                    Vec3d p = rod->getVertex( e );
                    Vec3d p1 = rod->getVertex( e + 1 );
                    p = ( p + p1 ) / 2.0;

                    glColor3d(1,0,0);
                    glVertex3d( p[0], p[1], p[2] );
                    glVertex3d( p[0] + m1[0], p[1] + m1[1], p[2] + m1[2] );

                    glColor3d(0,1,0);
                    glVertex3d( p[0], p[1], p[2] );
                    glVertex3d( p[0] + m2[0], p[1] + m2[1], p[2] + m2[2] );

                    glColor3d(0,0,1);
                    glVertex3d( p[0], p[1], p[2] );
                    glVertex3d( p[0] + m3[0], p[1] + m3[1], p[2] + m3[2] );
                }
                glEnd();
            }
        }
    }
	}*/

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
    MString i_shortName, MFnNumericData::Type i_type, double i_defaultValue, bool i_isInput,
    bool i_isArray )
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
    
    if ( i_isArray )
        nAttr.setArray( true );

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

    addNumericAttribute( ia_youngsModulus, "youngsModulus", "ymo", MFnNumericData::kDouble, 1000.0, true );
    stat = attributeAffects( ia_youngsModulus, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_youngsModulus->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_shearModulus, "shearModulus", "shm", MFnNumericData::kDouble, 340.0, true );
    stat = attributeAffects( ia_shearModulus, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_shearModulus->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_viscosity, "internalDamping", "ind", MFnNumericData::kDouble, 10.0, true );
    stat = attributeAffects( ia_viscosity, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_viscosity->ca_syncAttrs" ); return stat; }

	addNumericAttribute( ia_density, "density", "dns", MFnNumericData::kDouble, 1.3, true);
	stat = attributeAffects(ia_density, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_density->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_minorRadius, "minorRadius", "mir", MFnNumericData::kDouble, 0.05, true );
    stat = attributeAffects( ia_minorRadius, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_minorRadius->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_majorRadius, "majorRadius", "mar", MFnNumericData::kDouble, 0.05, true );
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
	if ( !stat ) { stat.perror( "attributeAffects ia_cacheFrame->oa_rodsChanged" ); return stat; }
    stat = attributeAffects( ia_cacheFrame, ca_simulationSync );
	if ( !stat ) { stat.perror( "attributeAffects ia_cacheFrame->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_readFromCache, "readFromCache", "rfc", MFnNumericData::kBoolean, false, true );
    stat = attributeAffects( ia_readFromCache, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_readFromCache->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_cacheFrame, ca_simulationSync );
    if ( !stat ) { stat.perror( "attributeAffects ia_cacheFrame->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_hairSprayScaleFactor, "hairSprayScaleFactor", "hsf", MFnNumericData::kDouble, 1.0, true );
    stat = attributeAffects( ia_hairSprayScaleFactor, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_hairSprayScaleFactor->oa_rodsChanged" ); return stat; }

    addNumericAttribute( ia_massDamping, "massDamping", "mda", MFnNumericData::kDouble, 10.0, true );
    stat = attributeAffects( ia_massDamping, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_massDamping->oa_rodsChanged" ); return stat; }

    addNumericAttribute( ia_percentageOfFozzieStrands, "percentageOfFozzieStrands", "pfs", MFnNumericData::kInt, 5, true );
    stat = attributeAffects( ia_percentageOfFozzieStrands, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_fozzieVertices->oa_rodsChanged" ); return stat; }

    addNumericAttribute( ia_drawMaterialFrames, "drawMaterialFrames", "dmf", MFnNumericData::kBoolean, false, true );
    stat = attributeAffects( ia_drawMaterialFrames, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_drawMaterialFrames->oa_rodsChanged" ); return stat; }

    addNumericAttribute( ia_lockFirstEdgeToInput, "lockFirstEdgeToInput", "lfe", MFnNumericData::kBoolean, true, true );
    stat = attributeAffects( ia_lockFirstEdgeToInput, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_lockFirstEdgeToInput->oa_rodsChanged" ); return stat; }

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


    {
        ia_hairSpray = MRampAttribute::createCurveRamp("hairSpray", "hs");
        stat = addAttribute( ia_hairSpray );
        CHECK_MSTATUS( stat );
    }
    stat = attributeAffects( ia_hairSpray, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_hairSpray->oa_rodsChanged" ); return stat; }

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

    {
        MFnTypedAttribute tAttr;
        oa_materialFrames = tAttr.create( "materialFrames", "maf",
                                           MFnData::kVectorArray, &stat );
        CHECK_MSTATUS( stat );
        CHECK_MSTATUS( tAttr.setReadable( true ) );
        CHECK_MSTATUS( tAttr.setWritable( false ) );
        stat = addAttribute( oa_materialFrames );
        CHECK_MSTATUS( stat );
    }
    {
        MFnTypedAttribute tAttr;
        oa_undeformedMaterialFrames = tAttr.create( "undeformedMaterialFrames", "umf",
                                           MFnData::kVectorArray, &stat );
        CHECK_MSTATUS( stat );
        CHECK_MSTATUS( tAttr.setReadable( true ) );
        CHECK_MSTATUS( tAttr.setWritable( false ) );
        stat = addAttribute( oa_undeformedMaterialFrames );
        CHECK_MSTATUS( stat );
    }
    {
        MFnTypedAttribute tAttr;
        ia_strandRootFrames = tAttr.create( "strandRootFrames", "srf",
                                           MFnData::kVectorArray, &stat );
        CHECK_MSTATUS( stat );
        CHECK_MSTATUS( tAttr.setReadable( false ) );
        CHECK_MSTATUS( tAttr.setWritable( true ) );
        stat = addAttribute( ia_strandRootFrames );
        CHECK_MSTATUS( stat );
    }

    stat = attributeAffects( ia_time, oa_verticesInEachRod );
	stat = attributeAffects( ia_time, oa_nonSimulatedVertices );
	stat = attributeAffects( ia_time, oa_simulatedVertices );
    stat = attributeAffects( ia_time, oa_materialFrames );
    stat = attributeAffects( ia_time, oa_undeformedMaterialFrames );
    stat = attributeAffects( ia_startTime, oa_undeformedMaterialFrames );
    stat = attributeAffects( ia_strandRootFrames, oa_undeformedMaterialFrames );
 
    addNumericAttribute( oa_numberOfRods, "numberOfRods", "nor", MFnNumericData::kInt, 0, false );
    stat = attributeAffects( ia_simStepTaken, oa_numberOfRods );
	if ( !stat ) { stat.perror( "attributeAffects ia_simStepTaken->oa_numberOfRods" ); return stat; }
 
    // Controlling and being controlled by external objects
    /*{
        MFnMatrixAttribute mAttr;
        ia_edgeTransforms = mAttr.create( "inEdgeTransforms", "iet", MFnMatrixAttribute::kDouble, &stat );
        if ( !stat ) 
        {
            stat.perror("create ia_edgeTransforms attribute");
            return stat;
        }
        CHECK_MSTATUS( mAttr.setWritable( true ) );
        CHECK_MSTATUS( mAttr.setReadable( false ) );
        CHECK_MSTATUS( mAttr.setArray( true ) );
        stat = addAttribute( ia_edgeTransforms );
        if ( !stat ) { stat.perror( "addAttribute ia_edgeTransforms" ); return stat; }
    }*/
    
    addNumericAttribute( ia_edgeTransforms, "edgeTransforms", "iet", MFnNumericData::kBoolean, false, true, true );
    stat = attributeAffects( ia_edgeTransforms, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_edgeTransforms->oa_rodsChanged" ); return stat; }
 
 /*   {
        MFnMatrixAttribute mAttr;
        oa_edgeTransforms = mAttr.create( "outEdgeTransforms", "oet", MFnMatrixAttribute::kDouble, &stat );
        if ( !stat ) 
        {
            stat.perror("create oa_edgeTransforms attribute");
            return stat;
        }
        CHECK_MSTATUS( mAttr.setWritable( false ) );
        CHECK_MSTATUS( mAttr.setReadable( true ) );
        CHECK_MSTATUS( mAttr.setArray( true ) );
        stat = addAttribute( oa_edgeTransforms );
        if ( !stat ) { stat.perror( "addAttribute oa_edgeTransforms" ); return stat; }
    }*/
    
    addNumericAttribute( oa_edgeTransforms, "outEdgeTransforms", "oet", MFnNumericData::kBoolean, false, false );
    stat = attributeAffects( ia_edgeTransforms, oa_edgeTransforms );
	if ( !stat ) { stat.perror( "attributeAffects ia_edgeTransforms->oa_edgeTransforms" ); return stat; }
    stat = attributeAffects( ia_simStepTaken, oa_edgeTransforms );
	if ( !stat ) { stat.perror( "attributeAffects ia_simStepTaken->oa_edgeTransforms" ); return stat; }
 
    addNumericAttribute( ia_userDefinedColors, "userDefinedColors", "udc", MFnNumericData::k3Double, 0, true, true );
    // This affects nothing as it is only cared about during drawing when we read it directly.

	return MS::kSuccess;
}
