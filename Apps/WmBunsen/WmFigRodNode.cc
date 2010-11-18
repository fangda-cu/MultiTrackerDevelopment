#include "WmFigRodNode.hh"
#include "WmFigConnectionNode.hh"

#include <maya/MFnMatrixAttribute.h>
#include <maya/MPlugArray.h>

using namespace BASim;


// Required by Maya to identify the node
/* static */ MTypeId WmFigRodNode::typeID ( 0x001135, 0x19 );
/* static */ MString WmFigRodNode::typeName( "wmFigRodNode" );

// Input attributes
/* static */ MObject WmFigRodNode::ia_time;
/* static */ MObject WmFigRodNode::ia_startTime;
/* static */ MObject WmFigRodNode::ia_nurbsCurves;
/* static */ MObject WmFigRodNode::oa_nurbsCurves;
/* static */ MObject WmFigRodNode::ia_barberShopVertices;
/* static */ MObject WmFigRodNode::ia_percentageOfBarberShopStrands;
/* static */ MObject WmFigRodNode::ia_simStepTaken;
/* static */ MObject WmFigRodNode::ia_strandRootFrames;

// Drawing
/* static */ MObject WmFigRodNode::ia_userDefinedColors;
/* static */ MObject WmFigRodNode::ia_draw3DRod;
/* static */ MObject WmFigRodNode::ia_drawScale;
/* static */ MObject WmFigRodNode::ca_drawDataChanged;

// Disk caching
/* static */ MObject WmFigRodNode::ia_cachePath;
/* static */ MObject WmFigRodNode::ia_cacheFrame;
/* static */ MObject WmFigRodNode::ia_readFromCache;
/* static */ MObject WmFigRodNode::ca_cachingHasBeenToggled;

// Output attributes
/* static */ MObject WmFigRodNode::oa_rodsChanged;
/* static */ MObject WmFigRodNode::oa_simulatedVertices;
/* static */ MObject WmFigRodNode::oa_nonSimulatedVertices;
/* static */ MObject WmFigRodNode::oa_verticesInEachRod;
/* static */ MObject WmFigRodNode::oa_materialFrames;
/* static */ MObject WmFigRodNode::oa_undeformedMaterialFrames;
/* static */ MObject WmFigRodNode::oa_numberOfRods;

// Cache attributes
/* static */ MObject WmFigRodNode::ca_simulationSync;
/* static */ MObject WmFigRodNode::ca_syncAttrs;

// User adjustable rod Options
/* static */ MObject WmFigRodNode::ia_solverType;
/* static */ MObject WmFigRodNode::ia_gravity;
/* static */ MObject WmFigRodNode::ia_cvsPerRod;
/* static */ MObject WmFigRodNode::ia_youngsModulus;
/* static */ MObject WmFigRodNode::ia_shearModulus;
/* static */ MObject WmFigRodNode::ia_viscosity;
/* static */ MObject WmFigRodNode::ia_density;
/* static */ MObject WmFigRodNode::ia_minorRadius;
/* static */ MObject WmFigRodNode::ia_majorRadius;
/* static */ MObject WmFigRodNode::ia_vertexSpacing;
/* static */ MObject WmFigRodNode::ia_minimumRodLength;
/* static */ MObject WmFigRodNode::ia_hairSpray;
/* static */ MObject WmFigRodNode::ia_hairSprayScaleFactor;
/* static */ MObject WmFigRodNode::ia_massDamping;
/* static */ MObject WmFigRodNode::ia_drawMaterialFrames;
/* static */ MObject WmFigRodNode::ia_lockFirstEdgeToInput;
/* static */ MObject WmFigRodNode::ia_simulationSet;

// For being controlled by external objects and for controlling external objects
/* static */ MObject WmFigRodNode::ia_edgeTransforms;
/* static */ MObject WmFigRodNode::oa_edgeTransforms;

WmFigRodNode::WmFigRodNode() : m_massDamping( 10 ), m_initialised( false ),
    /*mx_rodData( NULL ),*/ mx_world( NULL ), m_numberOfInputCurves( 0 ), 
    m_percentageOfBarberShopStrands( 100 ), m_verticesPerRod( -1 ), m_cacheFilename( "" ),
    m_pRodInput( NULL ), m_vertexSpacing( 0.0 ), m_minimumRodLength( 2.0 ),
    m_readFromCache( false ), m_writeToCache( false ), m_cachePath ( "" ),
    m_solverType( RodTimeStepper::IMPL_EULER ), m_gravity( 0, -980, 0 ), m_rodDisplayList(0),
    m_rodGeoChanged(true)
{
    m_rodOptions.YoungsModulus = 1000.0; /* megapascal */
    m_rodOptions.ShearModulus = 340.0;   /* megapascal */
    m_rodOptions.viscosity = 10.0;       /* poise */
    m_rodOptions.density = 1.3;          /* grams per cubic centimeter */
    m_rodOptions.radiusA = 0.05;         /* millimeter */
    m_rodOptions.radiusB = 0.05;         /* millimeter */
    m_rodOptions.refFrame = BASim::ElasticRod::TimeParallel;
    m_strandRootFrames.clear();
    m_controlledEdgeTransforms.clear();
    m_simulationSet.clear();
}


WmFigRodNode::~WmFigRodNode()
{
    if(m_rodDisplayList)
    {
        glDeleteLists(m_rodDisplayList, 1);
    }
}

MStatus WmFigRodNode::compute( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    MStatus stat;

    cerr << "WmFigRodNode::compute() with i_plug = " << i_plug.name() << endl;

    if ( i_plug == oa_numberOfRods )
    {
        // oa_numberOfRods only depends on time. This is because it is displayed in the AE
        // and so even if the node is hidden then it will evaluate everytime the AE is refreshed.
        // If it depends on the actual rod change inputs then the rods will do a lot of work.
        // This is not ideal...

        i_dataBlock.inputValue( ia_time, &stat ).asTime().value();
        CHECK_MSTATUS( stat );

        MDataHandle numRodsH = i_dataBlock.outputValue( oa_numberOfRods, &stat);
        CHECK_MSTATUS( stat );
    
        numRodsH.set( (int)m_rodGroup.numberOfRods() );
        
        stat = i_dataBlock.setClean( i_plug );
        if ( !stat )
        {
            stat.perror("WmFigRodNode::compute setClean");
            return MS::kFailure;
        }
    }
    else if ( i_plug == oa_rodsChanged )
    {
        // One of the inputs to the node has changed that should cause the simulation to 
        // either take a step or change some parameter.
        compute_oa_rodsChanged( i_plug, i_dataBlock );
    }
    else if ( i_plug == ca_simulationSync )
    {
        // The simulation has moved forward in time, we need to cache the result or update 
        // our OpenGL data for the rods.
        compute_ca_simulationSync( i_plug, i_dataBlock );
    }
    else if ( i_plug == oa_simulatedVertices )
    {
        // The Barbershop dynamic wire deformer is asking for the simulated rod vertices so it can
        // deform the hair appropriately.
        compute_oa_simulatedVertices( i_plug, i_dataBlock );   
    }
    else if ( i_plug == oa_nonSimulatedVertices )
    {
        // The Barbershop dynamic wire deformer is asking for the unsimulated rod vertices so it can
        // initialise its deformations correctly.
        compute_oa_nonSimulatedVertices( i_plug, i_dataBlock );
    }
    else if ( i_plug ==  oa_verticesInEachRod )
    {
        // The Barbershop dynamic wire deformer is asking for the number of vertices in each rod so
        // it can reconstruct the rods from the simulated and unsimulated vertices arrays.
        compute_oa_verticesInEachRod( i_plug, i_dataBlock );
    }
    else if ( i_plug == oa_materialFrames )
    {
        // Barbershop is asking for the material frames for each rod so that it knows how the
        // rods are rotating each frame.
        compute_oa_materialFrames( i_plug, i_dataBlock );
    }
    else if ( i_plug == oa_undeformedMaterialFrames )
    {
        // The Barbershop dynamic wire deformer is asking for the material frames for each rod when not
        // simulated. This is just the input curve to the rod stuck on the moving surface each frame.
        // This is used to compare the undeformed and deformed to get an idea of how twisted the rod is
        // each frame.
        compute_oa_materialFrames( i_plug, i_dataBlock );
    }
    else if ( i_plug == oa_edgeTransforms )
    {
        // The connection node wants the edge transform data, we're not going to give it as 
        // it's complicated but by having the connection it can follow it come get it when it wants
        // it and if we ever need to do some work before giving it we can do it here.
        compute_oa_EdgeTransforms( i_plug, i_dataBlock );
    }
    else if ( i_plug == ca_drawDataChanged )
    {
        // Update the list of colours we're holding for rods that the user has drawn in a different
        // colour. 
        compute_ca_drawDataChanged( i_plug, i_dataBlock );
    }
    else
    {
        return MS::kUnknownParameter;
    }

    return MS::kSuccess;
}

void WmFigRodNode::readCacheRelatedInputs( MDataBlock& i_dataBlock )
{
    MStatus stat;

    m_writeToCache = i_dataBlock.inputValue( ia_cacheFrame, &stat ).asBool();
    CHECK_MSTATUS( stat );

    m_cachePath = i_dataBlock.inputValue( ia_cachePath, &stat ).asString();
    CHECK_MSTATUS( stat );

    m_readFromCache = i_dataBlock.inputValue( ia_readFromCache, &stat ).asDouble();
}

void WmFigRodNode::writeCacheIfNeeded( MDataBlock& i_dataBlock )
{
    if ( m_writeToCache && !m_rodGroup.simulationNeedsReset() )
    {
        if ( !m_readFromCache )
        {
            getCacheFilename( i_dataBlock );
            WmFigRodFileIO::writeRodDataToCacheFile( m_cacheFilename, m_rodGroup );
        }
        else
        {
            MGlobal::displayWarning( "Ignoring request to cache frame as reading from cache!\n" );
        }
    }
}

/** @detail Initialises the RodData vector to the correct size and the initial
    rod data. It is called when we are about to recreate all the rods so time has just been
    set back to 'startTime'.

    @param i_rodData A pointer to a vector of pointers to rod data. This????
*/
void WmFigRodNode::initialiseRodData( MDataBlock& i_dataBlock )
{
    MStatus stat;

    readCacheRelatedInputs( i_dataBlock );

    // Remove all rods as we're about to build new ones from one of our inputs
    m_rodGroup.removeAllRods();

    // Delete the old rod input class as the user may have switched inputs and we're
    // recreating it no matter what.
    if ( m_pRodInput != NULL )
        delete m_pRodInput;

    //
    // We need to create an instance of the correct input class. Apart from the constructor
    // these all have the same interface so they provide a uniform way to get input and turn
    // it into rod data.
    //
    // If multiple inputs are connected or reading from the cache is enabled then there is
    // an ordering to what takes precidence. Only one input is currently used at a time,
    // although it may not be too hard to handle multiple inputs.
    //
   
    MPlug barberShopPlug( thisMObject(), ia_barberShopVertices );
   
    if ( m_readFromCache)
    {
        getCacheFilename( i_dataBlock );
        m_pRodInput = new WmFigRodFileInput( m_cacheFilename, m_rodGroup, m_rodOptions );
    }
    else if ( barberShopPlug.isConnected() )
    {
        m_pRodInput = new WmFigRodBarbInput( ia_barberShopVertices, ia_strandRootFrames,
                                             m_percentageOfBarberShopStrands, m_verticesPerRod,
                                             m_lockFirstEdgeToInput, m_vertexSpacing,
                                             m_minimumRodLength,
                                             m_rodOptions, m_massDamping, m_gravity, m_rodGroup,
                                             m_solverType, m_simulationSet );
    }
    else // Assume we have nurbs connected
    {
        m_pRodInput = new WmFigRodNurbsInput( ia_nurbsCurves, m_lockFirstEdgeToInput, m_rodGroup,
                           m_vertexSpacing, m_minimumRodLength, m_rodOptions, m_massDamping,
                           m_gravity, m_solverType, m_simulationSet );
    }
    
    m_pRodInput->initialiseRodDataFromInput( i_dataBlock );

    // The sim has been re-initialised so it no longer needs reset.
    m_rodGroup.setSimulationNeedsReset( false );

    // Write out this frame as it doesn't need simulated, we just need to store what was just
    // initialised.
    writeCacheIfNeeded( i_dataBlock );

    // We need to make sure we have the spline attr data for the rods since compute may not have been called yet
    // FIXME: These are broken so I'm removing them for just now.
    //updateHairsprayScales( dataBlock );
}

void WmFigRodNode::updateKinematicEdgesFromInput()
{

    for ( int r = 0; r < (int)m_rodGroup.numberOfRods(); ++r )
    {
        EdgeTransformRodMap::iterator rodIt = m_controlledEdgeTransforms.find( r );
        if ( rodIt != m_controlledEdgeTransforms.end() )
        {
            for ( EdgeTransformMap::iterator edgeIt = m_controlledEdgeTransforms[ r ].begin();
                    edgeIt != m_controlledEdgeTransforms[ r ].end();
                    edgeIt++ )
            {
                if ( m_currentTime == m_startTime )
                {
                    m_rodGroup.resetKinematicEdge( r, edgeIt->first, edgeIt->second.materialFrame );
                }
                else
                {
                    m_rodGroup.updateKinematicEdge( r, edgeIt->first, edgeIt->second.materialFrame );
                }
                
                // Set the next positions of these vertices to be wherever the input controller
                // says they should be.
                BASim::Vec3d position =  edgeIt->second.position;
                double length = m_rodGroup.elasticRod( r )->getEdge( edgeIt->first ).norm();
                BASim::Vec3d edge = edgeIt->second.materialFrame.m1;
                edge.normalize();
                BASim::Vec3d start = position - ( edge * length / 2.0 );
                BASim::Vec3d end = position + ( edge * length / 2.0 );
                m_rodGroup.setNextVertexPosition( r, edgeIt->first, start );
                m_rodGroup.setNextVertexPosition( r, edgeIt->first + 1, end );
            }
        }
    }
}

void WmFigRodNode::updateOrInitialiseRodDataFromInputs( MDataBlock& i_dataBlock )
{
    MStatus stat;

    readCacheRelatedInputs( i_dataBlock );

    if ( m_currentTime == m_startTime )
    {
        initialiseRodData( i_dataBlock );     
    }
    else
    {
        if ( m_pRodInput == NULL )
        {
            return;
        }

        if ( !m_rodGroup.simulationNeedsReset() )
        {
            if ( !m_readFromCache )
            {
                m_pRodInput->updateRodDataFromInput( i_dataBlock );
            }
            else
            {
                WmFigRodFileIO::updateRodDataFromCacheFile( m_cacheFilename, m_rodGroup );
            }
        }
        else
        {
            MGlobal::displayWarning( "Please rewind simulation to reset\n" );
        }

        updateKinematicEdgesFromInput();
    }

    m_rodGroup.setRodParameters( m_rodOptions.radiusA, m_rodOptions.radiusB,
                                 m_rodOptions.YoungsModulus,
                                 m_rodOptions.ShearModulus,
                                 m_rodOptions.viscosity,
                                 m_rodOptions.density );
}

/** @detail Returns the material frame matrix for a specific rod's edge.
    @param i_rod The index of the rod you wish the edge of (0 based indices).
    @param i_edge The index of the edge you are interested in (0 based indices).
    @return An MMatrix containing the matrial frame and the position of the edge. 
            | m1x, m1y, m1z, 0 |
            | m2x, m2y, m2z, 0 |
            | m3x, m3y, m3z, 0 |
            |   x,   y,   z, 1 |
            The position is defined as the mid point of the edge:

    @note If data does not exist for this rod, for example it hasn't yet been created
          or one of the indices is out of range then the identity matrix will be
          returned.
  */
MMatrix WmFigRodNode::getRodEdgeMatrix( int i_rod, int i_edge )
{
    MMatrix identMatrix;
    identMatrix.setToIdentity();
    
    // Check if the input parameters index a valid rod and edge, if not
    // return the identity matrix.

    if ( i_rod > m_rodGroup.numberOfRods() )
        return identMatrix;
    
    if ( m_rodGroup.elasticRod( i_rod ) == NULL )
        return identMatrix;
    
    if ( i_edge >= m_rodGroup.numberOfEdgesInRod( i_rod ) )
        return identMatrix;

    // If we got here then the input data indexes a valid rod.

    ElasticRod* rod = m_rodGroup.elasticRod( i_rod );

    // The position of the edge is defined as the mid point of the ege.   
    BASim::Vec3d edgePos = ( rod->getVertex( i_edge ) + 
                      rod->getVertex( i_edge + 1 ) ) / 2.0;
    

    // The material frame is taken directly from the rod and packaged up into a Maya MMatrix format.

    MMatrix edgeMatrix;
    edgeMatrix( 3, 0 ) = edgePos[ 0 ]; edgeMatrix( 3, 1 ) = edgePos[ 1 ]; edgeMatrix( 3, 2 ) = edgePos[ 2 ];
    
    BASim::Vec3d edge = rod->getEdge( i_edge );
    edgeMatrix( 0, 0 ) = edge[ 0 ]; edgeMatrix( 0, 1 ) = edge[ 1 ]; edgeMatrix( 0, 2 ) = edge[ 2 ];
    
    BASim::Vec3d material1 = rod->getMaterial1( i_edge );
    edgeMatrix( 1, 0 ) = material1[ 0 ]; edgeMatrix( 1, 1 ) = material1[ 1 ]; edgeMatrix( 1, 2 ) = material1[ 2 ];
    
    BASim::Vec3d material2 = rod->getMaterial2( i_edge );
    edgeMatrix( 2, 0 ) = material2[ 0 ]; edgeMatrix( 2, 1 ) = material2[ 1 ]; edgeMatrix( 2, 2 ) = material2[ 2 ];
    
    return edgeMatrix;
}

MString WmFigRodNode::getCacheFilename( MDataBlock& i_dataBlock )
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

void WmFigRodNode::updateHairsprayScales( MDataBlock& i_dataBlock )
{
    /*if ( mx_rodData == NULL )
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
        }
    }*/
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Computing oa_rodsChanged
// This means that one of the inputs to the node has changed that should cause the simulation to 
// either take a step or change some parameter.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

void WmFigRodNode::compute_oa_rodsChanged( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    MStatus stat;

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

    const double3 &gravity = i_dataBlock.inputValue( ia_gravity, &stat ).asDouble3();
    CHECK_MSTATUS( stat );
    m_gravity = BASim::Vec3d( gravity[0], gravity[1], gravity[2] );

    m_lockFirstEdgeToInput = i_dataBlock.inputValue( ia_lockFirstEdgeToInput, &stat ).asBool();
    CHECK_MSTATUS( stat );

    // This will get and store the cache filename for the current frame
    getCacheFilename( i_dataBlock );

    m_verticesPerRod = i_dataBlock.inputValue( ia_cvsPerRod, &stat ).asInt();
    CHECK_MSTATUS( stat );

    m_vertexSpacing = i_dataBlock.inputValue( ia_vertexSpacing, &stat ).asDouble();
    CHECK_MSTATUS( stat );

    m_minimumRodLength = i_dataBlock.inputValue( ia_minimumRodLength, &stat ).asDouble();
    CHECK_MSTATUS( stat );

    m_percentageOfBarberShopStrands = i_dataBlock.inputValue( ia_percentageOfBarberShopStrands, &stat ).asDouble();
    CHECK_MSTATUS( stat );
    
    m_solverType = ( RodTimeStepper::Method)  ( i_dataBlock.inputValue( ia_solverType, &stat ).asInt() );
    CHECK_MSTATUS( stat );
    
    MString simulationSet = i_dataBlock.inputValue( ia_simulationSet, &stat ).asString();
    CHECK_MSTATUS( stat );

    updateSimulationSet( simulationSet );

    updateControlledEdgeArrayFromInputs( i_dataBlock );

    MDataHandle verticesH = i_dataBlock.inputValue( ia_barberShopVertices, &stat );
    CHECK_MSTATUS( stat );    

    //////////////////////////////////////////////////////////////////////////////////////////
    
    // FIXME:
    // We should really just disable all rods at when this is toggled and make the user
    // go back to start time before they can see the cached rods?

    bool readFromCache = i_dataBlock.inputValue( ia_readFromCache, &stat ).asDouble();

    if ( m_readFromCache != readFromCache )
    {
        // We changed reading from cache status not at startTime so stop doing anything until
        // the user resets the sim to be sure all is safe.
        // Make sure that every rod is disabled before we read from the cache file
        
        cerr << "read status changed disabling all rods\n";
        m_rodGroup.setSimulationNeedsReset( true );
        
        m_readFromCache = readFromCache;
    }

    m_rodGroup.setIsReadingFromCache( m_readFromCache );
    updateOrInitialiseRodDataFromInputs( i_dataBlock );
        
    stat = i_dataBlock.setClean( i_plug );
    if ( !stat )
    {
        stat.perror("WmFigRodNode::compute setClean");
        return;
    }
}

void WmFigRodNode::updateSimulationSet( MString i_simulationSetString )
{
    m_simulationSet.clear();

    MStringArray subStrings;
    
    i_simulationSetString.split( ',', subStrings );

    for ( int c=0; c< subStrings.length(); ++c )
    {
        m_simulationSet.insert( subStrings[ c ].asInt() );
    }
}

void WmFigRodNode::updateControlledEdgeArrayFromInputs( MDataBlock& i_dataBlock )
{
    MStatus stat;

    MDataHandle inputCurveH;
    MObject inputCurveObj;    
    MArrayDataHandle inEdgeArrayH = i_dataBlock.inputArrayValue( ia_edgeTransforms, &stat );
    CHECK_MSTATUS(stat);
    
    if ( m_rodGroup.numberOfRods() == 0 )
    {
        return;
    }

    size_t numNodesConnected = inEdgeArrayH.elementCount();

    for ( unsigned int e = 0; e < numNodesConnected; e++ )
    {
        // Ask for the input on them all to make sure the connection nodes evaluate
        inEdgeArrayH.jumpToElement( e );
        inEdgeArrayH.inputValue( &stat );
        CHECK_MSTATUS( stat );
    }
    
    // First get rid of all the constraints as some of them may have been removed.
    // We could do this in connectionBroken but there is always weirdness when loading
    // and tracking stuff in connection made/broken.
    
    // FIXME: This is pretty darn innefficient, deleting then recreating
    // Although it does let us create and delete them at any frame in the simulation
    
    for ( EdgeTransformRodMap::iterator rodIt = m_controlledEdgeTransforms.begin();
        rodIt != m_controlledEdgeTransforms.end(); rodIt++ )
    {
        for ( EdgeTransformMap::iterator edgeIt = rodIt->second.begin();
                edgeIt != rodIt->second.end();
                edgeIt++ )
        {
            m_rodGroup.removeKinematicEdge( rodIt->first, edgeIt->first );
        }
        m_controlledEdgeTransforms[ rodIt->first ].clear();
    }
    m_controlledEdgeTransforms.clear();
    
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
                    unsigned int rodIndex, edgeIndex;
                    EdgeTransform edgeTransform;
                    connectionNode->getControlledRodInfo( rodIndex, edgeIndex, edgeTransform );
                    
                    m_controlledEdgeTransforms[ rodIndex ][ edgeIndex ] = edgeTransform;
                    
                    m_rodGroup.addKinematicEdge( rodIndex, edgeIndex, &edgeTransform.materialFrame );
                }
            }
        }
    }
    inEdgeArrayH.setClean();
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Compute ca_simulationSync and oa_numberOfRods
//
// These are recomputed whenever the simulation has moved forward in time, we need to cache
// the result or update our OpenGL data for the rods.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

void WmFigRodNode::compute_ca_simulationSync( const MPlug& i_plug, 
    MDataBlock i_dataBlock )
{
    MStatus stat;

    i_dataBlock.inputValue( ia_simStepTaken, &stat ).asBool();
    CHECK_MSTATUS( stat );

    readCacheRelatedInputs( i_dataBlock );

    writeCacheIfNeeded( i_dataBlock );
    
    MDataHandle numRodsH = i_dataBlock.outputValue( oa_numberOfRods, &stat);
    CHECK_MSTATUS( stat );

    numRodsH.set( (int)m_rodGroup.numberOfRods() );
    
    m_rodGeoChanged = true;

    stat = i_dataBlock.setClean( i_plug );
    if ( !stat )
    {
        stat.perror("WmFigRodNode::compute setClean");
        return;
    }
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Compute oa_simulatedVertices
//
// The Barbershop dynamic wire deformer is asking for the simulated rod vertices so it can
// deform the hair appropriately.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

void WmFigRodNode::compute_oa_simulatedVertices( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    MStatus stat;

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

    int numRods = m_rodGroup.numberOfRods();
    unsigned int idx = 0;
    for ( int r = 0; r < numRods; r++ )
    {
        unsigned int verticesInRod = m_rodGroup.numberOfVerticesInRod( r );
        simulatedVerticesArray.setLength( (unsigned int) ( simulatedVerticesArray.length() + verticesInRod ) );

        for ( unsigned int v = 0; v < verticesInRod; v++ )
        {
            BASim::Vec3d pos = m_rodGroup.elasticRod( r )->getVertex( v );
            simulatedVerticesArray[ idx ] = MVector( pos[0], pos[1], pos[2] );
            idx++;
        }
    }


    simulatedVerticesH.setClean();
    i_dataBlock.setClean( i_plug );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Compute oa_simulatedVertices
//
// The Barbershop dynamic wire deformer is asking for the simulated rod vertices so it can
// deform the hair appropriately.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

void WmFigRodNode::compute_oa_nonSimulatedVertices( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    MStatus stat;

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

    int numRods = m_rodGroup.numberOfRods();

    unsigned int idx = 0;

    for ( int r = 0; r < numRods; r++ )
    {
        unsigned int verticesInRod = m_rodGroup.numberOfVerticesInRod( r );
        nonSimulatedVerticesArray.setLength( (unsigned int) ( nonSimulatedVerticesArray.length() + verticesInRod ) );

        for ( unsigned int v = 0; v < verticesInRod; v++ )
        {
            BASim::Vec3d pos = m_rodGroup.nextVertexPosition( r, v );
            nonSimulatedVerticesArray[ idx ] = MVector( pos[0], pos[1], pos[2] );
            idx++;
        }
    }        

    nonSimulatedVerticesH.setClean();
    i_dataBlock.setClean( i_plug );
}


////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Compute oa_verticesInEachRod
//
// The Barbershop dynamic wire deformer is asking for the number of vertices in each rod so it can
// reconstruct the rods from the simulated and unsimulated vertices arrays.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

void WmFigRodNode::compute_oa_verticesInEachRod( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    MStatus stat;

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

    MIntArray verticesPerRodArray = verticesPerRodArrayData.array( &stat );
    int numRods = m_rodGroup.numberOfRods();
    verticesPerRodArray.setLength( (unsigned int)numRods );
    unsigned int idx = 0;

    for ( int r = 0; r < numRods; r++ )
    {
        unsigned int verticesInRod = m_rodGroup.numberOfVerticesInRod( r );
        verticesPerRodArray[ (int)r ] = verticesInRod;
    }

    verticesPerRodH.setClean();
    i_dataBlock.setClean( i_plug );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Compute oa_materialFrames
//
// The Barbershop dynamic wire deformer is asking for the material frames so that it can 
// see how much each rod rotates every frame.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

void WmFigRodNode::compute_oa_materialFrames( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    MStatus stat;

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

    int numRods = m_rodGroup.numberOfRods();
    unsigned int idx = 0;
    for ( int r = 0; r < numRods; r++ )
    {
        ElasticRod* rod = m_rodGroup.elasticRod( r );
        if ( rod != NULL )
        {
            unsigned int edgesInRod = rod->ne();
            materialFramesArray.setLength( (unsigned int) ( materialFramesArray.length() + edgesInRod*3 ) );

            for ( unsigned int e = 0; e < edgesInRod; e++ )
            {
                BASim::Vec3d m1 =  rod->getMaterial1( e );
                BASim::Vec3d m2 =  rod->getMaterial2( e );
                BASim::Vec3d m3 =  rod->getEdge( e );
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

    materialFramesH.setClean();
    i_dataBlock.setClean( i_plug );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Compute oa_undeformedMaterialFrames
//
// The Barbershop dynamic wire deformer is asking for the material frames for each rod when not
// simulated. This is just the input curve to the rod stuck on the moving surface each frame.
// This is used to compare the undeformed and deformed to get an idea of how twisted the rod is
// each frame.
//
////////////////////////////////////////////////////////////////////////////////////////////////////

void WmFigRodNode::compute_oa_undeformedMaterialFrames( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    MStatus stat;

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

        int numRods = m_rodGroup.numberOfRods();
        unsigned int idx = 0;
        for ( int r = 0; r < numRods; r++ )
        {
            ElasticRod* rod = m_rodGroup.elasticRod( r );
            if ( rod != NULL )
            {
                undeformedMaterialFramesArray.setLength( (unsigned int) ( undeformedMaterialFramesArray.length() + rod->ne()*3 ) );

                for ( int e=0; e<rod->ne(); e++ )
                {
                    // currently we only store the undeformed frame for the first vertex
                    BASim::Vec3d m1 =  rod->getMaterial1( e );
                    BASim::Vec3d m2 =  rod->getMaterial2( e );
                    BASim::Vec3d m3 =  rod->getEdge( e );
                    m3.normalize();

                    /*(*mx_rodData)[ r ]->undeformedMaterialFrame[ e ].m1 = m1;
                    (*mx_rodData)[ r ]->undeformedMaterialFrame[ e ].m2 = m2;
                    (*mx_rodData)[ r ]->undeformedMaterialFrame[ e ].m3 = m3;*/

                    m_rodGroup.setUndeformedMaterialFrame( r, e, m1, m2, m3 );

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
    else
    {
        // We're not at startTime so we must be simulating and need to work out what the
        // unsimulated material frames would be at this point based on the input curve.

        int numRods = m_rodGroup.numberOfRods();

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
            for ( int r = 0; r < numRods; r++ )
            {
                ElasticRod* rod = m_rodGroup.elasticRod( r );
                if ( rod != NULL )
                {
                /*  m_strandRootFrames is what we create the first matrix from then we create
                    the second matrix from the strand root frames we were given this frame.
                    The rod frames should not be in either matrix!!!!*/

                    BASim::Vec3d im1 = m_strandRootFrames[r].m1;
                    BASim::Vec3d im2 = m_strandRootFrames[r].m2;
                    BASim::Vec3d im3 = m_strandRootFrames[r].m3;

                    BASim::Vec3d cm1 = strandRootFrames[r].m1;
                    BASim::Vec3d cm2 = strandRootFrames[r].m2;
                    BASim::Vec3d cm3 = strandRootFrames[r].m3;

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

                    for ( int e=0; e<rod->ne(); e++ )
                    {

                        MaterialFrame materialFrame = m_rodGroup.undeformedMaterialFrame( r, e );
                        // currently we only store the undeformed frame for the first vertex
                        BASim::Vec3d m1 =  materialFrame.m1;
                        MVector mayaM1( m1[0], m1[1], m1[2] );
                        BASim::Vec3d m2 =  materialFrame.m2;
                        MVector mayaM2( m2[0], m2[1], m2[2] );
                        BASim::Vec3d m3 = materialFrame.m3;
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
    undeformedMaterialFramesH.setClean();
    i_dataBlock.setClean( i_plug );
}

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Compute oa_edgeTransforms
//
// The connection node wants the edge transform data, we're not going to give it as it's complicated
// but by having the connection it can follow it come get it when it wants it and if we ever need
// to do some work before giving it we can do it here.
// 
//
////////////////////////////////////////////////////////////////////////////////////////////////////

void WmFigRodNode::compute_oa_EdgeTransforms( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    MStatus stat;

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

////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Compute ca_drawDataChanged
//
// Update the list of colours we're holding for rods that the user has drawn in a different
// colour. 
//
////////////////////////////////////////////////////////////////////////////////////////////////////

void WmFigRodNode::compute_ca_drawDataChanged( const MPlug& i_plug, MDataBlock& i_dataBlock )
{
    MStatus stat;

    bool draw3DRod = i_dataBlock.inputValue( ia_draw3DRod, &stat ).asBool();
    CHECK_MSTATUS( stat );

    double drawScale = i_dataBlock.inputValue( ia_drawScale, &stat ).asDouble();
    CHECK_MSTATUS( stat );

    RodRenderer::DrawMode drawMode = RodRenderer::SIMPLE;

    if ( draw3DRod )
    {
        drawMode = RodRenderer::SMOOTH;
    }

    int numRods = m_rodGroup.numberOfRods();
    for ( int r = 0; r < numRods; r++ )
    {
        m_rodGroup.setDrawScale( drawScale );
        m_rodGroup.setDrawMode( drawMode );
    }
    
    MDataHandle inputColourHandle;
    MObject inputColourObj;
    MArrayDataHandle inArrayH = i_dataBlock.inputArrayValue( ia_userDefinedColors, &stat );
    CHECK_MSTATUS(stat);

    int numColoursSet = inArrayH.elementCount();

    // Get rid of all previously stored colours as the user may have removed some
    m_rodColourMap.clear();
    
    for ( unsigned int i = 0; i < numColoursSet; i++ )
    {
        inArrayH.jumpToArrayElement( i );
        
        // We're iterating over the actual indices 'i', as the array
        // is sparse we need to find out the logical index as seen
        // by the user in MEL.
        unsigned int elementIndex = inArrayH.elementIndex();
        
        const double3& colour = inArrayH.inputValue( &stat ).asDouble3();
        CHECK_MSTATUS( stat );

        // FIXME:
        // We use -1 to indicate not to colour this rod any more. Should really
        // remove the element from the array, but the APIs unclear.
        if ( colour[ 0 ] != -1 )
            m_rodColourMap[ elementIndex ] = BASim::Vec3d( colour[0], colour[1], colour[2] );
    }
    inArrayH.setClean();
    i_dataBlock.setClean( i_plug );

    m_rodGeoChanged = true;
}

void WmFigRodNode::getStrandRootFrames( MDataBlock& i_dataBlock, vector<MaterialFrame>& o_strandRootFrames )
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
    int idx = 0;
    for ( int rIdx=0; rIdx<o_strandRootFrames.size(); rIdx++ )
    {
        v = strandRootFrameVec[idx++];
        o_strandRootFrames[rIdx].m1 = BASim::Vec3d( v[0], v[1], v[2] );
        v = strandRootFrameVec[idx++];
        o_strandRootFrames[rIdx].m2 = BASim::Vec3d( v[0], v[1], v[2] );
        v = strandRootFrameVec[idx++];
        o_strandRootFrames[rIdx].m3 = BASim::Vec3d( v[0], v[1], v[2] );
    }
}


void WmFigRodNode::draw( M3dView& i_view, const MDagPath& i_path,
                            M3dView::DisplayStyle i_style,
                            M3dView::DisplayStatus i_status )
{
    
    MStatus stat;
    MObject thisNode = thisMObject();

    bool b;
    MPlug cachingToggledPlug( thisNode, ca_cachingHasBeenToggled );
    cachingToggledPlug.getValue( b );

    MPlug syncPlug( thisNode, ca_simulationSync );
    double d;
    stat = syncPlug.getValue( d );
    if ( !stat )
    {
        stat.perror( "WmFigRodNode::draw getting ca_simulationSync" );
        return;
    }
    
    // Pull on the draw data changed as this is the only function that ever will
    // as no one else cares.
    MPlug drawPlug( thisNode, ca_drawDataChanged );
    stat = drawPlug.getValue( d );
    if ( !stat )
    {
        stat.perror( "WmFigRodNode::draw getting ca_drawDataChanged" );
        return;
    }

    // FIXME: Debugging code, this also happens in compute_ca_DrawDataChanged
    // but it appears to happen after the rods are created so the render
    // is not getting update correctly.
    {
        bool draw3DRod;
        MPlug draw3DPlug( thisNode, ia_draw3DRod );
        stat = draw3DPlug.getValue( draw3DRod );
        CHECK_MSTATUS( stat );

        double drawScale;
        MPlug drawScalePlug( thisNode, ia_drawScale );
        stat = drawScalePlug.getValue( drawScale );
        CHECK_MSTATUS( stat );

        RodRenderer::DrawMode drawMode = RodRenderer::SIMPLE;

        if ( draw3DRod )
        {
            drawMode = RodRenderer::SMOOTH;
        }

        size_t numRods = m_rodGroup.numberOfRods();
        m_rodGroup.setDrawScale( drawScale );
        m_rodGroup.setDrawMode( drawMode );        
    }


    i_view.beginGL();
    //glPushAttrib( GL_CURRENT_BIT | GL_POINT_BIT | GL_LINE_BIT );
	glPushAttrib( GL_ALL_ATTRIB_BITS );

/*    if ( mx_rodData != NULL )
    {
        GLfloat currentColour[4];
        
        for ( size_t r=0; r<mx_rodData->size(); r++ )
        {
            bool colourOverride = false;
            if ( m_rodColourMap.find( r ) != m_rodColourMap.end() )
            {
                colourOverride = true;
                glGetFloatv( GL_CURRENT_COLOR, currentColour );
                
                BASim::Vec3d colour = m_rodColourMap[ r ];
                glColor3ub( colour[0], colour[1], colour[2] );
            }
            
            (*mx_rodData)[ r ]->rodRenderer->render();
            
            if ( colourOverride )
                glColor4fv( currentColour );
        }
        
    }*/

    if(m_rodDisplayList == 0)
    {
        m_rodDisplayList = glGenLists(1);

        glNewList(m_rodDisplayList, GL_COMPILE);
            m_rodGroup.render();
        glEndList();

        m_rodGeoChanged = false;
    }
    else if(m_rodGeoChanged)
    {
        glDeleteLists(m_rodDisplayList, 1);

        glNewList(m_rodDisplayList, GL_COMPILE);
            m_rodGroup.render();
        glEndList();

        m_rodGeoChanged = false;
    }

    glCallList(m_rodDisplayList);
    
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
                    BASim::Vec3d m1 =  rod->getMaterial1( e );
                    BASim::Vec3d m2 =  rod->getMaterial2( e );
                    BASim::Vec3d m3 =  rod->getEdge( e );
                    m3.normalize();

                    BASim::Vec3d p = rod->getVertex( e );
                    BASim::Vec3d p1 = rod->getVertex( e + 1 );
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

                    BASim::Vec3d m1 =  m_strandRootFrames[ r ].m1;
                    BASim::Vec3d m2 =  m_strandRootFrames[ r ].m2;
                    BASim::Vec3d m3 =  m_strandRootFrames[ r ].m3;

                    BASim::Vec3d p = rod->getVertex( 0 );

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
                    BASim::Vec3d m1 = (*mx_rodData)[ r ]->undeformedMaterialFrame[ e ].m1;
                    BASim::Vec3d m2 = (*mx_rodData)[ r ]->undeformedMaterialFrame[ e ].m2;
                    BASim::Vec3d m3 = (*mx_rodData)[ r ]->undeformedMaterialFrame[ e ].m3;

                    BASim::Vec3d p = rod->getVertex( e );
                    BASim::Vec3d p1 = rod->getVertex( e + 1 );
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
	//view.setDrawColor ( WmFigRodNode );

	glPopAttrib();
	i_view.endGL();
}

MStatus WmFigRodNode::connectionMade( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc )
{
    MStatus stat;
    MStatus retVal(MS::kUnknownParameter );

    if( plug == ia_nurbsCurves )
    {
        m_numberOfInputCurves++;
        
//        cerr << "Increasing number of curves (" << m_numberOfInputCurves << ")\n";
    }

    return retVal;
}

MStatus WmFigRodNode::connectionBroken( const  MPlug & plug, const  MPlug & otherPlug, bool asSrc )
{
    MStatus stat;
    MStatus retVal( MS::kUnknownParameter );

    if( plug == ia_nurbsCurves )
    {
        m_numberOfInputCurves--;
  //      cerr << "Decreasing number of curves (" << m_numberOfInputCurves << ")\n";
    }

    return retVal;
}


bool WmFigRodNode::isBounded() const
{
	return false;
}

void* WmFigRodNode::creator()
{
	return new WmFigRodNode();
}

/*static */ MStatus WmFigRodNode::addNumericAttribute( MObject& i_attribute, MString i_longName, 
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

/* static */ MStatus WmFigRodNode::initialize()
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

    addNumericAttribute( ia_vertexSpacing, "vertexSpacing", "vsp", MFnNumericData::kDouble, 1.0, true );
    stat = attributeAffects( ia_vertexSpacing, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_vertexSpacing->oa_rodsChanged" ); return stat; }

    addNumericAttribute( ia_minimumRodLength, "minimumRodLength", "mrl", MFnNumericData::kDouble, 2.0, true );
    stat = attributeAffects( ia_minimumRodLength, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_minimumRodLength->oa_rodsChanged" ); return stat; }

    addNumericAttribute( ia_cvsPerRod, "cvsPerRod", "cvr", MFnNumericData::kInt, -1, true );
    stat = attributeAffects( ia_cvsPerRod, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_cvsPerRod->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ia_cacheFrame, "cacheFrame", "caf", MFnNumericData::kBoolean, false, true );
    stat = attributeAffects( ia_cacheFrame, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_cacheFrame->oa_rodsChanged" ); return stat; }
    stat = attributeAffects( ia_cacheFrame, ca_simulationSync );
	if ( !stat ) { stat.perror( "attributeAffects ia_cacheFrame->ca_syncAttrs" ); return stat; }

    addNumericAttribute( ca_cachingHasBeenToggled, "cachingHasBeenToggled", "cht", MFnNumericData::kBoolean, false, false );

    addNumericAttribute( ia_readFromCache, "readFromCache", "rfc", MFnNumericData::kBoolean, false, true );
    stat = attributeAffects( ia_readFromCache, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_readFromCache->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_readFromCache, ca_simulationSync );
    if ( !stat ) { stat.perror( "attributeAffects ia_cacheFrame->ca_syncAttrs" ); return stat; }
    stat = attributeAffects( ia_readFromCache, ca_cachingHasBeenToggled );
    if ( !stat ) { stat.perror( "attributeAffects ia_cacheFrame->ca_cachingHasBeenToggled" ); return stat; }
    
    addNumericAttribute( ia_hairSprayScaleFactor, "hairSprayScaleFactor", "hsf", MFnNumericData::kDouble, 1.0, true );
    stat = attributeAffects( ia_hairSprayScaleFactor, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_hairSprayScaleFactor->oa_rodsChanged" ); return stat; }

    addNumericAttribute( ia_massDamping, "massDamping", "mda", MFnNumericData::kDouble, 10.0, true );
    stat = attributeAffects( ia_massDamping, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_massDamping->oa_rodsChanged" ); return stat; }

    addNumericAttribute( ia_percentageOfBarberShopStrands, "percentageOfFozzieStrands", "pfs", MFnNumericData::kDouble, 1.0, true );
    stat = attributeAffects( ia_percentageOfBarberShopStrands, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_percentageOfBarberShopStrands->oa_rodsChanged" ); return stat; }

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
        MFnTypedAttribute tAttr;
        MFnStringData fnStringData;
        MObject defaultString = fnStringData.create("");
        ia_simulationSet = tAttr.create( "simulationSet", "sse", MFnData::kString, defaultString, &stat );
        if ( !stat )
        {
            stat.perror( "create ia_simulationSet attribute" );
            return stat;
        }
        tAttr.setWritable( true );
        tAttr.setReadable( false );
        tAttr.setConnectable( true );
        stat = addAttribute( ia_simulationSet );
        if (!stat) { stat.perror( "addAttribute ia_simulationSet" ); return stat; }
    }
    stat = attributeAffects( ia_simulationSet, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_simulationSet->oa_rodsChanged" ); return stat; }

    
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
        ia_barberShopVertices = tAttr.create( "fozzieVertices", "fov",
                                          MFnData::kVectorArray, &stat );
        CHECK_MSTATUS( stat );
        CHECK_MSTATUS( tAttr.setReadable( false ) );
        CHECK_MSTATUS( tAttr.setWritable( true ) );
        stat = addAttribute( ia_barberShopVertices );
        CHECK_MSTATUS( stat );
    }
    stat = attributeAffects( ia_barberShopVertices, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_barberShopVertices->oa_rodsChanged" ); return stat; }


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
        oa_simulatedVertices = tAttr.create( "simulatedCurveCVs", "sve",
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
        oa_verticesInEachRod = tAttr.create( "curveCVCounts", "ver",
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
 
    addNumericAttribute( oa_numberOfRods, "curveCount", "nor", MFnNumericData::kInt, 0, false );
    stat = attributeAffects( ia_time, oa_numberOfRods );
	if ( !stat ) { stat.perror( "attributeAffects ia_time->oa_numberOfRods" ); return stat; }
 
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

    addNumericAttribute( ia_gravity, "gravity", "gra", MFnNumericData::k3Double, 0, true, false );
    stat = attributeAffects( ia_gravity, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_gravity->oa_rodsChanged" ); return stat; }
    
    addNumericAttribute( ia_edgeTransforms, "edgeTransforms", "iet", MFnNumericData::kBoolean, false, true, true );
    stat = attributeAffects( ia_edgeTransforms, oa_rodsChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_edgeTransforms->oa_rodsChanged" ); return stat; }
    
    addNumericAttribute( oa_edgeTransforms, "outEdgeTransforms", "oet", MFnNumericData::kBoolean, false, false );
    stat = attributeAffects( ia_edgeTransforms, oa_edgeTransforms );
	if ( !stat ) { stat.perror( "attributeAffects ia_edgeTransforms->oa_edgeTransforms" ); return stat; }
    stat = attributeAffects( ia_simStepTaken, oa_edgeTransforms );
	if ( !stat ) { stat.perror( "attributeAffects ia_simStepTaken->oa_edgeTransforms" ); return stat; }
 
    // This array holds the colours of any rods the user has decided to colour differently.
    addNumericAttribute( ia_userDefinedColors, "userDefinedColors", "udc", MFnNumericData::k3Double, 0, true, true );
    // This output attribute is purely so that we can get notified when the colour changes on the input
    // It will also be affected by the rods changing when we add the VBO code in here.
    addNumericAttribute( ca_drawDataChanged, "drawDataChanged", "ddc", MFnNumericData::kBoolean, true, false, false );
    stat = attributeAffects( ia_userDefinedColors, ca_drawDataChanged );
	if ( !stat ) { stat.perror( "attributeAffects ia_userDefinedColors->ca_drawDataChanged" ); return stat; }
    
    addNumericAttribute( ia_draw3DRod, "draw3DRod", "d3d", MFnNumericData::kBoolean, false, true );
    stat = attributeAffects( ia_draw3DRod, ca_drawDataChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_draw3DRod->ca_drawDataChanged" ); return stat; }

    addNumericAttribute( ia_drawScale, "drawScale", "dsc", MFnNumericData::kDouble, 10.0, true );
    stat = attributeAffects( ia_drawScale, ca_drawDataChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_drawScale->ca_drawDataChanged" ); return stat; }

    {
        MFnEnumAttribute enumAttrFn;
        ia_solverType = enumAttrFn.create( "solverType", "sot", (short) RodTimeStepper::IMPL_EULER, & stat );
        CHECK_MSTATUS( stat );
        enumAttrFn.addField( "Implicit Euler",   (short) RodTimeStepper::IMPL_EULER );
        enumAttrFn.addField( "Symplectic Euler",  (short) RodTimeStepper::SYMPL_EULER );
        enumAttrFn.addField( "Statics",   (short) RodTimeStepper::STATICS );
        enumAttrFn.setKeyable( false );
        enumAttrFn.setStorable( true );
        enumAttrFn.setWritable( true );
        enumAttrFn.setReadable( true );
        stat = addAttribute( ia_solverType );
        CHECK_MSTATUS( stat );
    }
    stat = attributeAffects( ia_solverType, oa_rodsChanged );
    if ( !stat ) { stat.perror( "attributeAffects ia_solverType->oa_rodsChanged" ); return stat; }

    return MS::kSuccess;

}
