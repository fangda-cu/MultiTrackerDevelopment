#include "WmFigRodBarbInput.hh"

#include <maya/MGlobal.h>

WmFigRodBarbInput::WmFigRodBarbInput( MObject& i_verticesAttribute,
    MObject& i_strandRootFramesAttribute, double i_percentageOfBarbStrands,
    size_t i_verticesPerRod, bool i_lockFirstEdgeToInput ) : 
    m_verticesAttribute( i_verticesAttribute ),
    m_strandRootFramesAttribute( i_strandRootFramesAttribute ),
    m_percentageOfBarbStrands( i_percentageOfBarbStrands ),
    m_verticesPerRod( i_verticesPerRod ), m_lockFirstEdgeToInput( i_lockFirstEdgeToInput )
{
    // we need to get pass the attribute here so that when we initialise data or
    // reload it we can just pull on the attribute
}

WmFigRodBarbInput::~WmFigRodBarbInput()
{
}

void WmFigRodBarbInput::initialiseRodDataFromInput( MDataBlock& i_dataBlock, vector<RodData*>* i_pRodData )
{
    MStatus stat;

    MVectorArray vertices;
    size_t numStrands;
    getStrandVertices( i_dataBlock, vertices, numStrands );

    // store the material frames coming from barbershop
    vector<MaterialFrame> strandRootFrames;
    getStrandRootFrames( i_dataBlock, strandRootFrames );

    size_t inputVertexIndex = 0;
    for ( size_t i = 0; i < numStrands; i++ )
    {
        // Make sure we have enough space to store the date for each CV.
        (*i_pRodData)[i]->allocateStorage( m_verticesPerRod );

        vector<Vec3d> inputStrandVertices;
        inputStrandVertices.resize( m_verticesPerRod );

        for ( size_t v=0; v<m_verticesPerRod; v++ )
        {
            MVector vertex = vertices[ (int)inputVertexIndex++ ];

            inputStrandVertices[ v ] = Vec3d( vertex[0], vertex[1], vertex[2] );                                                                    
        }
        (*i_pRodData)[ i ]->resetVertexPositions( inputStrandVertices );

        // We need to add edge data so that each from the first edge will be locked to the input curve
        if ( strandRootFrames.size() > i && m_lockFirstEdgeToInput )
        {
            // The strand root frames may not have been connected for some reason so don't
            // rely on this having data

            // rod is probably null at this point but rodData will deal with it nicely
            (*i_pRodData)[ i ]->addKinematicEdge( 0, (*i_pRodData)[ i ]->rod, &strandRootFrames[ i ] );
        }
        else
        {   // If there's no root frame data then just lock the edge and it will not rotate.
            (*i_pRodData)[ i ]->addKinematicEdge( 0 );
        }
    }
}

void WmFigRodBarbInput::updateRodDataFromInput( MDataBlock& i_dataBlock, 
    std::vector<RodData*>* i_pRodData )
{
 
    MVectorArray strandVertices;
    size_t numStrands;

    getStrandVertices( i_dataBlock, strandVertices, numStrands );

    // If there are no strands to use to update then don't update
    if ( numStrands == 0 )
        return;

    // Every Barbershop strand has the same number of vertices
    size_t numStrandVertices = strandVertices.length() / numStrands;

    // store the material frames coming from barbershop
    vector<MaterialFrame> strandRootFrames;
    getStrandRootFrames( i_dataBlock, strandRootFrames );

    if ( numStrands != i_pRodData->size() )
    {
        MGlobal::displayError( "Number of Barbershop strands does not equal number of rods!"
                                "Did you change the input? Rewind sim to reset." );
        return;
    }

    size_t inputVertexIndex = 0;
    for ( size_t i = 0; i < numStrands; i++ )
    {
        BASim::ElasticRod* pRod = (*i_pRodData)[ i ]->rod;
        if ( pRod != NULL )
        {
            vector<Vec3d> inputStrandVertices;
            inputStrandVertices.resize( m_verticesPerRod );

            for ( int c = 0; c < m_verticesPerRod ; ++c )
            {
                MVector vertex = strandVertices[ (int)inputVertexIndex++ ];

                inputStrandVertices[ c ] = Vec3d( vertex.x, vertex.y, vertex.z );
            }

            (*i_pRodData)[ i ]->updateNextRodVertexPositions( inputStrandVertices );

            // The strand root frames may not have been connected for some reason so don't
            // rely on having that data
            if ( strandRootFrames.size() > i && m_lockFirstEdgeToInput )
            {
                    (*i_pRodData)[ i ]->updateKinematicEdge( 0, strandRootFrames[ i ] );
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

size_t WmFigRodBarbInput::numberOfInputs( MDataBlock& i_dataBlock )
{
    MStatus stat;

    MDataHandle verticesH = i_dataBlock.inputValue( m_verticesAttribute, &stat );
    CHECK_MSTATUS( stat );
    MFnVectorArrayData verticesData( verticesH.data(), &stat );
    CHECK_MSTATUS( stat );
    MVectorArray vertices = verticesData.array( &stat );

    return ( vertices.length() / m_verticesPerRod ) * ( m_percentageOfBarbStrands / 100.0 );
}

void WmFigRodBarbInput::getStrandVertices( MDataBlock& i_dataBlock, MVectorArray& o_vertices, 
    size_t& o_numStrands )
{
    MStatus stat;

    MDataHandle verticesH = i_dataBlock.inputValue( m_verticesAttribute, &stat );
    CHECK_MSTATUS( stat );
    MFnVectorArrayData verticesData( verticesH.data(), &stat );
    CHECK_MSTATUS( stat );
    o_vertices = verticesData.array( &stat );

    o_numStrands = o_vertices.length() / m_verticesPerRod;

    // The user may have requested that we use less than the total number
    // of Fozzie strands, so scale by that %
    o_numStrands *= m_percentageOfBarbStrands / 100.0;
}

void WmFigRodBarbInput::getStrandRootFrames( MDataBlock& i_dataBlock, 
    vector<MaterialFrame>& o_strandRootFrames )
{
    MStatus stat;

    MDataHandle strandRootFramesH = i_dataBlock.inputValue( m_strandRootFramesAttribute, &stat );
    CHECK_MSTATUS( stat );
    MFnVectorArrayData strandRootFrameVecData( strandRootFramesH.data(), &stat );
    
    if ( !stat )    // Strand root frames are not connected so just leave
        return;
        
    CHECK_MSTATUS( stat );

    MVectorArray strandRootFrameVec;
    strandRootFrameVec = strandRootFrameVecData.array( &stat );

    o_strandRootFrames.resize( strandRootFrameVec.length()/3 );
    MVector v;
    unsigned int idx = 0;
    for ( unsigned int rIdx=0; rIdx < o_strandRootFrames.size(); rIdx++ )
    {
        v = strandRootFrameVec[ idx++ ];
        o_strandRootFrames[ rIdx ].m1 = Vec3d( v[0], v[1], v[2] );
        v = strandRootFrameVec[ idx++ ];
        o_strandRootFrames[ rIdx ].m2 = Vec3d( v[0], v[1], v[2] );
        v = strandRootFrameVec[ idx++ ];
        o_strandRootFrames[ rIdx ].m3 = Vec3d( v[0], v[1], v[2] );
    }
}
