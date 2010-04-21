#include "WmFigRodBarbInput.hh"

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

    MDataHandle verticesH = i_dataBlock.inputValue( m_verticesAttribute, &stat );
    CHECK_MSTATUS( stat );
    MFnVectorArrayData verticesData( verticesH.data(), &stat );
    CHECK_MSTATUS( stat );
    MVectorArray vertices = verticesData.array( &stat );

    size_t numStrands = vertices.length() / m_verticesPerRod;

    // The user may have requested that we use less than the total number
    // of Fozzie strands, so scale by that %
    numStrands *= m_percentageOfBarbStrands / 100.0;

    // store the material frames coming from barbershop
    vector<MaterialFrame> strandRootFrames;
    getStrandRootFrames( i_dataBlock, strandRootFrames );

    size_t inputVertexIndex = 0;
    for ( size_t i = 0; i < numStrands; i++ )
    {
        // Make sure we have enough space to store the date for each CV.
        (*i_pRodData)[i]->allocateStorage( m_verticesPerRod );

        vector<Vec3d> inputCurveVertices;
        inputCurveVertices.resize( m_verticesPerRod );

        for ( size_t v=0; v<m_verticesPerRod; v++ )
        {
            MVector vertex = vertices[ (int)inputVertexIndex++ ];

            inputCurveVertices[ v ] = Vec3d( vertex[0], vertex[1], vertex[2] );                                                                    
        }
        (*i_pRodData)[ i ]->resetVertexPositions( inputCurveVertices );

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
