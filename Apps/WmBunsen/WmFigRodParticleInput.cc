#include "WmFigRodParticleInput.hh"

#include <maya/MGlobal.h>


WmFigRodParticleInput::WmFigRodParticleInput( MObject& i_verticesAttribute, MObject& i_perRodParticleCountAttribute,
                        RodOptions& i_rodOptions,
                        double i_massDamping, BASim::Vec3d& i_gravity, WmFigRodGroup& i_rodGroup,
                        RodTimeStepper::Method i_solverType, std::set< int >& i_simulationSet )
 : 
    m_verticesAttribute( i_verticesAttribute ),
    m_perRodParticleCountAttribute( i_perRodParticleCountAttribute ),        
    m_rodGroup( i_rodGroup ), m_rodOptions( i_rodOptions ), m_massDamping( i_massDamping ),
    m_gravity( i_gravity ),
    m_solverType( i_solverType ), m_simulationSet( i_simulationSet )
{
}

WmFigRodParticleInput::~WmFigRodParticleInput()
{
}

void WmFigRodParticleInput::initialiseRodDataFromInput( MDataBlock& i_dataBlock )
{
    MStatus stat;

    MVectorArray vertices;
    MIntArray perRodParticleCount;
    getRodVertices( i_dataBlock, vertices, perRodParticleCount );

    m_rodGroup.removeAllRods();

    int numRods = perRodParticleCount.length();

    // The vertices array contains the vertices of *all* input strands
    // so we need to keep a seperate counter for it or it will always be giving
    // us the first strand.
    unsigned int inputStrandVertexIndex = 0;
    vector< vector<BASim::Vec3d> > inputStrandVertices;    
    
    inputStrandVertices.resize( numRods );
    
    for ( int r = 0; r < numRods; ++r )
    {               
        int verticesInRod = perRodParticleCount[ r ];

        inputStrandVertices[ r ].resize( verticesInRod );
        
        for ( int v = 0; v < verticesInRod; v++ )
        {
            MVector vertex = vertices[ inputStrandVertexIndex++ ];

            inputStrandVertices[ r ][ v ] = BASim::Vec3d( vertex[0], vertex[1], vertex[2] );                                          
        }
    }
    
    // Pretend the rods are coming from a cache as they are not going to be simmed
    // FIXME: Sending mass damping is pointless!
    m_rodGroup.addRodsFromCache( inputStrandVertices, m_rodOptions, 10.0 );    
}

void WmFigRodParticleInput::updateRodDataFromInput( MDataBlock& i_dataBlock )
{
    MVectorArray vertices;
    MIntArray perRodParticleCount;
    getRodVertices( i_dataBlock, vertices, perRodParticleCount );

    int numInputRods = perRodParticleCount.length();

    // If there are no particle curves to use for update then don't update
    if ( numInputRods == 0 )
    {
        cerr << "No input to update with\n";
        return;
    }
    
    if ( numInputRods != m_rodGroup.numberOfRods() )
    {
        MGlobal::displayError( "Number of particle curves does not equal number of rods!"
                                "Did you change the input? Rewind sim to reset." );
        return;
    }

    int inputStrandVertexIndex = 0;
    for ( int i = 0; i < numInputRods; i++ )
    {
        int numVerticesInRod = m_rodGroup.numberOfVerticesInRod( i );

        vector<BASim::Vec3d> inputStrandVertices;         
        inputStrandVertices.resize( numVerticesInRod );

        for ( int v = 0; v < numVerticesInRod; v++ )
        {
            MVector vertex = vertices[ inputStrandVertexIndex++ ];
            inputStrandVertices[ v ] = BASim::Vec3d( vertex[ 0 ], vertex[ 1 ], vertex[ 2 ] );
        }

        m_rodGroup.updateRodNextVertexPositions( i, inputStrandVertices );        
    }
}

int WmFigRodParticleInput::numberOfInputs( MDataBlock& i_dataBlock )
{
    MStatus stat;

    MDataHandle perRodParticleCountH = i_dataBlock.inputValue( m_perRodParticleCountAttribute, &stat );
    CHECK_MSTATUS( stat );
    MFnIntArrayData perRodParticleCountData( perRodParticleCountH.data(), &stat );
    CHECK_MSTATUS( stat );
    MIntArray perRodParticleCount = perRodParticleCountData.array( &stat );

    return perRodParticleCount.length();
}

void WmFigRodParticleInput::getRodVertices( MDataBlock& i_dataBlock, MVectorArray& o_vertices, 
                                             MIntArray& o_perRodParticleCount )

{
    MStatus stat;

    MDataHandle verticesH = i_dataBlock.inputValue( m_verticesAttribute, &stat );
    CHECK_MSTATUS( stat );
    MFnVectorArrayData verticesData( verticesH.data(), &stat );
    CHECK_MSTATUS( stat );
    o_vertices = verticesData.array( &stat );
    
    MDataHandle perRodParticleCountH = i_dataBlock.inputValue( m_perRodParticleCountAttribute, &stat );
    CHECK_MSTATUS( stat );
    MFnIntArrayData perRodParticleCountData( perRodParticleCountH.data(), &stat );
    CHECK_MSTATUS( stat );
    o_perRodParticleCount = perRodParticleCountData.array( &stat );    
}
