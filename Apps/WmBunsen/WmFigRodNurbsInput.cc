#include "WmFigRodNurbsInput.hh"

#include <maya/MGlobal.h>

WmFigRodNurbsInput::WmFigRodNurbsInput( MObject& i_nurbsAttribute, bool i_lockFirstEdgeToInput,
    WmFigRodGroup& i_rodGroup, double i_vertexSpacing, double i_minimumRodLength, RodOptions& i_rodOptions,
    double i_massDamping, BASim::Vec3d& i_gravity,  RodTimeStepper::Method i_solverType, 
    std::set< int >& i_simulationSet, const bool i_doReverseHairdo ):    
    m_inputNurbsAttribute( i_nurbsAttribute ), m_lockFirstEdgeToInput( i_lockFirstEdgeToInput ),
    m_rodGroup( i_rodGroup ), m_vertexSpacing( i_vertexSpacing ), m_minimumRodLength( i_minimumRodLength ),
    m_rodOptions( i_rodOptions ), m_massDamping( i_massDamping ), m_gravity( i_gravity ),  m_solverType( i_solverType ),
    m_simulationSet( i_simulationSet ), m_doReverseHairdo( i_doReverseHairdo )
{
    m_simulating = true;    
}

WmFigRodNurbsInput::~WmFigRodNurbsInput()
{
}

void WmFigRodNurbsInput::getAndResampleInputCurves( MDataBlock& i_dataBlock, vector< vector<BASim::Vec3d > >& o_inputCurveVertices )
{
    MStatus stat;

    MDataHandle inputCurveH;
    MObject inputCurveObj;
    MArrayDataHandle inArrayH = i_dataBlock.inputArrayValue( m_inputNurbsAttribute, &stat );
    CHECK_MSTATUS(stat);

    size_t numCurvesConnected = inArrayH.elementCount();

    o_inputCurveVertices.resize( numCurvesConnected );

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
        int numCVs = inCurveFn.numCVs();

        vector< MVector > curve;
        curve.resize( numCVs );

        // If the user has specificed a vertex spacing > 0 then they want to override the vertices
        // that came from Barbershop. Very sensible of them as the rods stability is much reduced 
        // with vertices that are very close together.
        if ( m_vertexSpacing > 0.0 && m_simulating)
        {
            // First find how many vertices we need based on the length of the strand
            
            // FIXME: This should be world, why isn't it?
            // stat = inCurveFn.getCV( c,cv,MSpace::kWorld );
            MPoint cv;
            double length = 0.0;

            for ( int v = 0; v < numCVs; v++ )
            {
                stat = inCurveFn.getCV( v, cv, MSpace::kObject );
                CHECK_MSTATUS( stat );

                curve[ v ] = cv;

                if ( v > 0 )
                {
                    length += ( curve[ v ] - curve[ v - 1 ] ).length();
                }
            }

            if ( length < m_minimumRodLength )
            {
                MGlobal::displayWarning( MString( "Input curve " ) + i + " is less than minimum length (" + m_minimumRodLength + "), not turning it into a rod." );

                // Remove all vertices from the curve so it doesn't get turned into a rod later
                curve.clear();
            }
            else
            {
                int numVerticesRequired = length / m_vertexSpacing;
    
                if ( numVerticesRequired < 2 )
                {
                    numVerticesRequired = 2;
                    MGlobal::displayWarning( MString( "Input strand " ) + i + " was going to have less than 2 vertices in the rod, setting it to 2." );
                }
    
                vector< MVector > resampledCurve;
                resampleCurve( numVerticesRequired, curve, resampledCurve );

                // Store the new curve back on curve as we use it below to create the actual rod
                curve = resampledCurve;
            }
        }
        else
        {
            for ( int c = 0; c < numCVs; ++c )
            {
                // FIXME: This should be world, why isn't it?
                // stat = inCurveFn.getCV( c,cv,MSpace::kWorld );
                stat = inCurveFn.getCV( c, cv, MSpace::kObject );
                CHECK_MSTATUS( stat );

                curve[ c ] = cv;
            }
        }

        // convert cvs to BASim::Vec3d
        o_inputCurveVertices[ i ].resize( curve.size() );
        for ( int c = 0; c < curve.size(); ++c )
        {
            o_inputCurveVertices[ i ][ c ] = BASim::Vec3d( curve[ c ].x, curve[ c ].y, curve[ c ].z );
        }
    }
}

void WmFigRodNurbsInput::initialiseRodDataFromInput( MDataBlock& i_dataBlock)
{
    vector< vector< BASim::Vec3d > > inputCurveVertices;

    getAndResampleInputCurves( i_dataBlock, inputCurveVertices );

    for ( int c=0; c< (int)inputCurveVertices.size(); ++c )
    {
        bool isPlaceHolderRod = false;

        if ( m_simulationSet.size() != 0 )
        {
            isPlaceHolderRod = ( m_simulationSet.count( c ) == 0 );
        }

        if ( inputCurveVertices[ c ].size() == 0 || isPlaceHolderRod)
        {
             // Add a placeholder rod, most likely too small to actually be a real rod.
             m_rodGroup.addRod();
             cerr << "Adding fake rod!\n";
        }
        else
        {
            RodOptions rodOptions = m_rodOptions;
            rodOptions.numVertices = (int)inputCurveVertices[ c ].size();

            // Mass damping should be in rod options, it's dumb to pass it seperately.
            int rodIndex = m_rodGroup.addRod( inputCurveVertices[ c ], rodOptions, m_massDamping, 
                                            m_gravity, m_solverType, false, m_doReverseHairdo );
                                            
            if ( m_lockFirstEdgeToInput && !m_rodGroup.isPlaceHolderRod( rodIndex ) )
            {
                m_rodGroup.addKinematicEdge( rodIndex, 0 );
            }
        }
    }
    
    matchRodToInputIfRequired( m_rodGroup, inputCurveVertices );
}

void WmFigRodNurbsInput::updateRodDataFromInput( MDataBlock& i_dataBlock )
{
    vector< vector< BASim::Vec3d > > inputCurveVertices;

    getAndResampleInputCurves( i_dataBlock, inputCurveVertices );

    for ( int c=0; c< (int)inputCurveVertices.size(); ++c )
    {
        if ( !m_rodGroup.isPlaceHolderRod( c ) )        
        {
            m_rodGroup.updateRodNextVertexPositions( c, inputCurveVertices[ c ] );
        }
    }
    
    matchRodToInputIfRequired( m_rodGroup, inputCurveVertices );
}

int WmFigRodNurbsInput::numberOfInputs( MDataBlock& i_dataBlock )
{
    MStatus stat;

    MArrayDataHandle inArrayH = i_dataBlock.inputArrayValue( m_inputNurbsAttribute, &stat );
    CHECK_MSTATUS(stat);
    return inArrayH.elementCount(); 
}
