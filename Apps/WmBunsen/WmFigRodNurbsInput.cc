#include "WmFigRodNurbsInput.hh"

#include <maya/MGlobal.h>

WmFigRodNurbsInput::WmFigRodNurbsInput( MObject& i_nurbsAttribute, bool i_lockFirstEdgeToInput,
    WmFigRodGroup& i_rodGroup, double i_vertexSpacing, double i_minimumRodLength, RodOptions& i_rodOptions,
    double i_massDamping ) : 
    m_inputNurbsAttribute( i_nurbsAttribute ), m_lockFirstEdgeToInput( i_lockFirstEdgeToInput ),
    m_rodGroup( i_rodGroup ), m_vertexSpacing( i_vertexSpacing ), m_minimumRodLength( i_minimumRodLength ),
    m_rodOptions( i_rodOptions ), m_massDamping( i_massDamping )
{
    // we need to get pass the attribute here so that when we initialise data or
    // reload it we can just pull on the attribute
}

WmFigRodNurbsInput::~WmFigRodNurbsInput()
{
}

void WmFigRodNurbsInput::initialiseRodDataFromInput( MDataBlock& i_dataBlock )
{
    MStatus stat;

    MDataHandle inputCurveH;
    MObject inputCurveObj;
    MArrayDataHandle inArrayH = i_dataBlock.inputArrayValue( m_inputNurbsAttribute, &stat );
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
        int numCVs = inCurveFn.numCVs();;

        vector< MVector > curve;
        curve.resize( numCVs );

        // If the user has specificed a vertex spacing > 0 then they want to override the vertices
        // that came from Barbershop. Very sensible of them as the rods stability is much reduced 
        // with vertices that are very close together.
        if ( m_vertexSpacing > 0.0 )
        {
            // First find how many vertices we need based on the length of the strand
            
            // FIXME: This should be world, why isn't it?
            // stat = inCurveFn.getCV( c,cv,MSpace::kWorld );
            MPoint cv;
            double length = 0.0;

            for ( size_t v = 0; v < numCVs; v++ )
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

        // convert cvs to Vec3d
        vector<Vec3d> inputCurveVertices;
        inputCurveVertices.resize( curve.size() );
        for ( int c = 0; c < curve.size(); ++c )
        {
            inputCurveVertices[ c ] = Vec3d( curve[ c ].x, curve[ c ].y, curve[ c ].z );
        }

        RodOptions rodOptions = m_rodOptions;
        rodOptions.numVertices = inputCurveVertices.size();

        // Mass damping should be in rod options, it's dumb to pass it seperately.
        size_t rodIndex = m_rodGroup.addRod( inputCurveVertices, rodOptions, m_massDamping );

        if ( m_lockFirstEdgeToInput )
        {
            m_rodGroup.addKinematicEdge( rodIndex, 0 );
        }
    }

    // FIXME:
    // The below code is what was used when this function was in WmBunsenRodNode to track
    // edges that were being controlled by attached Maya objects. It should work and should be
    // brought back to life as soon as the refactoring is done.

/*
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
  */  
}

void WmFigRodNurbsInput::updateRodDataFromInput( MDataBlock& i_dataBlock )
{
   /* MStatus stat;
    
    MArrayDataHandle inArrayH = i_dataBlock.inputArrayValue( m_inputNurbsAttribute, &stat );
    CHECK_MSTATUS(stat);
    size_t numCurvesConnected = inArrayH.elementCount(); 

    if ( i_pRodData->size() != numCurvesConnected )
    {
        MGlobal::displayError( "Number of rods does not equal number of input curves, rewind simulation to reset" );
        return;
    }

    for ( unsigned int i = 0; i < numCurvesConnected; i++ )
    {
        inArrayH.jumpToElement( i );
        MDataHandle inputCurveH = inArrayH.inputValue( &stat );
        CHECK_MSTATUS( stat );

        // Use asNurbsCurveTransformed to get the curve data as we
        // want it in world space.
        MObject inputCurveObj = inputCurveH.asNurbsCurveTransformed();
        MFnNurbsCurve inCurveFn( inputCurveObj );

        MPoint cv;
        int numCurveCVs = inCurveFn.numCVs();

        BASim::ElasticRod* pRod = (*i_pRodData)[ i ]->rod;
        if ( pRod != NULL )
        {
            int numRodVertices = pRod->nv();

            if ( numRodVertices != numCurveCVs )
            {
                MGlobal::displayError( "Number of vertices in rod does not equal number of CVs in input curve! Did you change the input? Rewind sim to reset." );
                return;
            }

            vector<Vec3d> inputCurveVertices;
            inputCurveVertices.resize( numCurveCVs );

            for ( int c = 0; c < numCurveCVs ; ++c )
            {
                MPoint cv;
                // stat = inCurveFn.getCV( c,cv,MSpace::kWorld );
                stat = inCurveFn.getCV( c,cv,MSpace::kObject );
                CHECK_MSTATUS( stat );

                inputCurveVertices[ c ] = Vec3d( cv.x, cv.y, cv.z );
            }

            (*i_pRodData)[ i ]->updateNextRodVertexPositions( inputCurveVertices );
            
            /*  FIXME: Again this should be resurrected as soon as the refactoring is complete.

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
            }*/
      /* }
    }*/
}

size_t WmFigRodNurbsInput::numberOfInputs( MDataBlock& i_dataBlock )
{
    MStatus stat;

    MArrayDataHandle inArrayH = i_dataBlock.inputArrayValue( m_inputNurbsAttribute, &stat );
    CHECK_MSTATUS(stat);
    return inArrayH.elementCount(); 
}