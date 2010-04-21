#include "WmFigRodNurbsInput.hh"

WmFigRodNurbsInput::WmFigRodNurbsInput( MObject& i_nurbsAttribute, bool i_lockFirstEdgeToInput ) : 
    m_inputNurbsAttribute( i_nurbsAttribute ), m_lockFirstEdgeToInput( i_lockFirstEdgeToInput )
{
    // we need to get pass the attribute here so that when we initialise data or
    // reload it we can just pull on the attribute
}

WmFigRodNurbsInput::~WmFigRodNurbsInput()
{
}

void WmFigRodNurbsInput::initialiseRodDataFromInput( MDataBlock& i_dataBlock, vector<RodData*>* i_pRodData )
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
        
        // Resize all the data vectors to make sure they're large enough for the number of cvs
        // we're getting from the input NURBS
        (*i_pRodData)[i]->allocateStorage( numCVs );

        vector<Vec3d> inputCurveVertices;
        inputCurveVertices.resize( numCVs );

        for ( int c = 0; c < numCVs; ++c )
        {
            MPoint cv;

            // FIXME: This should be world, why isn't it?
            // stat = inCurveFn.getCV( c,cv,MSpace::kWorld );
            stat = inCurveFn.getCV( c,cv,MSpace::kObject );
            CHECK_MSTATUS( stat );

            inputCurveVertices[ c ] = Vec3d( cv.x, cv.y, cv.z );            
        }
        (*i_pRodData)[ i ]->resetVertexPositions( inputCurveVertices );

        if ( m_lockFirstEdgeToInput )
        {
            // It doesnt matter that the edge data is garbage because we don't use the frames
            // from the nurb curve (since there are none!) We just add it so that Beaker
            // knows it isn't to be simulated. Since we didn't pass in a material frame
            // it will know not to use frame locking.
            (*i_pRodData)[ i ]->addKinematicEdge( 0 );
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