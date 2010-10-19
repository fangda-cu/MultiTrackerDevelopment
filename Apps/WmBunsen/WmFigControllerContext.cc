#include "WmFigControllerContext.hh"

//#include <weta/Wmaya.hh>

#include <GL/glu.h>

using namespace std;

WmFigControllerContext* WmFigControllerContext::activeContext = NULL;
MString WmFigControllerContext::typeName = "wmFigControllerContext";

WmFigControllerContext::WmFigControllerContext() : m_helpString( "Click on a rod to attach a constraint" )
{        
    setTitleString( "Figaro tool to add constraints to rods" );

    // Tell the context which XPM to use so the tool can properly
    // be a candidate for the 6th position on the mini-bar.
    setImage( "WmFigControllerContext.xpm", MPxContext::kImage1 );
}

void WmFigControllerContext::getClassName( MString& name ) const 
{
    name = typeName;
}

void WmFigControllerContext::toolOnSetup(MEvent &)
{
    MStatus stat;
    setHelpString( m_helpString );
}

void WmFigControllerContext::toolOffCleanup()
{
}

MStatus WmFigControllerContext::doPress( MEvent& i_event ) 
{    
    MStatus stat;

    i_event.getPosition( m_xStartMouse, m_yStartMouse );

    m_toolCommand = (WmFigControllerToolCommand*) newToolCommand();        

    return MS::kSuccess;
}

MStatus WmFigControllerContext::doDrag( MEvent& i_event ) 
{
    MStatus stat;
    
    i_event.getPosition( m_xMouse, m_yMouse );
    drawMarqueeSelectBox();

    MArgList args;
    args.addArg( m_xMouse );
    args.addArg( m_yMouse );

    m_toolCommand->doIt( args );
                   
    return MS::kSuccess;
}

MStatus WmFigControllerContext::doRelease( MEvent& i_event ) 
{
    MStatus stat;
    
    i_event.getPosition( m_xMouse, m_yMouse );
 
    // Work out which rods were selected
    vector<int> rodIndices;
    vector<int> edgeIndices;
    searchForRodsIn2DScreenRectangle( rodIndices, edgeIndices );

    M3dView view = M3dView::active3dView();
    view.refresh( true, true );
	
    // Tell the command which rods the user selected in the context
    m_toolCommand->setSelectedRods( rodIndices );
    
    m_toolCommand->finalize();

    
    if ( rodIndices.size() > 0 )
    {
        cerr << "rod " << rodIndices[ 0 ] << ", edge " << edgeIndices[ 0 ] << endl;

        MGlobal::executeCommand( MString( "$rodNode = `ls -sl`; \
                                $locator = `spaceLocator -p 0 0 0`; \
                                select -add $rodNode; \
                                wmFigaro -aeo -ro " ) + rodIndices[ 0 ] + " -ed " + edgeIndices[ 0 ] + ";" );
    }

    return MS::kSuccess;
}

bool WmFigControllerContext::searchForRodsIn2DScreenRectangle( vector<int>& o_rodIndices, vector<int>& o_edgeIndices )
{
    MStatus stat;
    
    //////////////////////////////////////////////////////
    //
    // At first, no strands are found.
    //
    //////////////////////////////////////////////////////

    o_rodIndices.clear();
    o_edgeIndices.clear();    

    short selectedWidth = abs( m_xMouse - m_xStartMouse );
    short selectedHeight = abs( m_yMouse - m_yStartMouse );
    short selectedCentreX = ( m_xMouse + m_xStartMouse ) / 2.0;
    short selectedCentreY = ( m_yMouse + m_yStartMouse ) / 2.0;
    
    if ( selectedWidth == 0 && selectedHeight == 0 )
        return false;
    
    //////////////////////////////////////////////////////
    //
    // Set projection and modelview matrices stored from 
    // the last draw call.
    //
    //////////////////////////////////////////////////////
	
    glMatrixMode( GL_PROJECTION );
    glPushMatrix();
    glLoadMatrixf( projectionMatrix() );
    glMatrixMode( GL_MODELVIEW );
    glPushMatrix();
    glLoadMatrixf( modelViewMatrix() );
    
    //////////////////////////////////////////////////////
    //
    // Find the selected rod node in which to find
    // rod indices.
    //
    //////////////////////////////////////////////////////
    
    MSelectionList selectionList;
    MGlobal::getActiveSelectionList( selectionList );
    
    MDagPath rodNodeDagPath;
    MObject rodNodeObj = MObject::kNullObj;
    bool foundRodNode = true;
    
    if ( selectionList.length() == 1 )
    {
        selectionList.getDagPath( 0, rodNodeDagPath, rodNodeObj );

        MDagPath childPath;
        
        // Look for a child as the user probably selected the transform
        MFnDagNode rodDagNodeFn( rodNodeDagPath.child( 0, &stat ), &stat );
        if ( stat )
        {
            stat = rodDagNodeFn.getPath( childPath );
            CHECK_MSTATUS( stat );
            childPath.extendToShape();
        
            MFnDependencyNode nodeFn( childPath.node( &stat ) );
            CHECK_MSTATUS( stat );
        
            if ( nodeFn.typeName() == WmFigRodNode::typeName )
            {
                rodNodeObj = childPath.node();
            }
            else
                foundRodNode = false;
        }
        else // Perhaps no child as the user selected the shape node directly
        {
            MFnDependencyNode nodeFn( rodNodeDagPath.node( &stat ) );
            CHECK_MSTATUS( stat );
        
            if ( nodeFn.typeName() == WmFigRodNode::typeName )
            {
                // rodNodeObj = nodeFn.node();
                rodNodeObj = rodNodeDagPath.node( &stat );
            }
            else
                foundRodNode = false;
        }
    }
    else
        foundRodNode = false;
    
    if ( !foundRodNode )
    {
        MGlobal::displayError( "Please select a single wmFigRodNode to find rods within." );
        return false;
    }
    
    MFnDependencyNode rodNodeDepFn( rodNodeObj );
    WmFigRodNode* rodNode = static_cast<WmFigRodNode*>( rodNodeDepFn.userNode() );
    
    if ( !rodNode )
        return false;
    
    //////////////////////////////////////////////////////
    //
    // Do GL selection.
    //
    //////////////////////////////////////////////////////
    
    vector<GLuint> selectionBuffer;
    GLint numHits = findRodsUsingOpenGLSelection( selectedCentreX, selectedCentreY, selectedWidth, 
                        selectedHeight, rodNode, selectionBuffer );
    
    
    //////////////////////////////////////////////////////
    //
    // Process hits. 
    //
    //////////////////////////////////////////////////////
    
    if ( numHits > 0 )
    {
        o_rodIndices.resize( numHits );
        o_edgeIndices.resize( numHits );
        
        ///////////////////////////////////////////////////
        // Skip through the results array picking out the 
        // indices of rods that were hit. We can safely
        // ignore most stuff as we know there is only
        // one name per rod due to the way we did selection
        // below.
        ///////////////////////////////////////////////////
        
        size_t index = 0;
        for ( size_t i=0; i<numHits; i++ )
        {
            index += 3;
            o_rodIndices[ i ] = selectionBuffer[ index ] >> 16;
            o_edgeIndices[ i ] = selectionBuffer[ index ] & 0xFFFF;
            index++;
        }
    }
    
    glMatrixMode( GL_MODELVIEW );
    glPopMatrix();
    glMatrixMode( GL_PROJECTION );
    glPopMatrix();
    
    return true;
}

GLint WmFigControllerContext::findRodsUsingOpenGLSelection( const double i_centreX, const double i_centreY,
    const double i_width, const double i_height, WmFigRodNode* i_rodNode,
    vector<GLuint>& o_selectionBuffer )
{
    WmFigRodGroup* rodGroup = i_rodNode->rodGroup();
    if ( rodGroup->numberOfRealRods() == 0 )
        return 0;
    
    GLint viewport[4];
    glGetIntegerv( GL_VIEWPORT, viewport );
    
    GLfloat projectionMatrix[16];
    glGetFloatv( GL_PROJECTION_MATRIX, projectionMatrix );
    
    glMatrixMode( GL_PROJECTION );
    glPushMatrix();
    glLoadIdentity();
    gluPickMatrix( i_centreX, i_centreY, i_width, i_height, viewport );
    glMultMatrixf( projectionMatrix );
    
    const int nRods = rodGroup->numberOfRods();

    // *4 because selection returns a bunch of stuff with each hit.
    o_selectionBuffer.resize( nRods * 4 );
    GLint numHits;

    glSelectBuffer( nRods * 4, &(o_selectionBuffer[0]) );
    glRenderMode( GL_SELECT );
    
    glInitNames();
	glPushName( 0 );
	
    for( int r = 0; r < nRods; ++r )
    {       
        //glLoadName( (GLuint) r << 16 );
 
     //   glBegin( GL_LINE_STRIP );        
        for( int v = 1, n = rodGroup->elasticRod( r )->nv(); v < n; ++v )
        {
            glLoadName( ( (GLuint) r << 16 ) + v );
            
            glBegin( GL_LINES );

            const Vec3d p = rodGroup->elasticRod( r )->getVertex( v - 1 );
            const Vec3d p2 = rodGroup->elasticRod( r )->getVertex( v );
            glVertex3f( p[0], p[1], p[2] );
            glVertex3f( p2[0], p2[1], p2[2] );

            glEnd();
        }

       // glEnd();
    }
    
    numHits = glRenderMode( GL_RENDER );
    
    glMatrixMode( GL_PROJECTION );
    glPopMatrix();

    return numHits;
}



MStatus WmFigControllerContext::doEnterRegion( MEvent & )
{
    return setHelpString( m_helpString );
}


MStatus WmFigControllerContext::drawMarqueeSelectBox() const
{
    M3dView view = M3dView::active3dView();
    view.refresh( true, true );
	
    view.beginGL();
	
    // Store all the previous draw settings so we can return them and dont break Maya
    GLboolean depthTest[ 1 ];
    glGetBooleanv( GL_DEPTH_TEST, depthTest );
    GLboolean colorLogicOp[ 1 ];
    glGetBooleanv( GL_COLOR_LOGIC_OP, colorLogicOp );
    GLboolean lineStipple[ 1 ];
    glGetBooleanv( GL_LINE_STIPPLE, lineStipple );

    glMatrixMode( GL_PROJECTION );
    glLoadIdentity();
    gluOrtho2D( 0.0, (GLdouble) view.portWidth(),
                0.0, (GLdouble) view.portHeight() );
    
    glMatrixMode( GL_MODELVIEW );
    glPushMatrix();
    glLoadIdentity();
    
    glIndexi( 2 );
    glLineWidth( 2.0 );
    glLineStipple( 1, 0x1111 );
    glEnable( GL_LINE_STIPPLE );
    glColor3f( 0.0, 0.0, 0.0 );	
    
    glBegin( GL_LINE_LOOP );
        glVertex2f( m_xStartMouse, m_yStartMouse );
        glVertex2f( m_xStartMouse, m_yMouse );
        glVertex2f( m_xMouse, m_yMouse );
        glVertex2f( m_xMouse, m_yStartMouse );
    glEnd();
    
    // Restore everything back the way it was so Maya is happy.
    glMatrixMode( GL_MODELVIEW );
    glPopMatrix();

    glLineWidth( 1.0f );

    glColor3f( 1.0f, 1.0f, 1.0f );
    glPointSize( 5.0f );
		
    glDisable( GL_COLOR_LOGIC_OP );
    glDisable( GL_BLEND );

    if ( colorLogicOp[0] ) glEnable( GL_COLOR_LOGIC_OP );
    else glDisable( GL_COLOR_LOGIC_OP );

    if ( depthTest[0] ) glEnable( GL_DEPTH_TEST );
    else glDisable( GL_DEPTH_TEST );

    if ( lineStipple[0] ) glEnable( GL_LINE_STIPPLE );
    else glDisable( GL_LINE_STIPPLE );

    glXSwapBuffers( view.display(), view.window() );

    view.endGL();

    return MStatus::kSuccess;
}

const float* WmFigControllerContext::projectionMatrix()
{
    MStatus status;
    
    // get the currently active M3dView.
    
    M3dView view = M3dView::active3dView( &status );
    
    if( status == MS::kSuccess )
    {
        // get the dag path of the current camera
        
        MDagPath cameraPath;
        status = view.getCamera( cameraPath );
        CHECK_MSTATUS( status );
        
        // create camera funcion set
        
        MFnCamera cameraFn( cameraPath, &status );
        CHECK_MSTATUS( status );
        
        // get the viewport width and height

        const int width = view.portWidth( &status );
        CHECK_MSTATUS( status );

        const int height = view.portHeight( &status );
        CHECK_MSTATUS( status );
        
        // window aspect

        const double aspect = double( width ) / double( height );
        
        // get the viewing frustum

        double dleft, dright, dbottom, dtop;
        status = cameraFn.getViewingFrustum( aspect, dleft, dright, dbottom, dtop, true, false );
        CHECK_MSTATUS( status );
        
        // get the near and far clipping planes

        const float near = cameraFn.nearClippingPlane( &status );
        CHECK_MSTATUS( status );

        const float far = cameraFn.farClippingPlane( &status );
        CHECK_MSTATUS( status );
        
        // is the camera orthographic
        
        const bool orthographic = cameraFn.isOrtho( &status );
        
        // stuff common to both a perspective and orthographic projection
        
        const float left   = static_cast<float>( dleft );
        const float right  = static_cast<float>( dright );
        const float bottom = static_cast<float>( dbottom );
        const float top    = static_cast<float>( dtop );
        
        const float rightMinusLeft = right - left;
        const float farMinusNear = far - near;
        const float topMinusBottom = top - bottom;
        
        m_projectionMatrix.matrix[0][1] = m_projectionMatrix.matrix[0][2] = m_projectionMatrix.matrix[0][3] = 0.f;
        m_projectionMatrix.matrix[1][0] = m_projectionMatrix.matrix[1][2] = m_projectionMatrix.matrix[1][3] = 0.f;
        
        if( orthographic )
        {
            // orthographic matrix
            
            m_projectionMatrix.matrix[2][0] = m_projectionMatrix.matrix[2][1] = m_projectionMatrix.matrix[2][3] = 0.f;
            m_projectionMatrix.matrix[3][3] = 1.f;
            
            m_projectionMatrix.matrix[0][0] = 2.f / rightMinusLeft;
            m_projectionMatrix.matrix[1][1] = 2.f / topMinusBottom;
            m_projectionMatrix.matrix[2][2] = -( 2.f / farMinusNear );
            
            m_projectionMatrix.matrix[3][0] = -( ( right + left ) / rightMinusLeft );
            m_projectionMatrix.matrix[3][1] = -( ( top + bottom ) / topMinusBottom );
            m_projectionMatrix.matrix[3][2] = -( ( far + near ) / farMinusNear );
        }
        else
        {
            // perspective matrix
            
            const float _2near = 2.f * near;
            
            m_projectionMatrix.matrix[3][0] = m_projectionMatrix.matrix[3][1] = m_projectionMatrix.matrix[3][3] = 0.f;
            m_projectionMatrix.matrix[2][3] = -1.f;
            
            m_projectionMatrix.matrix[0][0] = _2near / rightMinusLeft;
            m_projectionMatrix.matrix[1][1] = _2near / topMinusBottom;
            
            m_projectionMatrix.matrix[2][2] = -( ( far + near ) / farMinusNear );
            m_projectionMatrix.matrix[3][2] = -( ( _2near * far ) / farMinusNear );
            
            m_projectionMatrix.matrix[2][0] = ( right + left ) / rightMinusLeft;
            m_projectionMatrix.matrix[2][1] = ( top + bottom ) / topMinusBottom;
        }
    }
    else
    {
        MGlobal::displayError( "WmBarbFurSetNode::projectionMatrix() : failed to retrieve active M3dView." );
    }
    
    return &( m_projectionMatrix.matrix[0][0] );
}

const float* WmFigControllerContext::modelViewMatrix()
{
    MStatus status;
    
    // get the currently active M3dView.
    
    M3dView view = M3dView::active3dView( &status );
    
    if( status == MS::kSuccess )
    {
        // get the dag path of the current camera
        
        MDagPath cameraPath;
        status = view.getCamera( cameraPath );
        CHECK_MSTATUS( status );
        
        // create camera funcion set
        
        MFnCamera cameraFn( cameraPath, &status );
        CHECK_MSTATUS( status );
        
        // get the camera eye point, view direction and up vector
        
        MPoint eyePoint = cameraFn.eyePoint( MSpace::kWorld, &status );
        CHECK_MSTATUS( status );
        
        MVector viewDirection = cameraFn.viewDirection( MSpace::kWorld, &status );
        CHECK_MSTATUS( status );
        
        MVector upDirection = cameraFn.upDirection( MSpace::kWorld, &status );
        CHECK_MSTATUS( status );
        
        // build view matrix 
        // NOTE : the furSets are in world space we don't need transformation of the furSet.
        
        viewDirection.normalize();
        upDirection.normalize();
        
        const MVector side = viewDirection ^ upDirection;   // cross product
        const MVector up   = side          ^ viewDirection; // cross product
        
        m_modelViewMatrix.matrix[0][3] = m_modelViewMatrix.matrix[1][3] = m_modelViewMatrix.matrix[2][3] = 0.f;
        m_modelViewMatrix.matrix[3][0] = m_modelViewMatrix.matrix[3][1] = m_modelViewMatrix.matrix[3][2] = 0.f;
        m_modelViewMatrix.matrix[3][3] = 1.f;
        
        m_modelViewMatrix.matrix[0][0] = side[0];
        m_modelViewMatrix.matrix[1][0] = side[1];
        m_modelViewMatrix.matrix[2][0] = side[2];
        
        m_modelViewMatrix.matrix[0][1] = up[0];
        m_modelViewMatrix.matrix[1][1] = up[1];
        m_modelViewMatrix.matrix[2][1] = up[2];
        
        m_modelViewMatrix.matrix[0][2] = -( viewDirection[0] );
        m_modelViewMatrix.matrix[1][2] = -( viewDirection[1] );
        m_modelViewMatrix.matrix[2][2] = -( viewDirection[2] );
        
        MFloatMatrix eyeTranslation;
        eyeTranslation.matrix[3][0] = -( eyePoint[0] );
        eyeTranslation.matrix[3][1] = -( eyePoint[1] );
        eyeTranslation.matrix[3][2] = -( eyePoint[2] );
        
        m_modelViewMatrix = eyeTranslation * m_modelViewMatrix;
    }
    else
    {
        MGlobal::displayError( "WmBarbFurSetNode::modelViewMatrix() : failed to retrieve active M3dView." );
    }

    return &( m_modelViewMatrix.matrix[0][0] );
}


