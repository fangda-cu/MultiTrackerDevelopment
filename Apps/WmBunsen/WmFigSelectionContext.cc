#include "WmFigSelectionContext.hh"

#include <weta/Wmaya.hh>

#include <GL/glu.h>

using namespace std;

WmFigSelectionContext* WmFigSelectionContext::activeContext = NULL;
MString WmFigSelectionContext::typeName = "wmFigSelectionContext";

WmFigSelectionContext::WmFigSelectionContext() : m_helpString( "Click or drag to select rods" )
{        
    setTitleString( "Figaro rod selection tool" );

    // Tell the context which XPM to use so the tool can properly
    // be a candidate for the 6th position on the mini-bar.
    setImage( "WmFigSelectionContext.xpm", MPxContext::kImage1 );
}

void WmFigSelectionContext::getClassName( MString& name ) const 
{
    name = typeName;
}

void WmFigSelectionContext::toolOnSetup(MEvent &)
{
    MStatus stat;
    setHelpString( m_helpString );
}

void WmFigSelectionContext::toolOffCleanup()
{
}

MStatus WmFigSelectionContext::doPress( MEvent& i_event ) 
{    
    MStatus stat;

    i_event.getPosition( m_xStartMouse, m_yStartMouse );

    m_toolCommand = (WmFigSelectionToolCommand*) newToolCommand();        

    return MS::kSuccess;
}

MStatus WmFigSelectionContext::doDrag( MEvent& i_event ) 
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

MStatus WmFigSelectionContext::doRelease( MEvent& i_event ) 
{
    MStatus stat;

    m_toolCommand->finalize();

    return MS::kSuccess;		
}

MStatus WmFigSelectionContext::doEnterRegion( MEvent & )
{
    return setHelpString( m_helpString );
}


MStatus WmFigSelectionContext::drawMarqueeSelectBox() const
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

