#include "WmSwAddRodContext.hh"

using namespace std;

namespace sweeney {

WmSwAddRodContext* WmSwAddRodContext::activeContext = NULL;
MString WmSwAddRodContext::typeName = "wmSwAddRodContext";

WmSwAddRodContext::WmSwAddRodContext() : m_helpString( "Add guide rods to a wmSweeneyNode" ),
    m_sweeneyNode( NULL )
{        
    setTitleString( "Add guide rods" );

    // Tell the context which XPM to use so the tool can properly
    // be a candidate for the 6th position on the mini-bar.
    setImage( "WmSwAddRodContext.xpm", MPxContext::kImage1 );
}

void WmSwAddRodContext::getClassName( MString& name ) const 
{
    name = typeName;
}

void WmSwAddRodContext::toolOnSetup(MEvent &)
{
    MStatus stat;
    setHelpString( m_helpString );
}

void WmSwAddRodContext::toolOffCleanup()
{
}

MStatus WmSwAddRodContext::doPress( MEvent& i_event ) 
{    
    MStatus stat;

    i_event.getPosition( m_xStartMouse, m_yStartMouse );
    
    m_xMouse = m_xStartMouse;
    m_yMouse = m_yStartMouse;

    m_toolCommand = (WmSwAddRodToolCommand*) newToolCommand();    
    
    utils::findSelectedSweeneyNodeAndRodManager( m_sweeneyNode, m_rodManager );

    return MS::kSuccess;
}

MStatus WmSwAddRodContext::doDrag( MEvent& i_event ) 
{
    MStatus stat;
    
    i_event.getPosition( m_xMouse, m_yMouse );

    MArgList args;
    args.addArg( m_xMouse );
    args.addArg( m_yMouse );

    m_toolCommand->doIt( args );

    drawManipulator();
                   
    return MS::kSuccess;
}

MStatus WmSwAddRodContext::doRelease( MEvent& i_event ) 
{
    MStatus stat;
    
    i_event.getPosition( m_xMouse, m_yMouse );
 
    M3dView view = M3dView::active3dView();
    view.refresh( true, true );
	
    m_toolCommand->finalize();

    return MS::kSuccess;
}

MStatus WmSwAddRodContext::doEnterRegion( MEvent & )
{
    return setHelpString( m_helpString );
}

void WmSwAddRodContext::drawManipulatorGL()
{    
    glDisable(GL_DEPTH_TEST);

    // We don't draw any manipulators for this context tool but leave the function here 
    // as it gets called by the node as most context tools do need to draw a manipulator

//    m_GLArrow.draw( ??? );

    glEnable(GL_DEPTH_TEST);
}

void WmSwAddRodContext::drawManipulator()
{
    M3dView view = M3dView::active3dView();

    view.refresh( true, true );

    view.beginGL();

    // Done drawing
    
#ifdef _WIN32
	SwapBuffers( view.deviceContext() );
#elif defined (OSMac_)
	// On 64-bit OSX the GL context is an NSOpenGLContext ptr. To swap
	// buffers using that requires Cocoa code, which we're not yet set
	// up to do. So we currently only do the buffer swap on 32-bit OSX.
# if !defined(__x86_64__)
	::aglSwapBuffers( view.display());   
# endif
#else
	glXSwapBuffers( view.display(), view.window() );
#endif // _WIN32

    view.endGL();
}


}
