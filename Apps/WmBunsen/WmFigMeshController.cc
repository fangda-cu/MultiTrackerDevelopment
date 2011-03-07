#include "WmFigMeshController.hh"

using namespace BASim;

WmFigMeshController::WmFigMeshController( BASim::TriangleMesh* i_currentMesh, 
                                          BASim::TriangleMesh* i_previousMesh, 
                                          BASim::TriangleMesh* i_nextMesh, 
                                           double i_time, double i_dt ) 
                                         : m_currentMesh( i_currentMesh ),
                                           m_previousMesh( i_previousMesh ),
                                           m_nextMesh( i_nextMesh ),
                                           BASim::ScriptingController( i_time, i_dt ),
                                           m_startMeshTime( 0.0 ), m_endMeshTime( 0.0 )

// NOTE: We pass i_time and i_dt to BASim::ScriptingController but they are actually *ignored*
//       as BridsonStepper calls setDt() on the controller when it executes anyway.
//       We could probably just not bother getting them passed in here but we keep it this way
//       in case Breannan changes the way BridsonStepper works so we'll just copy his examples.

{   
    m_startTime = i_time;
    m_previousMayaTime = m_nextMayaTime = m_startTime;
}
 
bool WmFigMeshController::execute()
{
    // TODO: Assumes 24fps for Maya scenes, fix that

    //dt = fps / substeps
    //substeps = fps / dt

    //double interpolatedTime = getDt() / ( 1 / 24.0 );

    // TODO: Fix this hacked in time
    //double interpolatedTime = ( ( getTime() / getDt() ) / 10.0;

    std::cerr << "dt = " << getDt() << std::endl;
    std::cerr << "executing controller at time " << getTime() << "\n ";
    std::cerr << "executing controller at Maya time " << ( getTime() - m_startTime ) * 24.0 + m_startTime << "\n ";
    std::cerr << "current start/end mesh time = " << m_startMeshTime << " / " << m_endMeshTime << std::endl;
   
    std::cerr << "m_previousMayaTime = " << m_previousMayaTime << std::endl;
    std::cerr << "m_nextMayaTime = " << m_nextMayaTime << std::endl;

    double interpolation = ( getTime() - m_startTime ) * 24.0 + m_startTime - m_previousMayaTime;
    std::cerr << "interpolated mesh factor = " << interpolation << std::endl;
 
    for( TriangleMesh::vertex_iter previousMeshItr = m_previousMesh->vertices_begin(),
         nextMeshItr = m_nextMesh->vertices_begin(),
         currentMeshItr = m_currentMesh->vertices_begin(); 
         previousMeshItr != m_previousMesh->vertices_end(); 
         ++previousMeshItr, ++nextMeshItr, ++currentMeshItr )
    {
        BASim::Vec3d previousPosition = m_previousMesh->getVertex( *previousMeshItr );
        BASim::Vec3d nextPosition = m_nextMesh->getVertex( *nextMeshItr );
        
        m_currentMesh->getVertex( *currentMeshItr ) = ( 1 - interpolation ) * previousPosition + interpolation * nextPosition;
        //m_currentMesh->getVertex( *currentMeshItr ) = ( 1 - interpolation ) * previousPosition + interpolation ;

        //m_currentMesh->getVertex( *currentMeshItr ) = m_previousMesh->getVertex( *previousMeshItr ) * ( 1 - interpolation ) +
          //  m_nextMesh->getVertex( *nextMeshItr ) * interpolation; 
        
        //m_currentMesh->getVertex( *currentMeshItr ) = m_nextMesh->getVertex( *nextMeshItr );
        
    }

    return true;
}
