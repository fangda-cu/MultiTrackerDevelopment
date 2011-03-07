#ifndef WMFIGMESHCONTROLLER_HH_
#define WMFIGMESHCONTROLLER_HH_

#include <weta/Wfigaro/Core/TriangleMesh.hh>
#include <weta/Wfigaro/Core/ScriptingController.hh>


class WmFigMeshController : public BASim::ScriptingController
{
public:
    WmFigMeshController( BASim::TriangleMesh* i_currentMesh, 
                         BASim::TriangleMesh* i_previousMesh, 
                         BASim::TriangleMesh* i_nextMesh, 
                         double i_time, double i_dt );
    
    bool execute();

    void updateNextMayaTime( const double i_mayaTime )
    {
        m_previousMayaTime = m_nextMayaTime;
        m_nextMayaTime = i_mayaTime;
    }

private:
    
    // We have two meshes because Maya is only providing meshes on frame steps and the sim is 
    // stepping with smaller steps so we need to interpolate the mesh positions and 
    BASim::TriangleMesh* m_currentMesh;
    BASim::TriangleMesh* m_previousMesh;
    BASim::TriangleMesh* m_nextMesh;

    double m_startMeshTime;
    double m_endMeshTime;

    double m_startTime;

    double m_nextMayaTime;
    double m_previousMayaTime;
};


#endif