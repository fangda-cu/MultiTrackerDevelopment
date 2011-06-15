#ifndef WMFIGMESHCONTROLLER_HH_
#define WMFIGMESHCONTROLLER_HH_

#include <weta/Wfigaro/Core/TriangleMesh.hh>
#include <weta/Wfigaro/Core/ScriptingController.hh>
#include <weta/Wfigaro/Collisions/LevelSet.hh>

#include <ext/hash_map>

class WmFigMeshController : public BASim::ScriptingController
{
public:
    WmFigMeshController( BASim::TriangleMesh* i_currentMesh, 
                         BASim::TriangleMesh* i_previousMesh, 
                         BASim::TriangleMesh* i_nextMesh, 
                         double i_time, double i_dt );
    
    ~WmFigMeshController();
    
    bool execute();

    void setLevelSetCellSize( const float i_cellSize );
    void createLevelSet( const bool i_createLevelSet );    
    void drawLevelSet( const bool i_drawLevelSet );    

    void updateNextMayaTime( const double i_mayaTime );
    
    void setTriangleIndices( std::vector< unsigned int >& i_indices );
    
    void draw();
    
    BASim::LevelSet* currentLevelSet();
    
private:
    void calculateLevelSetSize( bridson::Vec3f &origin, Vec3ui &dims, Real &dx, Real length[3] );
    void buildLevelSet();
    
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
    
    // Level set data
    BASim::LevelSet *m_phiPrevious;
    BASim::LevelSet *m_phiCurrent;
    
    // Position and velocity of every vertex in the mesh
    std::vector<bridson::Vec3f> m_x;
    std::vector<bridson::Vec3f> m_v;
    
    Indices m_triIndices;
    Indices m_edgeIndices;
    
    Vec3Indices m_tri;
    
    float m_levelSetDX;
    bool m_createLevelSet;
    bool m_drawLevelSet;
};


#endif
