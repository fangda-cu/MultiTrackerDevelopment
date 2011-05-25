#include "WmFigMeshController.hh"

using namespace BASim;
using namespace std;

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
//       as BARodStepper calls setDt() on the controller when it executes anyway.
//       We could probably just not bother getting them passed in here but we keep it this way
//       in case Breannan changes the way BARodStepper works so we'll just copy his examples.

{   
    m_startTime = i_time;
    m_previousMayaTime = m_nextMayaTime = m_startTime;

    m_phiPrevious = new LevelSet;
    m_phiCurrent = new LevelSet;
}
 
bool WmFigMeshController::execute()
{
    // TODO: Assumes 24fps for Maya scenes, fix that
 
    buildLevelSet();
 
    double interpolation = ( getTime() - m_startTime ) * 24.0 + m_startTime - m_previousMayaTime;
    //  std::cerr << "interpolated mesh factor = " << interpolation << std::endl;
  
    for( TriangleMesh::vertex_iter previousMeshItr = m_previousMesh->vertices_begin(),
         nextMeshItr = m_nextMesh->vertices_begin(),
         currentMeshItr = m_currentMesh->vertices_begin(); 
         previousMeshItr != m_previousMesh->vertices_end(); 
         ++previousMeshItr, ++nextMeshItr, ++currentMeshItr )
    {
        BASim::Vec3d previousPosition = m_previousMesh->getVertex( *previousMeshItr );
        BASim::Vec3d nextPosition = m_nextMesh->getVertex( *nextMeshItr );
        
        m_currentMesh->getVertex( *currentMeshItr ) = ( 1 - interpolation ) * previousPosition + interpolation * nextPosition;
    }

    return true;
}

void WmFigMeshController::setTriangleIndices( std::vector< unsigned int >& i_indices )
{
    // Setup the data we need to build the level set
    m_triIndices.resize( m_currentMesh->nf() );

    // A copy of triangleIndices, but instead of a flat array
    // it is stored as a vector of Vec3ui, as needed by the level set code
    //
    m_tri.resize( m_currentMesh->nf() );

    for ( unsigned int i = 0; i < m_currentMesh->nf() ; ++i )
    {
        /*Vec3<Real> v1 = triangleVertex( 3 * i );
        Vec3<Real> v2 = triangleVertex( 3 * i + 1 );
        Vec3<Real> v3 = triangleVertex( 3 * i + 2 );

        centroids[i] = (v1 + v2 + v3) / 3.0; */
        m_triIndices[ i ] = i;

        m_tri[ i ][ 0 ] = i_indices[ 3 * i ];
        m_tri[ i ][ 1 ] = i_indices[ 3 * i + 1 ];
        m_tri[ i ][ 2 ] = i_indices[ 3 * i + 2 ];
    }	
}
/*
Vec3<Real>& WmFigMeshController::triangleVertex( uint faceID, uint i )
{
    return vertex( triangleIndices[ 3 * faceID + i ] );
}*/

void WmFigMeshController::buildLevelSet()
{   
    if ( m_tri.size() != m_currentMesh->nf() )
    {
        cerr << "WmFigMeshController::buildLevelSet() - triangle indices are not set, failed to build level set\n";
        return;
    }
    
    cerr << "WmFigMeshController::buildLevelSet() - building level set ...";
    
    // Make sure we have space for the per vertex data
    m_x.resize( m_currentMesh->nv() );
    m_v.resize( m_currentMesh->nv() );
    
    // First make the previous mesh equal the next mesh as we're about to update next mesh
    size_t vertexIndex = 0;
    for( TriangleMesh::vertex_iter meshItr = m_currentMesh->vertices_begin();
         meshItr != m_currentMesh->vertices_end(); ++meshItr )
    {
        BASim::Vec3d vertex = m_currentMesh->getVertex( *meshItr );
        
        for ( unsigned int c = 0; c < 3; ++c )
        {        
            m_x[ vertexIndex ][ c ] = vertex[ c ];
            
            // For now, assume a static mesh
            m_v[ vertexIndex ][ c ] = 0;
        }
        ++vertexIndex;
    }

    bridson::Vec3f origin;
    
    Real dx;
    Vec3ui dims;
    Real length[3];
    calculateLevelSetSize( origin, dims, dx, length );
    
    if ( m_phiCurrent->isInitialized() )
    {
        delete m_phiPrevious;
        m_phiPrevious = m_phiCurrent;
        m_phiCurrent = new LevelSet;
    }

    m_phiCurrent->buildLevelSet( m_tri, m_triIndices, m_x, m_v, origin, length, dx, dims[0], 
                                 dims[1], dims[2], m_currentMesh->nf() );
    
    cerr << "WmFigMeshController::buildLevelSet() - Complete!" << endl;
}

void WmFigMeshController::calculateLevelSetSize( bridson::Vec3f &origin, Vec3ui &dims, Real &dx, Real length[3] )
{
    // Work out where the centre of the level set should be, the dimensions and the cell size
    
    Vec3d xMin, xMax, dX;
    
    for( TriangleMesh::vertex_iter meshItr = m_currentMesh->vertices_begin();
         meshItr != m_currentMesh->vertices_end(); ++meshItr )
    {
        BASim::Vec3d vertex = m_currentMesh->getVertex( *meshItr );
        
        for ( unsigned int i = 0; i < 3; ++i )
        {
            if ( vertex[ i ] < xMin[ i ] )
            {
                xMin[ i ] = vertex[ i ];
            }
            if ( vertex[ i ] > xMax[ i ] )
            {
                xMax[ i ] = vertex[ i ];
            }
        }
    }
    
    // Hard code the cell size for now...
    dx = 1.0;
    
    // For now just pad it out a bit
    for ( unsigned int i = 0; i < 3; ++i )
    {
        xMin[ i ] -= 2 * dx;
        xMax[ i ] += 2 * dx;
    }
    
    for( uint i = 0; i < 3; ++i )
    {
        origin[ i ] = xMin[ i ];
        length[ i ] = xMax[ i ] - xMin[ i ];
    }
    
    
    for( uint i = 0; i < 3; i++ )
    {
        dims[ i ] = ( int )ceil( length[ i ] / dx );
    }
}

void WmFigMeshController::draw()
{
    if( m_phiCurrent->isInitialized() )
    {
        m_phiCurrent->draw();
    }
}
