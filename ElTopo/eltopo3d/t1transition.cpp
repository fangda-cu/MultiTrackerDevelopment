// ---------------------------------------------------------
//
//  t1transition.cpp
//  Fang Da 2013
//  
//  Functions handling T1 transitions (edge popping and vertex popping).
//
// ---------------------------------------------------------

#include <t1transition.h>
#include <broadphase.h>
#include <collisionqueries.h>
#include <runstats.h>
#include <subdivisionscheme.h>
#include <surftrack.h>
#include <trianglequality.h>


// ---------------------------------------------------------
//  Extern globals
// ---------------------------------------------------------

namespace ElTopo {

extern RunStats g_stats;

// ---------------------------------------------------------
// Member function definitions
// ---------------------------------------------------------

// ---------------------------------------------------------
///
/// Constructor.  Active SurfTrack object must be supplied.
///
// ---------------------------------------------------------

T1Transition::T1Transition(SurfTrack & surf, bool remesh_boundaries) :
    m_remesh_boundaries(remesh_boundaries),
    m_surf(surf)
{

}


// --------------------------------------------------------
///
/// Perform edge popping
///
// --------------------------------------------------------

bool T1Transition::pop_edges()
{
    if ( m_surf.m_verbose )
        std::cout << "---------------------- T1 Transition: Edge popping ----------------------" << std::endl;
    
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;
    bool pop_occurred = false;
    
    
    
    return pop_occurred;
}


// --------------------------------------------------------
///
/// Perform vertex popping
///
// --------------------------------------------------------

bool T1Transition::pop_vertices()
{
    if ( m_surf.m_verbose )
        std::cout << "---------------------- T1 Transition: Vertex popping ----------------------" << std::endl;
    
    NonDestructiveTriMesh & mesh = m_surf.m_mesh;
    bool pop_occurred = false;
    
    
    
    return pop_occurred;
}

}
