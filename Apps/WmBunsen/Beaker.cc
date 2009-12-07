/**
 * \file Beaker.cc
 *
 * \author acoull@wetafx.co.nz
 * \date 10/26/2009
 */

#include "Beaker.hh"
 
#include <BASim/BASim>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <sys/time.h>
#include <tclap/CmdLine.h>

using namespace BASim;
using namespace tr1;

Beaker::Beaker() : m_gravity( 0, -981.0, 0 )
{
    m_rodDataMap.clear();
    initialiseWorld();
    
    PetscTruth isInitialised;
    PetscInitialized( &isInitialised );
    if ( isInitialised != PETSC_TRUE )
    {
        cerr << "Initialising Petsc\n";
        int argc = 0; char **argv = NULL;
        PetscUtils::initializePetsc( &argc, &argv );
    }
    else
        cerr << "Skipping initalisation of Petsc\n";
    
    cerr << "Done with initialisation\n";
}

Beaker::~Beaker()
{
    cerr << "Finalising petsc\n";
///    PetscUtils::finalizePetsc();

    resetEverything();
}

void Beaker::initialiseWorld()
{
    m_world = new World();
    m_world->add_property( m_timeHandle, "time", 0.0 );
    m_world->add_property( m_dtHandle, "time-step", 0.01 );
    m_world->add_property( m_gravityHandle, "gravity", m_gravity );
    m_world->add_property( m_maxIterHandle, "maxIter", 100);
}

void Beaker::resetEverything()
{
    for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
    {
        vector<RodData*>& rodData = rdmItr->second;
        size_t numRods = rodData.size();
        for ( size_t r=0; r<numRods; r++ )
        {
            // We're safe to clear this vector as the individual destructors will safely delete the 
            // rod data in each vector element.
            delete rodData[ r ];
        }
    }
    m_rodDataMap.clear();
    
    delete m_world;
    
    initialiseWorld();
}

void Beaker::createSpaceForRods( size_t i_rodGroup, size_t i_numRods )
{
    //size_t numRods = m_rodDataMap[ i_rodGroup ].size();
    m_rodDataMap[ i_rodGroup ].resize( i_numRods );

    for ( size_t r=0; r<i_numRods; r++ )
    {
        m_rodDataMap[ i_rodGroup ][ r ] = new RodData();
    }
}

void Beaker::createRods( size_t i_rodGroup, ObjectControllerBase::SolverLibrary solverLibrary )
{
    vector<RodData*>& rodDataVector = m_rodDataMap[ i_rodGroup ];
    size_t numRods = rodDataVector.size();
    
    for ( size_t r=0; r<numRods; r++ )
    {
         // setupRod() is defined in ../BASim/src/Physics/ElasticRods/RodUtils.hh
        rodDataVector[ r ]->rod = setupRod( rodDataVector[ r ]->rodOptions, 
                                            rodDataVector[ r ]->initialVertexPositions, 
                                            rodDataVector[ r ]->undeformedVertexPositions );
    
        rodDataVector[ r ]->rod->fixVert( 0 );
        rodDataVector[ r ]->rod->fixVert( 1 );
        rodDataVector[ r ]->rod->fixEdge( 0 );
        
        rodDataVector[ r ]->stepper = setupRodTimeStepper( *(rodDataVector[ r ]->rod), solverLibrary );
        
        m_world->addObject( rodDataVector[ r ]->rod );
        m_world->addController( rodDataVector[ r ]->stepper );
        
        rodDataVector[ r ]->rodRenderer = new RodRenderer( *(rodDataVector[ r ]->rod) );
       // rodDataVector[ r ]->rodRenderer->setMode( RodRenderer::SIMPLE );
        
        rodDataVector[ r ]->shouldSimulate = true;
    }
    
    cerr << "Beaker - Created " << numRods << " rods\n";
}

void Beaker::takeTimeStep( int i_numberOfThreadsToUse, Scalar i_stepSize, 
  int i_subSteps, bool i_collisionsEnabled  )
{

  cerr << "i_collisionsEnabled = " << i_collisionsEnabled << endl;

    Scalar dt_save = getDt();
    Scalar startTime = getTime();
    Scalar currentTime = getTime();
    Scalar targetTime = currentTime + i_stepSize*i_subSteps;
    setDt( i_stepSize/i_subSteps );
    
    for ( int s=0; s<i_subSteps; s++ )
    {
        if ( (targetTime - currentTime) < getDt() + SMALL_NUMBER )
            setDt( targetTime - currentTime );

        // Update CollisionMeshData for this substep
        //
        Scalar interpolateFactor = ( (double)(s+1) / i_subSteps );
        if ( i_collisionsEnabled )
        {
            for ( CollisionMeshDataHashMapIterator cmItr = m_collisionMeshMap.begin();
                                                    cmItr != m_collisionMeshMap.end(); ++cmItr )
            {
                cmItr->second->interpolate( interpolateFactor ); 
            }
        }
        
        // interpolate fixed vertex positions and set timestep
        //
        for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
        {
            vector<RodData*>& rodData = rdmItr->second;
            size_t numRods = rodData.size();
            for ( size_t r=0; r<numRods; r++ )
            {
                dynamic_cast<RodTimeStepper*>(rodData[r]->stepper->getTimeStepper())->setTimeStep(getDt());
                BASim::ElasticRod* rod = rodData[r]->rod;
                for( int c = 0; c < rod->nv(); c++)
                {
                    if( rod->vertFixed( c ) )
                    {
                        rod->setVertex( c,interpolateFactor * rodData[r]->nextVertexPositions[c] + 
                               ( 1.0 - interpolateFactor ) * rodData[r]->prevVertexPositions[c] );
                    }
                }
            }
        }
        
#ifdef USING_INTEL_COMPILER
        m_world->executeInParallel( i_numberOfThreadsToUse );
#else

        if ( !i_collisionsEnabled )
        {
            m_world->execute();
        }
        else
        {
          // There is no way to pass in the collision meshes to world so I'm going to
          // iterate over its controllers myself.
          Controllers& controllers = m_world->getControllers();
          Controllers::iterator it;
          for (it = controllers.begin(); it != controllers.end(); ++it) 
          {
              cerr << "passing in dt for collisions of " << getDt() << endl;
              dynamic_cast<RodCollisionTimeStepper*>(*it)->execute(m_collisionMeshMap, getDt());
          }
        }
#endif
        setTime( currentTime + getDt() );
        currentTime = getTime();
    }
    
    // restore dt
    setDt( dt_save );
}

RodCollisionTimeStepper* Beaker::setupRodTimeStepper( BASim::ElasticRod& rod, 
    ObjectControllerBase::SolverLibrary solverLibrary )
{
    RodTimeStepper* stepper = new RodTimeStepper( rod );
    
    std::string integrator = "implicit";
    
    if (integrator == "symplectic") 
    {
        stepper->setDiffEqSolver( RodTimeStepper::SYMPL_EULER );
    } 
    else if (integrator == "implicit") 
    {
        stepper->setDiffEqSolver( RodTimeStepper::IMPL_EULER, solverLibrary );
    } 
    else 
    {
        std::cerr << "Unknown integrator. " << integrator
                  << "Using default instead." << std::endl;
    }
    
    stepper->setTimeStep(getDt());
    
    Scalar massDamping = 10.0;
    if (massDamping != 0) 
    {
        stepper->addExternalForce( new RodMassDamping( massDamping ) );
    }
    
    if (getGravity().norm() > 0) 
    {
        stepper->addExternalForce( new RodGravity( getGravity() ) );
    }
    
    // We use 1 iteration on purpose, it should work well. See 'Large steps in cloth simulation' paper
    //int iterations = 1;
    //stepper->setMaxIterations( iterations );

    RodCollisionTimeStepper* rodCollisionTimeStepper = new RodCollisionTimeStepper( stepper, &rod );
    
    return rodCollisionTimeStepper;
}

void Beaker::draw()
{
    for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
    {
        vector<RodData*>& rodData = rdmItr->second;
        size_t numRods = rodData.size();
        for ( size_t r=0; r<numRods; r++ )
        {
            rodData[ r ]->rodRenderer->render();
        }
    }

    for ( CollisionMeshDataHashMapIterator cmItr = m_collisionMeshMap.begin();
                                                   cmItr != m_collisionMeshMap.end(); ++cmItr )
    {
        cmItr->second->draw(); 
    }
}


bool Beaker::collisionMeshInitialised( const size_t id )
{
    BASim::CollisionMeshDataHashMapIterator itr = m_collisionMeshMap.find( id );
    if ( itr != m_collisionMeshMap.end() && itr->second )
        return itr->second->initialized();

    return false;
}

void Beaker::initialiseCollisionMesh( BASim::CollisionMeshData *collisionMeshData, size_t id )
{
    m_collisionMeshMap[ id ] = collisionMeshData;
    m_collisionMeshMap[ id ]->initialize();
}

void Beaker::removeCollisionMesh( const size_t id )
{
    std::cout << "Removing collision mesh with id " << id << std::endl;

    m_collisionMeshMap.erase(id);
}


