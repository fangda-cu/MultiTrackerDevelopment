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

Beaker::Beaker()
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
    m_world->add_property( m_time, "time", 0.0 );
    m_world->add_property( m_dt, "time-step", 0.1 );
    m_world->add_property( m_gravity, "gravity", Vec3d(0, -9.81, 0) );
}

void Beaker::resetEverything()
{
    cerr << "Beaker::resetEverything()\n";
    
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

void Beaker::createRods( size_t i_rodGroup )
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
        //rod->fixVert( rod->nv() - 2 );
        //rod->fixVert( rod->nv() - 1 );
        rodDataVector[ r ]->rod->fixEdge( 0 );
        //rod->fixEdge( rod->ne() - 1 );
        
        rodDataVector[ r ]->stepper = setupRodTimeStepper( *(rodDataVector[ r ]->rod) );
        
        m_world->addObject( rodDataVector[ r ]->rod );
        m_world->addController( rodDataVector[ r ]->stepper );
        
        rodDataVector[ r ]->rodRenderer = new RodRenderer( *(rodDataVector[ r ]->rod) );
    }
}

void Beaker::addRod( size_t i_rodGroup,
                     vector<Vec3d>& i_initialVertexPositions, 
                     vector<Vec3d>& i_undeformedVertexPositions,
                     RodOptions& i_options )
{
    // setupRod() is defined in ../BASim/src/Physics/ElasticRods/RodUtils.hh
    ElasticRod* rod = setupRod( i_options, i_initialVertexPositions, i_undeformedVertexPositions );

    rod->fixVert( 0 );
    rod->fixVert( 1 );
    rod->fixVert( rod->nv() - 2 );
    rod->fixVert( rod->nv() - 1 );
    rod->fixEdge( 0 );
    rod->fixEdge( rod->ne() - 1 );
    
    RodTimeStepper* stepper = setupRodTimeStepper( *rod );
    
    // FIXME: 
    // Why do we have to add rods and steppers seperately? It feels like the stepper should
    // be part of the rod?
    m_world->addObject( rod );
    m_world->addController( stepper );
    
    RodRenderer* rodRenderer = new RodRenderer( *rod );
    
    RodData* rodData = new RodData( rod, stepper, rodRenderer );
    m_rodDataMap[ i_rodGroup ].push_back( rodData );
}

void Beaker::takeTimeStep()
{
    m_world->execute();

    setTime( getTime() + getDt() );
}

RodTimeStepper* Beaker::setupRodTimeStepper( ElasticRod& rod )
{
    RodTimeStepper* stepper = new RodTimeStepper( rod );
    
    std::string integrator = "implicit";
    
    if (integrator == "symplectic") 
    {
        stepper->setDiffEqSolver( RodTimeStepper::SYMPL_EULER );
    } 
    else if (integrator == "implicit") 
    {
        stepper->setDiffEqSolver( RodTimeStepper::IMPL_EULER );
    } 
    else 
    {
        std::cerr << "Unknown integrator. " << integrator
                  << "Using default instead." << std::endl;
    }
    
    stepper->setTimeStep(0.1);
    
    Scalar massDamping = 0.0;
    if (massDamping != 0) 
    {
        stepper->addExternalForce( new RodMassDamping( massDamping ) );
    }
    
    if (getGravity().norm() > 0) 
    {
        stepper->addExternalForce( new RodGravity( getGravity() ) );
    }
    
    int iterations = 1;
    stepper->setMaxIterations( iterations );
    
    return stepper;
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
}

