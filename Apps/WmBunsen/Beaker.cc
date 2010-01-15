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
#include <BASim/src/Physics/ElasticRods/RodHairsprayForce.hh>

using namespace BASim;
using namespace tr1;

Beaker::Beaker() : m_plasticDeformations( false ), m_gravity( 0, -981.0, 0 ), m_timing( false )
{
    m_rodDataMap.clear();
    initialiseWorld();
    
    cerr << "Done with initialisation\n";
}

Beaker::~Beaker()
{
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

void Beaker::startTimer()
{
    if ( m_timing )
        gettimeofday( &m_timerStart, NULL );
}

double Beaker::stopTimer()
{
    if (!m_timing)
        return 0.0;

    timeval end;
    gettimeofday( &end, NULL );

    return (((         end.tv_sec * 1000000 +          end.tv_usec) -
             (m_timerStart.tv_sec * 1000000 + m_timerStart.tv_usec)) / 1000000.0);
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
    m_rods.clear();

    for ( size_t r=0; r<numRods; r++ )
    {
         // setupRod() is defined in ../BASim/src/Physics/ElasticRods/RodUtils.hh
        rodDataVector[ r ]->rod = setupRod( rodDataVector[ r ]->rodOptions, 
                                            rodDataVector[ r ]->initialVertexPositions, 
                                            rodDataVector[ r ]->undeformedVertexPositions );

        m_rods.push_back( rodDataVector[ r ]->rod );
    
        rodDataVector[ r ]->rod->fixVert( 0 );
        rodDataVector[ r ]->rod->fixVert( 1 );
        rodDataVector[ r ]->rod->fixEdge( 0 );
        
        rodDataVector[ r ]->stepper = setupRodTimeStepper( rodDataVector[ r ] );
        
        m_world->addObject( rodDataVector[ r ]->rod );
        m_world->addController( rodDataVector[ r ]->stepper );
        
        rodDataVector[ r ]->rodRenderer = new RodRenderer( *(rodDataVector[ r ]->rod) );
       // rodDataVector[ r ]->rodRenderer->setMode( RodRenderer::SIMPLE );
        
        rodDataVector[ r ]->shouldSimulate = true;
    }
    
    //cerr << "Beaker - Created " << numRods << " rods\n";
}

void Beaker::takeTimeStep( int i_numberOfThreadsToUse, Scalar i_stepSize, 
  int i_subSteps, bool i_collisionsEnabled,  bool i_selfCollisionPenaltyForcesEnabled,
  bool i_fullSelfCollisionsEnabled, int i_fullSelfCollisionIters,
  double i_selfCollisionCOR )
{
    // Check if anything has actually been initialised yet. We may still be being loaded by Maya.
    if ( m_rodDataMap.size() == 0 )
        return;

    Scalar dt_save = getDt();
    Scalar startTime = getTime();
    Scalar currentTime = getTime();
    Scalar targetTime = currentTime + i_stepSize*i_subSteps;
    setDt( i_stepSize/i_subSteps );
    
    double meshInterpolationTime = 0.0;
    double vertexInterpolationTime = 0.0;
    double executePart1Time = 0.0;
    double executePart2Time = 0.0;
    double executePart3Time = 0.0;
    double selfCollisionPenaltyForceTime = 0.0;
    double selfCollisionsRodResponseTime = 0.0;
    double integrationStepTime = 0.0;
    
    m_timing = true;
    
    if ( m_timing )
        cerr << "Taking step using " << i_subSteps << " substeps (dt=" << getDt() << ")\n";
    
    for ( int s=0; s<i_subSteps; s++ )
    {
        if ( (targetTime - currentTime) < getDt() + SMALL_NUMBER )
            setDt( targetTime - currentTime );

        // Update CollisionMeshData for this substep
        //
        startTimer();
        Scalar interpolateFactor = ( (double)(s+1) / i_subSteps );
        if ( i_collisionsEnabled )
        {
            for ( CollisionMeshDataHashMapIterator cmItr = m_collisionMeshMap.begin();
                                                    cmItr != m_collisionMeshMap.end(); ++cmItr )
            {
                cmItr->second->interpolate( interpolateFactor ); 
            }
        }
        meshInterpolationTime += stopTimer();
        
        // interpolate fixed vertex positions and set timestep
        //
        startTimer();
        for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
        {
            // FIXME: Should check if rod is enabled, if not then continue
            
            vector<RodData*>& rodData = rdmItr->second;
            size_t numRods = rodData.size();
            for ( size_t r=0; r<numRods; r++ )
            {
                dynamic_cast<RodTimeStepper*>(rodData[r]->stepper->getTimeStepper())->setTimeStep(getDt());
                RodBoundaryCondition* boundary = dynamic_cast<RodTimeStepper*>(rodData[r]->stepper->getTimeStepper())->getBoundaryCondition();
                BASim::ElasticRod* rod = rodData[r]->rod;

                
                for( int c = 0; c < rod->nv(); c++)
                {
                    // Calculate the position of the input curve for this rod at the current substep.
                    // This is used to set fixed vertices and also to calculate the hairspray force.
                    rodData[r]->currVertexPositions[c] = ( interpolateFactor * rodData[r]->nextVertexPositions[c] + 
                               ( 1.0 - interpolateFactor ) * rodData[r]->prevVertexPositions[c] );

                    // FIXME: THIS IS INSANE, This MUST go in RodCollisionTimeStepper but there are include issues there and I can't
                    // includ RodTimeStepper. FIX THESE AND MOVE THIS
                    
                    // First enforce the boundary conditions. Setting a vertex to fixed on the rod does
                    // not enforce anything any more. The fixed vertices must be explicitly added in here/                    
                    if ( rod->vertFixed(c) )
                    {  
                        boundary->setDesiredVertexPosition(c, rodData[r]->currVertexPositions[c]);
                        if (c>0)
                          if (rod->vertFixed(c-1))  
                            boundary->setDesiredEdgeAngle(c-1, rod->getTheta(c-1));
                    }
                }
                
                //
                // This is a bit weird, I'm setting stuff in the RodCollisionTimeStepper as well
                // as the RodTimeStepper it contains. It would work better as class based on
                // RodTimeStepper or if setting it in RodCollisionTimeStepper also set it in
                // the contained RodTimeStepper. I need to fix some #include conflicts before
                // I can do that.
                dynamic_cast<RodCollisionTimeStepper*>(rodData[r]->stepper)->doCollisions( i_collisionsEnabled );
                dynamic_cast<RodCollisionTimeStepper*>(rodData[r]->stepper)->setTimeStep( getDt() );
                dynamic_cast<RodCollisionTimeStepper*>(rodData[r]->stepper)->setCollisionMeshesMap( &m_collisionMeshMap );
        
            }
        }
        vertexInterpolationTime += stopTimer();

        startTimer();
        //m_world->executeInParallel( i_numberOfThreadsToUse );
        // Self collisions have to be done on all rods together so we split the time step into two steps.
        Controllers controllers = m_world->getControllers();
        int numControllers = (int)controllers.size();
    
        #pragma omp parallel for num_threads( i_numberOfThreadsToUse )
        for ( int i=0; i<numControllers; ++i )
        {
            int threadID = omp_get_thread_num();
            // only master thread prints the number of threads
            if ( threadID == 0 )
            {
                int nthreads = omp_get_num_threads();
                // printf("World::executeInParallel, using %d threads\n", nthreads );
            }
            // Dangerous cast...
            dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ])->executePart1();
        }
        executePart1Time += stopTimer();

        startTimer();
        // FIXME: Think this can be parallelised
        if(i_selfCollisionPenaltyForcesEnabled)
        {
            RodCollisionTimeStepper::getProximities(m_rods);
        }
        selfCollisionPenaltyForceTime += stopTimer();
  
        startTimer();
        #pragma omp parallel for num_threads( i_numberOfThreadsToUse )
        for ( int i=0; i<numControllers; ++i )
        {
            dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ])->execute();
        }
        integrationStepTime += stopTimer();
  
        //////////////////////////////////////////////
        //
        // Check forces to see if any rods need to be simulated at a slower pace
       // checkAllRodForces();

       startTimer();
        #pragma omp parallel for num_threads( i_numberOfThreadsToUse )
        for ( int i=0; i<numControllers; ++i )
        {
            dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ])->executePart2();
        }
        executePart2Time += stopTimer();

        startTimer();
        if (i_fullSelfCollisionsEnabled)
           RodCollisionTimeStepper::respondRodCollisions( m_rods, getDt(), i_fullSelfCollisionIters,
                                                          i_selfCollisionCOR );
        selfCollisionsRodResponseTime += stopTimer();
        
        startTimer();
        #pragma omp parallel for num_threads( i_numberOfThreadsToUse )
        for ( int i=0; i<numControllers; ++i )
        {
            dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ])->executePart3();
        }
        executePart3Time += stopTimer();


        // This sets the undeformed rod configuration to be the same as the current configuration. 
        // I'm not sure if this is actually what we want to do, is this just like pushing around line
        // segments. Maybe we should only run this on rods that have collided.

        /*if ( m_plasticDeformations )
        {
            for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
            {
                vector<RodData*>& rodData = rdmItr->second;
                size_t numRods = rodData.size();
                for ( size_t r=0; r<numRods; r++ )
                {
                    BASim::ElasticRod* rod = rodData[r]->rod;                         
                    rod->updateReferenceProperties();
                }
            }
        }*/

        setTime( currentTime + getDt() );
        currentTime = getTime();
    }
    
    cerr << "Finished all substeps, timings:\n";
    cerr << "meshInterpolationTime = " << meshInterpolationTime << endl;
    cerr << "vertexInterpolationTime = " << vertexInterpolationTime << endl;
    cerr << "executePart1Time = " << executePart1Time << endl;
    cerr << "executePart2Time = " << executePart2Time << endl;
    cerr << "executePart3Time = " << executePart3Time << endl;
    cerr << "selfCollisionPenaltyForceTime = " << selfCollisionPenaltyForceTime << endl;
    cerr << "selfCollisionsRodResponseTime = " << selfCollisionsRodResponseTime << endl;
    cerr << "integrationStepTime = " << integrationStepTime << endl;
    
    
    // restore dt
    setDt( dt_save );
}

void Beaker::checkAllRodForces()
{
    for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
    {
        vector<RodData*>& rodData = rdmItr->second;
        size_t numRods = rodData.size();
        for ( size_t r=0; r<numRods; r++ )
        {
            VecXd& forces = dynamic_cast<RodTimeStepper*>(rodData[r]->stepper->getTimeStepper())->getForcesAtLastStep();
            cerr << "forces at last step = \n" << forces << endl;
        }
    }
}

RodCollisionTimeStepper* Beaker::setupRodTimeStepper( RodData* i_rodData )
{

    ElasticRod& rod = (*i_rodData->rod);
    RodTimeStepper* stepper = new RodTimeStepper( rod );
    
    std::string integrator = "implicit";
    
    if (integrator == "symplectic") 
    {
        stepper->setDiffEqSolver( RodTimeStepper::SYMPL_EULER );
    } 
    else if (integrator == "implicit") 
    {
        //stepper->setDiffEqSolver( RodTimeStepper::IMPL_EULER, solverLibrary );
        stepper->setDiffEqSolver( RodTimeStepper::IMPL_EULER );
    } 
    else 
    {
        std::cerr << "Unknown integrator. " << integrator
                  << "Using default instead." << std::endl;
    }
    
    stepper->setTimeStep(getDt());
    
    Scalar massDamping = i_rodData->massDamping;
    if (massDamping != 0) 
    {
        stepper->addExternalForce( new RodMassDamping( massDamping ) );
    }
    
    if (getGravity().norm() > 0) 
    {
        stepper->addExternalForce( new RodGravity( getGravity() ) );
    }
  
    // Add the hairspray force to attract the rods back to their input curves    
    WsplineAttr  splineArr;
    splineArr.set( i_rodData->forceWeightMap );

    vector<Scalar> ks( rod.nv() );
    for ( int i=0; i<rod.nv(); ++i )
    {
        // Caclulate the hairspray force for this CV
        // FIXME: this should be how far along the rod the cv is rather than just cv/nCVs
        // cvs may not be evenly spaced (although that will make them more stable).
        double t = double( i ) / double( rod.nv() - 1  );
        ks[ i ] = splineArr.getValue( t );
        ks[ i ] *= i_rodData->hairSprayScaleFactor;
 
        // If the force is 100% then there is no point in using a force, just lock the 
        // vertex in place and it won't move.
        if ( ks[ i ] >= 1.0 )
        {
            rod.fixVert( i );
        }
    }
    
    //stepper->addExternalForce( new RodHairsprayForce( rod, ks, i_rodData->currVertexPositions ) );  
    
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

  /*  for ( CollisionMeshDataHashMapIterator cmItr = m_collisionMeshMap.begin();
                                                   cmItr != m_collisionMeshMap.end(); ++cmItr )
    {
        cmItr->second->draw(); 
    }*/

#if 0
    for ( RodDataMapIterator rdmItr = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
    {
        vector<RodData*>& rodData = rdmItr->second;
        size_t numRods = rodData.size();
        for ( size_t r=0; r<numRods; r++ )
        {
            BASim::ElasticRod* rod = rodData[r]->rod;
            if ( rodData[r]->ALLprevVertexPositions.size() != 0 )
            {       
                glColor3d(0,0,1);
                glBegin( GL_LINE_STRIP ); 
                size_t total = rodData[r]->ALLprevVertexPositions.size();
                for ( int c=0; c<total; c++ )
                {
                    if ( c>0 && c%rod->nv()==0 )
                    {
                        glEnd();
                        glBegin( GL_LINE_STRIP );
                    }
                    glVertex3d( rodData[r]->ALLprevVertexPositions[c][0],  
                                rodData[r]->ALLprevVertexPositions[c][1],  
                                rodData[r]->ALLprevVertexPositions[c][2] );
                }
                glEnd();
            }
            if ( rodData[r]->ALLnextVertexPositions.size() != 0 )
            {       
                glColor3d(1,0,0);
                glBegin( GL_LINE_STRIP ); 
                size_t total = rodData[r]->ALLnextVertexPositions.size();
                for ( int c=0; c<total; c++ )
                {
                    if ( c>0 && c%rod->nv()==0 )
                    {
                        glEnd();
                        glBegin( GL_LINE_STRIP );
                    }
                    glVertex3d( rodData[r]->ALLnextVertexPositions[c][0],  
                                rodData[r]->ALLnextVertexPositions[c][1],  
                                rodData[r]->ALLnextVertexPositions[c][2] );
                }
                glEnd();
            }
            if ( rodData[r]->ALLcurrVertexPositions.size() != 0 )
            {       
                glColor3d(0,1,0);
                glBegin( GL_LINE_STRIP ); 
                size_t total = rodData[r]->ALLcurrVertexPositions.size();
                for ( int c=0; c<total; c++ )
                {
                    if ( c>0 && c%rod->nv()==0 )
                    {
                        glEnd();
                        glBegin( GL_LINE_STRIP );
                    }
                    glVertex3d( rodData[r]->ALLcurrVertexPositions[c][0],  
                                rodData[r]->ALLcurrVertexPositions[c][1],  
                                rodData[r]->ALLcurrVertexPositions[c][2] );
                }
                glEnd();
            }
        }
    }
#endif
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


