/**
 * \file Beaker.cc
 *
 * \author acoull@wetafx.co.nz
 * \date 10/26/2009
 */

#include "Beaker.hh"

#ifdef WETA
#include <weta/Wfigaro/Physics/ElasticRods/RodHairsprayForce.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodMassDamping.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodGravity.hh>
#else
#include <BASim/src/Physics/ElasticRods/RodHairsprayForce.hh>
#include <BASim/src/Physics/ElasticRods/RodMassDamping.hh>
#include <BASim/src/Physics/ElasticRods/RodGravity.hh>
#endif

#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <sys/time.h>

#include <boost/format.hpp>
#include <boost/lexical_cast.hpp>

using namespace BASim;
using namespace tr1;
using namespace boost;

Beaker::Beaker() : m_plasticDeformations( false ), m_gravity( 0, -981.0, 0 ), 
    m_timingEnabled( false ), m_timingsFile( "" ),
    m_meshInterpolationTime( 0.0 ), m_vertexInterpolationTime( 0.0 ), m_objectCollisionForces( 0.0 ),
    m_objectCollisionResponse( 0.0 ), m_collisionStructuresTidyup( 0.0 ),
    m_selfCollisionPenaltyForceTime( 0.0 ), m_selfCollisionsResponseTime( 0.0 ),
    m_integrationStepTime( 0.0 ), m_slowestIntegrationTime( 0.0 ), m_fastestIntegrationTime( 99999999999999.0 ),
    m_fastestSelfCollisionPenaltyForceTime( 99999999999999.0 ), m_slowestSelfCollisionPenaltyForceTime( 0.0 ),
    m_fastestSelfCollisionsResponseTime( 999999999999.0 ), m_slowestSelfCollisionsResponseTime( 0.0 ),
    m_slowestCollisionForcesTime( 0.0 ), m_fastestCollisionForcesTime( 999999999.0 ),
    m_slowestCollisionResponseTime( 0.0 ), m_fastestCollisionResponseTime( 9999999999999.0 ),
    m_fastestFrameTime( 9999999999999999.0 ), m_slowestFrameTime( 0.0 ), m_totalSimTime( 0.0 ), 
    m_numberOfFramesSimulated( 0 ), m_numberofThreadsUsed( 0 ), m_numRods( 0 ),
    m_shouldDrawSubsteppedVertices( false ), m_isClumpingEnabled( false ), m_clumpingCoefficient( 0.3 ),
    m_isXMLLoggingEnabled( false ), m_sceneXML( NULL )
{
    m_rodDataMap.clear();
    m_initialRodConfigurations.clear();

    initialiseWorld();
}

Beaker::~Beaker()
{
    resetEverything();

    if ( m_timingsFP.is_open() )
        m_timingsFP.close();

    if ( m_sceneXML != NULL )
        delete m_sceneXML;
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
    /*for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
    {
        vector<RodData*>& rodData = rdmItr->second;
        size_t numRods = rodData.size();
        for ( size_t r=0; r<numRods; r++ )
        {
            // We're safe to clear this vector as the individual destructors will safely delete the
            // rod data in each vector element.
            delete rodData[ r ];
        }
    }*/
    m_rodDataMap.clear();

    m_rods.clear();

    delete m_world;

    initialiseWorld();
}

void Beaker::setTimingsFile( std::string i_fileName )
{
    m_timingsFile = i_fileName;
}

void Beaker::setTimingEnabled( bool i_timingsEnabled )
{
    if ( i_timingsEnabled == true && m_timingEnabled == false )
    {
        m_timingEnabled = i_timingsEnabled;
        resetTimers();

        // Open timings file for timer code to output to
        if ( m_timingsFP.is_open() )
            m_timingsFP.close();

        if ( m_timingsFile != "" )
        {
            m_timingsFP.open( m_timingsFile.c_str(), ios::out | ios::trunc );
        }
    }
    else if ( i_timingsEnabled == false && m_timingEnabled == true )
    {
        printTimingInfo();

        // Close timings file
        if ( m_timingsFP.is_open() )
            m_timingsFP.close();
    }
    else
    {
        // Do nothing as the state hasn't changed
    }

    m_timingEnabled = i_timingsEnabled;
}

void Beaker::resetTimers()
{
    cerr << "Resetting timers \n";

    m_meshInterpolationTime = 0.0;
    m_vertexInterpolationTime = 0.0;
    m_objectCollisionForces = 0.0;
    m_objectCollisionResponse = 0.0;
    m_collisionStructuresTidyup = 0.0;
    m_selfCollisionPenaltyForceTime = 0.0;
    m_selfCollisionsResponseTime = 0.0;
    m_fastestSelfCollisionPenaltyForceTime = 0.0;
    m_slowestSelfCollisionPenaltyForceTime = 0.0;
    m_fastestSelfCollisionsResponseTime = 0.0;
    m_slowestSelfCollisionsResponseTime = 0.0;
    m_integrationStepTime = 0.0;
    m_slowestIntegrationTime = 0.0;
    m_fastestIntegrationTime = 0.0;
    m_slowestCollisionForcesTime = 0.0;
    m_fastestCollisionForcesTime = 0.0;
    m_slowestCollisionResponseTime = 0.0;
    m_fastestCollisionResponseTime = 0.0;
    m_fastestFrameTime = 99999999999.999;
    m_slowestFrameTime = 0.0;
    m_numberOfFramesSimulated = 0;
    m_totalSimTime = 0.0;
}

std::string Beaker::makeString( double i_val )
{
    return str( boost::format( "%.2f" ) % i_val );
}

void Beaker::printTimingInfo()
{
    if ( m_numberOfFramesSimulated == 0 )
        return;

    // Calculate the timing results and output them to stderr and file.

    double totalSimTime = m_totalSimTime;

    /*double totalSimTime =  m_meshInterpolationTime + m_vertexInterpolationTime +
                                m_objectCollisionForces + m_objectCollisionResponse +
                                m_collisionStructuresTidyup + m_selfCollisionPenaltyForceTime +
                                m_selfCollisionsResponseTime + m_integrationStepTime;*/
    double averageFrameTime = totalSimTime / m_numberOfFramesSimulated;
    double averageIntgeration = m_integrationStepTime / m_numberOfFramesSimulated;
    double averageObjectCollisionForces = m_objectCollisionForces / m_numberOfFramesSimulated;
    double averageObjectCollisionResponse = m_objectCollisionResponse / m_numberOfFramesSimulated;
    double averageSelfCollisionPenalty = m_selfCollisionPenaltyForceTime / m_numberOfFramesSimulated;
    double averageSelfCollisionsResponseTime = m_selfCollisionsResponseTime / m_numberOfFramesSimulated;

    string resultString = "\nSummary of simulation (in seconds)\n============================\n";
    resultString += "Number of rods = " + makeString(m_numRods)+"\n";
    resultString += "Number of threads used = " + makeString(m_numberofThreadsUsed)+"\n";
    resultString += "Total frames simulated = " + makeString(m_numberOfFramesSimulated) + "\n";
    resultString += "Total simulation time = " + makeString(totalSimTime) + "\n";
    resultString += "Fastest frame = " + makeString(m_fastestFrameTime) + "\n";
    resultString += "Slowest frame = " + makeString(m_slowestFrameTime) + "\n";
    resultString += "Average frame = " + makeString(averageFrameTime)  + "\n";
    resultString += "\nDetailed breakdown (Total / Average):\n";
    resultString += "Mesh Interpolation (SINGLE THREADED) = " + makeString(m_meshInterpolationTime)
        + " / " + makeString(m_meshInterpolationTime / m_numberOfFramesSimulated) +"\n";
    resultString += "Vertex Interpolation (SINGLE THREADED) = " + makeString(m_vertexInterpolationTime) + " / " +
        makeString(m_vertexInterpolationTime / m_numberOfFramesSimulated) +"\n";
    resultString += "Object Collision Forces = " + makeString(m_objectCollisionForces) + " / "
        + makeString(m_objectCollisionForces / m_numberOfFramesSimulated) +"\n";
    resultString += "\tSlowest Frame = " + makeString(m_slowestCollisionForcesTime) + "\n";
    resultString += "\tFastest Frame = " + makeString(m_fastestCollisionForcesTime) + "\n";
    resultString += "Object Collision Response = " + makeString(m_objectCollisionResponse) + " / "
        + makeString(m_objectCollisionResponse / m_numberOfFramesSimulated) +"\n";
    resultString += "\tSlowest Frame = " + makeString(m_slowestCollisionResponseTime) + "\n";
    resultString += "\tFastest Frame = " + makeString(m_fastestCollisionResponseTime) + "\n";
    resultString += "Collision Structures Tidy up = " + makeString(m_collisionStructuresTidyup)
        + " / " + makeString(m_collisionStructuresTidyup / m_numberOfFramesSimulated) +"\n";
    resultString += "Self Collision Penalty Force (SINGLE THREADED) = " + makeString(m_selfCollisionPenaltyForceTime) + " / "
        + makeString(m_selfCollisionPenaltyForceTime / m_numberOfFramesSimulated) +"\n";
    resultString += "\tSlowest Frame = " + makeString(m_slowestSelfCollisionPenaltyForceTime) + "\n";
    resultString += "\tFastest Frame = " + makeString(m_fastestSelfCollisionPenaltyForceTime) + "\n";
    resultString += "Self Collision Response (SINGLE_THREADED) = " + makeString(m_selfCollisionsResponseTime) +
        " / " + makeString(m_selfCollisionsResponseTime / m_numberOfFramesSimulated) +"\n";
    resultString += "\tSlowest Frame = " + makeString(m_slowestSelfCollisionsResponseTime) + "\n";
    resultString += "\tFastest Frame = " + makeString(m_fastestSelfCollisionsResponseTime) + "\n";
    resultString += "Integration = " + makeString(m_integrationStepTime) + " / "
        + makeString(m_integrationStepTime / m_numberOfFramesSimulated) +"\n";
    resultString += "\tSlowest Frame = " + makeString(m_slowestIntegrationTime) + "\n";
    resultString += "\tFastest Frame = " + makeString(m_fastestIntegrationTime) + "\n";

    cerr << resultString;

    if ( m_timingsFP.is_open() )
        m_timingsFP << resultString;
}

void Beaker::startTimer( timeval& i_startTimer )
{
    if ( m_timingEnabled )
        gettimeofday( &i_startTimer, NULL );
        //gettimeofday( &m_timerStart, NULL );
}

double Beaker::stopTimer( timeval& i_startTimer )
{
    if ( !m_timingEnabled )
        return 0.0;

    timeval end;
    gettimeofday( &end, NULL );

    return (((         end.tv_sec * 1000000 +          end.tv_usec) -
             ( i_startTimer.tv_sec * 1000000 + i_startTimer.tv_usec)) / 1000000.0);
}


/*void Beaker::createSpaceForRods( size_t i_rodGroup, size_t i_numRods )
{
    //size_t numRods = m_rodDataMap[ i_rodGroup ].size();
    m_rodDataMap[ i_rodGroup ].resize( i_numRods );

    for ( size_t r=0; r<i_numRods; r++ )
    {
        m_rodDataMap[ i_rodGroup ][ r ] = new RodData();
    }
}*/

void Beaker::addRodsToWorld( size_t i_rodGroupIndex, WmFigRodGroup* i_rodGroup )
{
    m_rodDataMap[ i_rodGroupIndex ] = i_rodGroup;
    
    size_t numRods = m_rodDataMap[ i_rodGroupIndex ]->numberOfRods();
    
    m_initialRodConfigurations.clear();

    for ( size_t r=0; r<numRods; r++ )
    {
        if ( !m_rodDataMap[ i_rodGroupIndex ]->shouldSimulateRod( r ) )
            continue;

        // Store data so it can be written to an XML file later        
        InitialRodConfiguration initialRodConfiguration;

        initialRodConfiguration.rodOptions = m_rodDataMap[ i_rodGroupIndex ]->getRodOptions( r );
        initialRodConfiguration.gravity = m_rodDataMap[ i_rodGroupIndex ]->getGravity( r );
        initialRodConfiguration.massDamping = m_rodDataMap[ i_rodGroupIndex ]->getMassDamping( r );

        initialRodConfiguration.initialRodVertexPositions.resize( m_rodDataMap[ i_rodGroupIndex ]->elasticRod( r )->nv() );
        
        for ( size_t v=0;v<initialRodConfiguration.initialRodVertexPositions.size(); ++v )
        {
            initialRodConfiguration.initialRodVertexPositions[ v ] = 
                m_rodDataMap[ i_rodGroupIndex ]->elasticRod( r )->getVertex( v );
        }
    
        m_initialRodConfigurations.push_back( initialRodConfiguration );
      
        m_rods.push_back( m_rodDataMap[ i_rodGroupIndex ]->elasticRod( r ) );

        m_world->addObject( m_rodDataMap[ i_rodGroupIndex ]->elasticRod( r ) );
        m_world->addController( m_rodDataMap[ i_rodGroupIndex ]->collisionStepper( r ) );
    }
}

void Beaker::startXMLLogging( std::string& i_xmlFilePath, std::string& i_mayaSceneFilename )
{
    m_sceneXML = new SceneXML();

    double stepSize = 999.99;

    cerr << "Starting xml logging to file '" << i_xmlFilePath << "'\n";

    m_sceneXML->setInitialSceneState( i_xmlFilePath, m_initialRodConfigurations, i_mayaSceneFilename, stepSize );
}

void Beaker::writeXMLFileToDisk()
{
    cerr << "Writing xml data to disk\n";

    if ( m_sceneXML != NULL )
    {
        m_sceneXML->writeFile(); 
    
        delete m_sceneXML;
        m_sceneXML = NULL;
    }
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
    Scalar targetTime = currentTime + i_stepSize;
    setDt( i_stepSize/i_subSteps );

    double frameObjectCollisionForceTime = 0.0;
    double frameObjectCollisionResponse = 0.0;
    double frameSelfCollisionPenaltyForceTime = 0.0;
    double frameSelfCollisionsResponseTime = 0.0;
    double frameIntegrationStepTime = 0.0;

    double frameTime = 0.0;
    
    // Do a quick check, if no rods are enabled then do nothing.
    bool noRodsToSimulate = true;
    for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
    {
        WmFigRodGroup* pRodGroup = rdmItr->second;

        size_t numRods = pRodGroup->numberOfRods();
        for ( size_t r=0; r<numRods; ++r )
        {
            if ( pRodGroup->shouldSimulateRod( r ) )
            {
              //  cerr << "rod " << r << " is enabled " << endl;
                noRodsToSimulate = false;
                break;
            }   
            else         
            {
                //cerr << "rod " << r << " is disabled " << endl;
            }
        }
        
        if ( !noRodsToSimulate )
        {
            break;
        }
    }

    if ( noRodsToSimulate )
    {
        cerr << "Not doing anything as all rods are disabled!\n";
        return;
    }

    cerr << "Simulating\n";

    // Create space to track the target vertex positions of each rod as they substep towards 
    // their goal
    m_subSteppedVertexPositions.resize( i_subSteps );

    for ( int s=0; s<i_subSteps; s++ )
    {
        int rod_gid = 0;

        m_numRods = 0;

        if ( (targetTime - currentTime) < getDt() + SMALL_NUMBER )
            setDt( targetTime - currentTime );

        // Update CollisionMeshData for this substep
        //
        timeval timer;
        startTimer(timer);
        Scalar interpolateFactor = ( (double)(s+1) / i_subSteps );
        if ( i_collisionsEnabled )
        {
            for ( CollisionMeshDataHashMapIterator cmItr = m_collisionMeshMap.begin();
                                                    cmItr != m_collisionMeshMap.end(); ++cmItr )
            {
                cmItr->second->interpolate( interpolateFactor );
            }
        }
        double timeTaken = stopTimer(timer);
        frameTime += timeTaken;
        m_meshInterpolationTime += timeTaken;

        // interpolate fixed vertex positions and set timestep
        //
        startTimer(timer);
        for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
        {
            WmFigRodGroup* pRodGroup = rdmItr->second;
            size_t numRods = pRodGroup->numberOfRods();
            m_numRods += numRods;

            // for visualising
            m_rodRootMaterialFrame.resize( numRods );
            m_strandRootMaterialFrame.resize( numRods );
            m_rodRefMaterialFrame.resize( numRods );

            m_subSteppedVertexPositions[ s ].resize( numRods );

            for ( size_t r=0; r<numRods; r++ )
            {
                // Check if this is a rod or just a fake place holder as the input was too short
                if ( !pRodGroup->shouldSimulateRod( r ) )
                {
                    continue;
                }
                
                //RodCollisionTimeStepper* rodCollisionTimeStepper = dynamic_cast<RodCollisionTimeStepper*>(rodData[r]->stepper);
                RodCollisionTimeStepper* rodCollisionTimeStepper = pRodGroup->collisionStepper( r );

                // Setup Rod collision time stepper before we do any work
                rodCollisionTimeStepper->shouldDoCollisions( i_collisionsEnabled );
                rodCollisionTimeStepper->setTimeStep( getDt() );
                rodCollisionTimeStepper->setCollisionMeshesMap( &m_collisionMeshMap );
                rodCollisionTimeStepper->setClumping( m_isClumpingEnabled, m_clumpingCoefficient );
    
                BASim::ElasticRod* rod = pRodGroup->elasticRod( r );

                rod->global_rodID = rod_gid;
                rod_gid++;

                m_subSteppedVertexPositions[ s ][ r ].resize( rod->nv() );                
            }

            pRodGroup->updateCurrentVertexPositions( interpolateFactor );

            pRodGroup->updateAllBoundaryConditions();

            // Update the xml data with the rod positions for this time step
            // FIXME: The scene xml data only works if there is *ONLY* one rod group.
            if ( m_sceneXML )
            {
                vector< FrameData > frameData;
                frameData.resize( numRods );

                for ( size_t r=0; r<numRods; ++r )
                {
                    ElasticRod* rod = pRodGroup->elasticRod( r );
                    for ( int v=0; v<rod->nv(); ++v )
                    {
                        if ( rod->vertFixed( v ) )
                        {
                            frameData[ r ].fixedVertices[ v ] = rod->getVertex( v );
                        }
                    }
                }

                m_sceneXML->addFrameData( frameData );
            }
        }
        timeTaken = stopTimer( timer );
        frameTime += timeTaken;
        m_vertexInterpolationTime += timeTaken;

        // Now we have two code paths depending on if we're doing self collisions. If
        // we're not then we are safe to parallelise the entire thing. If we are then
        // we need to break it into blocks so we can run the self collisions in one thread.

        Controllers controllers = m_world->getControllers();
        int numControllers = (int)controllers.size();

        // Trying to use more threads than we have controllers would be dumb.
        // I think OpenMP will not let you as it doesn't make sense but I like
        // to be deliberate about things.
        int actualNumThreadsToUse = i_numberOfThreadsToUse;
        if ( i_numberOfThreadsToUse > numControllers )
            actualNumThreadsToUse = numControllers;

        m_numberofThreadsUsed = actualNumThreadsToUse;

        if ( !i_selfCollisionPenaltyForcesEnabled && !i_fullSelfCollisionsEnabled )
        {
            // Fantastic let's just run it all threaded!
            timeval threadFrameTimer;
            startTimer( threadFrameTimer );

            #pragma omp parallel for num_threads( actualNumThreadsToUse )
            for ( int i=0; i<numControllers; ++i )
            {                
                RodCollisionTimeStepper* rodCollisionTimeStepper = dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ]);

                rodCollisionTimeStepper->initialiseCollisionsAndApplyObjectCollisionForces();
            }

            // Clumping is like self collisions, has to all be done at the same time
            if ( m_isClumpingEnabled )
            {
                RodCollisionTimeStepper::getClumpingPairs( m_rods );
            }

            #pragma omp parallel for num_threads( actualNumThreadsToUse )
            for ( int i=0; i<numControllers; ++i )
            {
                RodCollisionTimeStepper* rodCollisionTimeStepper = dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ]);

                ElasticRod* elasticRod = rodCollisionTimeStepper->getRod();

                rodCollisionTimeStepper->execute();

                rodCollisionTimeStepper->respondToObjectCollisions();
                rodCollisionTimeStepper->tidyUpCollisionStructuresForNextStep();
            }
            
            frameTime += stopTimer( threadFrameTimer );
        }
        else
        {
            // Boo, the user wants to use self collisions so we need to split it all up into sections.
            startTimer(timer);
            #pragma omp parallel for num_threads( actualNumThreadsToUse )
            for ( int i=0; i<numControllers; ++i )
            {
                dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ])->initialiseCollisionsAndApplyObjectCollisionForces();
            }
            timeTaken = stopTimer(timer);
            frameObjectCollisionForceTime += timeTaken;
            m_objectCollisionForces += timeTaken;
            frameTime += timeTaken;
            
            if ( m_isClumpingEnabled )
            {
                RodCollisionTimeStepper::getClumpingPairs(m_rods);
            }

            startTimer(timer);
            // FIXME: Think this can be parallelised
            if( i_selfCollisionPenaltyForcesEnabled )
            {
                RodCollisionTimeStepper::getProximities(m_rods);
            }
            timeTaken = stopTimer(timer);
            frameSelfCollisionPenaltyForceTime += timeTaken;
            m_selfCollisionPenaltyForceTime += timeTaken;
            frameTime += timeTaken;

            startTimer(timer);
            #pragma omp parallel for num_threads( actualNumThreadsToUse )
            for ( int i=0; i<numControllers; ++i )
            {
                dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ])->execute();
            }
            timeTaken = stopTimer(timer);
            frameIntegrationStepTime += timeTaken;
            m_integrationStepTime += timeTaken;
            frameTime += timeTaken;

            //////////////////////////////////////////////
            //
            // Check forces to see if any rods need to be simulated at a slower pace
            // checkAllRodForces();
            //////////////////////////////////////////////

            startTimer(timer);
            #pragma omp parallel for num_threads( actualNumThreadsToUse )
            for ( int i=0; i<numControllers; ++i )
            {
                dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ])->respondToObjectCollisions();              
            }
            timeTaken = stopTimer(timer);
            frameObjectCollisionResponse += timeTaken;
            m_objectCollisionResponse += timeTaken;
            frameTime += timeTaken;

            startTimer(timer);
            if (i_fullSelfCollisionsEnabled)
               RodCollisionTimeStepper::respondRodCollisions( m_rods, getDt(), i_fullSelfCollisionIters,
                                                              i_selfCollisionCOR );
            timeTaken = stopTimer(timer);
            frameSelfCollisionsResponseTime += timeTaken;
            m_selfCollisionsResponseTime += timeTaken;
            frameTime += timeTaken;

            startTimer(timer);
            #pragma omp parallel for num_threads( actualNumThreadsToUse )
            for ( int i=0; i<numControllers; ++i )
            {
                dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ])->tidyUpCollisionStructuresForNextStep();
            }
            timeTaken = stopTimer(timer);
            m_collisionStructuresTidyup += timeTaken;
            frameTime += timeTaken;
        }
        // This sets the undeformed rod configuration to be the same as the current configuration.
        // I'm not sure if this is actually what we want to do, is this just like pushing around line
        // segments. Maybe we should only run this on rods that have collided.

        if ( m_plasticDeformations )
        {
            for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
            {
                WmFigRodGroup* pRodGroup = rdmItr->second;
                size_t numRods = pRodGroup->numberOfRods();

                for ( size_t r=0; r<numRods; r++ )
                {
                    BASim::ElasticRod* rod = pRodGroup->elasticRod( r );
                    rod->updateReferenceProperties();
                }
            }
        }

        setTime( currentTime + getDt() );
        currentTime = getTime();
    }

    if ( frameObjectCollisionForceTime > m_slowestCollisionForcesTime )
        m_slowestCollisionForcesTime = frameObjectCollisionForceTime;
    if ( frameObjectCollisionForceTime < m_fastestCollisionForcesTime )
        m_fastestCollisionForcesTime = frameObjectCollisionForceTime;

    if ( frameIntegrationStepTime > m_slowestIntegrationTime )
        m_slowestIntegrationTime = frameIntegrationStepTime;
    if ( frameIntegrationStepTime < m_fastestIntegrationTime )
        m_fastestIntegrationTime = frameIntegrationStepTime;

    if ( frameObjectCollisionResponse > m_slowestCollisionResponseTime )
        m_slowestCollisionResponseTime = frameObjectCollisionResponse;
    if ( frameObjectCollisionResponse < m_fastestCollisionResponseTime )
        m_fastestCollisionResponseTime = frameObjectCollisionResponse;

    if ( frameSelfCollisionPenaltyForceTime > m_slowestSelfCollisionPenaltyForceTime )
        m_slowestSelfCollisionPenaltyForceTime = frameSelfCollisionPenaltyForceTime;
    if ( frameSelfCollisionPenaltyForceTime < m_fastestSelfCollisionPenaltyForceTime )
        m_fastestSelfCollisionPenaltyForceTime = frameSelfCollisionPenaltyForceTime;

    if ( frameSelfCollisionsResponseTime > m_slowestSelfCollisionsResponseTime )
        m_slowestSelfCollisionsResponseTime = frameSelfCollisionsResponseTime;
    if ( frameSelfCollisionsResponseTime < m_fastestSelfCollisionsResponseTime )
        m_fastestSelfCollisionsResponseTime = frameSelfCollisionsResponseTime;

    m_numberOfFramesSimulated++;
    if ( frameTime < m_fastestFrameTime )
        m_fastestFrameTime = frameTime;
    if ( frameTime > m_slowestFrameTime )
        m_slowestFrameTime = frameTime;

    m_totalSimTime += frameTime;

    //storeMaterialFrames();

    // restore dt
    setDt( dt_save );

//    printTimingInfo();
}
/*
void Beaker::storeMaterialFrames()
{
    for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
    {
        vector<RodData*>& rodData = rdmItr->second;
        size_t numRods = rodData.size();
        for ( size_t r=0; r<numRods; r++ )
        {
            ElasticRod* rod = rodData[ r ]->rod;
            for ( size_t e=0; e<rod->ne(); e++ )
            {
                rodData[ r ]->materialFrame1[ e ] = rod->getMaterial1( e );
                rodData[ r ]->materialFrame2[ e ] = rod->getMaterial2( e );
                rodData[ r ]->materialFrame3[ e ] = rod->getEdge( e );
            }
        }

    }
}
*/

/*void Beaker::checkAllRodForces()
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
}*/

/*RodCollisionTimeStepper* Beaker::setupRodTimeStepper( RodData* i_rodData )
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
}*/

void Beaker::draw()
{   if ( m_shouldDrawSubsteppedVertices )
    {
        // Draw the onion skinned interpolated vertex positions
        for ( size_t s=0; s<m_subSteppedVertexPositions.size(); ++s )
        {
            for ( size_t r=0; r<m_subSteppedVertexPositions[ s ].size(); ++r )
            {
                float fraction = (float)(s)/m_subSteppedVertexPositions.size();
                glColor3f( fraction, fraction, fraction );
                glBegin( GL_LINE_STRIP );
                    
                for ( size_t c=0; c<m_subSteppedVertexPositions[ s ][ r ].size(); ++c )
                {
                    Vec3d p = m_subSteppedVertexPositions[ s ][ r ][ c ];
                    glVertex3d( p[ 0 ], p[ 1 ], p[ 2 ] );
                }
    
                glEnd();
            }
        }
    }

    /*for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
    {
        vector<RodData*>& rodData = rdmItr->second;
        size_t numRods = rodData.size();
        for ( size_t r=0; r<numRods; r++ )
        {
            rodData[ r ]->rodRenderer->render();
        }
    }*/

    /*glLineWidth(5.0);
    glBegin( GL_LINES );
    for ( size_t r=0; r<m_rodRootMaterialFrame.size(); r++ )
    {
        glColor3d(1,1,1);
        Vec3d p0 = Vec3d(0,0,0);
        Vec3d p1 = p0 + m_rodRootMaterialFrame[r].m1;
        glVertex3d( p0[0], p0[1], p0[2] );
        glVertex3d( p1[0], p1[1], p1[2] );
        p1 = p0 + m_rodRootMaterialFrame[r].m2;
        glVertex3d( p0[0], p0[1], p0[2] );
        glVertex3d( p1[0], p1[1], p1[2] );
        p1 = p0 + m_rodRootMaterialFrame[r].m3;
        glVertex3d( p0[0], p0[1], p0[2] );
        glVertex3d( p1[0], p1[1], p1[2] );

        glColor3d(0,0,1);
        p1 = p0 + m_strandRootMaterialFrame[r].m1;
        glVertex3d( p0[0], p0[1], p0[2] );
        glVertex3d( p1[0], p1[1], p1[2] );
        p1 = p0 + m_strandRootMaterialFrame[r].m2;
        glVertex3d( p0[0], p0[1], p0[2] );
        glVertex3d( p1[0], p1[1], p1[2] );
        p1 = p0 + m_strandRootMaterialFrame[r].m3;
        glVertex3d( p0[0], p0[1], p0[2] );
        glVertex3d( p1[0], p1[1], p1[2] );

        glColor3d(1,0,0);
        p1 = p0 + m_rodRefMaterialFrame[r].m1*2;
        glVertex3d( p0[0], p0[1], p0[2] );
        glVertex3d( p1[0], p1[1], p1[2] );
        p1 = p0 + m_rodRefMaterialFrame[r].m2*2;
        glVertex3d( p0[0], p0[1], p0[2] );
        glVertex3d( p1[0], p1[1], p1[2] );
        p1 = p0 + m_rodRefMaterialFrame[r].m3*2;
        glVertex3d( p0[0], p0[1], p0[2] );
        glVertex3d( p1[0], p1[1], p1[2] );
    }
    glEnd();
    glLineWidth(1.0);*/

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
