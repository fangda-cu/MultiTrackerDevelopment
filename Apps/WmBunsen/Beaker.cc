/**
 * \file Beaker.cc
 *
 * \author acoull@wetafx.co.nz
 * \date 10/26/2009
 */

#include "Beaker.hh"

#ifdef WETA
//#include <weta/Wfigaro/Physics/ElasticRods/RodHairsprayForce.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodMassDamping.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodGravity.hh>
#else
//#include <BASim/src/Physics/ElasticRods/RodHairsprayForce.hh>
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

Beaker::Beaker() :
    m_plasticDeformations(false), m_gravity(0, -981.0, 0), m_timingEnabled(false), m_timingsFile(""),
            m_meshInterpolationTime(0.0), m_vertexInterpolationTime(0.0), m_objectCollisionForces(0.0),
            m_objectCollisionResponse(0.0), m_collisionStructuresTidyup(0.0), m_selfCollisionPenaltyForceTime(0.0),
            m_selfCollisionsResponseTime(0.0), m_integrationStepTime(0.0), m_slowestIntegrationTime(0.0),
            m_fastestIntegrationTime(99999999999999.0), m_fastestSelfCollisionPenaltyForceTime(99999999999999.0),
            m_slowestSelfCollisionPenaltyForceTime(0.0), m_fastestSelfCollisionsResponseTime(999999999999.0),
            m_slowestSelfCollisionsResponseTime(0.0), m_slowestCollisionForcesTime(0.0),
            m_fastestCollisionForcesTime(999999999.0), m_slowestCollisionResponseTime(0.0),
            m_fastestCollisionResponseTime(9999999999999.0), m_fastestFrameTime(9999999999999999.0), m_slowestFrameTime(0.0),
            m_totalSimTime(0.0), m_numberOfFramesSimulated(0), m_numberofThreadsUsed(0), m_numRods(0),
            m_shouldDrawSubsteppedVertices(false), m_isClumpingEnabled(false), m_clumpingCoefficient(0.3),
            m_isXMLLoggingEnabled(false), /*m_sceneXML( NULL ), m_volumetricCollisions( NULL ),*/
            m_flip(), m_slip(), m_doVolumetricCollisions(false), m_targetEdgeDensity(100.0), m_volumetricRadius(1.0),
            m_gridDX(1.0), m_displayGrid(false), m_displayGridVelocitiesMultiplier(0.0), m_maxDisplayDensity(),
            m_displayCollisionBoundary(false), m_displayAirBoundary(false), m_stol(1.0e-6f * 0.01), m_atol(1.0e-6f * 0.01),
            m_rtol(1.0e-6f * 0.01), m_inftol(1.0e-6f * 0.01), m_BARodStepper(NULL), m_stopOnRodError(false)
{
    m_separationCondition[0] = m_separationCondition[1] = m_separationCondition[2] = -1.0;

    m_rodDataMap.clear();

    //m_initialRodConfigurations.clear();

    initialiseWorld(0.0, 0.01);
}

Beaker::~Beaker()
{
    //resetEverything();

    if (m_timingsFP.is_open())
        m_timingsFP.close();

    delete m_world;
    delete m_BARodStepper;

    // if ( m_sceneXML != NULL )
    // delete m_sceneXML;
}

void Beaker::initialiseWorld(const double i_time, const double i_dt)
{
    m_world = new World();
    cerr << "intialising world with time " << i_time << endl;
    m_world->add_property(m_timeHandle, "time", i_time);
    m_world->add_property(m_dtHandle, "time-step", i_dt);
    m_world->add_property(m_gravityHandle, "gravity", m_gravity);
    m_world->add_property(m_maxIterHandle, "maxIter", 100);
}

void Beaker::resetEverything(const double i_time, const double i_dt)
{
    /*for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
     {
     vector<RodData*>& rodData = rdmItr->second;
     int numRods = rodData.size();
     for ( int r=0; r<numRods; r++ )
     {
     // We're safe to clear this vector as the individual destructors will safely delete the
     // rod data in each vector element.
     delete rodData[ r ];
     }
     }*/
    m_rodDataMap.clear();

    m_rods.clear();

    m_collisionMeshDataHashMap.clear();

    delete m_world;

    initialiseWorld(i_time, i_dt);

    //delete m_volumetricCollisions;
    //m_volumetricCollisions = NULL;
}

void Beaker::setTimingsFile(std::string i_fileName)
{
    m_timingsFile = i_fileName;
}

void Beaker::setTimingEnabled(bool i_timingsEnabled)
{
    if (i_timingsEnabled == true && m_timingEnabled == false)
    {
        m_timingEnabled = i_timingsEnabled;
        resetTimers();

        // Open timings file for timer code to output to
        if (m_timingsFP.is_open())
            m_timingsFP.close();

        if (m_timingsFile != "")
        {
            m_timingsFP.open(m_timingsFile.c_str(), ios::out | ios::trunc);
        }
    }
    else if (i_timingsEnabled == false && m_timingEnabled == true)
    {
        printTimingInfo();

        // Close timings file
        if (m_timingsFP.is_open())
            m_timingsFP.close();
    }
    else
    {
        // Do nothing as the state hasn't changed
    }

    m_timingEnabled = i_timingsEnabled;
}

void Beaker::setStopOnRodError(bool i_stopOnRodError)
{
    m_stopOnRodError = i_stopOnRodError;
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

inline std::string makeString(int i_val)
{
    return str(boost::format("%d") % i_val);
}

inline std::string makeString(size_t i_val)
{
    return str(boost::format("%z") % i_val);
}

inline std::string makeString(double i_val)
{
    return str(boost::format("%.2f") % i_val);
}

void Beaker::printTimingInfo()
{
    if (m_numberOfFramesSimulated == 0)
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
    resultString += "Number of rods = " + makeString(m_numRods) + "\n";
    resultString += "Number of threads used = " + makeString(m_numberofThreadsUsed) + "\n";
    resultString += "Total frames simulated = " + makeString(m_numberOfFramesSimulated) + "\n";
    resultString += "Total simulation time = " + makeString(totalSimTime) + "\n";
    resultString += "Fastest frame = " + makeString(m_fastestFrameTime) + "\n";
    resultString += "Slowest frame = " + makeString(m_slowestFrameTime) + "\n";
    resultString += "Average frame = " + makeString(averageFrameTime) + "\n";
    resultString += "\nDetailed breakdown (Total / Average):\n";
    resultString += "Mesh Interpolation (SINGLE THREADED) = " + makeString(m_meshInterpolationTime) + " / " + makeString(
            m_meshInterpolationTime / m_numberOfFramesSimulated) + "\n";
    resultString += "Vertex Interpolation (SINGLE THREADED) = " + makeString(m_vertexInterpolationTime) + " / " + makeString(
            m_vertexInterpolationTime / m_numberOfFramesSimulated) + "\n";
    resultString += "Object Collision Forces = " + makeString(m_objectCollisionForces) + " / " + makeString(
            m_objectCollisionForces / m_numberOfFramesSimulated) + "\n";
    resultString += "\tSlowest Frame = " + makeString(m_slowestCollisionForcesTime) + "\n";
    resultString += "\tFastest Frame = " + makeString(m_fastestCollisionForcesTime) + "\n";
    resultString += "Object Collision Response = " + makeString(m_objectCollisionResponse) + " / " + makeString(
            m_objectCollisionResponse / m_numberOfFramesSimulated) + "\n";
    resultString += "\tSlowest Frame = " + makeString(m_slowestCollisionResponseTime) + "\n";
    resultString += "\tFastest Frame = " + makeString(m_fastestCollisionResponseTime) + "\n";
    resultString += "Collision Structures Tidy up = " + makeString(m_collisionStructuresTidyup) + " / " + makeString(
            m_collisionStructuresTidyup / m_numberOfFramesSimulated) + "\n";
    resultString += "Self Collision Penalty Force (SINGLE THREADED) = " + makeString(m_selfCollisionPenaltyForceTime) + " / "
            + makeString(m_selfCollisionPenaltyForceTime / m_numberOfFramesSimulated) + "\n";
    resultString += "\tSlowest Frame = " + makeString(m_slowestSelfCollisionPenaltyForceTime) + "\n";
    resultString += "\tFastest Frame = " + makeString(m_fastestSelfCollisionPenaltyForceTime) + "\n";
    resultString += "Self Collision Response (SINGLE_THREADED) = " + makeString(m_selfCollisionsResponseTime) + " / "
            + makeString(m_selfCollisionsResponseTime / m_numberOfFramesSimulated) + "\n";
    resultString += "\tSlowest Frame = " + makeString(m_slowestSelfCollisionsResponseTime) + "\n";
    resultString += "\tFastest Frame = " + makeString(m_fastestSelfCollisionsResponseTime) + "\n";
    resultString += "Integration = " + makeString(m_integrationStepTime) + " / " + makeString(
            m_integrationStepTime / m_numberOfFramesSimulated) + "\n";
    resultString += "\tSlowest Frame = " + makeString(m_slowestIntegrationTime) + "\n";
    resultString += "\tFastest Frame = " + makeString(m_fastestIntegrationTime) + "\n";

    cerr << resultString;

    if (m_timingsFP.is_open())
        m_timingsFP << resultString;
}

void Beaker::startTimer(timeval& i_startTimer)
{
    if (m_timingEnabled)
        gettimeofday(&i_startTimer, NULL);
    //gettimeofday( &m_timerStart, NULL );
}

double Beaker::stopTimer(timeval& i_startTimer)
{
    if (!m_timingEnabled)
        return 0.0;

    timeval end;
    gettimeofday(&end, NULL);

    return (((end.tv_sec * 1000000 + end.tv_usec) - (i_startTimer.tv_sec * 1000000 + i_startTimer.tv_usec)) / 1000000.0);
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

void Beaker::addRodsToWorld(int i_rodGroupIndex, WmFigRodGroup* i_rodGroup, double startTime, int numberOfThreads,
        PerformanceTuningParameters perf_param)
{
    std::cout << "Performance Tuning Parameters " << std::endl;
    std::cout << "Penalty Response " << perf_param.m_enable_penalty_response << std::endl;
    std::cout << "Implicit Thickenss " << perf_param.m_implicit_thickness << std::endl;
    std::cout << "Implicit Stiffness " << perf_param.m_implicit_stiffness << std::endl;
    std::cout << "Inextensibility Threshold " << perf_param.m_inextensibility_threshold << std::endl;
    std::cout << "Max Num Solver Iterations " << perf_param.m_solver.m_max_iterations << std::endl;
    std::cout << "Max Num Collision Iterations  " << perf_param.m_collision.m_max_iterations << std::endl;
    std::cout << "Explosion Detection  " << perf_param.m_enable_explosion_detection << std::endl;
    std::cout << "Explosion Dampening " << perf_param.m_explosion_damping << std::endl;
    std::cout << "Explosion Threshold " << perf_param.m_explosion_threshold << std::endl;
    std::cout << "Stretching threshold" << perf_param.m_stretching_threshold << std::endl;
    std::cout << "Solver failure " << perf_param.m_solver.m_in_case_of << std::endl;
    std::cout << "Max Number of Solver Substeps " << perf_param.m_solver.m_max_substeps << std::endl;
    std::cout << "Collison Failure  " << perf_param.m_collision.m_in_case_of << std::endl;
    std::cout << "Max Number of Collision Substeps " << perf_param.m_collision.m_max_substeps << std::endl;
    std::cout << "Explosion  Failure  " << perf_param.m_explosion.m_in_case_of << std::endl;
    std::cout << "Max Number of Explosion Substeps " << perf_param.m_explosion.m_max_substeps << std::endl;
    std::cout << "Stretching  Failure  " << perf_param.m_stretching.m_in_case_of << std::endl;
    std::cout << "Max Number of Stretching Substeps " << perf_param.m_stretching.m_max_substeps << std::endl;

    m_rodDataMap[i_rodGroupIndex] = i_rodGroup;

    int numRods = m_rodDataMap[i_rodGroupIndex]->numberOfRods();

    m_rods.clear();
    m_rodTimeSteppers.clear();
    m_triangleMeshes.clear();
    m_levelSets.clear();
    m_scriptingControllers.clear();

    //  cerr << "about to add rods to world, dt = " << getDt() << endl;

    bool areSimulatingAnyRods = false;
    for (int r = 0; r < numRods; r++)
    {
        if (!m_rodDataMap[i_rodGroupIndex]->shouldSimulateRod(r))
            continue;

        areSimulatingAnyRods = true;

        // cerr << "Adding rod " << r << " to world\n";
        m_rods.push_back(m_rodDataMap[i_rodGroupIndex]->elasticRod(r));
        m_world->addObject(m_rodDataMap[i_rodGroupIndex]->elasticRod(r));

        // cerr << "Adding rod time stepper " << r << " to world\n";
        m_rodTimeSteppers.push_back(m_rodDataMap[i_rodGroupIndex]->stepper(r));
    }

    if (!areSimulatingAnyRods)
    {
        cerr << "No rods being simulated so not setting up the simulation world.\n";
        return;
    }

    // Now add all the collision meshes and scripting controllers to the world
    for (CollisionMeshDataHashMap::iterator cmItr = m_collisionMeshDataHashMap.begin(); cmItr
            != m_collisionMeshDataHashMap.end(); ++cmItr)
    {
        m_scriptingControllers.push_back(cmItr->second->scriptingController());
        m_triangleMeshes.push_back(cmItr->second->triangleMesh());
        m_levelSets.push_back(cmItr->second->levelSet());

        cerr << "Added scripting controller\n";
        cerr << "Added triangle mesh\n";
        if (m_levelSets[m_levelSets.size() - 1] != NULL)
        {
            cerr << "Added level set\n";
        }
        else
        {
            cerr << "Mesh does not have level set data, adding null pointer\n";
        }
    }

    //cerr << "adding rods to world, dt = " << getDt() << endl;

    // FIXME: pass in timestep from Maya, it's ok to do this for test as the real timestep
    // is set at the beginning of takeTimeStep() but it's really sloppy to not bother setting it 
    // right to start with!
    m_BARodStepper = new BARodStepper(m_rods, m_triangleMeshes, m_scriptingControllers, m_rodTimeSteppers, 1.0 / 24.0,
            startTime, numberOfThreads, perf_param, &m_levelSets);
    m_world->addController(m_BARodStepper);

    // For convenience, give the RodData a pointer to bridsonStepper
    for (int r = 0; r < numRods; r++)
    {
        if (!m_rodDataMap[i_rodGroupIndex]->shouldSimulateRod(r))
            continue;
    }
    m_rodDataMap[i_rodGroupIndex]->setBridsonStepperOnAllRodData(m_BARodStepper);

    /* delete m_volumetricCollisions;

     if ( areSimulatingAnyRods )
     {
     m_volumetricCollisions = new VolumetricCollisionsCPU( m_rodDataMap );
     }*/
}

void Beaker::startXMLLogging(std::string& i_xmlFilePath, std::string& i_mayaSceneFilename)
{
    /*m_sceneXML = new SceneXML();

     double stepSize = 999.99;

     cerr << "Starting xml logging to file '" << i_xmlFilePath << "'\n";

     m_sceneXML->setInitialSceneState( i_xmlFilePath, m_initialRodConfigurations, i_mayaSceneFilename, stepSize );*/
}

void Beaker::writeXMLFileToDisk()
{
    /*cerr << "Writing xml data to disk\n";

     if ( m_sceneXML != NULL )
     {
     m_sceneXML->writeFile();

     delete m_sceneXML;
     m_sceneXML = NULL;
     }*/
}

bool Beaker::anyRodsActive()
{
    bool rodsToSimulate = false;
    for (RodDataMapIterator rdmItr = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr)
    {
        WmFigRodGroup* pRodGroup = rdmItr->second;

        int numRods = pRodGroup->numberOfRods();
        for (int r = 0; r < numRods; ++r)
        {
            if (pRodGroup->shouldSimulateRod(r))
            {
                rodsToSimulate = true;
                break;
            }
        }

        if (!rodsToSimulate)
        {
            break;
        }
    }

    return rodsToSimulate;
}
/*
 void Beaker::interpolateCollisionMeshesForSubstep( const float i_interpolateFactor, float* o_timeTaken )
 {
 // Update CollisionMeshData for this substep
 //
 timeval timer;
 startTimer( timer );

 for ( CollisionMeshDataHashMapIterator cmItr = m_collisionMeshMap.begin();
 cmItr != m_collisionMeshMap.end(); ++cmItr )
 {
 cmItr->second->interpolate( i_interpolateFactor );
 }

 *o_timeTaken = stopTimer( timer );
 }
 */

//void Beaker::setupRodTimeStepperForSubStep( WmFigRodGroup* i_pRodGroup, const int i_subStep,
//    const bool i_collisionsEnabled )
//{
//    int rod_gid = 0;
//    int numRods = i_pRodGroup->numberOfRods();
//
//    for ( int r=0; r<numRods; r++ )
//    {
//        // Check if this is a rod or just a fake place holder as the input was too short
//        if ( !i_pRodGroup->shouldSimulateRod( r ) )
//        {
//            continue;
//        }
//
//        // We don't need to do this any more as BARodStepper takes care
//        // of the setting the timestep in rodTimeStepper's as it basically
//        // owns them
//        /*
//        RodTimeStepper* rodTimeStepper = i_pRodGroup->stepper( r );
//
//        rodTimeStepper->setTimeStep( getDt() );
//        //rodTimeStepper->setCollisionMeshesMap( &m_collisionMeshMap );
//        //rodTimeStepper->setClumping( m_isClumpingEnabled, m_clumpingCoefficient );
//         * */
//
//        BASim::ElasticRod* rod = i_pRodGroup->elasticRod( r );
//
//        rod->global_rodID = rod_gid;
//        rod_gid++;
//
//      //  m_subSteppedVertexPositions[ i_subStep ][ r ].resize( rod->nv() );
//    }
//}

//int Beaker::calculateNumSubSteps(int numSubSteps, Scalar deltaT, double subDistMax)
//{
//     // Set an adaptive substep if velocity change is too great.  The subDistanceMax is a tunable parameter.
//     int steps=numSubSteps;
//
//     /*float biggestMax = 0.0f;
//     for ( CollisionMeshDataHashMapIterator cmItr = m_collisionMeshMap.begin();
//                                            cmItr != m_collisionMeshMap.end(); ++cmItr )
//     {
//         float maxVelMag = cmItr->second->m_maxVelocityMag;
//         if ( maxVelMag > biggestMax )
//             biggestMax = maxVelMag;
//     }
//
//     if ( ((biggestMax *deltaT)/numSubSteps) > subDistMax)
//     {
//         steps= (biggestMax * deltaT)/subDistMax;
//         if (steps < numSubSteps)
//         {
//             steps= numSubSteps;
//         }
//     }
//     std::cout<< "Max velocity = " << biggestMax <<" Substeps = " << steps << std::endl;
//*/
//     return steps;
//}

void Beaker::takeTimeStep(
        // int i_numberOfThreadsToUse,
        Scalar i_stepSize,
        // int i_subSteps,
        // double i_subDistanceMax,
        bool i_collisionsEnabled, bool i_selfCollisionPenaltyForcesEnabled, bool i_fullSelfCollisionsEnabled,
        int i_fullSelfCollisionIters, double i_selfCollisionCOR, FixedRodVertexMap* i_fixedVertices, bool i_zeroAllTwist,
        double i_constraintSrength, LockedRodVertexMap* i_lockedRodVertexMap)
{
    // i_fixedVertices is a thing I'm trying out where vertices can be fixed temporarily, only
    // for this time step. It's getting used in the rod shape where users can select vertices
    // to move around and those should be set as boundary conditions rather than simulated

    // Check if anything has actually been initialised yet. We may still be being loaded by Maya.
    if (m_rodDataMap.size() == 0)
        return;

    if (!anyRodsActive())
    {
        //cerr << "Not doing anything as all rods are disabled!\n";
        return;
    }

    //i_subSteps = calculateNumSubSteps( i_subSteps, i_stepSize, i_subDistanceMax);
    // std::cout << m_stopOnRodError << std::endl;
    m_BARodStepper->setStopOnRodError(m_stopOnRodError);
    Scalar dt_save = getDt();
    //Scalar startTime = getTime();
    //Scalar currentTime = getTime();
    Scalar targetTime = getTime() + i_stepSize;
    // setDt( i_stepSize / i_subSteps );

    // Initialise all the timers to 0

    double frameObjectCollisionForceTime = 0.0;
    double frameObjectCollisionResponse = 0.0;
    double frameSelfCollisionPenaltyForceTime = 0.0;
    double frameSelfCollisionsResponseTime = 0.0;
    double frameIntegrationStepTime = 0.0;

    double frameTime = 0.0;

    m_BARodStepper->skipRodRodCollisions(!i_fullSelfCollisionsEnabled);

    // Create space to track the target vertex positions of each rod as they substep towards 
    // their goal
    // m_subSteppedVertexPositions.resize( i_subSteps );

    // We should alter if self collisions are on or off here but they can only be changed
    // when we create bridsonStepper at time == startTime.

    // Interpolate fixed vertex positions and set timestep
    //
    //timeval timer;
    //startTimer(timer);
    for (RodDataMapIterator rdmItr = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr)
    {
        WmFigRodGroup* pRodGroup = rdmItr->second;
        int numRods = pRodGroup->numberOfRods();
        m_numRods += numRods;

        //pRodGroup->updateCurrentVertexPositions( interpolateFactor );
        pRodGroup->updateAllBoundaryConditions();
    }
    //timeTaken = stopTimer( timer );
    //frameTime += timeTaken;
    //m_vertexInterpolationTime += timeTaken;

    //  for ( int s=0; s<i_subSteps; s++ )
    {
        const int s = 0;
        //cout << "\nframe "<<24*currentTime<<", substep " << s << " / " << i_subSteps << endl;
        m_numRods = 0;

        // Make sure we don't step past the end of this frame, depending on
        // the size of dt, we may need to take a smaller step so we land exactly
        // on the frame boundary.
        //   if ( ( targetTime - getTime() ) < getDt() + SMALL_NUMBER )
        //    {
        //       setDt( targetTime - getTime() );
        //   }

        //    float timeTaken;
        //   Scalar interpolateFactor = ( (double)( s + 1 ) / i_subSteps );

        m_world->execute();
    }
    // restore dt
    //  setDt( dt_save );

    /*// Now we have two code paths depending on if we're doing self collisions. If
     // we're not then we are safe to parallelise the entire thing. If we are then
     // we need to break it into blocks so we can run the self collisions in one thread.

     World::Controllers controllers = m_world->getControllers();
     int numControllers = (int)controllers.size();

     // Trying to use more threads than we have controllers would be dumb.
     // I think OpenMP will not let you as it doesn't make sense but I like
     // to be deliberate about things.
     int actualNumThreadsToUse = i_numberOfThreadsToUse;
     if ( i_numberOfThreadsToUse > numControllers )
     actualNumThreadsToUse = numControllers;

     m_numberofThreadsUsed = actualNumThreadsToUse;

     //rodCollisionTimeStepper->updateVertexPositionConstraints();

     //  RodCollisionTimeStepper* collisionStepper = pRodGroup->collisionStepper(

     //if ( !i_selfCollisionPenaltyForcesEnabled && !i_fullSelfCollisionsEnabled && !m_doVolumetricCollisions )
     {
     // Fantastic let's just run it all threaded!
     timeval threadFrameTimer;
     startTimer( threadFrameTimer );

     /*#pragma omp parallel for num_threads( actualNumThreadsToUse )
     for ( int i=0; i<numControllers; ++i )
     {
     RodCollisionTimeStepper* rodCollisionTimeStepper = dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ]);

     rodCollisionTimeStepper->initialiseCollisionsAndApplyObjectCollisionForces();
     }

     // Clumping is like self collisions, has to all be done at the same time
     if ( m_isClumpingEnabled )
     {
     RodCollisionTimeStepper::getClumpingPairs( m_rods );
     }*/

    /*  #pragma omp parallel for num_threads( actualNumThreadsToUse )
     for ( int i=0; i<numControllers; ++i )
     {
     RodTimeStepper *stepper = dynamic_cast<RodTimeStepper*>(controllers[ i ]);
     /* ElasticRod* elasticRod = stepper->getRod();

     if ( stepper != NULL )
     {
     stepper->set_stol( m_stol );
     stepper->set_atol( m_atol );
     stepper->set_rtol( m_rtol );
     stepper->set_inftol( m_inftol );
     }*/

    /*    ElasticRod* rod = stepper->getRod();

     if ( stepper->isEnabled() )
     {
     if ( !stepper->execute() )
     {
     cerr << "rod " << i << " could not solve at minimum step size. DISABLING\n";
     stepper->setEnabled( false );
     continue;
     }
     else
     cerr << "Rod executed succesfully!\n";
     }
     /// ElasticRod* rod = stepper->getRod();

     // FIXME: Should I keep the rods vector with all actual active simulated rods in
     // it? Then just pass it each timestep? Would make it simpler than tracking
     // active rods in the volumetric code

     //if(m_volumetricCollisionsEnabled)
     //m_volumetricCollisions->respondVolumetricCollisions( m_rodDataMap, m_targetEdgeDensity,
     //           m_volumetricRadius, m_gridDx, m_separationCondition, m_collisionMeshMap );
     /*rodCollisionTimeStepper->collisionsBegin();
     rodCollisionTimeStepper->respondToObjectCollisions();
     rodCollisionTimeStepper->tidyUpCollisionStructuresForNextStep();*/

    /*}

     frameTime += stopTimer( threadFrameTimer );
     }
     /*else
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
     RodCollisionTimeStepper* rodCollisionTimeStepper = dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ]);
     RodTimeStepper *pStepper = dynamic_cast< RodTimeStepper *>( rodCollisionTimeStepper->getTimeStepper() );
     if ( pStepper != NULL )
     {
     pStepper->set_stol( m_stol );
     pStepper->set_atol( m_atol );
     pStepper->set_rtol( m_rtol );
     pStepper->set_inftol( m_inftol );
     }

     if ( rodCollisionTimeStepper->isEnabled() )
     {
     if( !rodCollisionTimeStepper->execute())
     {
     cerr << "Rod " << i << " could not solve at minimum step size. DISABLING\n";
     rodCollisionTimeStepper->setEnabled( false );
     continue;
     }
     }
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

     #pragma omp parallel for num_threads( actualNumThreadsToUse )
     for ( int i=0; i<numControllers; ++i )
     {
     dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ])->collisionsBegin();
     }

     /////////////// Volumetric
     if ( m_doVolumetricCollisions && m_volumetricCollisions != NULL )
     {
     m_volumetricCollisions->respondVolumetricCollisions( m_rodDataMap, m_targetEdgeDensity,
     m_volumetricRadius, m_gridDX, m_separationCondition, m_flip, m_slip,
     m_collisionMeshMap );
     }

     // We messed with the velocity so update the stored velocities on the rod before
     // doing anything else
     #pragma omp parallel for num_threads( actualNumThreadsToUse )
     for ( int i=0; i<numControllers; ++i )
     {
     dynamic_cast<RodCollisionTimeStepper*>(controllers[ i ])->updateEndPositions();
     }
     ///////////////////////////

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
     }*/
    // This sets the undeformed rod configuration to be the same as the current configuration.
    // I'm not sure if this is actually what we want to do, is this just like pushing around line
    // segments. Maybe we should only run this on rods that have collided.

    /// FIX ME, USE JUNGSEOCK's PLASTIC DEFORMATION STUFF

    /*      if ( m_plasticDeformations )
     {
     for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
     {
     WmFigRodGroup* pRodGroup = rdmItr->second;
     int numRods = pRodGroup->numberOfRods();

     for ( int r=0; r<numRods; r++ )
     {
     BASim::ElasticRod* rod = pRodGroup->elasticRod( r );
     rod->updateReferrenceProperties();
     }
     }
     }*/
    /*
     setTime( currentTime + getDt() );
     currentTime = getTime();
     //   cerr << "End of substep\n";
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

     //    printTimingInfo();*/
}
/*
 void Beaker::storeMaterialFrames()
 {
 for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
 {
 vector<RodData*>& rodData = rdmItr->second;
 int numRods = rodData.size();
 for ( int r=0; r<numRods; r++ )
 {
 ElasticRod* rod = rodData[ r ]->rod;
 for ( int e=0; e<rod->ne(); e++ )
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
 int numRods = rodData.size();
 for ( int r=0; r<numRods; r++ )
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
{
    /*if( m_doVolumetricCollisions && m_volumetricCollisions != NULL)
     m_volumetricCollisions->draw( m_displayGrid, m_displayGridVelocitiesMultiplier,
     m_maxDisplayDensity, m_displayCollisionBoundary ,m_displayAirBoundary );
     */

    if (m_shouldDrawSubsteppedVertices)
    {
        // Draw the onion skinned interpolated vertex positions
        for (int s = 0; s < (int) m_subSteppedVertexPositions.size(); ++s)
        {
            for (int r = 0; r < m_subSteppedVertexPositions[s].size(); ++r)
            {
                float fraction = (float) (s) / m_subSteppedVertexPositions.size();
                glColor3f(fraction, fraction, fraction);
                glBegin( GL_LINE_STRIP);

                for (int c = 0; c < m_subSteppedVertexPositions[s][r].size(); ++c)
                {
                    BASim::Vec3d p = m_subSteppedVertexPositions[s][r][c];
                    glVertex3d(p[0], p[1], p[2]);
                }

                glEnd();
            }
        }
    }

    /*for ( RodDataMapIterator rdmItr  = m_rodDataMap.begin(); rdmItr != m_rodDataMap.end(); ++rdmItr )
     {
     vector<RodData*>& rodData = rdmItr->second;
     int numRods = rodData.size();
     for ( int r=0; r<numRods; r++ )
     {
     rodData[ r ]->rodRenderer->render();
     }
     }*/

    /*glLineWidth(5.0);
     glBegin( GL_LINES );
     for ( int r=0; r<m_rodRootMaterialFrame.size(); r++ )
     {
     glColor3d(1,1,1);
     BASim::Vec3d p0 = BASim::Vec3d(0,0,0);
     BASim::Vec3d p1 = p0 + m_rodRootMaterialFrame[r].m1;
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
        int numRods = rodData.size();
        for ( int r=0; r<numRods; r++ )
        {
            BASim::ElasticRod* rod = rodData[r]->rod;
            if ( rodData[r]->ALLprevVertexPositions.size() != 0 )
            {
                glColor3d(0,0,1);
                glBegin( GL_LINE_STRIP );
                int total = rodData[r]->ALLprevVertexPositions.size();
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
                int total = rodData[r]->ALLnextVertexPositions.size();
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
                int total = rodData[r]->ALLcurrVertexPositions.size();
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

bool Beaker::collisionMeshInitialised(const int i_collisionMeshIndex)
{
    CollisionMeshDataHashMap::iterator itr = m_collisionMeshDataHashMap.find(i_collisionMeshIndex);

    if (itr != m_collisionMeshDataHashMap.end())
    {
        return true;
    }

    return false;
}

void Beaker::initialiseCollisionMesh(TriangleMesh* i_collisionMesh, LevelSet* i_levelSet,
        ScriptingController* i_scriptingController, const int i_collisionMeshIndex)
{
    cout << "Beaker: Initialising collision mesh " << i_collisionMeshIndex << endl;

    CollisionMeshData* collisionMeshData = new CollisionMeshData(i_collisionMesh, i_levelSet, i_scriptingController);
    m_collisionMeshDataHashMap[i_collisionMeshIndex] = collisionMeshData;
    //m_collisionMeshHashMap[ i_collisionMeshIndex ]->initialize();
}

void Beaker::removeCollisionMesh(const int i_collisionMeshIndex)
{
    std::cout << "Removing collision mesh with id " << i_collisionMeshIndex << std::endl;

    m_collisionMeshDataHashMap.erase(i_collisionMeshIndex);
}
