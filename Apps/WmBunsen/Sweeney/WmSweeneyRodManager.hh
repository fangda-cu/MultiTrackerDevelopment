#ifndef WMSWEENEYRODMANAGER_HH_
#define WMSWEENEYRODMANAGER_HH_

#ifdef WETA
#include <weta/Wfigaro/Core/EigenIncludes.hh>
#include <weta/Wfigaro/Physics/World.hh>
#include <weta/Wfigaro/Core/ObjectControllerBase.hh>
#include <weta/Wfigaro/Render/RodRenderer.hh>
#include <weta/Wfigaro/Core/TriangleMesh.hh>
#include <weta/Wfigaro/Core/ScriptingController.hh>
#include <weta/Wfigaro/Physics/ElasticRods/BridsonStepper.hh>
#else
#include <BASim/src/Core/EigenIncludes.hh>
#include <BASim/src/Core/ObjectControllerBase.hh>
#include <BASim/src/Physics/World.hh>
#include <BASim/src/Render/RodRenderer.hh>
#include <BASim/src/Physics/ElasticRods/BridsonStepper.hh>
#endif

/** \class WmSweeneyRodManager
 * \brief A class to create and control a simulation involving one of more rods.
 * 
 * This class has methods to allow the creation of a simulation world, populate it with rods
 * and then control the rods with various forces.
 */
class WmSweeneyRodManager
{
public:
    WmSweeneyRodManager();
    ~WmSweeneyRodManager();
private:
    
};

#endif
