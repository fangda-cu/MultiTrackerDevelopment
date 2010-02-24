#ifndef RODDATA_HH_
#define RODDATA_HH_

#include <weta/Wfigaro/Physics/ElasticRods/ElasticRod.hh>
#include <weta/Wfigaro/Physics/ElasticRods/AnisotropicRod.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodCollisionTimeStepper.hh>
#include <weta/Wfigaro/Physics/ElasticRods/RodUtils.hh>
#include <weta/Wfigaro/Render/RodRenderer.hh>
#include "SplineAttrEval.hh"
#include <tr1/unordered_map>

using namespace BASim;

class MaterialFrame
{
public:
    Vec3d m1;
    Vec3d m2;
    Vec3d m3;
};

class KinematicEdgeData
{
public:
    MaterialFrame materialFrame;
    
    // This returns the rotation angle between materialframe.m2 and the passed in vector v.
    // It can be used to know how much the rod should be rotated to line up with this edge
    // when defining boundary conditions.
    double getRotationAngle( Vec3d v )
    {
        return acos( materialFrame.m2.dot( v ) );
    }
};

typedef std::tr1::unordered_map< unsigned int, KinematicEdgeData > KinematicEdgeDataMap;
    

// This class holds all the info from Maya needed to simulate and render a single rod.
class RodData
{
public:
    RodData();
    RodData( ElasticRod* i_rod, RodCollisionTimeStepper* i_stepper, RodRenderer* i_rodRenderer );
    ~RodData();
    
    
    // If for some reason this rod shouldn't be simulated then set this flag to false. This
    // usually happens when the user has set the rod node to playback from cache.
    bool shouldSimulate;
    
    // Should we keep the ObjectHandle returned by World rather than the actual rod?
    ElasticRod* rod;
    
    RodCollisionTimeStepper* stepper;
    RodRenderer* rodRenderer;
    
    // These variables are updated directly by the WmBunsenRodNode each frame.
    std::vector<Vec3d> undeformedVertexPositions;
    std::vector<Vec3d> initialVertexPositions;
    RodOptions rodOptions;

    // These variables are for interpolating substeps between frames.
    // The position of the input curve vertices at the last frame
    std::vector<Vec3d> prevVertexPositions;
    // The position of the input curve vertices at the current substep
    std::vector<Vec3d> currVertexPositions;
    // The position of the input curve vertices at the next frame
    std::vector<Vec3d> nextVertexPositions;
    // whether this vertex is locked in place or not
   // vector<Bool> vertexFixed;

    double hairSprayScaleFactor;
    cvDataMap forceWeightMap;
    double massDamping;

    // debug info
    std::vector<Vec3d> ALLprevVertexPositions;
    std::vector<Vec3d> ALLnextVertexPositions;
    std::vector<Vec3d> ALLcurrVertexPositions;

    // Used to tell Beaker whether we are giving it root frames (from rods) or nothing (nurbs)
    bool rootFrameDefined;
    
    // The undeformed material frames for this rod at startTime
    std::vector<MaterialFrame> undeformedMaterialFrame;
    
    // A list of frames for this rod. Each entry corresponds to an edge and a frame for that
    // edge. This is used for the first edge when Barbershop is driving the rods and for
    // other edges when the user is controlling the rod edge by a Maya transform.
    KinematicEdgeDataMap kinematicEdgeDataMap;
};

#endif
